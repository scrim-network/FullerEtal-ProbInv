# Copyright 2016 Robert W. Fuller <hydrologiccycle@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# daisassim.R

source("assim.R")
source("roblib.R")
dynReload("dais", srcname=c("dais.c", "r.c"), extrasrc="r.h")


daisModel <- function(mp, assimctx)
{
    np     <- nrow(assimctx$frc)
    Rad    <- numeric(length=np)               # Radius of ice sheet
    Vais   <- numeric(length=np)               # Ice volume
    SLE    <- numeric(length=np)               # Sea-level equivalent [m]

    mp <- c(mp,
      Tf    = -1.8,             #Freezing temperature of sea water
      rho_w = 1.03*1000,        #Density of sea water [g/cm^3]
      rho_i = 0.917*1000,       #Density of ice water [g/cm^3]
      rho_m = 4.0*1000,         #Density of rock [g/cm^3]
      Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
      Rad0  = 1.8636e6          #Steady state AIS radius for present day Ta and SL [m]
      )

    .Call("daisOdeC", list(mp=mp, frc=assimctx$frc, out=list(Rad, Vais, SLE)), PACKAGE="dais")

    return (SLE)
}


iceflux <- function(iceflux, forcings)
{
    mp <- c(
      b0    = iceflux[10],
      slope = iceflux[11],
      mu    = iceflux[3],
      h0    = iceflux[8],
      c     = iceflux[9],
      P0    = iceflux[5],
      kappa = iceflux[6],
      nu    = iceflux[4],
      f0    = iceflux[7],
      Gamma = iceflux[1],
      alpha = iceflux[2],
      Tf    = -1.8,             #Freezing temperature of sea water
      rho_w = 1.03*1000,        #Density of sea water [g/cm^3]
      rho_i = 0.917*1000,       #Density of ice water [g/cm^3]
      rho_m = 4.0*1000,         #Density of rock [g/cm^3]
      Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
      Rad0  = 1.8636e6          #Steady state AIS radius for present day Ta and SL [m]
    )

    np     <- nrow(forcings)
    Rad    <- numeric(length=np)               # Radius of ice sheet
    Vais   <- numeric(length=np)               # Ice volume
    SLE    <- numeric(length=np)               # Sea-level equivalent [m]

    .Call("daisOdeC", list(mp=mp, frc=forcings, out=list(Rad, Vais, SLE)), PACKAGE="dais")

    return (SLE)
}


if (!exists("daisassimctx")) {
    daisassimctx <- env()
}


daisLogLik <- function(mp, sp, assimctx)
{
    y.mod <- assimctx$modelfn(mp, assimctx)
    sigma.y <- sp["sigma"]
  
    #get the residuals
    r1 <- median(assimctx$windows[1,]) - (y.mod[120000] - mean(y.mod[assimctx$SL.1961_1990]))
    r2 <- median(assimctx$windows[2,]) - (y.mod[220000] - mean(y.mod[assimctx$SL.1961_1990]))
    r3 <- median(assimctx$windows[3,]) - (y.mod[234000] - mean(y.mod[assimctx$SL.1961_1990]))
    r4 <- median(assimctx$windows[4,]) - (y.mod[240002] - mean(y.mod[assimctx$SL.1961_1990]))

    resid.y <- c(r1, r2, r3, r4)
    sterr.y <- assimctx$obs.errs #This makes the model heteroskedastic
  
    #Calculate the likelihood. The observations are not correlated. They are independent
    llik <- sum(dnorm(resid.y, sd=sqrt(sigma.y + sterr.y^2), log=TRUE))
  
    return (llik)
}


daisConfigAssim <- function(
    assimctx=daisassimctx
    )
{
    GSL <- scan("../../rucker_dais/Data/future_GSL.txt", what=numeric(), quiet=T)  #Time rate of change of sea-level
    TA  <- scan("../../rucker_dais/Data/future_TA.txt",  what=numeric(), quiet=T)  #Antarctic temp reduced to sea-level
    TO  <- scan("../../rucker_dais/Data/future_TO.txt",  what=numeric(), quiet=T)  #High latitude subsurface ocean temp
    SL  <- scan("../../rucker_dais/Data/future_SL.txt",  what=numeric(), quiet=T)  #Reconstructed sea-level

    project.forcings  <- matrix(c(TA,TO,GSL,SL), ncol=4, nrow=240300)
    hindcast.forcings <- matrix(c(TA[1:240010], TO[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

    # daisModel() uses this
    assimctx$frc <- hindcast.forcings

    # Best Case (Case #4) from Shaffer (2014)
    IP <- c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

    AIS_melt     <- iceflux(IP, hindcast.forcings)
    #Project_melt <- iceflux(IP, project.forcings)

    estimate.SLE.rate <- abs(-71/360)/1000
    time.years <- 2002-1992
    mid.cum.SLE_2002 <- estimate.SLE.rate*time.years

    estimate.SLE.error <- abs(-53/360)/1000     #1- sigma error
    SE2_2002 <- estimate.SLE.error*2            #2-sigma error

    positive_2SE <- mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
    negative_2SE <- mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

    upper.wind <- c(6.0, -6.9, -1.25, positive_2SE)
    lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
    assimctx$windows <- matrix(c(lower.wind, upper.wind), nrow=4, ncol=2)

    assimctx$obs.errs <- c(abs(median(assimctx$windows[1,])-assimctx$windows[1,1]),
                           abs(median(assimctx$windows[2,])-assimctx$windows[2,1]),
                           abs(median(assimctx$windows[3,])-assimctx$windows[3,1]),
                           SE2_2002)

    # Create a vector with each observation year
    #120kyr, 20Kyr, 6kyr, 2002
    obs.years <- c(120000, 220000, 234000, 240002)

    #Set up equation to find the residuals and then the prior sigma
    assimctx$SL.1961_1990 <- 239961:239990
    resid <- rep(NA, length(obs.years))         # Create a vector of the residuals
    for (i in 1:length(obs.years)) {
        resid[i] <- median(assimctx$windows[i,]
                  - (AIS_melt[obs.years[i]] - mean(AIS_melt[assimctx$SL.1961_1990])))
                    #/sd(assimctx$windows[i,])
	}
    sigma <- sd(resid)^2                        #calculate the variance (sigma^2)

    bound.lower <- IP - (IP*0.5)
    bound.upper <- IP + (IP*0.5)

    #Set bounds for gamma and alpha
    bound.lower[1:2] <- c(1/2,  0)
    bound.upper[1:2] <- c(17/4, 1)
    
    #Set bounds for bo and s
    bound.lower[10:11] <- c(725, 0.00045)
    bound.upper[10:11] <- c(825, 0.00075) 

    assimctx$modelfn <- daisModel
    assimctx$lbound  <- bound.lower
    assimctx$ubound  <- bound.upper

    # TODO:  set assimctx$units

    # read in optimized parameters
    raw     <- scan("../../rucker_dais/DAIS_matlab/OtimizedInitialParameters.txt", what=numeric(), quiet=T)
    init_p  <- matrix(raw, ncol=2, byrow=T)
    init_mp <- init_p[1:11, 2]
    init_sp <- numeric()

    # from matlab, minit(12)=[1.16], BUT upper bound is set to 1.0 on sigma!  so 1.16 is improbable
    init_sp["sigma"] <- 1.0
    #init_sp["sigma"] <- init_p[12, 2]

    print(init_mp)
    print(init_sp)

    names(assimctx$lbound) <- names(assimctx$ubound) <- names(init_mp) <-
        c("Gamma", "alpha", "mu", "nu", "P0", "kappa", "f0", "h0", "c", "b0", "slope")

    configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T)
}


daisRunAssim <- function(
    nbatch=1000,
    initial=is.null(assimctx$chain),
    assimctx=daisassimctx
    )
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (initial || ncol(assimctx$chain) != length(init_mp) + length(init_sp)) {
        print("using initial scale")
        scale <- abs(c(init_mp, init_sp) / 150)
    } else {
        print("using proposal matrix")
        scale <- assimProposalMatrix(assimctx$chain, mult=0.5)
    }

    runAssim(assimctx, nbatch, scale)

    results <<- assimctx$chain
}
