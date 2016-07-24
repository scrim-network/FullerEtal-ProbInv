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
# Reference:
# Ruckert, KL, Shaffer, G, Pollard, D, Forest, FE, and Keller, K.
# Neglecting cliff instability in Antarctic ice sheet models can reduce
# melting during warming periods, (in prep.).
#
# This program is distributed in the hope that it will be useful,
# but WITH NO WARRANTY (NEITHER EXPLICIT NOR IMPLICIT). We are not liable
# for the behavior of these codes in your own application. You are free
# to share this code so long as the authors(s) and version history remain
# intact. 
#
# Kelsey Ruckert, klr324@psu.edu
# Yawen Guan, yig5031@psu.edu
#
# daisassim.R

source("assim.R")


dynReload("../fortran/dais", makevars='PKG_FCFLAGS="-I../fortran -J../fortran"',
    srcname=paste("../fortran/src/", c("dais.f90", "run_dais.f90", "global.f90"), sep=""))
source("daisF.R")
Volo <- 2.4789e16

F_daisModel <- function(iceflux, assimctx)
{
    Volume_F <- daisF(
        tstep = 1,
        b0    = iceflux[10],
        slope = iceflux[11],
        mu    = iceflux[3],
        h0    = iceflux[8],
        c     = iceflux[9],
        P0    = iceflux[5],
        kappa = iceflux[6],
        nu    = iceflux[4],
        f0    = iceflux[7],
        gamma = iceflux[1],
        alpha = iceflux[2],
        Tf    = -1.8,             #Freezing temperature of sea water
        rho_w = 1030,             #Density of sea water [g/cm^3]
        rho_i = 917,              #Density of ice water [g/cm^3]
        rho_m = 4000,             #Density of rock [g/cm^3]
        Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
        Rad0  = 1.8636e6,         #Steady state AIS radius for present day Ta and SL [m]
        Ta     = assimctx$frc[, 1], 
        SL     = assimctx$frc[, 4],
        Toc    = assimctx$frc[, 2],
        dSL    = assimctx$frc[, 3]
        )

    return (Volume_F)    
}


# allocate globally for efficiency
SLE <- Vais <- Rad <- Flow <- Depth <- numeric()

C_daisModel <- function(mp, assimctx)
{
    np <- nrow(assimctx$frc)
    if (np != length(Rad)) {
        SLE   <<- numeric(length=np)            # Sea-level equivalent [m]
        Vais  <<- numeric(length=np)            # Ice volume
        Rad   <<- numeric(length=np)            # Radius of ice sheet
        Flow  <<- numeric(length=np)            # Ice flow
        Depth <<- numeric(length=np)            # Water depth
    }

    mp <- c(
        mp,
        Tf    = -1.8,             #Freezing temperature of sea water
        rho_w = 1030,             #Density of sea water [g/cm^3]
        rho_i =  917,             #Density of ice water [g/cm^3]
        rho_m = 4000,             #Density of rock [g/cm^3]
        Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
        Rad0  = 1.8636e6          #Steady state AIS radius for present day Ta and SL [m]
    )

    .Call("daisOdeC", list(mp=mp, frc=assimctx$frc, out=list(SLE, Vais, Rad, Flow, Depth)))

    return (SLE)
}


iceflux <- function(mp, forcings)
{
    assimctx     <- list()
    assimctx$frc <- forcings

    return (daisModel(mp, assimctx))
}


if (!exists("daisassimctx")) {
    daisassimctx <- env()
}


daisLogLik <- function(mp, sp, assimctx)
{
    y.mod <- assimctx$modelfn(mp, assimctx)
    sigma.y <- sp["sigma"]
  
    #get the residuals
    # could pre-calculate median(assimctx$windows[n,])
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


daisConfigAssim <- function(alex=T, fortran=F, assimctx=daisassimctx)
{
    if (fortran) {
        daisModel <<- F_daisModel
    } else {
        daisModel <<- C_daisModel
        dynUnload("dais_alex")
        dynUnload("dais_kelsey")
        if (alex) {
            dynLoad("dais_alex",   srcname=c("dais_alex.c",   "r.c"), extrasrc="r.h")
        } else {
            dynLoad("dais_kelsey", srcname=c("dais_kelsey.c", "r.c"), extrasrc="r.h")
        }
    }
    assimctx$alex    <- alex
    assimctx$fortran <- fortran

    GSL <- scan("../../ruckert_dais/Data/future_GSL.txt", what=numeric(), quiet=T)  #Time rate of change of sea-level
    TA  <- scan("../../ruckert_dais/Data/future_TA.txt",  what=numeric(), quiet=T)  #Antarctic temp reduced to sea-level
    TO  <- scan("../../ruckert_dais/Data/future_TO.txt",  what=numeric(), quiet=T)  #High latitude subsurface ocean temp
    SL  <- scan("../../ruckert_dais/Data/future_SL.txt",  what=numeric(), quiet=T)  #Reconstructed sea-level

    project.forcings  <- matrix(c(TA,TO,GSL,SL), ncol=4, nrow=240300)
    hindcast.forcings <- matrix(c(TA[1:240010], TO[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

    paramNames <- c("gamma", "alpha", "mu",    "nu",                "P0", "kappa", "f0", "h0", "c", "b0", "slope")
    assimctx$units <- c("", "",     "m^(0.5)", "m^(-0.5) yr^(-0.5)", "m", "1/K", "m/yr", "m", "m/K", "m", "")

    # daisModel() uses this
    assimctx$frc <- hindcast.forcings

    # Best Case (Case #4) from Shaffer (2014)
    IP <- c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)
    names(IP) <- paramNames

    AIS_melt     <- iceflux(IP, hindcast.forcings)
    #Project_melt <- iceflux(IP, project.forcings)

    estimate.SLE.rate <- abs(-71/360)/1000
    time.years        <- 2002-1992
    mid.cum.SLE_2002  <- estimate.SLE.rate*time.years

    estimate.SLE.error <- sqrt(time.years)*abs(-53/360)/1000    #1-sigma error
    SE2_2002 <- estimate.SLE.error*2                            #2-sigma error

    positive_2SE <- mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
    negative_2SE <- mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

    upper.wind <- c(6.0, -6.9, -1.25, positive_2SE)
    lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
    assimctx$windows <- matrix(c(lower.wind, upper.wind), nrow=4, ncol=2)

    assimctx$obs.errs <- (assimctx$windows[,2]-assimctx$windows[,1])*.5

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
    bound.lower[1:2] <- c( 1/2, 0)
    bound.upper[1:2] <- c(17/4, 1)
    
    #Set bounds for bo and s
    bound.lower[10:11] <- c(725, 0.00045)
    bound.upper[10:11] <- c(825, 0.00075)

    assimctx$modelfn <- daisModel
    assimctx$lbound  <- bound.lower
    assimctx$ubound  <- bound.upper

    if (0) {
        # read in optimized parameters
        raw     <- scan("../../ruckert_dais/DAIS_matlab/OtimizedInitialParameters.txt", what=numeric(), quiet=T)
        init_p  <- matrix(raw, ncol=2, byrow=T)
        init_p  <- init_p[, 2]
    } else {
        #init_p <- c(2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.6) # Kelsey Random guesses
        init_p <- IP
    }
    init_sp          <- numeric()
    init_mp          <- init_p[1:11]
    init_sp["sigma"] <- sigma

    names(assimctx$lbound) <- names(assimctx$ubound) <- names(init_mp) <- paramNames

   #configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T, sigma_max=Inf)
    configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T)
}


daisRunAssim <- function(nbatch=1e6, adapt=T, assimctx=daisassimctx)
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (adapt) {
        scale <- abs(c(init_mp, init_sp) / 25)
    } else {
        scale <- c(0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15) / 5
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale, adapt=adapt)

    results <<- assimctx$chain
}


daisRunFit <- function(assimctx=daisassimctx, useDE=F)
{
    init_p <- assimMaxLikelihood(assimctx, init_mp=assimctx$init_mp, init_sp=assimctx$init_sp, useDE=useDE)
    assimctx$init_mp <- init_p[  assimctx$mp_indices ]
    assimctx$init_sp <- init_p[ -assimctx$mp_indices ]

    return (init_p)
}
