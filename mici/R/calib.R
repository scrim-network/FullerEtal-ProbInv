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
# calib.R

source("assim.R")
source("ts.R")
source("Scripts/plot_PdfCdfSf.R")


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
        Ta    = assimctx$frc[, 1], 
        SL    = assimctx$frc[, 4],
        Toc   = assimctx$frc[, 2],
        dSL   = assimctx$frc[, 3]
        )

    return (Volume_F)
}


F_daisFastDynModel <- function(iceflux, assimctx)
{
    Volume_F <- dais_fastdynF(
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
        Tcrit = iceflux[12],
        lambda = iceflux[13],
        Tf    = -1.8,             #Freezing temperature of sea water
        rho_w = 1030,             #Density of sea water [g/cm^3]
        rho_i = 917,              #Density of ice water [g/cm^3]
        rho_m = 4000,             #Density of rock [g/cm^3]
        Toc_0 = 0.72,             #Present day high latitude ocean subsurface temperature [K]
        Rad0  = 1.8636e6,         #Steady state AIS radius for present day Ta and SL [m]
        Ta    = assimctx$frc[, 1],
        SL    = assimctx$frc[, 4],
        Toc   = assimctx$frc[, 2],
        dSL   = assimctx$frc[, 3]
        )

    return (Volume_F$Vais)
}


# allocate globally for efficiency
SLE <- Vais <- Rad <- Flow <- Depth <- numeric()


C_daisModel <- function(mp, assimctx)
{
    n <- nrow(assimctx$frc)
    if (n != length(Rad)) {
        SLE   <<- numeric(length=n)             # Sea-level equivalent [m]
        Vais  <<- numeric(length=n)             # Ice volume
        Rad   <<- numeric(length=n)             # Radius of ice sheet
        Flow  <<- numeric(length=n)             # Ice flow
        Depth <<- numeric(length=n)             # Water depth
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

    .Call(assimctx$daisCmodel, list(mp=mp, frc=assimctx$frc, out=list(SLE, Vais, Rad, Flow, Depth), sw=assimctx$sw))

    return (SLE)
}


iceflux <- function(mp, forcings, assimctx=daisassimctx)
{
    assimlst            <- list()
    assimlst$frc        <- forcings

    # TODO:  models.R needs to grab assimctx$sw from daisassimctx if it's available
    assimlst$sw         <- assimctx$sw
    assimlst$daisCmodel <- assimctx$daisCmodel

    return (assimctx$modelfn(mp, assimlst))
}


if (!exists("daisassimctx")) {
    daisassimctx <- env()
}


daisLogLik <- function(mp, sp, assimctx)
{
    llik  <- 0
    resid <- error <- numeric()

    y.mod <- assimctx$modelfn(mp, assimctx)

    # paleo constraints
    if (assimctx$paleo) {
        resid <- append(    resid, assimctx$obsonly[1:3] - (y.mod[assimctx$obs_ind[1:3]] - mean(y.mod[assimctx$SL.1961_1990])))
        error <- append(    error, assimctx$error  [1:3])
    }

    # instrumental contraint
    resid <- append(        resid, assimctx$obsonly[4]   - (y.mod[assimctx$obs_ind[4]]   -      y.mod[assimctx$SL.1992]))
    error <- append(        error, assimctx$error  [4])

    # future expert assessment constraint
    if (assimctx$expert) {
        if (exists("expert_prior", env=assimctx)) {
            llik <- llik + assimctx$expert_prior$dens(      y.mod[assimctx$obs_ind[5]]   -      y.mod[assimctx$SL.2010])
        } else {
            resid <- append(resid, assimctx$obsonly[5]   - (y.mod[assimctx$obs_ind[5]]   -      y.mod[assimctx$SL.2010]))
            error <- append(error, assimctx$error  [5])
        }
    }
  
    # Calculate the likelihood. The observations are not correlated. They are independent. This makes the model heteroskedastic.
    llik <- llik + sum(dnorm(resid, sd=sqrt(sp["var"] + error^2), log=TRUE))
  
    return (llik)
}


daisLoadModel <- function(cModel="rob")
{
    if (is.null(cModel)) {
        dynReload("../fortran/dais",         srcname=paste("../fortran/src/",
            c("global.f90", "dais.f90",         "run_dais.f90"),         sep=""))
        dynReload("../fortran/dais_fastdyn", srcname=paste("../fortran/src/",
            c("global.f90", "dais_fastdyn.f90", "run_dais_fastdyn.f90"), sep=""))
    } else {
        daisLib <- paste(cModel, "_dais", sep="")
        dynReload(daisLib, srcname=c(paste(daisLib, ".c", sep=""), "r.c"), extrasrc="r.h")
    }
}


daisLoadModel()
daisLoadModel("kelsey")
daisLoadModel(NULL)
source("daisF.R")
source("dais_fastdynF.R")


# cModel can be either rob, kelsey, or NULL right now.  NULL selects the Fortran model.
daisConfigAssim <- function(cModel="rob", fast_dyn=F, rob_dyn=F, paleo=T, pfeffer=F, unif=F, pollard=F, assimctx=daisassimctx)
{
    # configure model to run
    #
    if (is.null(cModel)) {
        if (fast_dyn) {
            assimctx$modelfn <- F_daisFastDynModel
        } else {
            assimctx$modelfn <- F_daisModel
        }
    } else {
        assimctx$modelfn     <- C_daisModel
        assimctx$daisCmodel <- paste("dais", toupper(substring(cModel, 1, 1)), substring(cModel, 2), "OdeC", sep="")
        daisLoadModel(cModel)
    }
    assimctx$cModel <- cModel


    # set up forcings
    #
    GSL <- scan("../../../ruckert_dais/Data/future_GSL.txt", what=numeric(), quiet=T)  #Time rate of change of sea-level
    TA  <- scan("../../../ruckert_dais/Data/future_TA.txt",  what=numeric(), quiet=T)  #Antarctic temp reduced to sea-level
    TO  <- scan("../../../ruckert_dais/Data/future_TO.txt",  what=numeric(), quiet=T)  #High latitude subsurface ocean temp
    SL  <- scan("../../../ruckert_dais/Data/future_SL.txt",  what=numeric(), quiet=T)  #Reconstructed sea-level
    assimctx$forcings <- cbind( time=(1L:length(SL) - 238000L), TA=TA, TO=TO, GSL=GSL, SL=SL )
    assimctx$expert   <- pfeffer | pollard
    assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=ifelse(assimctx$expert, 2100, 2010))
    assimctx$frc      <- assimctx$frc_ts[ , 2:ncol(assimctx$forcings) ]


    # set up observations and observation errors
    #


    # instrumental period from Shepherd et al. (2012)
    #

    estimate.SLE.rate <- abs(-71/360)/1000
    time.years        <- 2002-1992
    mid.cum.SLE_2002  <- estimate.SLE.rate*time.years

    # (*sqrt(10) because 10 years of potentially accumulated error:
    #  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
    #                = 10*year X error^2)
    estimate.SLE.error <- sqrt(time.years)*abs(-53/360)/1000    #1-sigma error
    SE2_2002           <- estimate.SLE.error*2                  #2-sigma error

    positive_2SE <- mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
    negative_2SE <- mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value


    # Precal windows 1-3:
    # from Shaffer (2014). modified by Kelsey
    #

    upper.wind <- c(6.0, -6.9, -1.25, positive_2SE)
    lower.wind <- c(1.8, -15.8, -4.0, negative_2SE)
    assimctx$windows  <- cbind(lower.wind, upper.wind)
    assimctx$obsonly  <- rowMean(assimctx$windows)

    # the windows are +-2 sigma.  normal functions in R expect 1 sigma.  divide by 2 to get 2 sigma
    # and again to get 1 sigma
    assimctx$error    <- (assimctx$windows[, 2] - assimctx$windows[, 1]) / 4

    # Create a vector with each observation year
    #120kyr, 20Kyr, 6kyr, 2002
    assimctx$obs_ind       <- tsGetIndices(       assimctx$frc_ts, c(-118000, -18000, -4000, 2002))
    assimctx$SL.1961_1990  <- tsGetIndicesByRange(assimctx$frc_ts, lower=1961, upper=1990)
    assimctx$SL.1992       <- tsGetIndices(       assimctx$frc_ts, 1992)
    assimctx$SL.2010       <- tsGetIndices(       assimctx$frc_ts, 2010)


    # Pfeffer et al. (2008)
    #
    if (pfeffer) {
        assimctx$windows <- rbind(assimctx$windows, c(146/1000, 619/1000))
       #assimctx$windows <- rbind(assimctx$windows, c(128/1000, 619/1000))
    }


    if (assimctx$expert) {
        assimctx$obs_ind    [5] <- tsGetIndices(assimctx$frc_ts, 2100)
        if (unif) {
            assimctx$expert_prior <- uniformPrior(min=assimctx$windows[5, 1], max=assimctx$windows[5, 2])
        } else {
            rmif(expert_prior, envir=assimctx)
            assimctx$obsonly[5] <- mean(assimctx$windows[5, ])
            assimctx$error  [5] <- (assimctx$obsonly[5] - assimctx$windows[5, 1]) / 2  # treating as 2-sigma
        }
    }


    # set up model parameters and priors
    #

    paramNames <- c("gamma", "alpha", "mu",    "nu",                "P0", "kappa", "f0", "h0", "c", "b0", "slope")
    assimctx$units <- c("",  "",    "m^(0.5)", "m^(-0.5) yr^(-0.5)", "m", "1/K", "m/yr", "m", "m/K", "m", "")

    #                  c('gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c'  , 'b0','slope')
    assimctx$lbound <- c( 0.50,  0,     4.35, 0.006, 0.175,  0.02,  0.6,  735.5,  47.5, 725, 0.00045)
    assimctx$ubound <- c( 4.25,  1,    13.05, 0.018, 0.525,  0.06,  1.8, 2206.5, 142.5, 825, 0.00075)

    # Best Case (Case #4) from Shaffer (2014)
    init_mp         <- c(  2.0, 0.35,   8.7,  0.012, 0.35,   0.04,  1.2, 1471,    95,   775, 0.0006)

    if (fast_dyn) {
        paramNames      <- c(paramNames,     "Tcrit", "lambda")
        assimctx$units  <- c(assimctx$units,     "K",   "m/yr")
        assimctx$lbound <- c(assimctx$lbound,  -20.0,    0.005)
        assimctx$ubound <- c(assimctx$ubound,  -10.0,    0.015)
        init_mp         <- c(init_mp,          -15.0,    0.010)
    }
    if (rob_dyn) {
        paramNames      <- c(paramNames,     "Hcrit",     "fa")
        assimctx$units  <- c(assimctx$units,     "m",    "1/m")
        assimctx$lbound <- c(assimctx$lbound,  200.0, gtzero())
        assimctx$ubound <- c(assimctx$ubound, 2000.0,     10.0)
        init_mp         <- c(init_mp,          400.0,      0.5)
    }

    assimctx$sw             <- logical()
    assimctx$sw["fast_dyn"] <- fast_dyn
    assimctx$sw["rob_dyn"]  <- rob_dyn
    assimctx$paleo          <- paleo

    names(init_mp) <- names(assimctx$lbound) <- names(assimctx$ubound) <- paramNames


    # calculate variance (sigma^2) from residuals
    #
    AIS_melt <- assimctx$modelfn(init_mp, assimctx)
    resid    <- assimctx$obsonly - (AIS_melt[assimctx$obs_ind] - mean(AIS_melt[assimctx$SL.1961_1990]))
    init_sp        <- numeric()
    init_sp["var"] <- sd(resid)^2


    # configure assimilation engine
   #configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T, var_max=Inf)
    configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T, var_max=100)
}


daisRunAssim <- function(nbatch=ifelse(adapt, 5e5, 4e6), adapt=T, assimctx=daisassimctx)
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (adapt) {
        scale <- abs(c(init_mp, init_sp) / 25)
    } else {
        if (assimctx$sw["fast_dyn"]) {
            scale <- c(0.1, 0.015, 0.2, 0.035, 0.1, 0.01, 0.1, 50, 10, 25, 0.0005,  0.5, 0.0005, 0.1) / 5
        } else if (assimctx$sw["rob_dyn"]) {
            scale <- c(0.1, 0.015, 0.2, 0.035, 0.1, 0.01, 0.1, 50, 10, 25, 0.0005, 25.0, 0.1,    0.1) / 5
        } else {
            scale <- c(0.1, 0.015, 0.2, 0.035, 0.1, 0.01, 0.1, 50, 10, 25, 0.0005,               0.1) / 5
        }
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale, adapt=adapt)

    results <<- assimctx$chain
}


daisRunFit <- function(useDE=F, assimctx=daisassimctx)
{
    init_p <- assimMaxLikelihood(assimctx, init_mp=assimctx$init_mp, init_sp=assimctx$init_sp, useDE=useDE)
    assimctx$init_mp <- init_p[  assimctx$mp_indices ]
    assimctx$init_sp <- init_p[ -assimctx$mp_indices ]

    return (init_p)
}


daisRunPredict <- function(nbatch=3500, endYear=2300, assimctx=daisassimctx)
{
    # Identify the burn-in period and subtract it from the chains.
    chain <- assimctx$chain[ burnedInd(assimctx$chain), ]
    DAIS_chains_burnin <<- chain

    # Find mean of the estimated parameters.
    mean.parameters <<- colMean(chain)
    print(mean.parameters)

    # Estimate mean bias.
    bias.mean <<- sqrt(mean.parameters["var"])

    #
    # Make projections.
    #
    frc   <- tsTrim(assimctx$forcings, endYear=endYear)[ , 2:ncol(assimctx$forcings) ]
    years <- nrow(frc)
    proj.mcmc.anomaly   <<- matrix(nrow=nbatch, ncol=years)
    proj.mcmc.1961_1990 <<- matrix(nrow=nbatch, ncol=years)

    # Sample from the chain.
    par.mcmc <- sampleChain(chain, nbatch)
    samples  <- 1:nbatch
    while (T) {
        for (i in samples) {

            # Run the model.
            sle   <- iceflux(par.mcmc[i, assimctx$mp_indices], frc, assimctx)

            # Standardize the anomaly.
            anom  <- sle - mean(sle[assimctx$SL.1961_1990])
            proj.mcmc.anomaly  [i, ] <<- anom

            # Add noise.
            proj.mcmc.1961_1990[i, ] <<- anom + rnorm(years, sd=sqrt(par.mcmc[i, "var"]))
        }

        # look for NaNs (non-finite)
        samples <- which(apply(proj.mcmc.1961_1990, MARGIN=1, FUN=function(x) { any(!is.finite(x)) }))
        if (!length(samples)) {
            break;
        }
        print(c("resampling", length(samples), "non-finite rows"))
        par.mcmc[samples, ] <- sampleChain(chain, length(samples))
    }

    #--------------------- Estimate PDFs, CDFs, and SFs in certain years --------------------------
    # Function to find SLE values in certain years 'fn.prob.proj'
    year.pcs <<- tsGetIndices(assimctx$forcings, c(-118000, -18000, -4000, 2002, 2050, 2100, 2300))

    mcmc.prob_proj <<- fn.prob.proj(proj.mcmc.1961_1990, year.pcs, nbatch, un.constr=T)

    # Calculate the pdf, cdf, and sf of AIS melt estimates in:
    LIG.sf.mcmc  <<- plot.sf(mcmc.prob_proj[,1], make.plot=F) # 120,000 BP (Last interglacial)
    LGM.sf.mcmc  <<- plot.sf(mcmc.prob_proj[,2], make.plot=F) # 20,000 BP (Last glacial maximum)
    MH.sf.mcmc   <<- plot.sf(mcmc.prob_proj[,3], make.plot=F) # 6,000 BP (Mid-holocene)
    sf.2002.mcmc <<- plot.sf(mcmc.prob_proj[,4], make.plot=F) # 2002 (Observed trend from 1993-2011)
    sf.2050.mcmc <<- plot.sf(mcmc.prob_proj[,5], make.plot=F) # 2050
    sf.2100.mcmc <<- plot.sf(mcmc.prob_proj[,6], make.plot=F) # 2100
    sf.2300.mcmc <<- plot.sf(mcmc.prob_proj[,7], make.plot=F) # 2300

    # Find out parameter relationships; set up a matrix
    d.pos_parameters <<- par.mcmc
    sub_chain        <<- par.mcmc

    # random stuff MCMC_plots.R uses
    enddate          <<- years
    mean.dais.par    <<- mean.parameters
    project.forcings <<- frc
    standards        <<- NULL
    windows          <<- assimctx$windows
    bound.lower      <<- assimctx$lbound
    bound.upper      <<- assimctx$ubound
    subset_length    <<- nbatch
    obs.years        <<- assimctx$obs_ind
}
