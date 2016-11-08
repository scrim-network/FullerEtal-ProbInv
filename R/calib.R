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
#source("Scripts/plot_PdfCdfSf.R")  # for Kelsey
#source("plot.R")  # pdfPlots()
loadLibrary("KernSmooth")  # bkde() in roblib.R


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
        tstep  = 1,
        b0     = iceflux[10],
        slope  = iceflux[11],
        mu     = iceflux[3],
        h0     = iceflux[8],
        c      = iceflux[9],
        P0     = iceflux[5],
        kappa  = iceflux[6],
        nu     = iceflux[4],
        f0     = iceflux[7],
        gamma  = iceflux[1],
        alpha  = iceflux[2],
        Tcrit  = iceflux[12],
        lambda = iceflux[13],
        Tf     = -1.8,             #Freezing temperature of sea water
        rho_w  = 1030,             #Density of sea water [g/cm^3]
        rho_i  = 917,              #Density of ice water [g/cm^3]
        rho_m  = 4000,             #Density of rock [g/cm^3]
        Toc_0  = 0.72,             #Present day high latitude ocean subsurface temperature [K]
        Rad0   = 1.8636e6,         #Steady state AIS radius for present day Ta and SL [m]
        Ta     = assimctx$frc[, 1],
        SL     = assimctx$frc[, 4],
        Toc    = assimctx$frc[, 2],
        dSL    = assimctx$frc[, 3]
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


iceflux <- function(mp, forcings, assimctx=daisctx)
{
    assimlst            <- list()
    assimlst$frc        <- forcings

    # TODO:  models.R needs to grab assimctx$sw from daisctx if it's available
    assimlst$sw         <- assimctx$sw
    assimlst$daisCmodel <- assimctx$daisCmodel

    return (assimctx$modelfn(mp, assimlst))
}


if (!exists("daisctx")) {
    daisctx <- env()
}


daisLogLik <- function(mp, sp, assimctx)
{
    llik  <- 0
    resid <- error <- numeric()

    y <- assimctx$modelfn(mp, assimctx)

    # paleo constraints
    if (assimctx$paleo) {
        resid <- append(resid, assimctx$obsonly[1:3] - (y[assimctx$obs_ind[1:3]]    - mean(y[assimctx$SL.1961_1990])))

        # constrain the LIG to 2-sigma
        lig <- resid[1]
        if (!is.finite(lig) || abs(lig) > 2*assimctx$error[1]) {
            return (-Inf)
        }

        error <- append(error, assimctx$error  [1:3])
    }

    # instrumental constraints
    if (assimctx$instrumental) {
        resid <- append(resid, assimctx$obsonly[4]   - (y[assimctx$obs_ind[4]]      -      y[assimctx$SL.1992]))
        error <- append(error, assimctx$error  [4])

        # IPCC rates
        resid <- append(resid, assimctx$trends_obs   - (y[assimctx$trends_ind[, 2]] -      y[assimctx$trends_ind[, 1]])
                                                     / (  assimctx$trends_ind[, 2]  -        assimctx$trends_ind[, 1]))
        error <- append(error, assimctx$trends_err)
    }

    # future expert assessment constraint
    if (assimctx$expert) {
        y_std <- y[assimctx$obs_ind[assimctx$expert_ind]] - y[assimctx$SL.expert]
        if (exists("expert_prior", env=assimctx)) {
            llik <- llik + assimctx$expert_prior$dens(                     y_std)
        } else {
            resid <- append(resid, assimctx$obsonly[assimctx$expert_ind] - y_std)
            error <- append(error, assimctx$error  [assimctx$expert_ind])
        }
        assimctx$y <- y_std
    }

    if (length(resid)) {
        var <- sp["var"]
        if (is.na(var)) {

            # no variance in chain, just use observation error
            llik <- llik + sum(dnorm(resid, sd=error, log=TRUE))
        } else {

            # Kelsey says, "The observations are not correlated. They are independent. This makes the model heteroskedastic."
            # TODO:  Tony does not include variance in IPCC rates
            llik <- llik + sum(dnorm(resid, sd=sqrt(sp["var"] + error^2), log=TRUE))
        }
    }
  
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
daisConfigAssim <- function(
    cModel="rob", fast_dyn=T, rob_dyn=F, instrumental=F, paleo=F, expert="pfeffer", prior="uniform", variance=F, assimctx=daisctx)
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
        assimctx$daisCmodel <- paste("dais", capitalize(cModel), "OdeC", sep="")
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
    assimctx$expert   <- !is.null(expert)

    # Kelsey runs the model to 2010, but her likelihood function only needs to be run through 2002;
    # using 2010 preserves the ability to run DAIScali_hetero_model_iid_mcmcmat.R;
    # nonetheless, now use 2011 in order to include the IPCC rates
    #
   #assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=ifelse(assimctx$expert, 2100, 2010))
    assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=ifelse(assimctx$expert, 2100, 2011))
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
   #upper.wind <- c(6.0,  -6.9, -1.25, positive_2SE)
   #lower.wind <- c(1.8, -15.8, -4.0,  negative_2SE)
    upper.wind <- c(7.4,  -6.9, -1.25, positive_2SE) # Windows 2-3 from Kelsey, Window 1 from DeConto and Pollard 2016
    lower.wind <- c(3.6, -15.8, -4.0,  negative_2SE)

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


    ## A fifth window is added to match IPCC AR5 Ch13 (page 1151) AIS SLR trend:
    ## 0.27 +/- 0.11 mm/year (convert to m/year here)
    ##
    ## Precal windows 5-?:
    ## Last "precalibration window" is 1993-2010 mean trend, from the IPCC AR5 Ch13
    ## (Page 1151), for AIS SLR contribution: 0.27 +- 0.11 mm/year
    ## Note that model output is in meters SLE and these trends are mm, so a
    ## conversion is necessary.
    assimctx$trends_obs <- c(0.27, 0.08,  0.40)  / 1000
    assimctx$trends_err <- c(0.11, 0.185, 0.205) / 1000
    assimctx$trends_ind <- matrix(ncol=2, nrow=length(assimctx$trends_obs), byrow=T, data=tsGetIndices(assimctx$frc_ts,
                           c(1993, 2010, 1992, 2001, 2002, 2011)))


    # Expert assessment
    #

    rmif(expert_prior, envir=assimctx)  # keep it clean

    if (assimctx$expert) {
        switch (expert,
            pfeffer={
                #
                # Pfeffer et al. (2008)
                #
               #assimctx$windows <- rbind(assimctx$windows, c(146/1000, 619/1000))
                assimctx$windows <- rbind(assimctx$windows, c(128/1000, 619/1000))
                assimctx$expert_std_yr <- 2010
            },
            pollard={
                stop("pollard unimplemented in daisConfigAssim()")
            }, {
                stop("unknown expert in daisConfigAssim()")
            })

        assimctx$expert_ind <- nrow(assimctx$windows)  # this might pollute the environment, but should be harmless
        assimctx$obsonly[assimctx$expert_ind] <- mean(assimctx$windows[assimctx$expert_ind, ])
        assimctx$error  [assimctx$expert_ind] <- (assimctx$obsonly[assimctx$expert_ind] - assimctx$windows[assimctx$expert_ind, 1]) / 2  # treating as 2-sigma
        assimctx$obs_ind[assimctx$expert_ind] <- tsGetIndices(assimctx$frc_ts, 2100)
        assimctx$SL.expert                    <- tsGetIndices(assimctx$frc_ts, assimctx$expert_std_yr)

        switch (prior,
            uniform={
                assimctx$expert_prior <- uniformPrior(min=assimctx$windows[assimctx$expert_ind, 1], max=assimctx$windows[assimctx$expert_ind, 2])
            },
            beta={
                # a=2, b=3 taken from Lempert, Sriver, and Keller (2012)
                assimctx$expert_prior <- betaPrior(   min=assimctx$windows[assimctx$expert_ind, 1], max=assimctx$windows[assimctx$expert_ind, 2], a=2, b=3)
            },
            normal={
                ;
            }, {
                stop("unknown prior in daisConfigAssim()")
            })
    }


    # set up model parameters and priors
    #

    paramNames <- c("gamma", "alpha", "mu",    "nu",                "P0", "kappa", "f0", "h0", "c", "b0", "slope")
    assimctx$units <- c("",  "",    "m^(0.5)", "m^(-0.5) yr^(-0.5)", "m", "1/K", "m/yr", "m", "m/K", "m", "")

    #                  c('gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c'  , 'b0','slope')
   #assimctx$lbound <- c( 0.50,  0,     4.35, 0.006, 0.175,  0.02,  0.6,  735.5,  47.5, 725, 0.00045)  # Kelsey priors
   #assimctx$ubound <- c( 4.25,  1,    13.05, 0.018, 0.525,  0.06,  1.8, 2206.5, 142.5, 825, 0.00075)
    assimctx$lbound <- c( 0.50 , 0,     7.05, 0.003, 0.026,  0.025, 0.6,  735.5,  47.5, 740, 0.00045)  # Tony priors
    assimctx$ubound <- c( 4.25 , 1,    13.65, 0.015, 1.5,    0.085, 1.8, 2206.5, 142.5, 820, 0.00075)


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
    assimctx$instrumental   <- instrumental
    assimctx$prior_name     <- prior
    assimctx$expert_name    <- expert

    names(init_mp) <- names(assimctx$lbound) <- names(assimctx$ubound) <- names(assimctx$units) <- paramNames


    # calculate variance (sigma^2) from residuals;  note that it's at best
    # approximate since not everything is standardized to 1961-1990
    #
    init_sp <- numeric()
    if (variance) {
        AIS_melt <- assimctx$modelfn(init_mp, assimctx)
        resid    <- assimctx$obsonly - (AIS_melt[assimctx$obs_ind] - mean(AIS_melt[assimctx$SL.1961_1990]))
        init_sp["var"] <- sd(resid)^2
        assimctx$units <- c(assimctx$units, "")  # add units for variance
    }


    # configure assimilation engine
   #configAssim(assimctx, init_mp, init_sp, ar=0, obserr=F, llikfn=daisLogLik, gamma_pri=T, var_max=Inf)
    configAssim(assimctx, init_mp, init_sp, ar=0, obserr=!variance, llikfn=daisLogLik, gamma_pri=T, var_max=100)
}


daisRunAssim <- function(nbatch=ifelse(adapt, 5e5, 4e6), adapt=T, assimctx=daisctx)
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
        if (!length(init_sp)) {
            scale <- scale[ 1:(length(scale) - 1) ]
        }
    }

    assimctx$ychain <- prmatrix(nbatch + ifelse(adapt, 0, 1), xvals=2100)

    runAssim(assimctx, nbatch=nbatch, scale=scale, adapt=adapt, extrafun=assimSaveY)

    #results <<- assimctx$chain
}


daisRunFit <- function(useDE=F, assimctx=daisctx)
{
    init_p <- assimMaxLikelihood(assimctx, init_mp=assimctx$init_mp, init_sp=assimctx$init_sp, useDE=useDE)
    assimctx$init_mp <- init_p[  assimctx$mp_indices ]
    assimctx$init_sp <- init_p[ -assimctx$mp_indices ]

    return (init_p)
}


if (!exists("prdaisctx")) {
    prdaisctx <- env()
}


# note that nbatch is not used here
daisRunPredict <- function(nbatch=3500, subsample=T, assimctx=daisctx, prctx=prdaisctx)
{
    prctx$assimctx <- assimctx

    prchain       <- assimFixOutput(assimctx, assimctx$ychain)
    burnIn        <- burnedInd(prchain)
    if (subsample) {
        samples   <- sample(burnIn, nbatch, replace=T)
    } else {
        samples   <- burnIn
    }
    prctx$prchain <- prmatrix(length(samples), xvals=attr(prchain, "xvals"))
    prctx$prchain[, 1:ncol(prchain)] <- prchain[samples, 1:ncol(prchain)]
}


daisRejSample <- function(prior=assimctx$expert_prior, assimctx=daisctx, prctx=prdaisctx)
{
    column <- as.character(2100)

    daisRunPredict(subsample=F, assimctx=assimctx, prctx=prctx)

    yvals    <- prctx$prchain
    burn_ind <- burnedInd(assimctx$chain)
    chain    <- assimctx$chain[ burn_ind, ]
    rej_ind  <- rejectSample(yvals[, column], tgt_dense_fn=assimctx$expert_prior$dens)
    len_ind  <- length(which(rej_ind))

    new_yvals                  <- prmatrix(len_ind, xvals=attr(yvals, "xvals"))
    new_yvals[, 1:ncol(yvals)] <- yvals[rej_ind, 1:ncol(yvals)]

    new_chain                  <- matrix(nrow=len_ind, ncol=ncol(chain))
    new_chain[, 1:ncol(chain)] <- chain[rej_ind, 1:ncol(chain)]
    colnames(new_chain)        <- colnames(chain)

    # save original chains
    prctx$prNoRejChain  <- assimFixOutput(assimctx, assimctx$ychain)
    assimctx$noRejChain <- assimctx$chain

    # save rejection sampled chains
    prctx$prchain       <- new_yvals
    assimctx$chain      <- new_chain

    print(paste("rejection sampling reduced rows from ", nrow(yvals), " to ", nrow(new_yvals), " (ratio=", format(nrow(yvals) / nrow(new_yvals), digits=3), ")", sep=""))

    #dev.new()
    #pdfPlots(new_yvals, column=column, burnin=F, col="black", lty="solid", legendloc=NULL)
}


daisRunPredictSlow <- function(nbatch=3500, endYear=2100, assimctx=daisctx, prctx=prdaisctx)
{
    prctx$assimctx <- assimctx

    print(colMean(assimctx$chain))

    frc_ts <- tsTrim(assimctx$forcings, endYear=endYear)
    frc    <- frc_ts[ , 2:ncol(assimctx$forcings) ]
   #prctx$prchain <- prmatrix(nbatch, xvals=frc_ts[, "time"])
    prctx$prchain <- prmatrix(nbatch, xvals= assimctx$expert_std_yr : endYear )

    years   <- nrow(frc)
    rows    <- 1:nbatch
    samples <- sample(burnedInd(assimctx$chain), nbatch, replace=T)
    while (T) {
        for (i in rows) {

            # Run the model.
            sle   <- iceflux(assimctx$chain[samples[i], assimctx$mp_indices], frc, assimctx)

            # Standardize the anomaly.
            anom  <- sle - sle[assimctx$SL.expert]

            # Add noise.  Or not.
           #prctx$prchain[i, ] <- anom # + rnorm(years, sd=sqrt(assimctx$chain[samples[i], "var"]))
            prctx$prchain[i, ] <- anom[ assimctx$SL.expert : length(anom) ]
        }

        # look for NaNs (non-finite)
        rows <- which(apply(prctx$prchain, MARGIN=1, FUN=function(x) { any(!is.finite(x)) }))
        if (!length(rows)) {
            break;
        }
        print(c("resampling", length(rows), "non-finite rows"))
        samples[rows] <- sample(burnedInd(assimctx$chain), length(rows), replace=T)
    }

    # reduce save file size
   #prTrimChains(prctx=prctx, lower=assimctx$expert_std_yr)
}


daisRunKelseyPredict <- function(nbatch=3500, endYear=2300, assimctx=daisctx)
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
