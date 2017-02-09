# Copyright (C) 2016, 2017 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# Parts of this code are also Copyright 2016 Tony Wong, twong@psu.edu
#
# calib.R

source("assim.R")
source("ts.R")
source("lhs.R")
#source("Scripts/plot_PdfCdfSf.R")  # for Kelsey
loadLibrary("KernSmooth")  # bkde() in roblib.R


F_daisModel <- function(iceflux, assimctx)
{
    iceflux <- c(assimctx$ep, iceflux)

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
    iceflux <- c(assimctx$ep, iceflux)

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

    .Call(assimctx$daisCmodel, list(mp=mp, ep=assimctx$ep, frc=assimctx$frc, out=list(SLE, Vais, Rad, Flow, Depth), sw=assimctx$sw))

    return (SLE)
}


if (!exists("daisctx")) {
    daisctx <- env()
}


daisLogPri <- function(mp, sp, assimctx)
{
    lpri <- 0

    # bounded uniform priors
    inBounds <- all(
        mp >= assimctx$lbound,    mp <= assimctx$ubound,
        sp >= assimctx$lbound_sp, sp <= assimctx$ubound_sp
        )
    if (!inBounds) {
        return (-Inf)  # zero prior probability
    }

    # inverse gamma prior for paleo variance
    var <- sp["var.paleo"]
    if (!is.na(var)) {
        lpri <- lpri + (-assimctx$alpha - 1) * log(var) + (-assimctx$beta / var)
    }

    # priors for non-uniform model parameters
    if (assimctx$gamma_pri) {
        lpri <- (lpri
              +  assimctx$lambda_prior$dens( mp["lambda"])
              +  assimctx$ Tcrit_prior$dens(-mp["Tcrit"]))  # negative requires multiplication of Tcrit by -1
    }

    return (lpri)
}


daisLogLik <- function(mp, sp, assimctx)
{
    llik  <- 0
    resid <- error <- numeric()

    y <- assimctx$modelfn(mp, assimctx)

    # for LHS sampling, always produce model output, even if the LF assigns zero probability
    if (assimctx$all_predict) {
        assimctx$y <- y - y[assimctx$SL.expert]
        y_std <- assimctx$y[assimctx$SL.2100]
    } else {
        y_std <- y[assimctx$SL.2100] - y[assimctx$SL.expert]
        assimctx$y <- y_std
    }

    # Heaviside likelihood function
    if (assimctx$heaviside) {
        resid <- assimctx$gmsl - (y[assimctx$SL.1900_201x] - mean(y[assimctx$SL.1880_1899]))
        if (!all(is.finite(resid)) || any(resid < 0)) {
            return (-Inf)
        }
    }

    # paleo constraints
    if (assimctx$paleo) {
        resid <- assimctx$obsonly[1:3] - (y[assimctx$obs_ind[1:3]] - mean(y[assimctx$SL.1961_1990]))

        # constrain the LIG to 2-sigma
        lig <- resid[1]
        if (!is.finite(lig) || abs(lig) > 2*assimctx$error[1]) {
            return (-Inf)
        }

        var <- sp["var.paleo"]
        if (is.na(var)) {
            llik <- llik + sum(dnorm(resid, sd=           assimctx$error[1:3]   , log=TRUE))
        } else {
            llik <- llik + sum(dnorm(resid, sd=sqrt(var + assimctx$error[1:3]^2), log=TRUE))
        }
    }

    # instrumental constraints
    if (assimctx$instrumental) {
        resid <- assimctx$obsonly[4] - (y[assimctx$obs_ind[4]] - y[assimctx$SL.1992])
        var <- sp["var.inst"]
        if (is.na(var)) {
            llik <- llik + dnorm(resid, sd=           assimctx$error[4]  ,  log=TRUE)
        } else {
            llik <- llik + dnorm(resid, sd=sqrt(var + assimctx$error[4]^2), log=TRUE)
        }

        # IPCC rates
        resid <- assimctx$trends_obs - (  (y[assimctx$trends_ind[, 2]] - y[assimctx$trends_ind[, 1]])
                                        / (  assimctx$trends_ind[, 2]  -   assimctx$trends_ind[, 1]))
        llik  <- llik + sum(dnorm(resid, sd=assimctx$trends_err, log=TRUE))
    }

    # future expert assessment constraint
    if (assimctx$expert) {
        llik <- llik + assimctx$expert_prior$dens(y_std)
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
    cModel="rob", fast_dyn=T, rob_dyn=F, fast_only=F, wide_prior=T,
    instrumental=F, paleo=F, expert="pfeffer", prior="uniform",
    gamma_pri=T, variance=T, all_predict=F, heaviside=F, assimctx=daisctx)
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
    GSL <- scan("../data/future_GSL.txt", what=numeric(), quiet=T)  #Time rate of change of sea-level
    TA  <- scan("../data/future_TA.txt",  what=numeric(), quiet=T)  #Antarctic temp reduced to sea-level
    TO  <- scan("../data/future_TO.txt",  what=numeric(), quiet=T)  #High latitude subsurface ocean temp
    SL  <- scan("../data/future_SL.txt",  what=numeric(), quiet=T)  #Reconstructed sea-level
    assimctx$forcings <- cbind( time=(1L:length(SL) - 238000L), TA=TA, TO=TO, GSL=GSL, SL=SL )
    assimctx$expert   <- !is.null(expert)

    # Kelsey runs the model to 2010, but her likelihood function only needs to be run through 2002;
    # using 2010 preserves the ability to run DAIScali_hetero_model_iid_mcmcmat.R;
    # running through 2011 is necessary to include the IPCC rates;
    # finally, just run it all the way to the prediction date since the assimilation now saves model output
    #
   #assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=ifelse(assimctx$expert, 2100, 2010))
   #assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=ifelse(assimctx$expert, 2100, 2011))
    assimctx$frc_ts   <- tsTrim(assimctx$forcings, endYear=2100)
    assimctx$frc      <- assimctx$frc_ts[ , 2:ncol(assimctx$forcings) ]


    # set up observations and observation errors
    #


    # for the Heaviside likelihood function
    #
    if (heaviside) {
        raw <- read.table("../data/GMSL_ChurchWhite2011_yr_2015.txt", fill=T)
        gmsl <- as.matrix(raw)
        colnames(gmsl) <- c("time", "sealvl", "error")
        gmsl[, "time"] <- gmsl[, "time"] - 0.5                         # time is in half years
        gmsl[, 2:ncol(gmsl)] <- gmsl[, 2:ncol(gmsl)] / 1000            # convert from mm to m
        gmsl                  <- tsBias(gmsl, lower=1880, upper=1899)  # standardize to first 20 years
        assimctx$SL.1880_1899 <- tsGetIndicesByRange(assimctx$frc_ts, lower=1880, upper=1899)
        assimctx$SL.1900_201x <- tsGetIndicesByRange(assimctx$frc_ts, lower=1900, upper=last(gmsl[, "time"]))
        assimctx$gmsl         <- tsTrim(gmsl, 1900)[, "sealvl"]
    }
    assimctx$heaviside <- heaviside


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
    assimctx$obsonly  <- rowMeans(assimctx$windows)

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

    # always set these--not just when there is an expert prior--in order to be able to run the LHS simulation
    assimctx$SL.expert <- tsGetIndices(assimctx$frc_ts, 2010)  # this is very Pfeffer-centric, using 2010 as the standardization year
    assimctx$SL.2100   <- tsGetIndices(assimctx$frc_ts, 2100)

    if (assimctx$expert) {
        switch (expert,
            pfeffer={
                #
                # Pfeffer et al. (2008)
                #
               #assimctx$expert_window <- c(146/1000, 619/1000)
                assimctx$expert_window <- c(128/1000, 619/1000)
            },
            pollard={
                stop("pollard unimplemented in daisConfigAssim()")
            }, {
                stop("unknown expert in daisConfigAssim()")
            })

        switch (prior,
            uniform={
                assimctx$expert_prior <- uniformPrior(min=assimctx$expert_window[1],    max=assimctx$expert_window[2])
            },
            beta={
                # a=2, b=3 taken from Lempert, Sriver, and Keller (2012)
                assimctx$expert_prior <- betaPrior(   min=assimctx$expert_window[1],    max=assimctx$expert_window[2], a=2, b=3)
            },
            normal={
                assimctx$expert_prior <- normPrior(mean=mean(assimctx$expert_window), upper=assimctx$expert_window[2])
            }, {
                stop("unknown prior in daisConfigAssim()")
            })
    }


    # set up model parameters and priors
    #

    paramNames      <- character()
    assimctx$units  <- character()
    assimctx$lbound <- numeric()
    assimctx$ubound <- numeric()
    init_mp         <- numeric()

    fixedParamNames <- character()
    assimctx$ep     <- numeric()

   #assimctx$fast_only <- fast_only
    if (fast_only) {
        # from a normal run
       #assimctx$ep <- c(assimctx$ep, 1.331147e+00, 4.472352e-01, 9.458727e+00, 1.336823e-02, 6.139461e-01, 4.597166e-02, 1.120297e+00, 1.289038e+03, 7.736420e+01, 7.819001e+02, 7.212625e-04)

        # from a beta run
       #assimctx$ep <- c(assimctx$ep, 2.535529e+00, 7.180800e-01, 1.049154e+01, 9.672927e-03, 7.780419e-01, 7.498188e-02, 1.006963e+00, 1.988770e+03, 1.145772e+02, 7.680061e+02, 7.227607e-04)

        # older runs
        assimctx$ep <- c(assimctx$ep, 2.000037e+00, 3.502423e-01, 8.699992e+00, 1.200846e-02, 3.500915e-01, 4.021515e-02, 1.200098e+00, 1.471000e+03, 9.500020e+01, 7.750002e+02, 5.896267e-04)
       #assimctx$ep <- c(assimctx$ep, 1.966633e+00, 1.568612e-01, 1.098878e+01, 1.176743e-02, 3.033654e-01, 6.382912e-02, 1.297880e+00, 1.861615e+03, 1.061198e+02, 7.793684e+02, 5.052874e-04)
        fixedParamNames <- c(fixedParamNames, "gamma", "alpha", "mu", "nu", "P0", "kappa", "f0", "h0", "c", "b0", "slope")
    } else {
        paramNames <- c(paramNames,     "gamma", "alpha", "mu",    "nu",                "P0", "kappa", "f0", "h0", "c", "b0", "slope")
        assimctx$units <- c(assimctx$units, "",  "",    "m^(0.5)", "m^(-0.5) yr^(-0.5)", "m", "1/C", "m/yr", "m", "m/C", "m", "")

        #                                 c('gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c'  , 'b0','slope')
       #assimctx$lbound <- c(assimctx$lbound, 0.50,  0,     4.35, 0.006, 0.175,  0.02,  0.6,  735.5,  47.5, 725, 0.00045)  # Kelsey priors
       #assimctx$ubound <- c(assimctx$ubound, 4.25,  1,    13.05, 0.018, 0.525,  0.06,  1.8, 2206.5, 142.5, 825, 0.00075)
        assimctx$lbound <- c(assimctx$lbound, 0.50 , 0,     7.05, 0.003, 0.026,  0.025, 0.6,  735.5,  47.5, 740, 0.00045)  # Tony priors
        assimctx$ubound <- c(assimctx$ubound, 4.25 , 1,    13.65, 0.015, 1.5,    0.085, 1.8, 2206.5, 142.5, 820, 0.00075)

        # Best Case (Case #4) from Shaffer (2014)
        init_mp         <- c(init_mp,         2.0,   0.35,  8.7,  0.012, 0.35,   0.04,  1.2, 1471,    95,   775, 0.0006)
    }

    # gamma prior for Tcrit and lambda
    rmif( Tcrit_prior, envir=assimctx)  # keep it clean
    rmif(lambda_prior, envir=assimctx)  # keep it clean
    assimctx$gamma_pri <- gamma_pri
    if (assimctx$gamma_pri) {
        shape.lambda <- 8.1               # gives 5% quantile at lambda=0.005 and
        rate.lambda  <- 100*shape.lambda  # gives mean at 0.01 m/yr, DeConto and Pollard (2016)
        rate.Tcrit   <- 1.37              # gives 5% quantile at Tcrit = -10 deg C
        shape.Tcrit  <- 15*rate.Tcrit     # gives mean at -15 deg C (negative requires multiplication of Tcrit by -1)

        assimctx$Tcrit_prior  <- gammaPrior(shape=shape.Tcrit,  rate=rate.Tcrit)
        assimctx$lambda_prior <- gammaPrior(shape=shape.lambda, rate=rate.lambda)
    }

    if (fast_dyn) {
        paramNames          <- c(paramNames,     "Tcrit", "lambda")
        assimctx$units      <- c(assimctx$units, "deg C",   "m/yr")
        if (wide_prior) {
            assimctx$lbound <- c(assimctx$lbound,  -30.0,   -0.005)
            assimctx$ubound <- c(assimctx$ubound,    0.0,    0.025)
        } else {
            assimctx$lbound <- c(assimctx$lbound,  -20.0,    0.005)
            assimctx$ubound <- c(assimctx$ubound,  -10.0,    0.015)
        }
        init_mp             <- c(init_mp,          -15.0,    0.010)
    }
    if (rob_dyn) {
        paramNames          <- c(paramNames,     "Hcrit",     "fa")
        assimctx$units      <- c(assimctx$units,     "m",    "1/m")
        assimctx$lbound     <- c(assimctx$lbound,  200.0, gtzero())
        assimctx$ubound     <- c(assimctx$ubound, 2000.0,     10.0)
        init_mp             <- c(init_mp,          400.0,      0.5)
    }

    names(init_mp) <- names(assimctx$lbound) <- names(assimctx$ubound) <- names(assimctx$units) <- paramNames
    names(assimctx$ep) <- fixedParamNames

    assimctx$sw             <- logical()
    assimctx$sw["fast_dyn"] <- fast_dyn
    assimctx$sw["rob_dyn"]  <- rob_dyn

    assimctx$paleo          <- paleo
    assimctx$instrumental   <- instrumental
    assimctx$prior_name     <- prior
    assimctx$expert_name    <- expert
    assimctx$all_predict    <- all_predict
    assimctx$wide_prior     <- wide_prior

    init_sp            <- numeric()
    assimctx$lbound_sp <- numeric()
    assimctx$ubound_sp <- numeric()
    if (assimctx$paleo && variance) {
        #
        # calculate paleo variance (sigma^2) from residuals
        #
        AIS_melt <- assimctx$modelfn(init_mp, assimctx)
        resid    <- assimctx$obsonly[1:3] - (AIS_melt[assimctx$obs_ind[1:3]] - mean(AIS_melt[assimctx$SL.1961_1990]))
        init_sp["var.paleo"] <- sd(resid)^2
        assimctx$units <- c(assimctx$units, "")  # add units for variance

        assimctx$alpha <- 2
        assimctx$beta  <- 1
        assimctx$lbound_sp["var.paleo"] <- gtzero()
        assimctx$ubound_sp["var.paleo"] <- 100
    }
    if (assimctx$instrumental && variance) {
        #
        # calculate instrumental variance (sigma^2) from residual
        #
        AIS_melt <- assimctx$modelfn(init_mp, assimctx)
        resid    <- assimctx$obsonly[4] - (AIS_melt[assimctx$obs_ind[4]] - AIS_melt[assimctx$SL.1992])
       #init_sp["var.inst"] <- sd(resid)^2  # can't calculate standard deviation from only ONE residual
        init_sp["var.inst"] <- resid^2
        assimctx$units <- c(assimctx$units, "")  # add units for variance

        assimctx$lbound_sp["var.inst"] <- gtzero()
        assimctx$ubound_sp["var.inst"] <- 4e-4
    }

    initAssim(assimctx, init_mp, init_sp, daisLogPri, daisLogLik)
}


daisRunAssim <- function(nbatch=ifelse(adapt, 5e5, 4e6), adapt=T, n.chain=1, assimctx=daisctx)
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

    if (assimctx$all_predict) {
        # could trim to -150K or so
        assimctx$ychain <- prmatrix(nbatch + ifelse(adapt, 0, 1), xvals=assimctx$frc_ts[, 1])
    } else {
        assimctx$ychain <- prmatrix(nbatch + ifelse(adapt, 0, 1), xvals=2100)
    }

    runAssim(assimctx, nbatch=nbatch, n.chain=n.chain, dyn.libs=dllName(c("../fortran/dais", "../fortran/dais_fastdyn", "kelsey_dais", "rob_dais")), scale=scale, adapt=adapt, extrafun=assimSaveY)

    #results <<- assimctx$chain
}


daisRunLhs <- function(nbatch1=1e3, nbatch2=2*nbatch1, assimctx=daisctx)
{
    # this seed ensures we get fixed parameters that will allow us to wander
    # in to the bifurcated parameter space
    #
   #set.seed(1)  # shows disconnected parameter space with small peninsual
   #set.seed(4)  # shows disconnected parameter space with one random dot in peninsula
   #set.seed(5)  # plausible and non-plausible lambda connect here
    set.seed(7)  # shows disconnected parameter space with bigger peninsula
   #set.seed(8)  # shows disconnected parameter space with smaller peninsula
   #set.seed(9)  # shows disconnected parameter space with small peninsula
   #set.seed(10) # see above

    # this first run is just to get a good estimate for fixed parameters;
    # uniform prior is not a good choice since likelihood is equal for all samples
    #
    daisConfigAssim(fast_only=F, wide_prior=T,
                    instrumental=T, paleo=T, expert="pfeffer", prior="normal", gamma_pri=T, heaviside=T, assimctx=assimctx, variance=F)
    assimctx$lhs_ychain <- prmatrix(nbatch1, xvals=2100)
    lhs_est <- assimRunLhs(assimctx=assimctx, nbatch=nbatch1, extrafun=assimLhsSaveY)

    # TODO:  this is hackish, knowing last two parameters are Tcrit and lambda
    fixed_indices <- 1:(length(assimctx$init_mp) - 2)
    daisConfigAssim(fast_only=T, wide_prior=T,
                    instrumental=F, paleo=F, expert=NULL,      prior=NULL,     gamma_pri=F, heaviside=F, assimctx=assimctx)  # the LF is ignored
    assimctx$ep <- lhs_est$maxLikParam[fixed_indices]
    print(assimctx$ep)

    assimctx$lhs_ychain <- prmatrix(nbatch2, xvals=2100)
    assimctx$lhs <- assimRunLhs(assimctx=assimctx, nbatch=nbatch2, extrafun=assimLhsSaveY)
    assimctx$lhs$ychain <- assimctx$lhs_ychain
    rmif(lhs_ychain, envir=assimctx)  # keep it clean
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


daisRunPredict <- function(nbatch=( nrow(assimctx$chain) / 5 ), subsample=T, assimctx=daisctx, prctx=prdaisctx)
{
    prctx$assimctx <- assimctx

    prchain       <- assimFixOutput(assimctx, assimctx$ychain)
    burnIn        <- burnedInd(prchain)
    if (subsample) {
        samples   <- sample(burnIn, nbatch, replace=T)
    } else {
        # note that nbatch is not used here
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


daisSlrLab <- function(year=2100)
{
   #return (paste("Projected AIS Volume Loss in", year, "[SLE m]"))
    return (paste("Projected AIS SLR in", year, "(m)"))
}


daisTcritLab <- function()
{
    # .(tau) substitutes the value of tau
    #bquote(paste(tau[1]==.(tau), " a")

   #return (    bquote(T[crit]~(degree*C)))
    return (expression(T[crit]~(degree*C)))
}


daisLambdaLab <- function()
{
   #return (    bquote(lambda~(m~y^-1)))
    return (expression(lambda~(m~y^-1)))
}


daisInitRegress <- function(assimctx=daisctx)
{
    if (!is.null(assimctx$slope.Ta2Tg)) {
        return ()
    }

    if (is.null(assimctx$forcings)) {
        daisConfigAssim(assimctx=assimctx)
    }

    # linear regression of Antarctic surface temperature to GMST
    #

    file <- "../data/HadCRUT.4.4.0.0.annual_ns_avg.txt"
   #raw  <- scan(file, what=numeric(), quiet=T)
   #gmst <- matrix(raw, ncol=12, byrow=T)
    raw  <- read.table(file, fill=T)
    gmst <- as.matrix(raw)

    gmst <- gmst[ , 1:2 ]
    colnames(gmst) <- c("time", "temp")
    gmst <- tsBias(gmst, lower=1850, upper=1869)

    gmst_lm <- tsTrim(gmst, 1850, 1997)[, 2]
    ta_lm   <- tsTrim(assimctx$forcings, 1850, 1997)[, 2]
    fit <- lm(gmst_lm ~ ta_lm)
    assimctx$intercept.Ta2Tg <- coef(fit)[1]
    assimctx$slope.Ta2Tg     <- coef(fit)[2]
}


daisGmst <- function(x, assimctx=daisctx)
{
    daisInitRegress(assimctx=assimctx)

    y <- assimctx$slope.Ta2Tg * x + assimctx$intercept.Ta2Tg

    return (y)
}


daisTa <- function(y, assimctx=daisctx)
{
    daisInitRegress(assimctx=assimctx)

    x <- (y - assimctx$intercept.Ta2Tg) / assimctx$slope.Ta2Tg

    return (x)
}


daisStats <- function(assimctx=daisctx)
{
    burn_ind <- burnedInd(assimctx$chain)
    quants   <- quantile(assimctx$chain[burn_ind, "Tcrit"], probs=c(0.05, 0.50, 0.95))
    gmst     <- daisGmst(quants)
    rounded  <- signif(gmst, digits=2)

    return (rounded)
}
