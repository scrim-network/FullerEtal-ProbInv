# Copyright 2009, 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# allgrgis.R
# written by Robert W. Fuller on 091009

#source("assim.R")
source("grinsted.R")
#source("grinassim.R")
source("alley.R")

loadSlrModel()


allgrgisOde <- function(t, y, parms)
{
    mp       <- parms$mp
    gmst     <- parms$frc[[1]]
    gmstCol  <- parms$spl

    #s_total <- y[1]
    s_other  <- y[2]
    s_gis    <- min(y[3], mp["max_sle"])

    temp     <- tsFindByDate(gmst, t, gmstCol)

    if (parms$sw["old_ref"]) {
        masstemp <- tsFindByDate(parms$frc[[2]], t, gmstCol) * mp["gis_scl"]
    } else {
        masstemp <- (temp + mp["gis_temp"]) * mp["gis_scl"]
    }

    ds_other     <- (mp["a"] * temp     + mp["b"]              - s_other) / mp["tau"]
    if (masstemp <= 0) {
        ds_gis   <- (mp["c"] * masstemp                        - s_gis)   / mp["tau2"]
    } else {
        ds_gis   <- (mp["c"] * masstemp + mp["d"] * masstemp^2 - s_gis)   / mp["tau2"]
    }
    ds_total <- ds_other + ds_gis

    # save this for the GIS likelihood
    names(ds_gis)   <- "ds_gis"
    names(ds_total) <- "ds_total"

    return (list(c(ds_total, ds_other, ds_gis), ds_gis, ds_total))
}


allgrgis <- function(times, parms)
{
    if (length(parms$ep)) {
        parms$mp <- replace(parms$mp, names(parms$ep), parms$ep)
        parms$ep <- NULL
    }

    y <- c(rep(parms$mp["s0"], 2), parms$mp["gis_s0"])
    names(y) <- c("sealvl", "other", "gis")

    if (parms$sw["log_tau"]) {
        parms$mp["tau"]  <- exp(parms$mp["tau"])
    }
    if (parms$sw["log_tau2"]) {
        parms$mp["tau2"] <- exp(parms$mp["tau2"])
    }
    if (parms$sw["old_ref"]) {
        parms$frc[[2]] <- tsBias(parms$frc[[1]], 1, lower=1850, upper=1879)
    }

    # ode( y=y, times=times, func=allgrgisOde,    parms=parms, method="euler" )
    euler( y=y, times=times, func="allgrgisOdeC", parms=parms, dllname="slrmodel", initfunc="grinstedOdeInit", nout=2, outnames=c("ds_gis", "ds_total") )
}


# because allgrgisLogLik() overrides logLik(), this is called
# from exactly one place to initialize sp:  configAssim()
#
allgrgisModel <- function(mp, assimctx)
{
    slr <- allgrgis(assimctx$times, list(mp=mp, frc=list(assimctx$gmst), spl=2, ep=assimctx$ep, sw=assimctx$sw))

    return (slr[, "sealvl"])
}


allgrgisFitOde <- function(t, y, parms)
{
    mp       <- parms$mp
    gis_gmst <- parms$frc[[1]]
    gmstCol  <- parms$spl

    s <- min(y[1], mp["max_sle"])

    if (parms$sw["old_ref"]) {
        # could get the same effect be setting gis_temp to zero
        temp <-  tsFindByDate(gis_gmst, t, gmstCol) * mp["gis_scl"]
    } else {
        temp <- (tsFindByDate(gis_gmst, t, gmstCol) + mp["gis_temp"]) * mp["gis_scl"]
    }
    if (temp <= 0) {
        ds <- (mp["c"] * temp                    - s) / mp["tau2"]
    } else {
        ds <- (mp["c"] * temp + mp["d"] * temp^2 - s) / mp["tau2"]
    }

    names(ds) <- "ds_gis"

    return (list(ds, ds))
}


allgrgisFit <- function(times, mp, gis_gmst, ep, sw, gmstCol=2)
{
    if (length(ep)) {
        mp <- replace(mp, names(ep), ep)
    }

    y <- mp["gis_s0"]
    names(y) <- "sealvl"

    if (sw["log_tau2"]) {
        mp["tau2"] <- exp(mp["tau2"])
    }
    if (sw["old_ref"]) {
        gis_gmst <- tsBias(gis_gmst, 1, lower=1850, upper=1879)
    }

    parms <- list(mp=mp, frc=list(gis_gmst), spl=gmstCol, sw=sw)
    # ode( y=y, times=times, func=allgrgisFitOde,    parms=parms, method="euler" )
    euler( y=y, times=times, func="allgrgisFitOdeC", parms=parms,
           dllname="slrmodel", initfunc="grinstedOdeInit", nout=1, outnames="ds_gis" )
}


rignotPlotFit <- function(assimctx=allgrgisctx, lines=T, outfiles=F)
{
    newDev("rignotfit", outfiles)

    massBal <- loadRignot(pure=T)
    slr <- allgrgisFit(assimctx$gis_times, mp=assimctx$init_mp, assimctx$gmst, assimctx$ep, assimctx$sw)

    emptyPlot(
        xlim=c(massBal[1, "time"], last(massBal[, "time"])),
        #ylim=range(slr[assimctx$gis_ind, "ds_gis"]),
        ylim=c(0, 1e-3),
        xlab="Year",
        ylab=slrGreenlandLab()
        )

    ts <- cbind(assimctx$gis_times[assimctx$gis_ind], assimctx$gis_obs, assimctx$gis_err$error)
    colnames(ts) <- c("time", "SLE", "error")
    tsErrorBars(ts, shade=F, col="purple", xbeam=T, ibeam=T)

    if (lines) {
        lines(slr[, c("time", "ds_gis")])
    } else {
        points(slr[, c("time", "ds_gis")])
    }

    legend("topleft", c("Model", "Observations"), lty=c("solid", NA), col=c("black", "purple"), pch=c(NA, 3))

    mtext("Figure 1. Best fit of GIS model to Rignot mass balance", outer=TRUE, side=1)

    if (outfiles) { finDev() }
}


allgrgisMinFn <- function(mp, assimctx)
{
    slr <- allgrgisFit(assimctx$gis_times, mp, assimctx$gmst, assimctx$ep, assimctx$sw)
    err <- sse(assimctx$gis_obs, slr[assimctx$gis_ind, "ds_gis"])
    if (is.na(err)) {
        return (Inf)
    }

    return (err)
}


allgrgisConfig <- function(assimctx, startYear, log_tau2, old_ref)
{
    massBal <- loadRignot(pure=T)
    assimctx$gis_obs   <- massBal[, "SLE"]
    assimctx$gis_err   <- setErrorAssim(env(), massBal[, "error"])
    assimctx$gis_ind   <- massBal[, "time"] - startYear + 1
    assimctx$gis_times <- startYear:(massBal[nrow(massBal), "time"])

    if (old_ref) {
        gmsl     <- tsBias(loadJev()$annual, 1, lower=1850, upper=1879)
        s0center <- tsFindByDate(gmsl, startYear)
        stderror <- tsFindByDate(gmsl, startYear, "error")

        assimctx$gis_lbound <- c(s0center - stderror, 0.0, 0.0,   50, 1.5, 0.0)
        assimctx$gis_ubound <- c(s0center + stderror, 2.0, 2.0, 5000, 2.0, 1.5)
    } else {
        # temp is offset from AR4
        # scl is polar amplification
        # s0 is how far out of equilibrium in 1850 (at end of LIA)
        #                         s0    c    d  tau2  scl temp
        assimctx$gis_lbound <- c(0.0, 0.0, 0.0,   50, 1.5, 0.0)

        # so this is what old_ref was really about...
        #assimctx$gis_ubound <- c(0.5, 1.0, 0.5, 5000, 2.0, 1.5)

        # these do not reproduce
        #assimctx$gis_ubound <- c(0.5, 2.0, 1.0, 5000, 2.0, 1.5)
        #assimctx$gis_ubound <- c(0.5, 10.0, 2.5, 10000, 2.0, 1.5)

        # this reproduces:  tide gage has slightly higher mode
        #
        #assimctx$gis_ubound <- c(0.5, 2.0, 2.0, 5000, 2.0, 1.5)
        assimctx$gis_ubound <- c(0.5, 5.0, 2.0, 5000, 2.0, 1.5)
    }

    names(assimctx$gis_lbound) <- names(assimctx$gis_ubound) <- c("gis_s0", "c", "d", "tau2", "gis_scl", "gis_temp")

    if (log_tau2) {
        assimctx$gis_lbound["tau2"] <- log(assimctx$gis_lbound["tau2"])
        assimctx$gis_ubound["tau2"] <- log(assimctx$gis_ubound["tau2"])
    }

    # switches are created in various call sequences between the
    # two config functions;  don't inadvertently delete a switch
    #
    if (is.null(assimctx$sw)) {
        assimctx$sw <- logical()
    }
    assimctx$sw["log_tau2"] <- log_tau2

    assimctx$ep <- c(max_sle=7.5)
}


allgrgisRunFit <- function(assimctx=allgrgisctx)
{
    control <- list(CR=1.0, itermax=500)

    assimctx$fit <- named_DEoptim(
        FUN=allgrgisMinFn,
        lower=assimctx$gis_lbound, upper=assimctx$gis_ubound,
        control=control,
        assimctx=assimctx
        )
}


allgrgisLogLik <- function(mp, sp, assimctx)
{
    llik <- 0

    if (length(assimctx$ep)) {
        mp <- replace(mp, names(assimctx$ep), assimctx$ep)
    }

    # integrate for a few additional years so that the GIS likelihood
    # can be drawn from all eight GIS observations
    #
    slr <- allgrgis(assimctx$gis_times, list(mp=mp, frc=list(assimctx$gmst), spl=2, sw=assimctx$sw))

    # dmvnorm is sensitive to the order here when used purely with residuals
    # hence, res = obs - model
    #

    if (assimctx$sw["rignot_ts"]) {
        res_gis <- assimctx$gis_obs - slr[assimctx$gis_ind, "ds_gis"]
        llik    <- llik + llik_obs(res_gis, sp, assimctx$gis_err)
    }

    if (assimctx$sw["tide_ts"]) {
        res_total <- assimctx$obsonly - slr[1:length(assimctx$times), "sealvl"]
        llik      <- llik + assimctx$logLik(res_total, sp, assimctx)
    }

    if (assimctx$sw["tau2_prior"]) {
        llik <- llik + assimctx$tau2_prior$dens(mp["tau2"])
    }

    if (assimctx$sw["alley_prior"]) {
        llik <- llik + logLik(mp, sp, assimctx$alleyctx)
    }

    return (llik)
}


if (!exists("allgrgisctx")) {
    allgrgisctx <- env()
}


allgrgisConfigAssim <- function(
    ar=1, obserr=T,
    old_ref=F,
    tide_ts=T, rignot_ts=T, alley_prior=T, tau2_prior=T,
    log_tau=F, log_tau2=F,
    useSSE=(!tide_ts) | (ar==2 & obserr),
    fp=NULL,
    startYear, endYear,
    assimctx=allgrgisctx
    )
{
    if (useSSE) {
        grctx <- assimctx$grctx
        if (is.null(grctx$fit$optim$bestmem) || grctx$sw["log_tau"] != log_tau)
        {
            grctx <- assimctx$grctx <- env()

            print("running fit for initial Grinsted parameters")

            # this is the one case where grinstedConfig() is not needed below
            grinstedRunFit(grctx, startYear, endYear, log_tau)
        }
    } else {
        grctx <- env()
    }

    # startYear or endYear may have changed, so always run this
    grinstedConfig(grctx, startYear, endYear, log_tau)

    if (missing(startYear)) {
        startYear <- grctx$times[1]
    }
    if (missing(endYear)) {
        endYear <- last(grctx$times)
    }

    allgrgisConfig(assimctx=assimctx, startYear=startYear, log_tau2=log_tau2, old_ref=old_ref)

    # could call grinAssimConfig() here and override values

    # allgrgisModel() uses these
    assimctx$times <- startYear:endYear
    assimctx$gmst  <- loadHadcrut()$annual

    #print(assimctx$gis_gmst[1, "temp"] - assimctx$gmst[1, "temp"])

    # plotting functions use this
    assimctx$units <- c("m", "m/C", "m", "Years", "m", "m/C", "m/C^2", "Years", "", "C", "", "", "")

    # trim the time series AFTER applying ar4 bias
    gmsl <- tsTrim(loadJev()$annual, startYear, endYear)

    # assim.R uses these
    assimctx$obsonly <- gmsl[, "sealvl"]
    assimctx$modelfn <- allgrgisModel
    assimctx$lbound  <- c(grctx$lbound, assimctx$gis_lbound)
    assimctx$ubound  <- c(grctx$ubound, assimctx$gis_ubound)
    setErrorAssim(assimctx, gmsl[, "error"])

    # likelihood function uses these
    assimctx$sw["tide_ts"]      <- tide_ts
    assimctx$sw["rignot_ts"]    <- rignot_ts
    assimctx$sw["tau2_prior"]   <- tau2_prior
    assimctx$sw["alley_prior"]  <- alley_prior
    assimctx$sw["log_tau"]      <- log_tau

    # model function uses this
    assimctx$sw["old_ref"]      <- old_ref

    if (tau2_prior) {
        #
        # log normal prior on tau2 with 2 standard deviations between 4500 and 1500
        # log(4500) - log(1500) = log(4500/1500) = log(3)
        #
        if (log_tau2) {
            assimctx$tau2_prior <- llnormPrior(mean=1500, upper=4500)
        } else {
            assimctx$tau2_prior <- lnormPrior( mean=1500, upper=4500)
        }
    } else {
        rmif(tau2_prior, envir=assimctx)
    }

    if (alley_prior) {
        if (is.null(assimctx$alleyctx)) {
            assimctx$alleyctx <- env()
            cat("Alley prior assimilation: ")
        }
        alleyConfigAssim(assimctx=assimctx$alleyctx, herr=F, useSSE=useSSE)
    }

    if (useSSE) {

        # get initial conditions from best fit model
        #

        if (   is.null(assimctx$fit$optim$bestmem)
            || assimctx$fit$optim$bestmem["tau2"] < assimctx$gis_lbound["tau2"]
            || assimctx$fit$optim$bestmem["tau2"] > assimctx$gis_ubound["tau2"])
        {
            print("running fit for initial model parameters")

            assimctx$gis_lbound <- removeFixedParms(assimctx$gis_lbound, fp)
            assimctx$gis_ubound <- removeFixedParms(assimctx$gis_ubound, fp)
            assimctx$ep         <- replace(assimctx$ep, names(fp), fp)

            allgrgisRunFit(assimctx=assimctx)
        }

        init_mp <- c(grctx$fit$optim$bestmem, assimctx$fit$optim$bestmem)

        # use priors as initial conditions
        #

        if (tau2_prior) {
            init_mp["tau2"] <- assimctx$tau2_prior$mode()
        }

        if (alley_prior) {
            init_mp["c"] <- assimctx$alleyctx$init_mp["c"]
            init_mp["d"] <- assimctx$alleyctx$init_mp["d"]
        }

        if (old_ref) {
            fp <- c(fp, gis_s0=0)
        }

        init_mp <- removeFixedParms(init_mp, fp)

    } else {
        init_mp <- NULL
    }

    configFixedParms(assimctx, fp)

    configAssim(assimctx, init_mp, ar=ar, obserr=obserr, llikfn=allgrgisLogLik, itermax=1000)
    if (!tide_ts) {

        # avoid numerical instability problem with Jeffreys prior for variance
        # on random sigma;  otherwise, Jeffreys prior dominates and acceptance
        # rate becomes unacceptably low
        #
        assimctx$logPri <- lpri_bounds

        # does noise sampling from random sigma make sense?
        assimctx$noise     <- noise_zeros
        assimctx$init_sp   <- numeric()
        assimctx$lbound_sp <- numeric()
        assimctx$ubound_sp <- numeric()
    }
}


allgrgisRunAssim <- function(
    nbatch=1000000, initial=is.null(assimctx$chain), mult=NULL,
    assimctx=allgrgisctx
    )
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp
    ar      <- assimctx$ar
    obserr  <- assimctx$obserr

    if (initial || any(colnames(assimctx$chain) != names(c(init_mp, init_sp)))) {

       if (ar == 0) {
            if (!obserr) {
                scale <- c(rep(50, length(init_mp)), rep(100, length(init_sp)))
            } else {
                scale <- 100
            }
        } else if (ar == 1) {

            if (assimctx$sw["tide_ts"]) {
                scale <- 37

            } else if (assimctx$sw["rignot_ts"]) {
                scale <- 12

            } else if (assimctx$sw["alley_prior"]) {
                scale <- 10

            } else if (assimctx$sw["tau2_prior"]) {
                scale <- 2
            }

        } else {
            scale <- c(rep(25, length(init_mp)), 100, 100, rep(50, length(init_sp) - 2))
        }

        scale <- abs(c(init_mp, init_sp)) / scale

        if (assimctx$sw["log_tau"]) {
            scale["tau"]  <- scale["tau"]  / 10
        }
        if (assimctx$sw["log_tau2"]) {
            scale["tau2"] <- scale["tau2"] / 10
        }

        if (assimctx$sw["old_ref"]) {
            print("using initial scale:")
            print(scale)
        } else {
            scale <- NULL
        }

    } else {

        print("using proposal matrix")

        if (is.null(mult)) {
            if (ar == 0) {
                # inverse relation between mult and acceptance
                mult <- ifelse(obserr, 0.25, 0.305)
            } else if (ar == 1) {
                mult <- 0.25
            } else {
                mult <- ifelse(obserr, 0.40, 0.20)
            }
        }
        #print(paste("mult is", mult))
        scale <- assimProposalMatrix(assimctx$chain, mult=mult)
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale)
}


if (!exists("prallgrgisctx")) {
    prallgrgisctx <- env()
}


noise_gis <- function(sp, N, assimctx)
{
    nobs <- length(assimctx$gis_err$error)
    sparse <- noise_obs(sp, nobs, assimctx$gis_err)
    noise <- rep(0, N)
    noise[ assimctx$gis_ind ] <- sparse

    #print(noise);stop("foo")

    # save noise for noise_copy
    assimctx$noise_realization <- noise
}


allgrgisRunPredict <- function(
    ...,
    assimctx=allgrgisctx, prctx=prallgrgisctx,
    calibration=F
    )
{
    gmst <- loadUrban()

    if (calibration) {
        forcings <- list(assimctx$gmst)
    } else {
        forcings <- list(gmst)
    }

    runPredict(
        assimctx=assimctx, prctx=prctx,
        forcings=forcings,
        modelfn=allgrgis,

        # can use assimctx$noise for last item to check noise generation
        noisefn=list(noise_save, noise_copy, noise_gis, noise_copy, noise_zeros),

        outnames=c("prchain", "otherchain", "gischain", "ds_gis", "ds_total"),
        ep=assimctx$ep,
        sw=assimctx$sw,
        ...
        )
}


allgrgisRunUniformPredict <- function(
    ...,
    assimctx=allgrgisctx, prctx=prallgrgisctx
    )
{
    gmst <- loadUrban()

    runUniformPredict(
        assimctx=assimctx, prctx=prctx,
        forcings=list(gmst), modelfn=allgrgis,
        outnames=c("prunichain", "uni_otherchain", "uni_gischain", "uni_ds_gis", "uni_ds_total"),
        ep=assimctx$ep,
        sw=assimctx$sw,
        ...
        )
}


allgrgisRunSeqPredict <- function(
    nbatch=100000,
    startTemp = -2, endTemp = 6,
    assimctx=allgrgisctx, prctx=prallgrgisctx
    )
{
    prctx$assimctx <- assimctx

    mcmcChain <- assimctx$chain

    nChain <- sample(burnedInd(mcmcChain), nbatch, replace=T)

    temp <- seq(startTemp, endTemp, length.out=50)

    seq_gis   <- prmatrix(nbatch, temp)
    seq_other <- prmatrix(nbatch, temp)

    for (i in 1:nbatch) {
        p <- mcmcChain[ nChain[i], ]

        seq_other[i, ] <- p["a"] * temp + p["b"]

        # TODO:  gis_temp could be a fixed parameter
        masstemp <- temp + p["gis_temp"]

        ind <- (masstemp <= 0)
        seq_gis[i, ind] <- p["c"] * masstemp[ind]
        ind <- !ind
        seq_gis[i, ind] <- p["c"] * masstemp[ind] + p["d"] * (masstemp[ind])^2
    }

    # TODO:  this will have to change if we sample the uncertainty in max_sle
    max_sle <- assimctx$ep["max_sle"]
    seq_gis[ seq_gis > max_sle ] <- max_sle

    prctx$seq_gischain   <- seq_gis
    prctx$seq_otherchain <- seq_other

    if (nbatch <= 20) {
        return (seq_gis)
    }
}
