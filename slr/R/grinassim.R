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
# grinassim.R
# written by Robert W. Fuller on 0906xx

source("assim.R")
source("grinsted.R")


grinModel <- function(mp, assimctx)
{
    slr <- grinsted(assimctx$times, list(mp=mp, frc=list(assimctx$frc), spl=2, ep=assimctx$ep, sw=assimctx$sw))

    return (slr[assimctx$mod_ind, "sealvl"])
}


if (!exists("grinassimctx")) {
    grinassimctx <- env()
}


grinConfigAssim <- function(
    ar=1, obserr=T,
    startYear, endYear,
    log_tau=F, paleo=F, historical=F, ipcc=F, fp=NULL, ipccFile="miub_echo_g.txt",
    useSSE=obserr, fixrho=F, rholimit=0.99,
    assimctx=grinassimctx
    )
{
    assimctx$ipccFile <- ipccFile

    if (useSSE) { # && (is.null(assimctx$fit$optim$bestmem) || assimctx$sw["log_tau"] != log_tau)) {
        grinstedRunFit(assimctx, startYear, endYear, log_tau, paleo, historical, ipcc, fp, ipccFile)
    } else {
        grinstedConfig(assimctx, startYear, endYear, log_tau, paleo, historical, ipcc, fp, ipccFile)
        rmif(fit, envir=assimctx)
    }

    # assim.R uses these
    assimctx$modelfn <- grinModel

    # get initial conditions from best fit model
    # TODO:  is 0.85 a good choice for ipcc which doesn't seem to fit AR(p)?
    #configAssim(assimctx, assimctx$fit$optim$bestmem, ar=ar, obserr=obserr, rholimit=ifelse(ipcc, 0.85, 0.99))
    configAssim(assimctx, assimctx$fit$optim$bestmem, ar=ar, obserr=obserr, fixrho=fixrho, rholimit=rholimit)
}


grinConfigPerfectModel <- function(
    truth="mode", uniform=T, useUrban=F, year=ifelse(useUrban, 2100, last(assimctx$times)), spdiv=1, assimctx=grinassimctx
    )
{
    # TODO:  copy assimctx before destructively overwriting?
    # need deep_copy_env()?

    if (useUrban) {
        assimctx$frc <- loadUrbanAvg()
    }

    # avoid problems where sigma is fixed
    # TODO:  allow sp to be fixed in assim.R?
    if (spdiv != 1) {
        assimctx$chain[, "sigma"] <- assimctx$chain[, "sigma"] / spdiv
    }
    configPerfectModel(assimctx=assimctx, truth=truth, uniform=uniform)

    # TODO:  need to grep for colon and eradicate these
    # fgrep ":" *R | fgrep -i year
    # fgrep ":" *R | fgrep -i time
    #
    times <- assimctx$times[1]:year
    slr <- grinsted(times, list(mp=assimctx$init_mp, frc=list(assimctx$frc), spl=2, ep=assimctx$ep, sw=assimctx$sw))
    assimctx$true_predict <- last(slr[, "sealvl"])
    print(paste("true prediction:", assimctx$true_predict))
}


grinRunAssim <- function(
    nbatch=1000000, initial=is.null(assimctx$chain), mult=NULL,
    assimctx=grinassimctx
    )
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp
    ar      <- assimctx$ar
    obserr  <- assimctx$obserr

    if (initial || any(colnames(assimctx$chain) != names(c(init_mp, init_sp)))) {
        #print("using initial scale")

        if (ar == 0) {
            if (!obserr) {
                scale <- c(rep(50, length(init_mp)), rep(100, length(init_sp)))
            } else {
                scale <- c(rep(100, length(init_mp)))
            }
        } else if (ar == 1) {
            scale <- c(rep(25, length(init_mp) + length(init_sp)))
        } else {
            scale <- c(rep(25, length(init_mp)), 100, 100, rep(50, length(init_sp) - 2))
        }
        scale <- abs(c(init_mp, init_sp)) / scale

        if (assimctx$sw["log_tau"] && !is.na(scale["tau"])) {
            scale["tau"] <- scale["tau"] / 10
        }

        scale <- NULL

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


if (!exists("prgrinctx")) {
    prgrinctx <- env()
}


grinRunPredict <- function(
    ...,
    assimctx=grinassimctx, prctx=prgrinctx,
    calibration=F
    )
{
    gmst <- loadUrban()

    if (calibration) {
        forcings <- list(assimctx$frc)
    } else {
        forcings <- list(gmst)
    }

    runPredict(
        assimctx=assimctx, prctx=prctx,
        forcings=forcings, modelfn=grinsted,
        ep=assimctx$ep, sw=assimctx$sw,
        ...
        )
}


grinRunMobergHindcast <- function(..., assimctx=grinassimctx, prctx=prgrinctx)
{
    runPredict(
        assimctx=assimctx, prctx=prctx,
        year=last(assimctx$times),
        forcings=list(assimctx$frc), modelfn=grinsted,
        ep=assimctx$ep, sw=assimctx$sw,
        ...
        )
}


grinRunUniformPredict <- function(..., assimctx=grinassimctx, prctx=prgrinctx)
{
    gmst <- loadUrban()
    runUniformPredict(
        assimctx=assimctx, prctx=prctx,
        forcings=list(gmst), modelfn=grinsted,
        ep=assimctx$ep, sw=assimctx$sw,
        ...
        )
}


grinRunNoise <- function(
    nbatch=100000, year=2200,
    assimctx=grinassimctx
    )
{
    times <- assimctx$times[1]:year
    N     <- length(times)

    # could use prmatrix, except that names are set here instead of in assim.R
    noise           <- sampleNoiseAssim(assimctx, N, nbatch)
    colnames(noise) <- as.character(times)

    assimctx$prnoise <- noise

    if (nbatch <= 20) {
        return (noise)
    }
}
