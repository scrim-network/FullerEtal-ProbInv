# Copyright 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# pola.R

source("assim.R")
source("plot.R")
source("ts.R")


loadPolaSlr <- function(file="../data/kkeller/annual_PSMSLdata.ascii")
{
    raw <- scan(file=file, what=numeric(), quiet=T)
    slr <- matrix(raw, ncol=2, byrow=T)
    colnames(slr) <- c("time", "sealvl")

    bad <- which(slr[, "sealvl"] == -99999)
    for (row in bad) {
        slr[row, "sealvl"] <- (slr[row - 1, "sealvl"] + slr[row + 1, "sealvl"]) / 2
    }

    slr <- tsBias(slr, 1)

    return (slr)
}


pola <- function(times, parms)
{
    mp <- parms$mp
    if (length(parms$ep)) {
        mp <- replace(mp, names(parms$ep), parms$ep)
    }

    t   <- times - times[1]
    slr <- mp["a"] * t^2 + mp["b"] * t + mp["c"]

    return (cbind(times, slr))
}


polaModel <- function(mp, assimctx)
{
    slr <- pola(assimctx$times, list(mp=mp, ep=assimctx$ep))

    return (slr[, 2])
}


if (!exists("polassimctx")) {
    polassimctx <- env()
}


if (!exists("prpolactx")) {
    prpolactx <- env()
}


polaConfigAssim <- function(
    assimctx=polassimctx,
    fixc=F
    )
{
    slr <- loadPolaSlr()


    # assim.R uses these
    #

    assimctx$obsonly <- slr[, "sealvl"]
    assimctx$times   <- slr[, "time"]
    assimctx$modelfn <- polaModel

    c        <- slr[1, "sealvl"]
    stderror <- 100

    #                    a, b, c
    assimctx$lbound <- c(0, 0, c - 4 * stderror)
    assimctx$ubound <- c(5, 5, c + 4 * stderror)
    names(assimctx$lbound) <- names(assimctx$ubound) <- c("a", "b", "c")


    # plotting functions use these
    assimctx$units <- c("mm/a", "mm/a^2", "mm", "", "")
    assimctx$obs   <- slr

    # polaModel() uses these
    assimctx$ep <- numeric()

    if (fixc) {
        fp <- numeric()
        fp["c"] <- c
        configFixedParms(assimctx, fp)
    }

    # get initial conditions from best fit model
    configAssim(assimctx, ar=1, obserr=F, sigma_max=stderror)
}


polaRunAssim <- function(
    nbatch=100000,
    initial=is.null(assimctx$chain),
    assimctx=polassimctx
    )
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (initial || ncol(assimctx$chain) != length(init_mp) + length(init_sp)) {
        #print("using initial scale")

        scale <- NULL
    } else {
        print("using proposal matrix")

        mult  <- 0.5
        scale <- assimProposalMatrix(assimctx$chain, mult=mult)
    }

    runAssim(assimctx, nbatch, scale)
}


polaPlotFit <- function(assimctx=polassimctx, outfiles=F)
{
    newDev("polafit", outfile=outfiles)

    emptyPlot(xlim=range(assimctx$times), ylim=range(assimctx$obsonly), xlab="Year", ylab="Sea-level anomaly (mm)")
    points(assimctx$obs)

    if (T) {
        slr <- polaModel(assimctx$init_mp, assimctx)
        lines(assimctx$times, slr, lwd=2)
    } else {
        slr <- pola(assimctx$times, parms=list(mp=assimctx$init_mp, ep=assimctx$ep))
        lines(slr, lwd=2)
    }
}


polaRunPredict <- function(..., nbatch=10000, year=2100, assimctx=polassimctx, prctx=prpolactx, replace=F)
{
    runPredict(
        nbatch=nbatch,
        year=year,
        assimctx=assimctx, prctx=prctx,
        forcings=NULL,
        modelfn=pola,
        replace=replace,
        ep=assimctx$ep,
        ...
        )
}


polaPlotPredict <- function(prctx=prpolactx, outfiles=F)
{
    prqplot(
        prctx=prctx,
        ycaption="Sea-level anomaly (mm)",
        caption="Figure 1. Posterior predictive for sea-level anomaly",
        outfiles=outfiles,
        shade=F
        )
}
