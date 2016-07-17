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
# mass.R
# written by Robert W. Fuller on 0907xx

source("roblib.R")
source("data.R")
source("plot.R")
source("assim.R")
loadLibrary("DEoptim")


massBal <- function(mp, gmst)
{
    mass <- matrix(nrow=nrow(gmst), ncol=2)
    colnames(mass) <- c("time", "mass")

    temp <- gmst[, "temp"]
    mass[, "time"] <- gmst[, "time"]
    mass[, "mass"] <- mp["c"] * temp + mp["d"] * temp^2

    return (mass)
}


massMinFn <- function(mp, gmst, obs)
{
    mass <- mp["c"] * gmst + mp["d"] * gmst^2
    err <- sse(mass, obs)
    #return (ifelse(is.finite(err), err, Inf))
}


massModel <- function(mp, assimctx)
{
    mass <- mp["c"] * assimctx$gmst + mp["d"] * assimctx$gmst^2
}


if (!exists("massctx")) {
    massctx <- env()
}


massRunFit <- function(assimctx=massctx)
{
    massConfig(assimctx=assimctx)

    #control <- list(CR=1.0, itermax=500)
    control <- list()

    assimctx$fit <- named_DEoptim(
        FUN=massMinFn,
        lower=assimctx$lbound, upper=assimctx$ubound,
        control=control,
        gmst=assimctx$gmst, obs=assimctx$obsonly
        )
}


# Nathan did not like this
massPlotFit <- function(outfiles=F, startYear=1850, endYear=2007, assimctx=massctx)
{
    newDev("massfit", outfiles)

    gmst <- loadHadcrut(lower=1850, upper=1879)$annual
    gmst <- tsTrim(gmst, startYear, endYear)
    obs <- loadRignot(pure=F)

    model <- massBal(assimctx$fit$optim$bestmem, gmst)

    ylim <- range(model[, "mass"], obs[, "mass"])
    xlim <- c(startYear, endYear)
    emptyPlot(xlim, ylim, "Year", "GIS mass balance (Gton)")

    cex=c(0.25, 1.5)
    #pch=c(20, 10)
    pch=c(1, 10)

    points(model, pch=pch[1], cex=cex[1])
    points(x=obs[, "time"], y=obs[, "mass"], pch=pch[2], cex=cex[2])

    legend(
        "bottomleft",
        legend=c("Model", "Observations"),
        pch=pch, pt.cex=cex
        )

    parms <- paste(
        sep="",
        "Fitted parameters: ",
        "c = ",   format(assimctx$fit$optim$bestmem["c"], digits=4, nsmall=1), ", ",
        "d = ",   format(assimctx$fit$optim$bestmem["d"], digits=4, nsmall=1)
        )
    mtext(parms, outer=T, side=1)

    if (outfiles) { finDev() }
}


summation <- function(x)
{
    N <- length(x)
    #y <- vector(mode="numeric", length=N)
    y <- numeric(length=N)
    total <- 0
    for (i in safefor(1:N)) {
        total <- total + x[i]
        y[i] <- total
    }

    return (y)
}


massToSlr <- function(ts)
{
    # -360 Gton per mm SLE (sea level equivalent)
    #
    # 1000 mm per meter
    #
    y <- summation(ts[, "mass"]) / -360 / 1000
    ts[, "mass"] <- y
    colnames(ts)[2] <- "sealvl"

    return (ts)
}


# 2 on list of what Nathan wants
massPlotPred <- function(outfiles=F, assimctx=massctx)
{
    newDev("masspred", outfiles)

    gmst  <- loadUrbanAvg(lower=1850, upper=1879)
    model <- massBal(assimctx$fit$optim$bestmem, gmst)
    slr   <- massToSlr(model)

    plot(slr, type="l", xlab="Year", ylab="GIS attributed sea-level anomaly (m)")

    mtext("Figure 1. Predicted GIS contribution to sea level from fitted quadratic", outer=TRUE, side=1)

    if (outfiles) { finDev() }
}


massConfig <- function(assimctx=massctx)
{
    # ltzero is used where the constraint is < 0
    ltzero <- -1e-16

    #                         c       d
    assimctx$lbound <- c( -1000,  -1000)
    assimctx$ubound <- c(ltzero, ltzero)

    names(assimctx$lbound) <- names(assimctx$ubound) <- c("c", "d")

    obs  <- loadRignot(pure=F)
    gmst <- loadHadcrut(lower=1850, upper=1879)$annual
    gmst <- tsTrimForcing(gmst, obs)

    assimctx$obsonly <- obs[, "mass"]
    assimctx$gmst <- gmst[, "temp"]

    # TODO:  cannot eradicate this without specifying column in configAssim()
    setErrorAssim(assimctx, obs[, "mass error"])
}


massConfigAssim <- function(ar=0, obserr=T, useSSE=F, assimctx=massctx)
{
    massConfig(assimctx=assimctx)

    # plotting functions use this
    assimctx$units   <- c("Gton/C", "Gton/C^2")

    # assim.R uses these
    assimctx$modelfn <- massModel

    # get initial conditions from best fit model
    if (useSSE && is.null(assimctx$fit$optim$bestmem)) {
        massRunFit(assimctx=assimctx)
    }

    configAssim(assimctx, assimctx$fit$optim$bestmem, ar=ar, obserr=obserr)
}


massRunAssim <- function(nbatch=1000000, initial=T, assimctx=massctx)
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (initial || ncol(assimctx$chain) != length(init_mp) + length(init_sp)) {
        print("using initial scale")

        scale <- c(rep(2, length(init_mp) + length(init_sp)))
        scale <- abs(c(init_mp, init_sp)) / scale

    } else {
        print("using proposal matrix")
        scale <- assimProposalMatrix(assimctx$chain, mult=0.25)
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale)
}
