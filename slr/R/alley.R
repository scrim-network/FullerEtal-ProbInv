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
# alley.R
# written by Robert W. Fuller on 091026
# assimilate the Alley Prior, Figure 6.13, CCSP chapter 6, 2009

source("data.R")
source("assim.R")
source("predict.R")
loadLibrary("DEoptim")


if (!exists("alleyctx")) {
    alleyctx <- env()
}


alleyPlotPrior <- function(temp_bias=0)
{
    fig613 <- loadAlley(temp_bias)

    # pch=16, 19, 20 (tried 3 and 10 as well);  tried 2-3 for cex as well
    points(
        x=fig613[, "temp"], y=fig613[, "SLE"],
        col="seagreen", lwd=2, pch=16, cex=1.5
        )
    rect(
        fig613[, "temp_lbound"], fig613[, "SLE_lbound"],
        fig613[, "temp_ubound"], fig613[, "SLE_ubound"],
        border="seagreen", lwd=2
        )
}


alleyPlotPredict <- function(
    prctx=alleyctx, prchain=prctx$prchain,  
    outfiles=F, newdev=T, temp_bias=0,
    ...
    )
{
    if (newdev) {
        newDev("seqpred", outfiles)
    }

    prqplotsub(
        prchain,
        legends=c("Mean", "95% CI", "Alley"),
        col=c(rep("mediumblue", 2), "seagreen"),
        lty=c("solid", "dashed", "solid"),
        prctx=prctx,
        obs=F,
        shade=F,
        xcaption="Scaled Greenland temperature anomaly (C)",
        ycaption="Equilibrium sea-level anomaly (m)",
        caption="Figure 1. Posterior predictive for GIS equilibrium sea level",
        newdev=F,
        outfiles=outfiles,
        ...
        )
    alleyPlotPrior(temp_bias)

    if (newdev && outfiles) { finDev() }
}


alleyPlotSeqGis <- function(prctx=prallgrgisctx, outfiles=F, newdev=T, ...)
{
    mcmcChain <- prctx$assimctx$chain

    # for consistency, use mean
    # TODO:  gis_temp could be a fixed parameter
    bias <- mean(mcmcChain[ burnedInd(mcmcChain), "gis_temp"])

    alleyPlotPredict(
        prctx=prctx, prchain=prctx$seq_gischain, temp_bias=-bias,
        outfiles=outfiles, newdev=newdev,
        ...
        )
}


alleyPlotFit <- function(assimctx=alleyctx)
{
    if (is.null(assimctx$fit$optim$bestmem)) {
        print("running fit for initial model parameters")
        alleyRunFit(assimctx=assimctx)
    }

    emptyPlot(
        xlim=c(-11.75, 6.25), ylim=c(-4.75, 8),
        xlab="Temperature anomaly (C)", ylab=slrGreenlandLab()
        )
    alleyPlotPrior()
    curve(alleyModel(x, mp=assimctx$fit$optim$bestmem), add=T, lwd=2)
}


alleyModel <- function(temp, mp)
{
    sle  <- numeric(length=length(temp))
    ind  <- (temp <= 0)
    sle[ind] <- mp["c"] * temp[ind]
    ind  <- !ind
    sle[ind] <- mp["c"] * temp[ind] + mp["d"] * (temp[ind])^2

    return (sle)
}


alleyMinFn <- function(mp, assimctx)
{
    if (assimctx$herr) {
        xvals <- pindex("xtemp", mp)
        err <- sse(assimctx$obsonly, alleyModel(xvals, mp)) + sse(xvals, assimctx$xobs)
    } else {
        err <- sse(assimctx$obsonly, alleyModel(assimctx$xobs, mp))
    }

    if (is.na(err)) {
        return (Inf)
    }

    return (err)
}


alleyConfig <- function(assimctx=alleyctx, omitneg=F, stderror=4, herr=T)
{
    fig613 <- loadAlley()

    # cannot have a zero on the diagonal or matrix will be singular in dmvnorm()
    if (omitneg) {
        fig613 <- fig613[ (fig613[, "temp"] > 0), ]
    } else {
        fig613 <- fig613[ (fig613[, "temp"] != 0), ]
    }
    
    assimctx$obsonly <- fig613[, "SLE"]
    assimctx$xobs <- fig613[, "temp"]

    # treat height of box as 4 standard errors per Nathan
    setErrorAssim(assimctx, (fig613[, "SLE_ubound"] - fig613[, "SLE_lbound"]) / stderror)

    #                    c  d
    assimctx$lbound <- c(0, 0) 
    assimctx$ubound <- c(2, 2)
    names(assimctx$lbound) <- names(assimctx$ubound) <- c("c", "d")

    if (herr) {
        names <- names(assimctx$lbound)
        assimctx$lbound    <- append(assimctx$lbound, fig613[, "temp_lbound"])
        assimctx$ubound    <- append(assimctx$ubound, fig613[, "temp_ubound"])
        assimctx$xtemp_err <- setErrorAssim(env(),   (fig613[, "temp_ubound"] - fig613[, "temp_lbound"]) / stderror)
        names(assimctx$lbound) <- names(assimctx$ubound) <- append(names, paste("xtemp", 1:length(assimctx$xobs), sep=""))
    }

    assimctx$herr <- herr
}


alleyRunFit <- function(assimctx=alleyctx, ...)
{
    alleyConfig(assimctx=assimctx, ...)

    itermax <- ifelse(assimctx$herr, 500, 100)
    control <- list(CR=1.0, itermax=itermax)

    assimctx$fit <- named_DEoptim(
        FUN=alleyMinFn,
        lower=assimctx$lbound, upper=assimctx$ubound,
        control=control,
        assimctx=assimctx
        )
}


alleyLogLik <- function(mp, sp, assimctx)
{
    llik1 <- logLik(mp, sp, assimctx)

    xvals <- pindex("xtemp", mp)
    res   <- assimctx$xobs - xvals
    llik2 <- llik_obs(res, sp, assimctx$xtemp_err)

    #print(paste("llik1 is", llik1, "llik2 is", llik2))

    return (llik1 + llik2)
}


alleyConfigAssim <- function(assimctx=alleyctx, herr=T, useSSE=F, ...)
{
    # get initial conditions from best fit model
    bestmem <- assimctx$fit$optim$bestmem
    if (useSSE && (is.null(bestmem) || length(bestmem) != length(assimctx$lbound))) {
        print("running fit for initial model parameters")
        alleyRunFit(assimctx=assimctx, herr=herr, ...)
    } else {
        alleyConfig(assimctx=assimctx, herr=herr, ...)
        rmif(fit, envir=assimctx)
    }

    assimctx$modelfn <- function(mp, assimctx) alleyModel(assimctx$xobs, mp)
    assimctx$units   <- c("m/C", "m/C^2", rep("C", length(assimctx$xobs)))

    if (assimctx$herr) {
        #assimctx$units <- append(assimctx$units, rep("C", length(assimctx$xobs)))
        configAssim(assimctx, assimctx$fit$optim$bestmem, ar=0, obserr=T, llikfn=alleyLogLik)
    } else {
        configAssim(assimctx, assimctx$fit$optim$bestmem, ar=0, obserr=T)
    }
}


alleyRunAssim <- function(nbatch=1000000, initial=is.null(assimctx$chain), assimctx=alleyctx)
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    if (initial || ncol(assimctx$chain) != length(init_mp) + length(init_sp)) {
        print("using initial scale")

        # avoid a c parameter of 0 if omitneg=T
        init_mp[ iszero(init_mp) ] <- 0.2

        if (assimctx$herr) {
            init_mp[ 3:length(init_mp) ] <- assimctx$xtemp_err$error
        }

        # acceptance rate tuning
        # 10 ~ .5
        #  9 ~ .47
        #  8 ~ .42
        #  7 ~ .37
        #  6 ~ .33
        #  5 ~ .28
        scale <- c(rep(8, length(init_mp) + length(init_sp)))
        scale <- abs(c(init_mp, init_sp)) / scale

    } else {
        print("using proposal matrix")

        # accept=0.75723 for 2x 100K at 0.25
        scale <- assimProposalMatrix(assimctx$chain, mult=0.90)
    }

    runAssim(assimctx, nbatch, scale)
}


alleyRunPredict <- function(
    ...,
    assimctx=alleyctx, prctx=alleyctx
    )
{
    runPredict(
        assimctx=assimctx, prctx=prctx,
        forcings=NULL, modelfn=function(temp, parms) cbind(temp, alleyModel(temp, parms$mp)),
        noisefn=list(noise_zeros),
        outnames=c("prchain"),
        xvals=seq(-12, 7, by=.25),
        ...
        )
}
