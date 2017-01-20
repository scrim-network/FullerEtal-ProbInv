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
# linassim.R

source("assim.R")
source("plot.R")


linModel <- function(mp, assimctx)
{
    y <- mp["a"] * assimctx$x + mp["b"]

    return (y)
}


linLogLik <- function(mp, sp, assimctx)
{
    y     <- assimctx$modelfn(mp, assimctx)
    llik  <- assimctx$expert_prior$dens(y[assimctx$obs_ind])
  
    return (llik)
}


if (!exists("linctx")) {
    linctx <- env()
}


linConfigAssim <- function(assimctx=linctx, prior="uniform")
{
    assimctx$modelfn <- linModel

    # linModel() uses this
    assimctx$x <- 1:10

    assimctx$lbound <- c(-10, -10)
    assimctx$ubound <- c( 10,  10)
    names(assimctx$lbound) <- names(assimctx$ubound) <- c("a", "b")

    assimctx$obs_ind <- 1
    min  <- -0.5
    max  <-  1.5
    mean <- (min + max) / 2
    switch (prior,
        uniform={
            assimctx$expert_prior <- uniformPrior(min, max)
        },
        beta={
            # a=2, b=3 taken from Lempert, Sriver, and Keller (2012)
            assimctx$expert_prior <- betaPrior(min, max, a=2, b=3)
        },
        normal={
            assimctx$expert_prior <- normPrior(mean, max)
        }, {
            stop("unknown prior in linConfigAssim()")
        })

    # get initial conditions from best fit model
    configAssim(assimctx, ar=0, obserr=F, llikfn=linLogLik)
}


linRunAssim <- function(nbatch=5e6, adapt=T, n.chain=1, initial=is.null(assimctx$chain), assimctx=linctx)
{
    scale <- NULL

    if (!adapt && !initial) {
        print("using proposal matrix")

        mult  <- 1.5
        scale <- assimProposalMatrix(assimctx$chain, mult=mult)
    }

    runAssim(assimctx, nbatch=nbatch, n.chain=n.chain, scale=scale, adapt=adapt)
}


if (!exists("prlinctx")) {
    prlinctx <- env()
}


linRunPredict <- function(nbatch=1e5, assimctx=linctx, prctx=prlinctx)
{
    prctx$assimctx <- assimctx

    print(colMean(assimctx$chain))

    prctx$prchain <- prmatrix(nbatch, xvals=assimctx$x)

    rows    <- 1:nbatch
    samples <- sample(burnedInd(assimctx$chain), nbatch, replace=T)
    
    for (i in rows) {
        prctx$prchain[i, ] <- linModel(assimctx$chain[samples[i], assimctx$mp_indices], assimctx)
    }
}


linPlotPredict <- function(prctx=prlinctx, outfiles=T, filetype="pdf")
{
    newDev("cmp_prior_linear", outfile=outfiles, width=8.5, height=11/2, filetype=filetype)

    x <- as.character(prctx$assimctx$obs_ind)
    xlim <- 1.25 * range(prctx$prchain[, x])
    priorPdfPlot(prctx$prchain, column=x, prior=prctx$assimctx$expert_prior, xlab="y", smoothing=0.5, xlim=xlim)

    caption <- paste("Figure n. PDF of y at x =", x)
    mtext(caption, outer=F, line=4, side=1, font=2)

    finDev()
}
