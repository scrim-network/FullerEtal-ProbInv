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


linConfigAssim <- function(assimctx=linctx)
{
    assimctx$modelfn <- linModel

    # linModel() uses this
    assimctx$x <- 1:10

    assimctx$lbound <- c(-10, -10)
    assimctx$ubound <- c( 10,  10)
    names(assimctx$lbound) <- names(assimctx$ubound) <- c("a", "b")

    assimctx$obs_ind      <- 1
    assimctx$expert_prior <- uniformPrior(-0.5, 1.5)

    # get initial conditions from best fit model
    configAssim(assimctx, ar=0, obserr=T, llikfn=linLogLik)
}


linRunAssim <- function(nbatch=1e7, adapt=T, initial=is.null(assimctx$chain), assimctx=linctx)
{
    scale <- NULL

    if (!adapt && !initial) {
        print("using proposal matrix")

        mult  <- 1.5
        scale <- assimProposalMatrix(assimctx$chain, mult=mult)
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale, adapt=adapt)
}


if (!exists("prlinctx")) {
    prlinctx <- env()
}


linRunPredict <- function(nbatch=1e6, assimctx=linctx, prctx=prlinctx)
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


linPlotPredict <- function(prctx=prlinctx, outfiles=T, filetype="png")
{
    newDev("cmp_prior_linear", outfile=outfiles, width=8.5, height=11/2, filetype=filetype)

    lwd <- 2
    lty <- c("solid")
    x   <- prctx$assimctx$obs_ind

    pdfPlots(
        chains=list(prctx$prchain),
        column=as.character(x),
        lty=lty,
        legendloc=NULL,
        #legends=fname,
        #col="black",
        col="red",
        burnin=F,
        #xlim=c(0, max(cictx$range)),
        #xlim=c(-0.2, 1.1),
        xlab="x",
        #yline=2,
        lwd=lwd
        )

    shadecol <- rgb(255, 0, 0, alpha=32, maxColorValue=255)
    pr       <- prctx$assimctx$expert_prior
    priorPlot(pr, shade=T, border=NA, col=shadecol)

    legend(
        "topright",
        legend=c("inversion", "prior"),
        col=c("red", shadecol),
        lty=c(lty, NA),
        lwd=c(lwd, NA),
        pch=c(NA, 15),
        bg=c(NA, shadecol)
        )

    caption <- paste("Figure n. PDF of y at x =", x)
    mtext(caption, outer=F, line=4, side=1, font=2)

    finDev()
}
