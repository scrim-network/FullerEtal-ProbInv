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
# pdfplot.R
# written by Robert W. Fuller on 090407

source("plot.R")


diagpdf <- function(
    outfiles=F,
    mcmcChain=assimctx$chain, units=assimctx$units,
    chainCols=1:ncol(mcmcChain),
    chartRows=2, chartCols=2,
    maxSteps=2000,
    assimctx=gisassimctx
    )
{
    rows <- nrow(mcmcChain)
    if (rows > maxSteps) {
        rows <- seq(from=1, to=rows, length.out=maxSteps)
    } else {
        rows <- 1:rows
    }

    for (i in 1:length(chainCols)) {

        newpage <- !((i - 1) %% (chartRows * chartCols))
        col <- chainCols[i]

        if (newpage) {
            pageno <- ceiling( i / (chartRows * chartCols))
            fname <- paste(sep="", "chains", format(pageno))
            newDev(fname, outfiles)
            par(mfrow=c(chartRows, chartCols), omi=c(0.5, 0, 0, 0))
        }

        ylab <- paste(sep="", colnames(mcmcChain)[col],
            ifelse(nchar(units[col]), paste(sep="", " (", units[col], ")"), ""))
        plot(x=rows, y=(mcmcChain[rows, col]), type="l", xlab="MCMC step", ylab=ylab)

        if (newpage) {
            mtext(paste(sep="", "Figure ", format(pageno), ". MCMC chains"),
                outer=T,
                side=1)
        }
    }

    if (outfiles) { finDev() }
}


pairplot <- function(
    outfiles=F,
    mcmcChain=assimctx$chain,
    maxPairs=1250,
    assimctx=gisassimctx
    )
{
    newDev("pairs", outfiles)

    if (nrow(mcmcChain) > maxPairs) {
        byPairs <- floor(nrow(mcmcChain) / maxPairs)
    } else {
        byPairs = 1
    }

    pairs( mcmcChain[ seq(1, nrow(mcmcChain), by=byPairs), ] )

    if (outfiles) { finDev() }
}


pdfplot <- function(
    outfiles=F,
    mcmcChain=assimctx$chain, units=assimctx$units,
    chainCols=1:ncol(mcmcChain),
    chartRows=ifelse(outfiles, 4, 3), chartCols=3,
    burnin=T,
    caption="Posterior probability density functions",
    figno=1,
    na.rm=F,
    smooth=T,
    prior=F,
    newdev=T,
    xlim=NULL,
    ylim=NULL,
    bound=F,
    truth=F,
    assimctx=allgrgisctx,
    ...
    )
{
    missing_xlim <- missing(xlim)
    missing_ylim <- missing(ylim)

    if (burnin) {
        burnin <- -burninInd(mcmcChain)
    } else {
        burnin <- 1:nrow(mcmcChain)
    }

    if (nrow(mcmcChain) > 100000) {
        breaks <- 50
    } else {
        breaks <- "Scott"
    }

    if (truth) {
        true_param <- c(assimctx$init_mp, assimctx$init_sp)
    }

    for (i in 1:length(chainCols)) {

        newpage <- !((i - 1) %% (chartRows * chartCols))
        col <- chainCols[i]

        if (newdev && newpage) {
            pageno <- ceiling( i / (chartRows * chartCols))
            fname <- paste(sep="", "pdfs", format(pageno))
            newDev(fname, outfiles, ...)
            #par(mfrow=c(chartRows, chartCols), omi=c(0.5, 0, 0, 0))
            par(mfrow=c(chartRows, chartCols))

            # bott, left, top, right
            par(mar=c(3, 3, 1, 1))
        }

        xlab=paste(sep="", colnames(mcmcChain)[col],
            ifelse(nchar(units[col]), paste(sep="", " (", units[col], ")"), ""))

        if (smooth) {
            epdf <- density(mcmcChain[burnin, col], na.rm=na.rm)
            if (missing_xlim) {
                xlim <- range(epdf$x)
            }
            if (missing_ylim) {
                ylim <- range(epdf$y)
            }
            if (bound) {
                if (!is.na(assimctx$lbound[col])) {
                    xlim <- c(assimctx$lbound[col], assimctx$ubound[col])
                } else {
                    if (!is.na(pmatch("rho", colnames(assimctx$chain)[col]))) {
                        xlim <- c(0, 1)
                    }
                }
            }

            emptyPlot(xlim, ylim, xlab, "Probability density", rhs=F, top=T)
            #lines(x=xlim, y=c(0, 0))
            #lines(epdf, lwd=1, lty="solid", col="black")

            # this code is broken. x[0] is NA or integer(0)
            #polygon(c(epdf$x, epdf$x[0]), c(epdf$y, epdf$y[0]), col="gray", border=NA)
            polygon(c(epdf$x, epdf$x[1]), c(epdf$y, epdf$y[1]), col="gray", border=NA)            

            if (prior) {
                prname <- paste(sep="", colnames(mcmcChain)[col], "_prior")
                if (exists(prname, envir=assimctx)) {
                    pr <- get(prname, envir=assimctx)
                    curve(exp(pr$dens(x)), add=T, lty="dashed")
                }
            }
        } else {
            hist(mcmcChain[burnin, col], breaks=breaks, probability=T, main="", xlab=xlab)

            box()

            # lwd=1.5 has inconsistent device-dependent results here
            lines(density(mcmcChain[burnin, col], na.rm=na.rm), lwd=1)
        }

        if (truth) {
            abline(v=true_param[i], lty="dotted", lwd=1)
        }

        if (newdev && newpage) {
            mtext(paste(sep="", "Figure ", format(pageno + figno - 1), ". ", caption),
                outer=T, side=1)
        }
    }

    if (newdev && outfiles) { finDev() }
}


# TODO:  generalize to replace pdfplot()?

# chainload("~/runs/prassim", oldnames=c("grinassimctx"), newnames=c("assim"))
# chainload("/tmp/prassim", oldnames=c("grinassimctx"), newnames=c("assim"))

assimPdfPlot <- function(
    basename="assim",
#    chartRows=ifelse(outfiles, 4, 3), chartCols=3,
    color=c("black", "red", "blue", "green", "purple", "yellow", "gray"),
    chartRows=3, chartCols=2,
    subcaption="Perfect model",
    caption=paste("Posterior probability densities for", subcaption),
    truth=T,
    newdev=T, outfiles=F
    )
{
    if (newdev) {
        newDev("cmpassim", outfiles)
        par(mfrow=c(chartRows, chartCols))
    }

    refassim   <- get(paste(sep="", basename, 1), envir=as.environment(".GlobalEnv"))
    burnin     <- -burninInd(refassim$chain)
    units      <- refassim$units
    #true_param <- c(refassim$init_mp, refassim$init_sp)

    for (col in 1:ncol(refassim$chain)) {
        densities <- list()
        truths    <- list()

        xlab=paste(sep="", colnames(refassim$chain)[col],
            ifelse(nchar(units[col]), paste(sep="", " (", units[col], ")"), ""))

        xlim <- ylim <- numeric()

        i <- 1
        while (T) {
            assimname <- paste(basename, i, sep="")
            if (!exists(assimname, envir=as.environment(".GlobalEnv"))) {
                break;
            }
            assim <- get(assimname, envir=as.environment(".GlobalEnv"))

            mcmcChain <- assim$chain[burnin, col]
            densities[[i]] <- density(mcmcChain, na.rm=na.rm)
            truths[[i]]    <- c(assim$init_mp, assim$init_sp)[col]

            xlim <- range(xlim, densities[[i]]$x)
            ylim <- range(ylim, densities[[i]]$y)

            i <- i + 1
        }

        if (!is.na(refassim$lbound[col])) {
            xlim <- c(refassim$lbound[col], refassim$ubound[col])
        } else {
            if (!is.na(pmatch("rho", colnames(refassim$chain)[col]))) {
                xlim <- c(0, 1)
            }
        }

        emptyPlot(xlim=xlim, ylim=ylim, xlab=xlab, ylab="Probability density")
        for (i in 1:length(densities)) {
            lines(densities[[i]], lwd=1, col=color[i]) # , lty=lty[i], col=col[i])
        }

        if (truth) {
            for (i in 1:length(truths)) {
                abline(v=truths[[i]], lty="dotted", lwd=1)
            }
            #abline(v=true_param[col], lty="dotted", lwd=1)
        }
    }    

    mtext(caption, outer=TRUE, side=1)

    if (newdev && outfiles) { finDev() }
}
