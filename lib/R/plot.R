# Copyright 2009, 2010, 2016 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# plot.R
# written by Robert W. Fuller on 090809

source("roblib.R")  # prmatrix(), safefor(), colMeans()


newDev <- function(fname, outfile, height=11, width=8.5, filetype="pdf", horiz=F)
{
    fname <- paste("../figures/", fname, sep="")
    if (outfile) {
        switch(filetype, 
            png={
                fname <- paste(sep="", fname, ".png")
                png(fname, height=height, width=width, units="in", res=300)
            },
            pdf={
                fname <- paste(sep="", fname, ".pdf")
                # a4 and a4r are options for Europe
                pdf(fname, paper=ifelse(horiz, "USr", "letter"), onefile=F, height=height, width=width)
            },
            eps={
                fname <- paste(sep="", fname, ".eps")
                postscript(fname, horizontal=horiz, onefile=F, height=height, width=width)
            }, {
                stop("unimplemented filetype in newDev()")
            })
    } else {
        dev.new()
    }

    # there are problems copying to postscript and pdf devices;  hence, copy from
    # postscript and pdf devices;  this has to be enabled
    #
    if (!dev.interactive() && interactive()) {
        dev.control("enable")
    }

    # bottom, left, top, right (margins in inches)
    par(omi=c(0.25, 0.25, 0.25, 0.25))
}


# finalize devices
finDev <- function()
{
    devs <- dev.list()
    for (src in devs) {
        dev.set(which=src)

        # if printing to the device, close it
        if (!dev.interactive()) {

            # if there is a user, display what is printed
            if (interactive()) {
                dev.copy(device=dev.new)
            }

            dev.off(which=src)
        }
    }
}


lmPlot <- function(x, y, ...)
{
    plot(x, y, ...)
    abline(coef(lm(y ~ x)), lwd=2)
}


normHist <- function(x, ...)
{
    out <- hist(x, ...)
    box()

    # this is the standard, but ugly solution
    #mu <- mean(x)
    #sigma <- sd(x)
    #curve(dnorm(x, mean=mu, sd=sigma), add=T)

    xnorm <- seq(min(x), max(x), length.out=1000)
    lines(xnorm, dnorm(xnorm, mean=mean(x), sd=sd(x)))

    return (out)
}


tsColLines <- function(x, ...)
{
    for (col in safefor(2:ncol(x))) {
        lines(x[, 1], x[, col], ...)
    }
}


xshade <- function(xleft, xright, col=rgb(0.8, 0.8, 0.8))
{
    plotxy  <- par()$usr
    ybottom <- plotxy[3]
    ytop    <- plotxy[4]

    if (missing(xleft)) {
        xleft <- plotxy[1]
    }
    if (missing(xright)) {
        xright <- plotxy[2]
    }

    rect(xleft, ybottom, xright, ytop, col=col, border=F)
}


# avoid the problem of having to switch between low and high level
# plotting functions as objects are added and removed from a plot
#
emptyPlot <- function(xlim, ylim, xlab=NULL, ylab=NULL, rhs=T, box=T, top=F, xpad=T, line=2)
{
    plot.new()
    plot.window(xlim, ylim, xaxs=ifelse(xpad, "r", "i"))
    axis(1)
    axis(2)
    if (rhs) {
        axis(4, labels=F, tcl=-0.25)
    }
    if (top) {
        # positive values for tcl put the tickmarks inside the plot
        axis(3, labels=F, tcl=-0.10)
    }
    title(xlab=xlab, ylab=ylab, line=line)
    if (box) {
        box()
    }
}


labelPlot <- function(letter, line=3, where="topleft")
{
    switch(where,
        topleft={
            at <- par("usr")[4]
        },
        left={
            at <- NA
        }, {
            stop("unknown location in labelPlot()")
        })

    # left:  side=2
    # horizontal:  las=1
    # above axis label:  line=3
    # bold:  font=2
    # locate at y-axis:  at=par("usr")[4]
    #
    # font=2
    mtext(letter, side=2, las=1, line=line, at=at, cex=1.25, font=2)
}


tsErrorBars <- function(ts, shade=T, lines=!shade, xbeam=F, ibeam=!xbeam, obscol=2, lwd=1, col="gray", tick=0.5)
{
    x     <- ts[, "time"]
    upper <- ts[, obscol] + ts[, "error"]
    lower <- ts[, obscol] - ts[, "error"]

    if (lines) {
        segments(x, lower, x, upper, col=col, lwd=lwd)
        if (ibeam) {
            segments(x - tick, c(lower, upper), x + tick, c(lower, upper), col=col, lwd=lwd)
            #segments(x - tick, upper, x + tick, upper, col=col, lwd=lwd)
        }
        if (xbeam) {
            segments(x - tick, ts[, obscol], x + tick, ts[, obscol], col=col, lwd=lwd)
        }
    }

    if (shade) {
        polygon(c(x, rev(x), x[1]), c(lower, rev(upper), lower[1]), col=col, border=NA)
    }
}


plotInterval <- function(pred, mean=rep(0, length(noise)), ci=.68, ...)
{
    lower <- (1 - ci) / 2
    probs <- c(lower, 1 - lower)

    xvals <- as.numeric(colnames(pred))
    interval <- colQuantile(pred, probs)
    interval[, 1] <- interval[, 1] + mean
    interval[, 2] <- interval[, 2] + mean
    interval <- cbind(xvals, interval)
    tsColLines(interval, ...)
}


colQuantile <- function(x, probs=c(0.025, 0.975), ...)
{
    cols <- ncol(x)
    q    <- prmatrix(cols, probs)
    for (col in safefor(1:cols)) {
        q[col, ] <- quantile(x[, col], probs=probs, ...)
    }

    return (q)
}


prProbExceed <- function(chain=prallgrgisctx$prchain, threshold=(48*2.54/100))
{
    cols <- ncol(chain)
    p_exceed <- numeric(cols)

    for (i in 1:cols) {
        cdf <- ecdf(chain[, i])
        p_exceed[i] <- 1 - cdf(threshold)
    }

    names(p_exceed) <- colnames(chain)

    return (p_exceed)
}


# subtracting one year from the prior year, as prChainRate() does, subtracts
# some of the auto-correlated noise;  hence, prChainRate() produces
# tighter credible intervals
#
prChainRate <- function(chain=prallgrgisctx$prchain)
{
    rows <- nrow(chain)
    cols <- ncol(chain)

    xvals <- attr(chain, "xvals")
    xvals <- xvals[ -cols ]
    rates <- prmatrix(rows, xvals)

    # apply may make more sense;  see ~/src/mici/R/calib.R
    for (i in 1:rows) {
        rates[i, ] <- chain[i, 2:cols] - chain[i, 1:(cols-1)]
    }

    return (rates)
}


ciCalc <- function(..., xvals=attr(chains[[1]], "xvals"), probs=c(0.025, 0.975), chains=list(...))
{
    cictx <- env()

    cols <- which(attr(chains[[1]], "xvals") %in% xvals)
    cictx$cols  <- cols
    cictx$xvals <- xvals

    cictx$means <- list()
    cictx$cis   <- list()
    cictx$range <- numeric()

    for (i in 1:length(chains)) {
        if (1 == length(cols)) {
            cictx$means[[i]] <- mean(chains[[i]][, cols])
            cictx$cis[[i]]   <- quantile(chains[[i]][, cols], probs=probs)
        } else {
            cictx$means[[i]] <- colMeans(chains[[i]][, cols])
            cictx$cis[[i]]   <- colQuantile(chains[[i]][, cols], probs=probs)
        }
        cictx$range <- range(cictx$range, cictx$cis[[i]])
    }

    return (cictx)
}


ciPlot <- function(cictx, col="red", lty=c("solid", "dashed"), lwd=2, plotmeans=T, plotci=T)
{
    k <- 1
    for (i in 1:length(cictx$means)) {
        if (plotmeans) {
            lines(cictx$xvals, cictx$means[[i]], col=col[i], lty=lty[k], lwd=lwd)
            k <- k + 1
        }

        if (plotci) {
            for (j in 1:ncol(cictx$cis[[i]])) {
                lines(cictx$xvals, cictx$cis[[i]][, j], col=col[i], lty=lty[k], lwd=lwd)
            }
            k <- k + 1
        }
    }
}


cdfCalc <- function(..., column=NULL, chains=list(...))
{
    xlim <- numeric()
    cdfs <- list()

    for (i in 1:length(chains)) {
        if (is.null(column)) {
            chain <- chains[[i]]
        } else {
            chain <- chains[[i]][, column]
        }
        cdf  <- ecdf(chain)
        xlim <- range(chain, xlim)
        cdfs[[i]] <- cdf
    }

    return (env(xlim=xlim, cdfs=cdfs))
}


cdfPlot <- function(cdfctx, col, lty, lwd=2)
{
    for (i in 1:length(cdfctx$cdfs)) {

        # cannot pass the expression cdfs[[i]] directly to curve()
        # because it is not an expression of x or a SIMPLE function name
        #
        cdf <- cdfctx$cdfs[[i]]
        curve(cdf, add=T, col=col[i], lty=lty[i], lwd=lwd)
    }
}


cdfPlotWindow <- function(cdfctx,
    col, lty, lwd=2,
    xlab=NULL, ylab="Cumulative density",
    xlim=NULL)
{
    if (is.null(xlim)) {
        xlim <- cdfctx$xlim
    }

    ylim <- c(0, 1)

    plot.new()
    plot.window(xlim, ylim, xaxs="i")

    axis(1) # bottom
    axis(2) # left

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    # right
    axis(4, labels=F, tcl=-0.25)

    title(xlab=xlab, ylab=ylab, line=2)
    box()

    cdfPlot(cdfctx, col, lty, lwd)
}


cdfPlots <- function(
    ..., col, lty,
    lwd=2,
    xlab=gmslLab(column),
    ylab="Cumulative density",
    xlim=NULL,
    column=NULL,
    chains=list(...)
    )
{
    cdf <- cdfCalc(chains=chains, column=column)
    cdfPlotWindow(cdf, col, lty, lwd, xlab, ylab, xlim)
}


pdfCalc <- function(..., column=NULL, burnin=T, na.rm=F, chains=list(...))
{
    xlim <- ylim <- numeric()
    densities <- list()
    means     <- list()

    for (i in 1:length(chains)) {

        mcmcChain <- chains[[i]]
        if (burnin) {
            ind <- -burninInd(mcmcChain)
        } else {
            ind <- 1:nrow(mcmcChain)
        }
        if (is.null(column)) {
            mcmcChain <- mcmcChain[ind, ]
        } else {
            mcmcChain <- mcmcChain[ind, column]
        }

        densities[[i]] <- density(mcmcChain, na.rm=na.rm)
        means[[i]]     <-    mean(mcmcChain, na.rm=na.rm)

        xlim <- range(xlim, densities[[i]]$x)
        ylim <- range(ylim, densities[[i]]$y)
    }

    return (env(xlim=xlim, ylim=ylim, densities=densities, means=means))
}


pdfPlot <- function(pdfctx, col, lty, lwd, plotmeans=F)
{
    for (i in 1:length(pdfctx$densities)) {
        lines(pdfctx$densities[[i]], lty=lty[i], lwd=lwd, col=col[i])
        if (plotmeans) {
            abline(v=pdfctx$means[[i]], lty="dotted", lwd=lwd, col=col[i])
        }
    }
}


pdfPlotWindow <- function(pdfctx,
    col, lty, lwd=2,
    xlab=NULL, ylab="Probability density",
    xlim=NULL, ylim=NULL, plotmeans=F, yline=1)
{
    if (is.null(xlim)) {
        xlim <- pdfctx$xlim
    }
    if (is.null(ylim)) {
        ylim <- pdfctx$ylim
    }

    plot.new()
    plot.window(xlim, ylim, xaxs="i")

    # bottom
    axis(1)

    # left
    ticks <- axTicks(2)
    ticks <- c(0, last(ticks))
    axis(2, at=ticks)

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    # put the y-axis label between the two tick marks, by default
    title(ylab=ylab, line=yline)
    title(xlab=xlab, line=2)
    box()

    pdfPlot(pdfctx, col, lty, lwd, plotmeans)
}


pdfPlots <- function(
    ..., legends, col, lty,
    column=as.character(2100),
    chains=list(...),
    lwd=2,
    xlab=gmslLab(column),
    ylab="Probability density",
    burnin=T, na.rm=F,
    plotmeans=F,
    xlim=NULL, ylim=NULL,
    legendloc="topright",
    truth=F,
    trueval=grinassimctx$true_predict,
    yline=1
    )
{
    pdfctx <- pdfCalc(chains=chains, column=column, burnin=burnin, na.rm=na.rm)

    pdfPlotWindow(pdfctx, col, lty, lwd, xlab, ylab, xlim, ylim, plotmeans, yline)

    if (truth) {
        if (plotmeans) {
            abline(v=trueval, lty="solid",  lwd=lwd+0.5)
        } else {
            abline(v=trueval, lty="dotted", lwd=1)
        }
        #abline(v=-0.4371877, lty="dotted", lwd=1)
        #abline(v=2.627464, lty="dotted", lwd=1)
        #abline(v=0.592957455497411, lty="dotted", lwd=1)
        #abline(v=0.593337748954618, lty="dotted", lwd=1)
    }

    if (!is.null(legendloc)) {
        legend(
            legendloc,
            legend=legends,
            lty=lty,
            col=col,
            lwd=rep(lwd, length(lty))
            )
    }
}
