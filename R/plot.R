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
loadLibrary("KernSmooth")
loadLibrary("RColorBrewer")


newDev <- function(fname, outfile, single=T, height=9.7, width=ifelse(single, 3.5, 7.2), units="in", filetype="pdf", horiz=F, mar=c(3, 3, 0.25, 0.25))
{
    fname <- paste("../figures/", fname, sep="")
    if (outfile) {
        switch(filetype, 
            png={
                fname <- paste(sep="", fname, ".png")
                png(fname, units=units, res=ifelse(units=="in", 300, NA),   height=height, width=width, pointsize=7)
            },
            pdf={
                fname <- paste(sep="", fname, ".pdf")
                # a4 and a4r are options for Europe
                pdf(fname, paper=ifelse(horiz, "USr", "letter"), onefile=F, height=height, width=width, pointsize=7)
            },
            eps={
                fname <- paste(sep="", fname, ".eps")
                postscript(fname, horizontal=horiz,              onefile=F, height=height, width=width, pointsize=7)
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

    # bottom, left, top, right (margins in lines of text)
    par(oma=rep(0, 4))
    par(mar=mar)
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
                dev.copy(device=dev.new, width=par("din")[1], height=par("din")[2], units="in", pointsize=par("ps"))
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


labelPlot <- function(letter, line=3, where="topleft", cex=8/7)
{
    switch(where,
        topleft={
            at <- par("usr")[4]  # locate at y-axis
        },
        log={
            at <- 1
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
    # cex:  8/7 is the multiplier from 7 point to 8 point for Nature
    #
    mtext(letter, side=2, las=1, line=line, at=at, cex=cex, font=2)
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


plotUnits <- function(inches, horiz=F)
{
    plt_inches <- par("pin")[ ifelse(horiz, 1, 2) ]  # pin is (width, height)
    usr        <- par("usr")  # "usr" is (x1, x2, y1, y2)
    plt_units  <- ifelse(horiz, usr[2] - usr[1], usr[4] - usr[3])
    units      <- inches * plt_units / plt_inches
   #print(paste(inches, plt_units, plt_inches, units))
    return (units)
}


plotArrowX <- function(xlim, label, y=plotUnits(length/2), length=0.10, code=3, col="black", lty="solid", lwd=1, offset=0.50, font=1)
{
    arrows(xlim[1], y, xlim[2], y, length=length, code=code, col=col, lty=lty, lwd=lwd)
    text(mean(xlim), y + offset * plotUnits(length), labels=label, pos=3, col=col, font=font)
}


# TODO:  this could use apply
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


# recommend the qualitative palettes that are color blind friendly:  Dark2, Paired or Set2
plotGetColors <- function(n=3, alpha=255, pal="Set2")
{
    # display.brewer.all(type="qual", colorblindFriendly=T)
    pal <- brewer.pal(n=n, name=pal)

    # append the alpha value
    pal <- paste(pal, toupper(as.hexmode(alpha)), sep="")

    return (pal)
}


plotGetAlphaColors <- function(col, n, min=0, max=255)
{
    alphas <- floor(seq(min, max, length.out=n))
    pal    <- paste(ifelse(nchar(col) > 7, substring(col, 1, 7), col), toupper(as.hexmode(alphas)), sep="")

    return (pal)
}


priorPlot <- function(pr, col="gray", lty="dotted", lwd=2, xlim=par("usr")[1:2], shade=F, n=1001L, border=NA)
{
    x <- seq(from=xlim[1], to=xlim[2], length.out=n)
    y <- pr$dens(x, log=F)
    if (shade) {
        polygon(c(x, x[1]), c(y, y[1]), col=col, lty=lty, lwd=lwd, border=border)
    } else {
        lines(x=x, y=y, col=col, lty=lty, lwd=lwd)
    }
}


priorPdfPlot <- function(chain, column, prior, xlim=NULL, ylim=NULL, xlab, lty="solid", lwd=2, legends=c("prior", "inversion"), col="#FF0000", shadecol="#FF000020", smoothing=1, new=F)
{
    pdfPlots(
        chains=list(chain),
        column=column,
        lty=lty,
        legendloc=NULL,
        col=col,
        burnin=F,
        xlab=xlab,
        lwd=lwd,
        smoothing=smoothing,
        xlim=xlim,
        ylim=ylim,
        new=new
        )

    priorPlot(prior, shade=T, border=NA, col=shadecol)

    if (!is.null(legends)) {
        legend(
            "topright",
            legend=legends,
            col=c(shadecol, col),
            lty=c(NA, lty),
            lwd=c(NA, lwd),
            pch=c(15, NA),
            pt.cex=2,
            bg="white"
            )
    }
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


cdfCalc <- function(..., column=NULL, chains=list(...), survival=F)
{
    xlim   <- numeric()
    cdfs   <- list()
    yticks <- numeric()

    for (i in 1:length(chains)) {
        if (is.null(column)) {
            chain <- chains[[i]]
        } else {
            chain <- chains[[i]][, column]
        }
        xlim   <- range(chain, xlim)
        x      <- sort(chain, decreasing=!survival)
        n      <- length(chain)
       #y      <- seq(1, 1/n, by= -1/n)
        y      <- seq(1, 1/n, length.out=n)
       #yticks <- max(yticks, floor(log10(n)))
        yticks <- 3

        cdfs[[i]] <- list(x=x, y=y)
    }

    return (env(xlim=xlim, cdfs=cdfs, survival=survival, yticks=yticks))
}


cdfPlot <- function(cdfctx, col, lty, lwd=2)
{
    for (i in 1:length(cdfctx$cdfs)) {
        lines(cdfctx$cdfs[[i]], col=col[i], lty=lty[i], lwd=lwd)  # , type="s")
    }
}


cdfPlotWindow <- function(cdfctx,
    col, lty, lwd=2,
    xlab=NULL, ylab="Cumulative density",
    xlim=NULL, ylim=NULL, log=F, yline=ifelse(log, 2.5, 2),
    new=F)
{
    if (is.null(xlim)) {
        xlim <- cdfctx$xlim
    }

    if (is.null(ylim)) {
        ylim <- c(ifelse(log, 10^(-cdfctx$yticks), 0), 1)
    }

    if (!new) {
        plot.new()
        plot.window(xlim, ylim, xaxs="i", log=ifelse(log, "y", ""))
    }

    # bottom
    axis(1)

    # left
    if (log) {
        yticks <- (-cdfctx$yticks):0
        axis(2, at=10^yticks, label=parse(text=paste("10^", yticks, sep="")), las=1)
    } else {
        axis(2)
    }

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    # right
    if (log) {
        axis(4, at=10^yticks, labels=F, tcl=-0.25)
    } else {
        axis(4, labels=F, tcl=-0.25)
    }

    title(xlab=xlab, line=2)
    title(ylab=ylab, line=yline)
    box()

    cdfPlot(cdfctx, col, lty, lwd)
}


cdfPlots <- function(
    ..., col, lty,
    lwd=2,
    xlab=gmslLab(column),
    ylab=ifelse(survival, "Survival (1-CDF)", "Cumulative density"),
    xlim=NULL, ylim=NULL, column=NULL,
    survival=F,
    log=F, yline=ifelse(log, 2.5, 2),
    new=F,
    chains=list(...)
    )
{
    cdf <- cdfCalc(chains=chains, column=column, survival=survival)
    cdfPlotWindow(cdf, col, lty, lwd, xlab, ylab, xlim, ylim, log, yline, new)
}


pdfDensity <- function(x, kernel="box", smoothing=1)
{
    bw  <- dpik(x, kernel=kernel) * smoothing
    d   <- bkde(x, kernel=kernel, bandwidth=bw, truncate=T, range.x=range(x))
    d$x <- c(min(x), d$x, max(x))
    d$y <- c(     0, d$y, 0     )

    return (d)
}


pdfCalc <- function(..., column=NULL, burnin=T, na.rm=F, smoothing=rep(1, length(chains)), chains=list(...))
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

       #densities[[i]] <- density(mcmcChain, na.rm=na.rm)
        densities[[i]] <- pdfDensity(mcmcChain, smoothing=smoothing[i])
        means[[i]]     <-       mean(mcmcChain, na.rm=na.rm)

        # accumulating xlim/ylim, so keep feeding back in
        xlim <- range(xlim, densities[[i]]$x)
        ylim <- range(ylim, densities[[i]]$y)
    }

    return (env(xlim=xlim, ylim=ylim, densities=densities, means=means))
}


pdfPlot <- function(pdfctx, col, lty, lwd, plotmeans=F, reverse=F)
{
    for (i in 1:length(pdfctx$densities)) {
        if (reverse) {
            lines(pdfctx$densities[[i]]$y, pdfctx$densities[[i]]$x, lty=lty[i], lwd=lwd, col=col[i])
        } else {
            lines(pdfctx$densities[[i]], lty=lty[i], lwd=lwd, col=col[i])
        }
        if (plotmeans) {
            abline(v=pdfctx$means[[i]], lty="dotted", lwd=lwd, col=col[i])
        }
    }
}


plotDensityAxis <- function(side=2, ...)  # left
{
    ticks <- axTicks(side)
    ticks <- c(0, last(ticks))
    axis(side=side, at=ticks, ...)
}


pdfPlotWindow <- function(pdfctx,
    col, lty, lwd=2,
    xlab=NULL, ylab="Probability density",
    xlim=NULL, ylim=NULL, plotmeans=F, yline=1, new=F)
{
    if (is.null(xlim)) {
        xlim <- pdfctx$xlim
    }
    if (is.null(ylim)) {
        ylim <- pdfctx$ylim
    }

    if (!new) {
        plot.new()
        plot.window(xlim, ylim, xaxs="i")
    }

    axis(1)            # bottom
    plotDensityAxis()  # left

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
    smoothing=rep(1, length(chains)),
    yline=1,
    new=F
    )
{
    pdfctx <- pdfCalc(chains=chains, column=column, burnin=burnin, na.rm=na.rm, smoothing=smoothing)

    pdfPlotWindow(pdfctx, col, lty, lwd, xlab, ylab, xlim, ylim, plotmeans, yline, new)

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


plotLayout <- function(...)
{
    cex <- par("cex")
    ps  <- par("ps")
    n   <- layout(...)
    par(cex=cex, ps=ps)

    return (n)
}


pairPlot <- function(..., units=NULL, topColumn=NULL, sideColumn=NULL, legends=NULL, title="Prior", label=NULL,
    col=plotGetColors(length(chains)), shadecol=plotGetColors(length(chains), 48), ccol, pdfcol=col, lwd=2,
    burnin=T,
    xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, xline=2, yline=2,
    points=25000,
    smoothing=rep(1, length(chains)),
    layout=T, mar=c(3, ifelse(is.null(label), 3, 4)), method="points", plotfn=NULL,
    chains=list(...)
    )
{
    lty      <- rep("solid", length(chains))
    topPdf   <- pdfCalc(chains=chains, column=topColumn,  burnin=burnin, smoothing=smoothing)
    sidePdf  <- pdfCalc(chains=chains, column=sideColumn, burnin=burnin, smoothing=smoothing)

    bottom   <- mar[1]
    left     <- mar[2]

    if (is.null(xlim)) {
        xlim <- topPdf$xlim
    }
    if (is.null(ylim)) {
        ylim <- sidePdf$xlim
    }

    if (layout) {
        plotLayout(matrix(1:4, nrow = 2, byrow = T), widths = c(10, 3), heights = c(3, 10))
    }


    # top PDF
    #

    # bottom=0.25 gives space between the PDF and the main plot
    #    top=0.50 gives space for the upper y-axis tick label and the panel label
    par(mar=c(0.25, left, 0.50, 0))  # bottom, left, top, right
    plot.new()
    plot.window(xlim=xlim, ylim=topPdf$ylim, xaxs="i")
   #plotDensityAxis(labels=F)
    axis(side=2, at=c(0, topPdf$ylim), labels=F)
    title(ylab="PDF",          line=0)
   #title(ylab="Density",      line=0)
   #title(ylab="Prob density", line=2)
    pdfPlot(topPdf, col=pdfcol, lty=lty, lwd=lwd)

    if (!is.null(label)) {
        labelPlot(label)
    }


    # legend
    #

    # bottom, left, top, right
    #    top=0.50 gives space for the label
    # should not matter because legend should be centered
   #par(mar=c(0, 0, 0.50, 0))
    par(mar=rep(0, 4))
    plot.new()
    # cex=5/7 would reduce font to minimum size for Nature (assuming pointsize is 7)
    # y.intersp=0.5 mashes the symbols together
    # title.adj=0.1 left justifies the legend title with a little gap between the border and the title
    legend(
        "center",
        legend=legends,
        title=title,
        title.adj=0.1,
       #bty="n",
        fill=col,
        border=col
        )


    # main plot
    #

    # bottom, left, top, right
    par(mar=c(bottom, left, 0, 0))
    plot.new()
    plot.window(xlim=xlim, ylim=ylim, xaxs="i", yaxs="i")

    # the last tick mark's label can overlap with the probability density's first tick mark label
    ticks <- axTicks(1)
   #ticks <- ticks[ -length(ticks) ]
    axis(1, at=ticks)  # bottom
    ticks <- axTicks(2)
   #ticks <- ticks[ -length(ticks) ]
    axis(2, at=ticks)  # left
    axis(3, labels=F, tck=-0.01)  # top
    axis(4, labels=F, tck=-0.01)  # right

    for (i in 1:length(chains)) {
        samples <- sample(burnedInd(chains[[i]]), points, replace=T)
        x <- chains[[i]][samples, topColumn]
        y <- chains[[i]][samples, sideColumn]
        switch (method,
            contours={
                z <- cbind(x, y)
                d <- bkde2D(z, bandwidth=( c(dpik(x), dpik(y)) * smoothing[i] ))  #, range.x=list(range(x), range(y)))
                levels <- seq(0.05, 0.95, length.out=32) * max(d$fhat)
               #pal    <- plotGetAlphaColors(col=col[i], n=( length(levels) - 1 ), min=24, max=128)
                .filled.contour(d$x1, d$x2, d$fhat, col=ccol[[i]], levels=levels)
               #image(          d$x1, d$x2, d$fhat, col=ccol[[i]], breaks=levels, add=T)
            },
            points={
                points(x, y, pch=20, col=shadecol[i], cex=0.5)
            },
            plotfn={
                plotfn(samples=chains[[i]][samples, ], i=i, topColumn=topColumn, sideColumn=sideColumn, col=col, shadecol=shadecol, ccol=ccol)
            },
            outline={
                z <- cbind(x, y)
                d <- bkde2D(z, bandwidth=( c(dpik(x), dpik(y)) * smoothing[i] ))  #, range.x=list(range(x), range(y)))
               #levels <- c(0.10, 0.99) * max(d$fhat)
                levels <- c(0.10) * max(d$fhat)
                l      <- contourLines(d$x1, d$x2, d$fhat, levels=levels)
                for(j in 1:length(l)) {
                    lines(l[[j]]$x, l[[j]]$y, lwd=lwd, col=col[i])
                }
                max_d    <- which.max(d$fhat)
                max_ind  <- arrayInd(max_d, .dim=dim(d$fhat))
                shadecol <- plotGetColors(length(chains), 128)
                points(d$x1[max_ind[1]], d$x2[max_ind[2]], col=shadecol[i], cex=3, pch=19)
            }, {
                stop("unknown method in pairPlot()")
            })
    }

    box()

    if (is.null(xlab) || is.null(ylab)) {
        var_unit <- paste(names(units), " (", units, ") ", sep="")
        names(var_unit) <- names(units)
        if (is.null(xlab)) {
            xlab <- var_unit[topColumn]
        }
        if (is.null(ylab)) {
            ylab <- var_unit[sideColumn]
        }
    }

    # the line is an offset from the axis, so it is not the same as mar
    title(xlab=xlab, line=xline)
    title(ylab=ylab, line=yline)

    # these are the same as the above
   #mtext(xlab=xlab, side=1, line=2)
   #mtext(ylab=ylab, side=2, line=2)


    # right PDF
    #

    # bottom, left, top, right
    #  left=0.25 gives space between the PDF and the main plot
    # right=0.25 gives space for the lower x-axis tick mark (need 1.00 for label)
    par(mar=c(bottom, 0.25, 0, 0.25))
    plot.new()
    plot.window(ylim=ylim, xlim=sidePdf$ylim, yaxs="i")
   #plotDensityAxis(1, labels=F)  # bottom axis
    axis(side=1, at=c(0, sidePdf$ylim), labels=F)
    title(xlab="PDF",          line=0)
   #title(xlab="Density",      line=0)
   #title(xlab="Prob density", line=2)
    pdfPlot(sidePdf, col=pdfcol, lty=lty, lwd=lwd, reverse=T)
}


pdfCdfPlots <- function(...,
    legends,
    col=plotGetColors(length(chains)), lty=rep("solid", length(chains)),
    column=as.character(2100), chains=list(...),
    lwd=2,
    xlab=gmslLab(column),
    ylab_pdf="Probability density",
    ylab_cdf=ifelse(survival, "Survival (1-CDF)", "Cumulative density"),
    burnin=T,
    legendloc="topright",
    survival=F, log=F, yline=2,
    smoothing=rep(1, length(chains)),
    layout=T,
    vlines=NULL,
    xlim=NULL,
    labels=c("a", "b")
    )
{
    if (layout) {
        # reserve lines to use outer=T for the lower axis and label;
        # allows using 0.5 in par(fig) and getting equally sized plots;
        # one line is 0.2 inches:  par("csi") or par("cin")[2]
        #
        #par(oma=c(4, 1, 1, 1))
        # bottom, left, top, right
        par(omi=c(1.00, 0.25, 0.25, 0.25))
    }


    # top of figure (PDFs)
    #

    if (layout) {
        par(fig=c(0, 1, 0.55, 1))
    }

    par(mar=c(0, 4, 0.25, 1))
    plot.new()

    pdfctx <- pdfCalc(chains=chains, column=column, burnin=burnin, smoothing=smoothing)

    if (is.null(xlim)) {
        xlim <- pdfctx$xlim
    }

    #    xlim=c(0, max(cictx$range)),
    plot.window(xlim=xlim, ylim=pdfctx$ylim, xaxs="i")

    axis(1, labels=F, tcl=-0.10)  # bottom
    plotDensityAxis()             # left

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    plotDensityAxis(4, labels=F, tcl=-0.25)  # right

    title(ylab=ylab_pdf, line=yline)
    box()
    if (!is.null(vlines)) {
        abline(v=vlines, lty=last(lty), lwd=1.5, col=last(col))
    }
    pdfPlot(pdfctx, col=col, lty=lty, lwd=lwd)
    legend(
        legendloc,
        legend=legends,
        col=col,
        lty=lty,
       #lwd=c(rep(lwd, length(col)), 1.5),
        lwd=1.5,
        bg="white"
        )
    labelPlot(labels[1])


    # bottom of figure (CDFs)
    #

    if (layout) {
        par(fig=c(0, 1, 0, 0.55), new=T)
    }

    par(mar=c(3, 4, 0, 1))
    cdfPlots(
        chains=chains,
        column=column,
        xlim=xlim,
        ylab=ylab_cdf,
        xlab=xlab,
        lwd=lwd, col=col, lty=lty,
        log=log, survival=survival
        )
    if (!is.null(vlines)) {
        abline(v=vlines, lty=last(lty), lwd=1.5, col=last(col))
    }
    labelPlot(labels[2], where=ifelse(log, "log", "topleft"))
}


plotLinearFit <- function(x, y, col="black", lty="solid", lwd=2)
{
    fit   <- lm(y ~ x)
    c     <- coef(fit)
    pts_x <- range(x)
   #pts_y <- predict(fit, newdata=data.frame(x=pts_x))
    pts_y <- pts_x * c[2] + c[1]
    pts_y <- constrain(pts_y, range(y))
    pts_x <- (pts_y - c[1]) / c[2]
    lines(pts_x, pts_y, col=col, lty=lty, lwd=lwd)

    return (fit)
}
