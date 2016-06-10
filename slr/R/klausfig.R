# Copyright 2011 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# klausfig.R

source("newfig.R")


figPredict <- function(
    prctx=prallgrgisctx,
    #caption="Figure 1. Posterior predictive PDFs for global mean sea level",
    caption=NULL,
    outfiles=F
    )
{
    newDev("predict", width=4, height=6, outfiles)

    # reserve lines to use outer=T for the lower axis and label;
    # allows using 0.5 in par(fig) and getting equally sized plots;
    # one line is 0.2 inches:  par("csi") or par("cin")[2]
    #
    #par(oma=c(4, 1, 1, 1))
    par(omi=c(0.75, 0.25, 0.25, 0.25))


    # parameters common to CDF/PDF
    #

    taus <- c(pr1$assimctx$ep["tau"], pr2$assimctx$ep["tau"], pr3$assimctx$ep["tau"])
    legends <- expression()
    for (tau in taus) {
        legends <- append(legends, bquote(paste(tau[1]==.(tau), " a")))
    }

    col     <- c("black", "blue", "red")
    lty     <- c("solid", "dotdash", "longdash")
    lwd     <- 2

    column <- as.character(2100)

    pdflim <- c(0, 1.1)


    # top of figure (PDFs)
    #

    par(fig=c(0, 1, 0.5, 1))
    par(mar=c(0, 3, 0, 1))
    plot.new()

    pdfctx <- pdfCalc(pr1$prchain, pr2$prchain, pr3$prchain, column=column, burnin=F)

    plot.window(xlim=pdflim, ylim=pdfctx$ylim, xaxs="i")

    # bottom
    axis(1, labels=F, tcl=-0.10) # bottom

    # left
    ticks <- axTicks(2)
    ticks <- c(0, last(ticks))
    axis(2, at=ticks)

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    axis(4, at=ticks, labels=F, tcl=-0.25) # right

    title(ylab="Probability density", line=2)
    box()

    pdfPlot(pdfctx, col=col, lty=lty, lwd=lwd)

    legend(
        "topright",
        legend=legends,
        lty=lty,
        col=col,
        lwd=rep(lwd, length(lty))
        )

    labelPlot("a")


    # bottom of figure (CDFs)
    #

    par(fig=c(0, 1, 0, 0.5), new=T)
    par(mar=c(0, 3, 0, 1))

    cdfPlots(
        pr1$prchain, pr2$prchain, pr3$prchain,
        column=column,
        xlim=pdflim,
        xlab=NULL,
        col=col, lty=lty, lwd=lwd
        )
    labelPlot("b")

    title(xlab=gmslLab(year=2100), line=2, outer=T)


    # figure title
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# figure with forcing and probabilistic hindcasts and forecasts, per Klaus
figProbForce <- function(prctx=prallgrgisctx, year=2100, outfiles=F)
{
    newDev("probforce", height=8.5, width=11, outfiles, horiz=T)

    # reserve lines to use outer=T for the lower axis and label;
    # one line is 0.2 inches:  par("csi") or par("cin")[2]
    #
    #par(oma=c(4, 1, 1, 1))

    # bottom, left, top, right
    par(omi=c(0.75, 0.25, 0.25, 0.25))

    # stacked plots
    ngraphs=4
    layout(matrix(1:ngraphs, nrow=ngraphs, ncol=1))

    # layout() sets cex=0.66!
    # TODO:  what else does it set??
    #par(cex=1)
    #par(ps=7)

    #layout.show(len)


    times <- attr(prctx$prchain, "xvals")[1]:year
    #times <- prctx$assimctx$times[1]:year

    legends <- c("Mean", "95% credible interval")
    lwd <- 1


    # 1st figure:  lines on bottom, left, top, right
    # only need to set this once, apparently
    par(mar=c(0, 3, 0, 1))


    # temperature forcing
    #

    ar4bias=T

    gmst <- loadUrban(ar4bias=ar4bias)

    # get this into chain form for calculating interval
    xvals <- gmst[, 1]  # save years
    gmst <- gmst[, 2:ncol(gmst)]  # remove years
    gmst <- t(gmst)
    attr(gmst, "xvals") <- xvals  # restore years

    cictx <- ciCalc(gmst, xvals=times)

    plot.new()
    plot.window(xlim=range(times), ylim=cictx$range, yaxs="i", xaxs="i")
    axis(2) # left
    axis(1, labels=F, tcl=-0.10) # bottom
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right
    title(ylab="Temperature (C)", line=2)
    box()

    ciPlot(cictx, lwd=lwd)
    legend(
        "topleft",
        legend=c(legends, "Brohan et al. 2005"),
        lty=c("solid", "dashed", NA),
        lwd=c(rep(lwd, 2), NA),
        col=c(rep("red", 2), "purple"),
        pch=c(rep(NA, 2), 3)
        )

    #brohan <- loadHadcrut(ar4bias=ar4bias, dataset="global")$annual
    brohan <- loadHadcrut(ar4bias=ar4bias)$annual
    points(brohan, pch=3, cex=0.5, col="purple")

    labelPlot("a", where="left")


    # other:  mean, and confidence interval
    #

    cictx <- ciCalc(prctx$otherchain, xvals=times)

    plot.new()
    plot.window(xlim=range(times), ylim=cictx$range, yaxs="i", xaxs="i")
    axis(2) # left
    axis(1, labels=F, tcl=-0.10) # bottom
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right
    title(ylab="Non-Greenland sea-level rise (m)", line=2)
    box()

    ciPlot(cictx, lwd=lwd)
    legend(
        "topleft",
        legend=c(legends),
        lty=c("solid", "dashed"),
        lwd=c(rep(lwd, 2)),
        col=c(rep("red", 2)),
        pch=c(rep(NA, 2))
        )

    labelPlot("b", where="left")


    # Rignot fit:  mean, confidence interval, and observations
    #

    cictx <- ciCalc(prctx$ds_gis, xvals=times)

    plot.new()
    plot.window(xlim=range(times), ylim=cictx$range, yaxs="i", xaxs="i")
    axis(2) # left
    axis(1, labels=F, tcl=-0.10) # bottom
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right
    title(ylab=slrGreenlandLab(), line=2)
    box()

    legend(
        "topleft",
        legend=c(legends, "Rignot et al. 2008"),
        lty=c("solid", "dashed", NA),
        lwd=c(rep(lwd, 2), NA),
        col=c(rep("red", 2), "purple"),
        pch=c(rep(NA, 2), 3)
        )

    #ciPlot(cictx, col=col, lty=lty, plotci=T)
    ciPlot(cictx, lwd=lwd)

    # could use xbeam=T, ibeam=T below, but would need to match to size of cex?
    #massBal <- loadRignot(pure=T)
    #points(massBal[, "time"], massBal[, "SLE"], pch=3, cex=0.5, col="purple")

    ts <- cbind(prctx$assimctx$gis_times[prctx$assimctx$gis_ind], prctx$assimctx$gis_obs, prctx$assimctx$gis_err$error)
    colnames(ts) <- c("time", "SLE", "error")
    tsErrorBars(ts, shade=F, ibeam=T, xbeam=T, col="purple")

    labelPlot("c", where="left")


    # tide gage:  mean, confidence interval, and observations
    #

    cictx <- ciCalc(prctx$prchain, xvals=times)

    plot.new()
    plot.window(xlim=range(times), ylim=cictx$range, yaxs="i", xaxs="i")
    axis(2) # left
    #axis(1, labels=F, tcl=-0.10) # bottom
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right
    title(ylab=gmslLab(), line=2)
    box()

    ts <- cbind(prctx$assimctx$times, prctx$assimctx$obsonly, prctx$assimctx$error)
    colnames(ts) <- c("time", "y", "error")
    tsErrorBars(ts, col="purple", shade=F)
    #ciPlot(cictx, col=col, lty=lty)
    ciPlot(cictx, lwd=lwd)
    legend(
        "topleft",
        legend=c(legends, "Jevrejeva et al. 2006"),
        lty=c("solid", "dashed", NA),
        lwd=c(rep(lwd, 2), NA),
        col=c(rep("red", 2), "purple"),
        pch=c(rep(NA, 2), 3)
        )

    labelPlot("d", where="left")


    # x-axis and x-axis title
    axis(1, outer=T)
    title(xlab="Year", outer=T, line=2)

    if (outfiles) { finDev() }
}


figHindcast <- function(
    prctx=prallgrgisctx,
    #caption="Figure 1. Hindcasts for global mean sea level",
    caption=NULL,
    outfiles=F
    )
{
    newDev("hindcast", width=8.5, height=5, outfiles)

    par(mfcol=c(1, 2))

    # bott, left, top, right
    par(mar=c(4, 3, 0, 2))

    taus <- c(pr1$assimctx$ep["tau"], pr2$assimctx$ep["tau"], pr3$assimctx$ep["tau"])
    legends <- expression()
    for (tau in taus) {
        legends <- append(legends, bquote(paste(tau[1]==.(tau), " a")))
    }

    col     <- c("black", "blue", "red")
    lty     <- c("solid", "dotdash", "longdash")


    # hindcast
    #

    times <- prctx$assimctx$times
    cictx <- ciCalc(pr1$prchain, pr2$prchain, pr3$prchain, xvals=times)
    emptyPlot(xlim=range(times), ylim=cictx$range, xlab="Year", ylab=gmslLab())
    ts <- cbind(times, prctx$assimctx$obsonly, prctx$assimctx$error)
    colnames(ts) <- c("time", "y", "error")
    tsErrorBars(ts, col="purple", shade=F)
    ciPlot(cictx, col=col, lty=lty, plotci=F)
    legend(
        "topleft",
        legend=c(legends, "Jevrejeva et al. 2006"),
        lty=c(lty, NA),
        lwd=c(rep(2, 3), NA),
        col=c(col, "purple"),
        pch=c(rep(NA, 3), 3)
        )
    labelPlot("a")


    # Rignot fit
    #

    massBal <- loadRignot(pure=T)
    times <- massBal[, "time"]
    times <- (times[1] - 1):(last(times) + 1)
    cictx <- ciCalc(pr1$ds_gis, pr2$ds_gis, pr3$ds_gis, xvals=times)
    emptyPlot(xlim=range(times), ylim=cictx$range, xlab="Year", ylab=slrGreenlandLab())
    legend(
        "topleft",
        legend=c(legends, "Rignot et al. 2008"),
        lty=c(lty, NA),
        lwd=c(rep(2, 3), NA),
        col=c(col, "purple"),
        pch=c(rep(NA, 3), 3)
        )

    # could use xbeam=T, ibeam=T below, but would need to match to size of cex?
    points(massBal[, "time"], massBal[, "SLE"], pch=3, cex=0.5, col="purple")

    ts <- cbind(prctx$assimctx$gis_times[prctx$assimctx$gis_ind], prctx$assimctx$gis_obs, prctx$assimctx$gis_err$error)
    colnames(ts) <- c("time", "SLE", "error")
    tsErrorBars(ts, shade=F, col="purple")
    ciPlot(cictx, col=col, lty=lty, plotci=F)

    labelPlot("b")


    # figure title
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}
