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
# ipcc.R

source('data.R')
source('plot.R')
source('roblib.R') # notDir()
source('grinassim.R') # ipccPlotFit()


ipccPlotControl <- function(xaxis=c(1850, 2000), outfiles=F)
{
    ipccPlotSeaLevel(path="../data/ipcc/control", caption="Figure 1. Comparison of IPCC pre-industrial control runs", domingues=F, xaxis=xaxis, outfiles=outfiles)
}


ipccPlotPredict <- function(xaxis=c(2001, 2100), outfiles=F)
{
    ipccPlotSeaLevel(path="../data/ipcc/predict", caption="Figure 1. Comparison of IPCC predictions (SRESA1B)", domingues=F, xaxis=xaxis, ar4bias=F, legendloc="topleft", outfiles=outfiles)
}


ipccPlotDriftCorrect <- function(xaxis=c(1850, 2100), outfiles=F)
{
    ipccPlotSeaLevel(path="../data/ipcc/control", caption="Figure 1. Comparison of drift-corrected predictions", domingues=F, xaxis=xaxis, legendloc="topleft", correct=T, outfiles=outfiles)
}


ipccPlotSeaLevel <- function(path="../data/ipcc/hindcast", color=c("black", "red", "blue", "green", "purple", "yellow", "gray", "orange"), caption="Figure 1. Comparison of IPCC predictions", domingues=T, newdev=T, xaxis=NULL, yaxis=NULL, ar4bias=T, legendloc="bottomright", correct=F, outfiles=F)
{
    if (newdev) {
        newDev("ipcc_pred", outfiles)
    }

    preds <- list()
    xrange <- numeric()
    yrange <- numeric()

    filenames <- notDir(list.files(path, full.names=T))
    for (i in 1:length(filenames)) {
        filename <- filenames[i]

        if (correct) {

            filename <- basename(filenames[i])
            ts <- loadIpccSeaLevelDriftCorrect(filename, truncate=T)$obs

        } else {
            ts <- loadIpccSeaLevel(filename, ar4bias=ar4bias)$annual
            if (!ar4bias) {
                ts[, "sealvl"] <- ts[, "sealvl"] - ts[1, "sealvl"]
            }
        }

        preds[[i]] <- ts

        xrange <- range(xrange, ts[, "time"])
        yrange <- range(yrange, ts[, "sealvl"])
    }

    if (domingues) {
        ts <- loadDomingues(ar4bias=ar4bias)
        if (!ar4bias) {
            ts[, "sealvl"] <- ts[, "sealvl"] - ts[1, "sealvl"]
        }

        preds[[i + 1]] <- ts
        filenames[[i + 1]] <- "Domingues et al."
        xrange <- range(xrange, ts[, "time"])
        yrange <- range(yrange, ts[, "sealvl"])
    }

    if (is.null(xaxis)) {
        xaxis <- xrange
    }
    if (is.null(yaxis)) {
        yaxis <- yrange
    }

    emptyPlot(xlim=xaxis, ylim=yaxis,
        xlab="Year", ylab="Thermosteric global mean sea-level anomaly (m)")

    for (i in 1:length(preds)) {
        ts <- preds[[i]]
        lines(ts[, "time"], ts[, "sealvl"], col=color[i])
    }

    modelnames <- sub("\\..*", "", basename(filenames))
    legend(
        legendloc,
        legend=modelnames,
        lty="solid",
        col=color,
        lwd=rep(1, length(preds))
        )

    mtext(caption, outer=TRUE, side=1)

    if (newdev && outfiles) { finDev() }
}


ipccPlotFits <- function(assimlist=list(grinassimctx), xlim=c(1850, 2100), endYear=2100, outfiles=F)
{
    newDev("gcm_fit2", width=4, height=8, outfile=outfiles)

    len       <- length(assimlist)
    assimYear <- last(assimlist[[1]]$times)

    # reserve lines to use outer=T for the lower axis and label;
    # one line is 0.2 inches:  par("csi") or par("cin")[2]
    #
    #par(oma=c(4, 1, 1, 1))
    par(omi=c(0.75, 0.75, 0.25, 0.25))

    # stacked plots
    layout(matrix(1:len, nrow=len, ncol=1))

    # layout() sets cex=0.66!
    # TODO:  what else does it set??
    #par(cex=1)
    #par(ps=7)

    #layout.show(len)

    pch <- 1
    label <- utf8ToInt("a") - 1

    ylim <- numeric()
    obslist <- list()
    slrlist <- list()
    for (i in 1:len) {
        assimctx <- assimlist[[i]]

        # ipcc observations
        idx <- which(assimctx$obs[, "time"] <= endYear)
        obs <- assimctx$obs[idx, ]

        # best-fit model
        times <- obs[1, "time"]:endYear
        slr <- grinsted(times, list(mp=assimctx$init_mp, frc=list(assimctx$frc), spl=2, ep=assimctx$ep, sw=assimctx$sw))

        ylim <- range(ylim, slr[, "sealvl"], obs[, "sealvl"])

        obslist[[i]] <- obs
        slrlist[[i]] <- slr
    }

    for (i in 1:len) {
        assimctx <- assimlist[[i]]
        obs      <-   obslist[[i]]
        slr      <-   slrlist[[i]]

        #ylim <- range(slr[, "sealvl"], obs[, "sealvl"])

        par(mar=c(0, 1, 0, 1))

        plot.new()
        plot.window(xlim=xlim, ylim=ylim, xaxs="i")

        axis(2) # left
        axis(1, labels=F, tcl=-0.10) # bottom
        axis(3, labels=F, tcl=-0.10) # top
        axis(4, labels=F, tcl=-0.25) # right

        box()

        points(obs, col="blue", pch=pch)
        lines(slr, col="black", lwd=2)
        legend <- "Simple model (max. likelihood)"
        legend <- append(legend, paste("AOGCM:", sub("\\..*", "", basename(assimlist[[i]]$ipccFile))))

        legend(
            "topleft",
            #legend=c("Simple model (maximum likelihood)", "AOGCM: MIUB ECHO-G"),

            legend=legend,
            lty=c("solid", NA),
            lwd=c(2, NA),
            col=c("black", "blue"),
            cex=0.75,
            pch=c(NA, pch)
            )

        if (ceiling(len/2) == i) {
            arrowx <- (ylim[2] - ylim[1]) * 0.6 + ylim[1]

            # Klaus didn't like italics
            font <- 1
            #xoff <- 5
            xoff <- 0.5 * par("cxy")[1]
            arrowlen <- 1.2 * strwidth("Hindcast", font=font) # + par("cxy")[1]

            arrows(assimYear + xoff, arrowx, assimYear + arrowlen + xoff, lwd=2, length=par("ps")/72/2)
            arrows(assimYear - xoff, arrowx, assimYear - arrowlen - xoff, lwd=2, length=par("ps")/72/2)

            text(assimYear, arrowx + strheight("Forecast"), labels="Forecast", pos=4, font=font)
            text(assimYear, arrowx + strheight("Forecast"), labels="Hindcast", pos=2, font=font)
        }

        abline(v=assimYear, lwd=1)

        # left:  side=2
        # horizontal:  las=1
        # above axis label:  line=3
        # bold:  font=2
        # locate at y-axis:  at=par("usr")[4]
        #
        # font=2
        # at=par("usr")[4], 
        #mtext(intToUtf8(label + i), side=2, las=1, line=3, cex=1.25, font=2)
        labelPlot(intToUtf8(label + i), where="left")

        #labelPlot(labels[i])
    }

    axis(1, outer=T)
    title(xlab="Year", outer=T, line=2)
    title(ylab=gmslLab(), outer=T, line=1)

    if (outfiles) { finDev() }
}
