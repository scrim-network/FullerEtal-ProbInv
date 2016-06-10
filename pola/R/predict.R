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
# predict.R
# written by Robert W. Fuller on 090725

source('plot.R')


prqplot <- function(prctx=prallgrgisctx, ...)
{
    prqplotsub(
        prctx$prchain,
        legends=c("Mean", "95% credible interval"),
        col=c("blue", "red"),
        lty=c("solid", "dashed"),
        prctx=prctx,
        ...
        )
}


prqplotsub <- function(
    ..., legends, col, lty,
    prctx=NULL,
    args=list(...),
    N=ifelse(plotci, 2 * length(args), length(args)),
    lwd   =rep( 2, N),
    pch   =rep(NA, N),
    pt.bg =rep(NA, N),
    pt.lwd=rep(NA, N),
    pt.cex=rep( 1, N),
    xmin=as.numeric(colnames(args[[1]])[1]),
    xmax=as.numeric(last(colnames(args[[1]]))),
    obs=T, shade=T, invert=F, obserr=F,
    caption="Figure 1. Posterior predictive for global mean sea level",
    outfiles=F,
    xvals= as.numeric(colnames(args[[1]]) [
           as.numeric(colnames(args[[1]])) >= xmin
         & as.numeric(colnames(args[[1]])) <= xmax ]),
    xcaption="Year",
    ycaption=gmslLab(),
    plotci=T,
    ymin=Inf,
    ymax=-Inf,
    height=6,
    width=6,
    newdev=T
    )
{
    if (newdev) {
        newDev("ppred", outfiles, height=height, width=width)
    }

    missing_ymin <- missing(ymin)
    missing_ymax <- missing(ymax)

    prmean <- list()
    ci     <- list()

    for (i in 1:length(args)) {

        chain <- args[[ i ]]
        if (!missing(xmin) || !missing(xmax)) {
            chain <- chain[, as.character(xvals)]
        }

        prmean[[i]] <- cbind(xvals, colMeans(chain))
        #prmean[[i]] <- cbind(xvals, colMode(chain))
        ci[[i]]     <- cbind(xvals, colQuantile(chain))
        if (missing_ymin) {
            ymin <- c(min(ymin, ci[[i]][, "0.025"]))
        }
        if (missing_ymax) {
            ymax <- c(max(ymax, ci[[i]][, "0.975"]))
        }
    }

    xlim <- c(xvals[1], last(xvals))
    ylim <- c(ymin, ymax)
    emptyPlot(xlim, ylim, xcaption, ycaption)

    if (obs) {
        legends <- append(legends, "Observations")
        lty     <- append(lty,     NA)
        col     <- append(col,     "purple")
        lwd     <- append(lwd,     NA)
        pch     <- append(pch,      3)
        pt.bg   <- append(pt.bg,   NA)
        pt.lwd  <- append(pt.lwd,  NA)
        pt.cex  <- append(pt.cex, 1.0)
    }

    # assimilation range as gray shaded area;
    # draw this first so that it is in the back in terms of Z-order
    #
    if (shade) {
        shadecol <- rgb(0.8, 0.8, 0.8)

        if (invert) {
            if (xmin < prctx$assimctx$times[1]) {
                xshade(xmin, prctx$assimctx$times[1], col=shadecol)
            }
            if (xmax > last(prctx$assimctx$times)) {
                xshade(last(prctx$assimctx$times), xmax, col=shadecol)
            }
            legends <- append(legends, "Prediction range")
        } else {
            xshade(prctx$assimctx$times[1], last(prctx$assimctx$times), col=shadecol)
            legends <- append(legends, "Assimilation range")
        }

        # a sneaky way of doing boxes is to use pch=22 with a border color (col="black"),
        # a fill color (pt.bg), a larger character size (pt.cex=2.0), and a smaller border (pt.lwd=0.5)
        lty     <- append(lty,     NA)
        col     <- append(col,     "black")
        lwd     <- append(lwd,     NA)
        pch     <- append(pch,     22)
        pt.bg   <- append(pt.bg,   shadecol)
        pt.lwd  <- append(pt.lwd, 0.5)
        pt.cex  <- append(pt.cex, 2.0)
    }

    # draw borders AFTER xshade() so that shading doesn't overwrite borders
    box()

    # observations as plus signs
    if (obs) {
        times <- prctx$assimctx$obstime
        if (is.null(times)) {
            times <- prctx$assimctx$times
        }

        if (obserr) {
            ts <- cbind(times, prctx$assimctx$obsonly, prctx$assimctx$error)
            colnames(ts) <- c("time", "y", "error")
            tsErrorBars(ts, col="purple", shade=F)
        }

        cex <- ifelse(any(par("mfrow") > 1, length(xvals) > 200), 0.5, 1)

        # a temporary hack for showing all observations when using prqplot()
        # for GCM assimilation;  not needed anymore with ipccPlotFits()
        #
        #points(prctx$assimctx$obs, pch=3, cex=cex, col="purple")
        points(times, prctx$assimctx$obsonly, pch=3, cex=cex, col="purple")
    }

    j <- 1
    for (i in 1:length(args)) {

        # mean of runs as a line
        lines(prmean[[i]], col=col[j], lty=lty[j], lwd=lwd[j])
        j <- j + 1

        if (plotci) {
            # credible interval as dashed lines
            tsColLines(ci[[i]], col=col[j], lty=lty[j], lwd=lwd[j])
            j <- j + 1
        }
    }

    legend(
        "topleft",
        legend=legends,
        lty=lty,
        lwd=lwd,
        col=col,

        # this forces the Z-order
        bg="white",

        pch=pch,
        pt.bg=pt.bg,
        pt.lwd=pt.lwd,
        pt.cex=pt.cex
        )

    if (newdev) {
        mtext(caption, outer=TRUE, side=1)
    }

    if (newdev && outfiles) { finDev() }
}
