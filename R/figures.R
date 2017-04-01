# Copyright (C) 2016, 2017 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# figures.R

source("plot.R")
source("calib.R")  # daisRejSample()
loadLibrary("fields")  # image.plot()


figLhs <- function(assimctx=daisctx, outfiles=T, filetype="pdf", display=T)
{
    ais2100        <- assimctx$lhs$ychain
    parameters.lhs <- assimctx$lhs$chain
    obs.pfeffer    <- assimctx$expert_window
    lo.Tcrit       <- assimctx$lbound["Tcrit"]
    hi.Tcrit       <- assimctx$ubound["Tcrit"]
    lo.lambda      <- assimctx$lbound["lambda"]
    hi.lambda      <- assimctx$ubound["lambda"]

    igood    <- which(!is.na(ais2100))
    ais2100  <- ais2100[igood]
    lambda   <- parameters.lhs[igood, "lambda"]
    Tcrit    <- parameters.lhs[igood, "Tcrit"]
    ipfeffer <- which(ais2100 >= obs.pfeffer[1] & ais2100 <= obs.pfeffer[2])

    lims     <- c(min(ais2100[ which(!is.na(ais2100)) ]), max(ais2100[ which(!is.na(ais2100)) ]) )
    ncols    <- 50

    binwidth <- diff(lims) / ncols
    bounds   <- seq(lims[1], lims[2], length.out = ncols - 4)  # -5 for extra bins, +1 for bounding end point
    bounds   <- c(lims[1] - 2 * binwidth,  # no 3*bindwidth here bc interval closed on left,
                  lims[1] -     binwidth,  # open on right by .bincode()
                  bounds,
                  lims[2] +     binwidth,
                  lims[2] + 2 * binwidth,
                  lims[2] + 3 * binwidth)
    col.bin  <- .bincode(ais2100, bounds, right=F)
   #print(c(length(bounds), range(col.bin)))

    # rainbow palette is no good for the color blind
   #cols <- colorRampPalette(c("red", "orange", "yellow", "green", "blue"), space="Lab")(ncols)

    # for some color blind friendly schemes:  display.brewer.all(type="seq", colorblindFriendly=T)
   #cols <- colorRampPalette(c("yellow", "orange", "red"), space="Lab")(ncols)
   #cols <- colorRampPalette(c("blue", "green", "yellow"), space="Lab")(ncols)

    midpt <- (mean(obs.pfeffer) - bounds[1]) / (last(bounds) - bounds[1])
   #cols  <- designer.colors(n=ncols, col=c("blue", "white", "red"),       x=c(0, midpt, 1))
   #cols  <- designer.colors(n=ncols, col=c("purple3", "white", "red"),    x=c(0, midpt, 1))
    cols  <- designer.colors(n=ncols, col=c("darkviolet", "white", "red"), x=c(0, midpt, 1))


    newDev("fig_lhs", outfile=outfiles, width=3.5, height=9.7 / 2, filetype=filetype, mar=rep(0, 4))
    plotLayout(matrix(1:4, nrow=2, byrow=T), widths=c(8.25, 1.75))
   #mar <- c(3.5, 4, 1.5, 0)
    mar <- c(3.5, 4, 1.0, 0)


    par(mar=mar)
   #plot(Tcrit, lambda, pch=16,                       col=cols[col.bin], cex=0.75, xlim=c(lo.Tcrit, hi.Tcrit), ylim=c(lo.lambda, hi.lambda), ann=F, xaxs="i", yaxs="i")
    plot(Tcrit, lambda, pch=21, col="black", lwd=0.25, bg=cols[col.bin], cex=0.75, xlim=c(lo.Tcrit, hi.Tcrit), ylim=c(lo.lambda, hi.lambda), ann=F)

    mtext(daisTcritLab(),  side=1, line=2.25)
    mtext(daisLambdaLab(), side=2, line=2.0)
   #mtext(daisSlrLab(),    side=3, line=0.7, adj=1.4)
    labelPlot("a")

    plotColorBar(limits=c(min(bounds), max(bounds)), col=cols, mar=c(mar[1] + 2, 0.5, mar[3] + 2, 4))
    mtext(daisSlrLab(),    side=4, line=2.5)


    par(mar=mar)
   #plot(Tcrit[ipfeffer], lambda[ipfeffer], pch=16,                       col=cols[col.bin[ipfeffer]], cex=0.75, xlim=c(lo.Tcrit,hi.Tcrit), ylim=c(lo.lambda,hi.lambda), ann=F, xaxs="i", yaxs="i")
    plot(Tcrit[ipfeffer], lambda[ipfeffer], pch=21, col="black", lwd=0.25, bg=cols[col.bin[ipfeffer]], cex=0.75, xlim=c(lo.Tcrit,hi.Tcrit), ylim=c(lo.lambda,hi.lambda), ann=F)

    mtext(daisTcritLab(),  side=1, line=2.25)
    mtext(daisLambdaLab(), side=2, line=2.0)
   #mtext(daisSlrLab(),    side=3, line=0.7, adj=1.4)
    labelPlot("b")

    plotColorBar(limits=c(min(bounds), max(bounds)), col=cols, mar=c(mar[1] + 2, 0.5, mar[3] + 2, 4))
    mtext(daisSlrLab(),    side=4, line=2.5)


    if (outfiles) { finDev(display=display) }
}


legendfn <- function(legends, title, lwd, col, shadecol, ccol, ...)
{
   #lnames <- rep(expression('' <= -17.5*degree*C, '' > -17.5*degree*C), 2)  # , '' <= 'fit', '' > 'fit')
    lnames <- expression('' <= '-'*17.5*degree*C, '' > '-'*17.5*degree*C)
    lwd    <- 2
    legend(
        "center",
        legend=lnames,
        title=title,
       #title.adj=0.1,
        col=col,
        lwd=c(NA, NA), # , lwd, lwd),
        pch=c(15, 15), #, NA,  NA),
       #x.intersp=0.25,
        x.intersp=0.125,
        y.intersp=0.80,
        text.width=0.925*strwidth(lnames)[1],
       #pt.cex=1.67
        seg.len=0.50,
        pt.cex=1.25
        )
}


plotfn <- function(samples, i, topColumn, sideColumn, col, shadecol, ccol, assimctx, ...)
{
    cold <- which(samples[, "Tcrit"] <= -17.5)
    warm <- which(samples[, "Tcrit"] >  -17.5)

    coldCol <- col[3]
    warmCol <- col[4]
    if (!is.null(assimctx$wide_prior) && assimctx$wide_prior) {
        lwd <- 1.50
        cex <- 0.25
    } else {
        lwd <- 2.00
        cex <- 0.50
    }
    inches  <- 0.05
   #f       <- 1/3

    x <- samples[, topColumn]
    y <- samples[, sideColumn]

    points(      x[warm], y[warm], pch=20, col=shadecol[2], cex=cex)
    points(      x[cold], y[cold], pch=20, col=shadecol[1], cex=cex)
   #lines(lowess(x[warm], y[warm], f=f), lwd=lwd, col=warmCol)
   #lines(lowess(x[cold], y[cold], f=f), lwd=lwd, col=coldCol)

    # top equations
    if (length(warm)) {
        fit <- plotLmFit(x[warm], y[warm], lwd=lwd, col=warmCol)
        plotLmText(fit, col=warmCol, ...)
    }

    # bottom equations
    if (length(cold)) {
        fit <- plotLmFit(x[cold], y[cold], lwd=lwd, col=coldCol)
        plotLmText(fit, col=coldCol, where="bottomleft", ...)
    }
}


figDiagFast <- function(assimctx=daisctx, prctx=prdaisctx, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_diag_fast", outfile=outfiles, height=(2/3)*9.7, filetype=filetype, mar=rep(0, 4))

    nfig <- 3

   #plotLayout(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T), widths=c(10, 3), heights=rep(c(3, 10), nfig))
   #plotLayout(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T), widths=c(10, 2), heights=rep(c(3, 10), nfig))
   plotLayout(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T), widths=c(10, 1.9), heights=rep(c(3, 10), nfig))

    # limits for SLE
   #xlim   <- c(0.075, 0.675)
    xlim   <- c(0.035, 0.700)
    title  <- daisTcritLab()
   #title  <- expression(T[crit])
    lnames <- expression('' <= '-'*17.5*degree*C, '' > '-'*17.5*degree*C)
    lfn    <- NULL
    lfn    <- legendfn
    slrCol <- as.character(2100)

    # see file junk and function figLambda for the rejection sampling version
    if (is.null(assimctx$diagChain)) {
        daisRunPredict(subsample=F, assimctx=assimctx, prctx=prctx)
        burned_ind <- burnedInd(assimctx$chain)
        assimctx$diagChain <- cbind(assimctx$chain[burned_ind, ], prctx$prchain)
    }

    if (!is.null(assimctx$wide_prior) && assimctx$wide_prior) {
        points <- 3e3
    } else {
        points <- 6e3
    }
    method <- "plotfn"
    col    <- c(plotGetColors(3)[1:2], plotGetColors(4, pal="Dark2")[3:4])
    pdfcol <- "black"
    lwd    <- 1
    left   <- 4
    bottom <- 4.25
    xline  <- 2.25
    smooth <- rep(1.5, 3)

    limitTcrit   <- c(assimctx$lbound["Tcrit"]  - 0.5,   assimctx$ubound["Tcrit"]  + 0.5)
    limitLambda  <- c(assimctx$lbound["lambda"] - 0.001, assimctx$ubound["lambda"] + 0.001)
    limitLambda2 <- c(assimctx$lbound["lambda"] - 0.001, assimctx$ubound["lambda"] + 0.004)

    pairPlot(chains=list(assimctx$diagChain), layout=F, legends=lnames, points=points, method=method,
        pdfcol=pdfcol, lwd=lwd, col=col, smoothing=smooth, label="a", title=title, mar=c(bottom, left),
        ylim=xlim, xlim=limitTcrit,  xname="italic(T[crit])", yname="italic(SLR)", ylab=daisSlrLab(), xlab=daisTcritLab(),
        assimctx=assimctx, sideColumn=slrCol,   topColumn="Tcrit",  xline=xline, legendfn=lfn)

    pairPlot(chains=list(assimctx$diagChain), layout=F, legends=lnames, points=points, method=method,
        pdfcol=pdfcol, lwd=lwd, col=col, smoothing=smooth, label="b", title=title, mar=c(bottom, left),
        ylim=xlim, xlim=limitLambda, xname="lambda",          yname="italic(SLR)", ylab=daisSlrLab(), xlab=daisLambdaLab(),
        assimctx=assimctx, topColumn="lambda", sideColumn=slrCol,   xline=xline, legendfn=lfn)

    pairPlot(chains=list(assimctx$diagChain), layout=F, legends=lnames, points=points, method=method,
        pdfcol=pdfcol, lwd=lwd, col=col, smoothing=smooth, label="c", title=title, mar=c(3.5, left),
        xlim=limitTcrit, ylim=limitLambda2, xname="italic(T[crit])", yname="lambda", xlab=daisTcritLab(),
        ylab=daisLambdaLab(),
        assimctx=assimctx, topColumn="Tcrit",  sideColumn="lambda", xline=xline, legendfn=lfn)

    if (outfiles) { finDev(display=display) }
}
