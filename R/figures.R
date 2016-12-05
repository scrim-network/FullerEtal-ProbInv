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
# figures.R

source("plot.R")
source("calib.R")  # daisRejSample()


figLambda <- function(assimctx=daisctx, prctx=prdaisctx, outline=F, lambda=T, outfiles=T, filetype="pdf")
{
    newDev(ifelse(lambda, "fig3_lambda", "fig3_Tcrit"), outfile=outfiles, width=4, height=8, filetype=filetype)

   #layout(cbind(matrix(1:4, nrow=2, byrow=T), matrix(5:8, nrow=2, byrow=T)), widths = rep(c(10, 3), 2), heights = c(3, 10))
    layout(cbind(matrix(1:8, nrow=4, byrow=T)), widths = c(10, 3), heights = rep(c(3, 10), 2))

    # limits for SLE
    xlim <- c(0, 0.8)

    if (lambda) {
        # limits for lambda
       #ylim <- c(-.016, .016)
        ylim <- c(.004, .016)
        sideColumn <- "lambda"
    } else {
       #ylim <- c(-41, -9)
        ylim <- c(-21, -9)
        sideColumn <- "Tcrit"
    }

    cnames <- "Uniform"
    units <- assimctx$units
    units <- append(units, "m")
    names(units)[length(units)] <- "2100"

    if (is.null(assimctx$noRejChain)) {
        daisRejSample(assimctx=assimctx, prctx=prctx)
    }

    pre_chain  <- cbind(assimctx$noRejChain, prdaisctx$prNoRejChain)
    burned_ind <- burnedInd(assimctx$noRejChain)
    pre_chain  <- pre_chain[burned_ind, ]
    post_chain <- cbind(assimctx$chain, prdaisctx$prchain)

    points <- ifelse(outline, 1e5, 6e3)
    method <- ifelse(outline, "outline", "points")
    col    <- plotGetColors(3)

    pairPlot(pre_chain,  layout=F, units=units, xlim=xlim, ylim=ylim, method=method,
        topColumn="2100", sideColumn=sideColumn, legends=cnames, points=points, label="a", col=col, smoothing=rep(1.5, 3))

    pairPlot(post_chain, layout=F, units=units, xlim=xlim, ylim=ylim, method=method,
        topColumn="2100", sideColumn=sideColumn, legends=cnames, points=points, label="b", col=col, smoothing=rep(1.5, 3))

    caption <- paste("Figure n. Diagnosing Uniform Inversion; (a) Before rejection sampling, (b) After rejection sampling")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


plotfn <- function(samples, i, topColumn, sideColumn, col, shadecol, ccol)
{
    cold <- which(samples[, "Tcrit"] <= -17.5)
    warm <- which(samples[, "Tcrit"] >  -17.5)
    
    #x <- samples[, topColumn]
    #y <- samples[, sideColumn]
    points(samples[cold, topColumn], samples[cold, sideColumn], pch=20, col=shadecol[1], cex=0.5)
    points(samples[warm, topColumn], samples[warm, sideColumn], pch=20, col=shadecol[2], cex=0.5)
}


figTony <- function(assimctx=daisctx, prctx=prdaisctx, outline=F, outfiles=T, filetype="pdf")
{

}


figTony2 <- function(assimctx=daisctx, prctx=prdaisctx, outline=F, outfiles=T, filetype="pdf")
{
    newDev("fig3_tony", outfile=outfiles, width=5, height=10, filetype=filetype)

    nfig <- 3

    layout(cbind(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T)), widths = c(10, 3), heights = rep(c(3, 10), nfig))

    # limits for SLE
    xlim <- c(0.1, 0.65)

    cnames <- c("<= -17.5", ">  -17.5")
    units <- assimctx$units
    units <- append(units, "m")
    names(units)[length(units)] <- "2100"

    if (is.null(assimctx$diagChain)) {
        daisRunPredict(subsample=F, assimctx=assimctx, prctx=prctx)
        burned_ind <- burnedInd(assimctx$chain)
        assimctx$diagChain <- cbind(assimctx$chain[burned_ind, ], prctx$prchain)
    }

    points <- 6e3
    method <- "plotfn"
    col    <- plotGetColors(3)
    pdfcol <- "black"
    lwd    <- 1

    limitTcrit  <- c(assimctx$lbound["Tcrit"],  assimctx$ubound["Tcrit"])
    limitLambda <- c(assimctx$lbound["lambda"], assimctx$ubound["lambda"])

    limitTcrit  <- c(-20.5, -9.5)
    limitLambda <- c(.004, .016)

    pairPlot(assimctx$diagChain, layout=F, units=units, legends=cnames, points=points, method=method,
        title="Tcrit", col=col, smoothing=rep(1.5, 3), pdfcol=pdfcol, lwd=lwd,
        topColumn="2100",  sideColumn="Tcrit",  label="a", xlim=xlim,       ylim=limitTcrit)

    pairPlot(assimctx$diagChain, layout=F, units=units, legends=cnames, points=points, method=method,
        title="Tcrit", col=col, smoothing=rep(1.5, 3), pdfcol=pdfcol, lwd=lwd,
        topColumn="2100",  sideColumn="lambda", label="b", xlim=xlim,       ylim=limitLambda)

    pairPlot(assimctx$diagChain, layout=F, units=units, legends=cnames, points=points, method=method,
        title="Tcrit", col=col, smoothing=rep(1.5, 3), pdfcol=pdfcol, lwd=lwd,
        topColumn="Tcrit", sideColumn="lambda", label="c", xlim=limitTcrit, ylim=limitLambda)

   #caption <- paste("Figure n. Diagnosing Uniform Inversion; (a) Before rejection sampling, (b) After rejection sampling")
   #mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}