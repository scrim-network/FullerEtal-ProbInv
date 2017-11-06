# Copyright (C) 2017 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# deconto.R

if (file.exists(       'ep="b";n=5e6.RData')) {
    iter <- "5e6"
} else if (file.exists('ep="b";n=2e6.RData')) {
    iter <- "2e6"
} else if (file.exists('ep="b";n=5e5.RData')) {
    iter <- "5e5"
}
#iter <- "5e6"
#iter <- "2e6"
#iter <- "5e5"


source('plot.R')
source('calib.R')  # daisSlrLab()


if (!exists("pr1")) {

    loadChains(paste('0SAVE_ME/ep="n";n=', iter, ".RData", sep=""))
    loadChains(paste('0SAVE_ME/ip="n";n=', iter, ".RData", sep=""), newnames=c("ias", "ipr"))
    loadChains(paste(         'ep="n";n=', iter, ".RData", sep=""), newnames=c( "das",  "dpr"))
    loadChains(paste(         'ip="n";n=', iter, ".RData", sep=""), newnames=c("dias", "dipr"))

      as1_scaled <- daisScaleChain(  as1$chain, assimctx=as1)
     ias1_scaled <- daisScaleChain( ias1$chain, assimctx=as1)

     das1_scaled <- daisScaleChain( das1$chain, assimctx=das1)
    dias1_scaled <- daisScaleChain(dias1$chain, assimctx=das1)

    Tcrit  <- -as1$ Tcrit_prior$rand(5e6)
    lambda <-  as1$lambda_prior$rand(5e6)
    priorChain <- cbind(Tcrit, lambda)
    colnames(priorChain) <- c("Tcrit", "lambda")
    prior_scaled <- daisScaleChain(priorChain, assimctx=as1)
    rm(Tcrit, lambda, priorChain)  # keep it clean
}


# "Deep uncertainty"
figDecontoPredict <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_deconto_predict", outfile=outfiles, height=9.7/3, filetype=filetype, mar=rep(0, 4))

    chains <- list( ipr1$prchain,  pr1$prchain,
                   dipr1$prchain, dpr1$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.0005, 0.9995))
    xlim <- cictx$cis[[4]]
    ylim <- c(0.00, 4.0*1.04)

    col <- plotGetColors(3)[1:2]
    col <- c(rep(col[1], 2), rep(col[2], 2))
    lty <- rep(c("solid", "dotted"), 2)
    lwd <- 2

    par(mar=c(3, 3, 0.25, 0.75))
    plot.new()
    plot.window(xlim, ylim, xaxs="i", yaxs="i")

    pdfPlots(
        chains=chains,
        column=as.character(2100),
        lty=lty,
        col=col,
        burnin=F,
        xlab=daisSlrLab(),
        lwd=lwd,
        legendloc=NULL,
        yline=1,
        new=T
        )

    legend=paste(c(rep("Pfeffer", 2), rep("DeConto & Pollard", 2)),
                 c("expert assessment", "+paleo+obs+IPCC"))
    legend(
        "topright",
        legend=legend,
        bg="white",
        col=col,
        lty=c("dotted", "solid"),
        lwd=lwd
        )

    if (outfiles) { finDev(display=display) }
}


figLegend <- function(legends, title, lwd, col, shadecol, ccol, ...)
{
    legend(
        "center",
        legend=legends,
        title=title,
       #title.adj=0.1,
       #bty="n",
        fill=  c(col),
        border=c(col),
        )

   #box()
   #lnames <- rep(expression('' <= -17.5*degree*C, '' > -17.5*degree*C), 2)  # , '' <= 'fit', '' > 'fit')
   #lnames <- expression('' <= '-'*17.5*degree*C, '' > '-'*17.5*degree*C)
   #lwd    <- 2
#    legend(
#        "center",
#        legend=lnames,
#        title=title,
#       #title.adj=0.1,
#        col=col,
#        lwd=c(NA, NA), # , lwd, lwd),
#        pch=c(15, 15), #, NA,  NA),
#       #x.intersp=0.25,
#        x.intersp=0.125,
#        y.intersp=0.80,
#        text.width=0.925*strwidth(lnames)[1],
#       #pt.cex=1.67
#        seg.len=0.50,
#        pt.cex=1.25
#        )
}


# Inferred prior probability:  paleo and instrumental observations sharpen inference
figDecontoInfer <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_deconto_infer", outfile=outfiles, height=(2/3)*9.7, filetype=filetype, mar=rep(0,4))

    nfig <- 2
    plotLayout(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T), widths = c(10, 4), heights = rep(c(4, 10), nfig))

    xlim    <- c(-6, 11)
    ylim    <- c( 2, 23)
    points  <- c(1e6, rep(min(nrow(assimctx$chain), 1e5), 3))
    method  <- "outline"
    col     <- c("gray", plotGetColors(3)[1:2], "white")
    lty     <- c("dotted", rep("solid", 3))
   #title   <- "Interpretation"
    title   <- NULL
    smooth  <- c(3, rep(2, 3))
    xline   <- 2.25
    lwd     <- 2.0
    cex     <- 2.0
    legends <- c("Wide prior", "Pfeffer", "DeConto", "& Pollard")

    pairPlot(chains=list(prior_scaled,  as1_scaled,  das1_scaled),
             layout=F, xlim=xlim, ylim=ylim, cex=cex, method=method, legends=legends, legendfn=figLegend,
             points=points, col=col, lwd=lwd, lty=lty, title=title, smoothing=smooth, xline=xline,
             xlab=daisTcritLab(), ylab=daisLambdaMmLab(), topColumn="Tcrit", sideColumn="lambda",
             label="a", mar=c(5, 4)) # mar=c(4.25, 4))

    pairPlot(chains=list(prior_scaled, ias1_scaled, dias1_scaled),
             layout=F, xlim=xlim, ylim=ylim, cex=cex, method=method, legends=legends, legendfn=figLegend,
             points=points, col=col, lwd=lwd, lty=lty, title=title, smoothing=smooth, xline=xline,
             xlab=daisTcritLab(), ylab=daisLambdaMmLab(), topColumn="Tcrit", sideColumn="lambda",
             label="b", mar=c(3.50, 4))

    if (outfiles) { finDev(display=display) }
}


 figDecontoPredict(outfiles=T, filetype="qpng", display=T)
 figDecontoInfer(  outfiles=T, filetype="qpng", display=T)
