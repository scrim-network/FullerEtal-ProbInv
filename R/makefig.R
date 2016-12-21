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
# makefig.R

outfiles <- T
year <- 2100
#filetype <- "png"
filetype <- "pdf"

if (file.exists(       'ep="b";n=5e6.RData')) {
    iter <- "5e6"
} else if (file.exists('ep="b";n=2e6.RData')) {
    iter <- "2e6"
} else if (file.exists('ep="b";n=5e5.RData')) {
    iter <- "5e5"
}
#iter <- "5e+6"
#iter <- "2e+6"
#iter <- "5e+5"

source('plot.R')


fnames <- c("uniform", "beta", "normal")
cnames <- capitalize(fnames)


if (!exists("pr1")) {
    loadChains(paste('ep="', substr(fnames, 1, 1), '";n=', iter, ".RData", sep=""))
    loadChains(paste('ip="', substr(fnames, 1, 1), '";n=', iter, ".RData", sep=""), newnames=c("ias", "ipr"))
   #load('p="u";n=5e5.RData')
   #load('prior="uniform";nbatch=5e6.RData')
   #as1 <-   daisctx
   #pr1 <- prdaisctx

    # rejection sample uniform prior
    source('calib.R')
    daisRejSample(assimctx=as1, prctx=pr1)
   #daisRunPredict(nbatch=5000, assimctx=as1, prctx=pr1)
}


checkSamples <- function(assimctx=as1, prctx=pr1)
{
    assimctx$daisCmodel <- "daisRobOdeC"
    bar <- txtProgressBar(min=1, max=nrow(assimctx$chain), style=3)
    for (i in safefor(1:nrow(assimctx$chain))) {

        # TODO:  can skip running models and comparison when mp don't change (MCMC reject);
        #        would give a four-fold speed up assuming accept rate around 0.25
       #y     <- assimctx$modelfn(assimctx$chain[i, ], assimctx)
        y1 <- C_daisModel(       assimctx$chain[i, ], assimctx)
        y2 <- F_daisFastDynModel(assimctx$chain[i, ], assimctx)

        y1_std <- y1[assimctx$obs_ind[assimctx$expert_ind]] - y1[assimctx$SL.expert]
        y2_std <- y2[assimctx$obs_ind[assimctx$expert_ind]] - y2[assimctx$SL.expert]
        if (!isTRUE(all.equal(y1_std, y2_std, prctx$prchain[i, ], check.names=F, check.attributes=F))) {
            print(paste(   i, y1_std, y2_std, prctx$prchain[i, ]))
        }
        setTxtProgressBar(bar, i)
    }
    close(bar)
}


xlab <- paste("Projected AIS Volume Loss in", year, "[SLE m]")


plotBounds <- function(assimctx=as1, lwd=1)
{
   #abline(v=assimctx$windows[assimctx$expert_ind, ], lty="dotted", lwd=lwd)
    abline(v=assimctx$windows[assimctx$expert_ind, ], lty="solid",  lwd=lwd)
}



# "Deep uncertainty"
figAisPriors <- function(assimctx=as1)
{
    newDev("fig_deep", outfile=outfiles, height=9.7/3, filetype=filetype, mar=rep(0, 4))

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.0005, 0.9995))
    xlim   <- cictx$cis[[3]]
    ylim   <- c(0, 4)

    nfig <- 1
    plotLayout(cbind(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T)), heights = c(1, rep(10, nfig)))

    par(mar=c(0, 3, 0.25, 0.25))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$windows[assimctx$expert_ind, ], label="Range by Pfeffer et al. (2008)", offset=0)

    col    <- plotGetColors(3)
    lty    <- rep("solid", 3)
    lwd    <- 2

    par(mar=c(3, 3, 0.25, 0.25))
    plot.new()
    plot.window(xlim, ylim, xaxs="i")
    plotBounds()
    pdfPlots(
        chains=chains,
        column=as.character(2100),
        lty=lty,
        col=col,
        burnin=F,
        xlab=xlab,
        lwd=lwd,
        legendloc=NULL,
        yline=1,
        new=T
        )

    legend(
        "topright",
        legend=c("Pfeffer et al. (2008)", paste(cnames, "interpretation")),
        bg="white",
        col=c("black", col),
        lty=c("solid", lty),
        lwd=c(1, rep(lwd, 3))
        )

    if (outfiles) { finDev() }
}


# "Probabilistic inversion works"
figCmpPriors <- function(assimctx=as1)
{
    newDev("fig_invert", outfile=outfiles, height=9.7 / 2, filetype=filetype, mar=rep(0, 4))

    nfig <- 3
    plotLayout(cbind(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T)), heights = c(1.5, rep(10, nfig)))

    xlim <- c(0, 0.8)
    ylim <- c(0, 4.0)

    par(mar=c(0, 3.5, 0.25, 0.75))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$windows[assimctx$expert_ind, ], label="Range by Pfeffer et al. (2008)", offset=0)

    labels   <- c("a", "b", "c")
    col      <- plotGetColors(3)
    shadecol <- plotGetColors(3, 96)

    prctxs <- list(pr1, pr2, pr3)
    for (i in 1:length(prctxs)) {
        par(mar=c(ifelse(i==length(prctxs), 3, 2), 3.5, 0.25, 0.75))
        plot.new()
        plot.window(xlim, ylim, xaxs="i")
        plotBounds()

        prctx    <- prctxs[[i]]
        assimctx <- prctx$assimctx
        pr       <- assimctx$expert_prior
        if (is.null(pr)) {
            pr   <- normPrior(assimctx$obsonly[assimctx$expert_ind], assimctx$windows[assimctx$expert_ind, 2])
        }

        priorPdfPlot(prctx$prchain, "2100", prior=pr, xlim=xlim, ylim=ylim, xlab=ifelse(i==length(prctxs), xlab, ""), col=col[i], shadecol=shadecol[i], legends=NULL, new=T)
        labelPlot(labels[i], line=2.5)

        if (i==2) {
            legend(
                "topright",
                legend=c("Pfeffer et al. (2008)", paste(cnames, "inversion"), paste(cnames, "prior")),
                col=c("black", col, shadecol),
                lty=c("solid", rep("solid", 3), rep(NA, 3)),
                lwd=c(1,  rep(2,  3), rep(NA, 3)),
                pch=c(NA, rep(NA, 3), rep(15, 3)),
                pt.cex=2.25,
                bg="white"
                )
        }
    }

    if (outfiles) { finDev() }
}


# Inferred prior probability:  paleo and instrumental observations sharpen inference
figInfer <- function(assimctx=as1, outline=T)
{
    newDev("fig_infer", outfile=outfiles, height=(2/3)*9.7, filetype=filetype, mar=rep(0,4))

    nfig <- 2
    plotLayout(cbind(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T)), widths = c(10, 3), heights = rep(c(3, 10), nfig))

   #xlim <- c(assimctx$lbound["Tcrit"],  assimctx$ubound["Tcrit"])
   #ylim <- c(assimctx$lbound["lambda"], assimctx$ubound["lambda"])
    xlim <- c(-21, -9)
    ylim <- c(.004, .016)
    points <- ifelse(outline, min(nrow(assimctx$chain), 1e5), 6e3)
    method <- ifelse(outline, "outline", "points")
    col    <- plotGetColors(3)
    title  <- "Interpretation"
    smooth <- rep(2, 3)

    pairPlot( as1$chain,  as2$chain,  as3$chain, layout=F, units=assimctx$units, xlim=xlim, ylim=ylim,
             method=method, legends=cnames, points=points, col=col, title=title, smoothing=smooth,
             topColumn="Tcrit", sideColumn="lambda", label="a", mar=c(4.25, 4))

    pairPlot(ias1$chain, ias2$chain, ias3$chain, layout=F, units=assimctx$units, xlim=xlim, ylim=ylim,
             method=method, legends=cnames, points=points, col=col, title=title, smoothing=smooth,
             topColumn="Tcrit", sideColumn="lambda", label="b")

    if (outfiles) { finDev() }
}


figPdfCdf <- function(chains, col, lty, ylim=c(0, 4), xval=2100, column=as.character(xval), labels=c("a", "b"), assimctx=as1)
{
    cictx  <- ciCalc(chains=chains, xvals=xval, probs=c(0.0005, 0.9995))
    xlim   <- cictx$cis[[3]]

    par(mar=c(0, 4, 0.25, 1))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$windows[assimctx$expert_ind, ], label="Range by Pfeffer et al. (2008)", offset=0)

    par(mar=c(2, 4, 0.25, 1))
    plot.new()
    plot.window(xlim, ylim=ylim, xaxs="i")
    plotBounds()
    pdfPlots(
        chains=chains,
        column=column,
        col=col,
        lty=lty,
        burnin=F,
        xlab=xlab,
        legendloc=NULL,
        yline=2,
        new=T
        )
    labelPlot(labels[1])

    par(mar=c(3, 4, 0.25, 1))
    plot.new()
    ylim_cdf <- c(1e-3, 1)
    plot.window(xlim, ylim_cdf, xaxs="i", log="y")
    plotBounds()
    cdfPlots(
        chains=chains,
        column=column,
        xlim=xlim,
        ylim=ylim_cdf,
        xlab=xlab,
        col=col,
        lty=lty,
        log=T, survival=T,
        new=T
        )
    labelPlot(labels[2], where="log")
}


# Predicted AIS volume loss in 2100 with all observations.
figPredict <- function(assimctx=as1)
{
    newDev("fig_predict", outfile=outfiles, height=9.7/3, filetype=filetype, mar=rep(0, 4))

    nfig <- 2
    plotLayout(cbind(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T)), heights = c(1.5, rep(10, nfig)))

    chains <- list(ipr1$prchain, ipr2$prchain, ipr3$prchain)
    col    <- c(plotGetColors(3), "black")
    lty    <- c(rep("solid", 3), "solid")

    figPdfCdf(chains=chains, col=col, lty=lty)

    legend(
        "bottomleft",
        legend=c("Pfeffer et al. (2008)", paste(cnames, "interpretation")),
        bg="white",
        lty=c("solid", lty),
        col=c("black", col),
        lwd=c(1, rep(2, 3))
        )

    if (outfiles) { finDev() }
}


figCmpPredict <- function(assimctx=as1)
{
    newDev("fig_cdf", outfile=outfiles, width=8.5, height=11, filetype=filetype)

   #par(omi=c(0.25, 0, 0.25, 0))
    plotLayout(rbind(matrix(1:4, nrow=2), matrix(5:8, nrow=2)))

    legends <- c("exp only", "all data")
    lty  <- c("dotted", "solid", "dotted")
    col  <- plotGetColors(3)
    xlim <- c(0, 0.8)

    pdfCdfPlots(
        layout=F,
        legends=c(paste(legends, fnames[1]), "Pfeffer"),
        legendloc="topleft",
        col=c(rep(col[1], 2), "black"),
        lty=lty,
        column=as.character(2100),
        chains=list(pr1$prchain, ipr1$prchain),
        xlab=xlab,
        burnin=F,
        log=T, survival=T,
        vlines=assimctx$windows[assimctx$expert_ind, ],
        smoothing=c(rep(0.5, 2)),
        xlim=xlim,
        labels=c("a", "b")
        )

    pdfCdfPlots(
        layout=F,
        legends=c(paste(legends, fnames[2]), "Pfeffer"),
        legendloc="topright",
        col=c(rep(col[2], 2), "black"),
        lty=lty,
        column=as.character(2100),
        chains=list(pr2$prchain, ipr2$prchain),
        xlab=xlab,
        burnin=F,
        log=T, survival=T,
        vlines=assimctx$windows[assimctx$expert_ind, ],
        smoothing=c(rep(1.25, 2)),
        xlim=xlim,
        labels=c("c", "d")
        )

    pdfCdfPlots(
        layout=F,
        legends=c(paste(legends, fnames[3]), "Pfeffer"),
        legendloc="topleft",
        col=c(rep(col[3], 2), "black"),
        lty=lty,
        column=as.character(2100),
        chains=list(pr3$prchain, ipr3$prchain),
        xlab=xlab,
        burnin=F,
        log=T, survival=T,
        vlines=assimctx$windows[assimctx$expert_ind, ],
        smoothing=c(rep(1.25, 2)),
        xlim=xlim,
        labels=c("e", "f")
        )

   #caption <- paste("Figure 5. Expert assessment vs. all data")
   #mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# compare PDFs with/without all observations
figCmpInst <- function()
{
    newDev("fig_inst", outfile=outfiles, width=8.5, height=7, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain, ipr1$prchain, ipr2$prchain, ipr3$prchain)
    lty <- c(rep("dotted", 3), rep("solid", 3))
    col <- rep(plotGetColors(3), 2)
    lwd <- 1.5

    pdfPlots(
        chains=chains,
        column=as.character(2100),
        burnin=F,
        col=col,
        lty=lty,
        #xlim=c(0, max(cictx$range)),
        #xlim=c(-0.2, 1.1),
        xlab=xlab,
        smoothing=rep(c(0.50, rep(1.5, 2)), 2),
        legendloc=NULL,
        #yline=2,
        lwd=lwd
        )
    plotBounds()
    legend(
        "topright",
        legend=c(paste("exp only", fnames), paste("all data", fnames), "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(rep(lwd, 6), 1.5)
        )

    caption <- paste("Figure n. Expert assessment vs. all data")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


if (outfiles) { finDev() }
