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
}


checkSamples <- function(assimctx=as1, prctx=pr1)
{
    assimctx$daisCmodel <- "daisRobOdeC"
    bar <- txtProgressBar(min=1, max=nrow(assimctx$chain), style=3)
    for (i in safefor(1:nrow(assimctx$chain))) {
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


plotBounds <- function(assimctx=as1, lwd=1.5)
{
    abline(v=assimctx$windows[assimctx$expert_ind, ], lty="dotted", lwd=lwd)
}


# "Deep uncertainty"
figAisPriors <- function()
{
    newDev("fig1", outfile=outfiles, width=8.5, height=11/2, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    col <- plotGetColors(3)
    #lty <- c("solid", "dashed", "dotted")  # "dotdash"
    lty <- rep("solid", 3)
    lwd <- 2
    pdfPlots(
        chains=chains,
        column=as.character(2100),
        lty=lty,
        col=col,
        burnin=F,
    #    xlim=c(0, max(cictx$range)),
    #    xlim=c(-0.2, 1.1),
        xlab=xlab,
        lwd=lwd,
        legendloc=NULL,
        smoothing=c(0.50, rep(1.25, 2)),
        yline=2
        )
    plotBounds()
    legend(
        "topright",
        legend=c(cnames, "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(rep(lwd, 3), 1.5),
        cex=0.75
        )

    caption <- paste("Figure 1. Probabilistic inversion of expert assessment with different priors")
    mtext(caption, outer=F, line=4, side=1, font=2)

    if (outfiles) { finDev() }
}


# "Probabilistic inversion works"
figCmpPriors <- function()
{
    newDev("fig2", width=8.5, height=7, outfile=outfiles, filetype=filetype)

    par(mfrow=c(2, 2))
    par(mar=c(4, 3, 0, 3))

    labels   <- c("a", "b", "c")
    col      <- plotGetColors(3)
    shadecol <- plotGetColors(3, 48)

    prctxs <- list(pr1, pr2, pr3)
    for (i in 1:length(prctxs)) {
        prctx <- prctxs[[i]]
        fname <-  fnames[i]

        assimctx <- prctx$assimctx
        pr       <- assimctx$expert_prior
        if (is.null(pr)) {
            pr      <- normPrior(assimctx$obsonly[assimctx$expert_ind], assimctx$windows[assimctx$expert_ind, 2])
            xlim    <- NULL
        } else {
            xlim    <- assimctx$windows[assimctx$expert_ind, ]
            margin  <- (xlim[2] - xlim[1]) * .05
            xlim[1] <- xlim[1] - margin
            xlim[2] <- xlim[2] + margin
        }
        xlim <- c(0, 0.8)
        if (assimctx$prior_name == "uniform") {
            smoothing <- 0.50
        } else {
            smoothing <- 1.25
        }
        priorPdfPlot(prctx$prchain, "2100", prior=pr, xlim=xlim, xlab=xlab, col=col[i], shadecol=shadecol[i], legends=c(fnames[i], "inversion"), smoothing=smoothing)
        labelPlot(labels[i])
    }

    # figure title
    caption <- paste("Figure 2. Probabilistic inversion of expert assessment by prior")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# Inferred prior probability:  paleo and instrumental observations sharpen inference
figUber <- function(assimctx=as1, outline=T)
{
    newDev("fig3", outfile=outfiles, width=8.5, height=4.25, filetype=filetype)

    layout(cbind(matrix(1:4, nrow=2, byrow=T), matrix(5:8, nrow=2, byrow=T)), widths = rep(c(10, 3), 2), heights = c(3, 10))

#    xlim   <- c(assimctx$lbound["Tcrit"],  assimctx$ubound["Tcrit"])
#    ylim   <- c(assimctx$lbound["lambda"], assimctx$ubound["lambda"])
    xlim <- c(-21, -9)
    ylim <- c(.004, .016)
    points <- ifelse(outline, 1e5, 6e3)
    method <- ifelse(outline, "outline", "points")
    col    <- plotGetColors(3)

    pairPlot( as1$chain,  as2$chain,  as3$chain, layout=F, units=assimctx$units, xlim=xlim, ylim=ylim, method=method,
        topColumn="Tcrit", sideColumn="lambda", legends=cnames, points=points, label="a", col=col, smoothing=rep(2, 3))

    pairPlot(ias1$chain, ias2$chain, ias3$chain, layout=F, units=assimctx$units, xlim=xlim, ylim=ylim, method=method,
        topColumn="Tcrit", sideColumn="lambda", legends=cnames, points=points, label="b", col=col, smoothing=rep(2, 3))

    caption <- paste("Figure 3. Inferred prior probability; (a) Expert assessment only, (b) All data")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


figUber2 <- function(assimctx=as1)
{
    newDev("fig3_2", outfile=outfiles, width=8.5, height=4.25, filetype=filetype)

    layout(cbind(matrix(1:4, nrow=2, byrow=T), matrix(5:8, nrow=2, byrow=T)), widths = rep(c(10, 3), 2), heights = c(3, 10))

    points <- 1e5

    xlim <- c(assimctx$lbound["Tcrit"],  assimctx$ubound["Tcrit"])
    ylim <- c(assimctx$lbound["lambda"], assimctx$ubound["lambda"])

    col   <- c("#D00000", "#0000D0")
    n     <- 31
    reds  <- paste("#",     toupper(as.hexmode(floor(seq(128, 255, length.out=n)))), "0000", as.hexmode(200), sep="")
    blues <- paste("#0000", toupper(as.hexmode(floor(seq(128, 255, length.out=n)))),         as.hexmode(200), sep="")
    ccol  <- list(reds, blues)

    #n <- 31; pie(rep(1, n), col=paste("#", toupper(as.hexmode(floor(seq(128, 255, length.out=n)))), "0000", as.hexmode(200), sep=""))

    pairPlot( as1$chain,  as2$chain, topColumn="Tcrit", sideColumn="lambda",
        layout=F, units=assimctx$units, xlim=xlim, ylim=ylim, method="contours",
        legends=cnames[1:2], points=points, label="a", col=col, ccol=ccol)

    pairPlot(ias1$chain, ias2$chain, topColumn="Tcrit", sideColumn="lambda",
        layout=F, units=assimctx$units, xlim=xlim, ylim=ylim, method="contours",
        legends=cnames[1:2], points=points, label="b", col=col, ccol=ccol)

    caption <- paste("Figure 3. Inferred prior probability; (a) Expert assessment only, (b) All data")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# Predicted AIS volume loss in 2100 with all observations.
figPredict <- function(assimctx=as1)
{
    newDev("fig4", outfile=outfiles, width=8.5, height=8, filetype=filetype)

    chains <- list(ipr1$prchain, ipr2$prchain, ipr3$prchain)

    pdfCdfPlots(
        legends=c(cnames, "Pfeffer"),
        legendloc="topleft",
        col=c(plotGetColors(3), "black"),
        lty=c(rep("solid", 3), "dotted"),
        column=as.character(2100),
        chains=chains,
        xlab=xlab,
        burnin=F,
        log=T, survival=T,
        vlines=assimctx$windows[assimctx$expert_ind, ],
        smoothing=c(0.50, rep(1.25, 2))
        )

    # figure title
    caption <- paste("Figure 4. Probabilistic inversion with paleo and instrumental observations")
    mtext(caption, outer=F, line=4, side=1, font=2)

    if (outfiles) { finDev() }
}


figCmpPredict <- function(assimctx=as1)
{
    newDev("fig5", outfile=outfiles, width=8.5, height=11, filetype=filetype)

    par(omi=c(0.25, 0, 0.25, 0))
    layout(rbind(matrix(1:4, nrow=2), matrix(5:8, nrow=2)))

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

    caption <- paste("Figure 5. Expert assessment vs. all data")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# compare PDFs with/without all observations
figCmpInst <- function()
{
    newDev("fig5_2", outfile=outfiles, width=8.5, height=7, filetype=filetype)

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
        lwd=c(rep(lwd, 6), 1.5),
        cex=0.75
        )

    caption <- paste("Figure n. Expert assessment vs. all data")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


if (outfiles) { finDev() }
