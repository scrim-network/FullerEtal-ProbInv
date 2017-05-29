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
# makefig.R

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
source('calib.R')        # daisRejSample()
loadLibrary("MCMCpack")  # figPlotMarginals(), invGammaPrior()


fnames <- c("uniform", "beta", "normal")
cnames <- capitalize(fnames)


if (!exists("pr1")) {
    loadChains(paste('ep="', substr(fnames, 1, 1), '";n=', iter, ".RData", sep=""))
    loadChains(paste('ip="', substr(fnames, 1, 1), '";n=', iter, ".RData", sep=""), newnames=c("ias", "ipr"))
   #load('p="u";n=5e5.RData')
   #as1 <-   daisctx
   #pr1 <- prdaisctx

    # rejection sample uniform prior
    daisRejSample(assimctx=as1, prctx=pr1)
   #daisRunPredict(nbatch=5000, assimctx=as1, prctx=pr1)

     as1_scaled <- daisScaleChain( as1$chain, assimctx=as1)
     as2_scaled <- daisScaleChain( as2$chain, assimctx=as1)
     as3_scaled <- daisScaleChain( as3$chain, assimctx=as1)
    ias1_scaled <- daisScaleChain(ias1$chain, assimctx=as1)
    ias2_scaled <- daisScaleChain(ias2$chain, assimctx=as1)
    ias3_scaled <- daisScaleChain(ias3$chain, assimctx=as1)

    Tcrit  <- -as1$ Tcrit_prior$rand(5e6)
    lambda <-  as1$lambda_prior$rand(5e6)
    priorChain <- cbind(Tcrit, lambda)
    colnames(priorChain) <- c("Tcrit", "lambda")
    prior_scaled <- daisScaleChain(priorChain, assimctx=as1)
    rm(Tcrit, lambda, priorChain)  # keep it clean
}


checkRejSamples <- function(assimctx=as1, prctx=pr1)
{
    assimctx$daisCmodel <- "daisRobOdeC"
    bar <- txtProgressBar(min=1, max=nrow(assimctx$chain), style=3)
    for (i in safefor(1:nrow(assimctx$chain))) {

        # TODO:  can skip running models and comparison when mp don't change (MCMC reject);
        #        would give a four-fold speed up assuming accept rate around 0.25
       #y     <- assimctx$modelfn(assimctx$chain[i, ], assimctx)
        y1 <- C_daisModel(       assimctx$chain[i, ], assimctx)
        y2 <- F_daisFastDynModel(assimctx$chain[i, ], assimctx)

        y1_std <- y1[assimctx$SL.2100] - y1[assimctx$SL.expert]
        y2_std <- y2[assimctx$SL.2100] - y2[assimctx$SL.expert]
        if (!isTRUE(all.equal(y1_std, y2_std, prctx$prchain[i, ], check.names=F, check.attributes=F))) {
            print(paste(   i, y1_std, y2_std, prctx$prchain[i, ]))
        }
        setTxtProgressBar(bar, i)
    }
    close(bar)
}


figPlotBounds <- function(assimctx=as1, lwd=1)
{
   #abline(v=assimctx$expert_window, lty="dotted", lwd=lwd)
    abline(v=assimctx$expert_window, lty="solid",  lwd=lwd)
}


# "Deep uncertainty"
figAisPriors <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_deep", outfile=outfiles, height=9.7/3, filetype=filetype, mar=rep(0, 4))

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.0005, 0.9995))
    xlim   <- cictx$cis[[3]]
    ylim   <- c(0, 4*1.04)

    nfig <- 1
    plotLayout(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T), heights = c(1, rep(10, nfig)))

    par(mar=c(0, 3, 0.25, 0.25))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$expert_window, label="Expert assessment range", offset=0)

    col    <- plotGetColors(3)
    lty    <- rep("solid", 3)
    lwd    <- 2

    par(mar=c(3, 3, 0.25, 0.25))
    plot.new()
    plot.window(xlim, ylim, xaxs="i", yaxs="i")
    figPlotBounds()
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

    legend(
        "topright",
        legend=c("Expert assessment", paste(cnames, "interpretation")),
        bg="white",
        col=c("black", col),
        lty=c("solid", lty),
        lwd=c(1, rep(lwd, 3))
        )

    if (outfiles) { finDev(display=display) }
}


# "Probabilistic inversion works"
figCmpPriors <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_invert", outfile=outfiles, height=9.7 / 2, filetype=filetype, mar=rep(0, 4))

    nfig <- 3
    plotLayout(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T), heights = c(1.5, rep(10, nfig)))

    xlim <- c(0, 0.8)
    ylim <- c(0, 4.0*1.04)

    par(mar=c(0, 3.5, 0.25, 0.75))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$expert_window, label="Expert assessment range", offset=0)

    labels   <- c("a", "b", "c")
    col      <- plotGetColors(3)
    shadecol <- plotGetColors(3, 96)

    prctxs <- list(pr1, pr2, pr3)
    for (i in 1:length(prctxs)) {
        par(mar=c(ifelse(i==length(prctxs), 3, 2), 3.5, 0.25, 0.75))
        plot.new()
        plot.window(xlim, ylim, xaxs="i", yaxs="i")
        figPlotBounds()

        prctx <- prctxs[[i]]
        priorPdfPlot(prctx$prchain, as.character(2100), prior=prctx$assimctx$expert_prior, xlim=xlim, ylim=ylim, xlab=ifelse(i==length(prctxs), daisSlrLab(), ""), col=col[i], shadecol=shadecol[i], legends=NULL, new=T)
        labelPlot(labels[i], line=2.5)

        if (i==2) {
            legend(
                "topright",
                legend=c("Expert assessment", paste(cnames, "inversion"), paste(cnames, "prior")),
                col=c("black", col, shadecol),
                lty=c("solid", rep("solid", 3), rep(NA, 3)),
                lwd=c(1,  rep(2,  3), rep(NA, 3)),
                seg.len=1,
                pch=c(NA, rep(NA, 3), rep(15, 3)),
                pt.cex=2.25,
                bg="white"
                )
        }
    }

    if (outfiles) { finDev(display=display) }
}


figLegend <- function(legends, title, lwd, col, shadecol, ccol, ...)
{
    legends <- c("Wide prior", paste(cnames, "interp."))

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
figInfer <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_infer", outfile=outfiles, height=(2/3)*9.7, filetype=filetype, mar=rep(0,4))

    nfig <- 2
    plotLayout(matrix(1:(4*nfig), nrow=(2*nfig), byrow=T), widths = c(10, 4), heights = rep(c(4, 10), nfig))

    xlim    <- c( -6, 11)
    ylim    <- c(  2, 19)
    points  <- c(1e6, rep(min(nrow(assimctx$chain), 1e5), 3))
    method  <- "outline"
    col     <- c("gray", plotGetColors(3))
    lty     <- c("dotted", rep("solid", 3))
   #title   <- "Interpretation"
    title   <- NULL
    smooth  <- c(3, rep(2, 3))
    xline   <- 2.25
    lwd     <- 2.0
    cex     <- 2.0
    legends <- c("Wide prior", paste(cnames, "Interp."))

    pairPlot(chains=list(prior_scaled,  as1_scaled,  as2_scaled,  as3_scaled),
             layout=F, xlim=xlim, ylim=ylim, cex=cex, method=method, legends=legends, legendfn=figLegend,
             points=points, col=col, lwd=lwd, lty=lty, title=title, smoothing=smooth, xline=xline,
             xlab=daisTcritLab(), ylab=daisLambdaMmLab(), topColumn="Tcrit", sideColumn="lambda",
             label="a", mar=c(5, 4)) # mar=c(4.25, 4))

    pairPlot(chains=list(prior_scaled, ias1_scaled, ias2_scaled, ias3_scaled),
             layout=F, xlim=xlim, ylim=ylim, cex=cex, method=method, legends=legends, legendfn=figLegend,
             points=points, col=col, lwd=lwd, lty=lty, title=title, smoothing=smooth, xline=xline,
             xlab=daisTcritLab(), ylab=daisLambdaMmLab(), topColumn="Tcrit", sideColumn="lambda",
             label="b", mar=c(3.50, 4))

    if (outfiles) { finDev(display=display) }
}


# helper function
figPdfCdf <- function(chains, col, lty, xlim, ylim=c(0, 4.25), labels=c("a", "b"), column=as.character(2100), bottom=c(2, 3), assimctx=as1)
{
    par(mar=c(bottom[1], 4, 0.25, 1))
    plot.new()
    plot.window(xlim, ylim=ylim, xaxs="i", yaxs="i")
    figPlotBounds()
    pdfPlots(
        chains=chains,
        column=column,
        col=col,
        lty=lty,
        burnin=F,
        xlab=daisSlrLab(),
        legendloc=NULL,
        yline=2,
        new=T
        )
    labelPlot(labels[1])

    par(mar=c(bottom[2], 4, 0.25, 1))
    plot.new()
    ylim_cdf <- c(1e-3, 1)
    plot.window(xlim, ylim_cdf, xaxs="i", log="y")
    figPlotBounds()
    cdfPlots(
        chains=chains,
        column=column,
        xlim=xlim,
        ylim=ylim_cdf,
        xlab=daisSlrLab(),
        col=col,
        lty=lty,
        log=T, survival=T,
        new=T
        )
    labelPlot(labels[2], where="log")
}


figPlotMarginals <- function(fname, chains, assimctx=ias1, outfiles=T, filetype="pdf", display=T)
{
    # originally 6x9
    newDev(fname, outfile=outfiles, width=6, height=7, filetype=filetype, mar=rep(0, 4))
    plotLayout(matrix(1:16, ncol=4, byrow=T))

    xlabs <- list(expression(gamma), expression(alpha), expression(mu), expression(nu), "P0", expression(kappa), "f0", "h0", "c", "b0", "slope", expression(T[crit]), expression(lambda), expression(sigma[P]^{2}), expression(sigma[I]^{2}))

    lty <- c("solid", "dashed")
    col <- c("black", "red")
    lwd <- 1

    par(mar=c(4.5, 4, 1.5, 0.50))

    lbound <- c(assimctx$lbound, assimctx$lbound_sp)
    ubound <- c(assimctx$ubound, assimctx$ubound_sp)
    ubound["var.paleo"] <-   4
    lbound["lambda"]    <-   0
    lbound["Tcrit"]     <- -25
    ubound["Tcrit"]     <-  -5

    for (i in safefor(1:ncol(chains[[1]]))) {
        column <- colnames(assimctx$chain)[i]

        if (i > ncol(chains[[2]])) {
            pdfctx  <- pdfCalc(chains=list(chains[[1]]), column=column, smoothing=rep(2, 2))
        } else {
            pdfctx  <- pdfCalc(chains=chains,            column=column, smoothing=rep(2, 2))
        }
        ylim    <- pdfctx$ylim
        ylim[2] <- pdfctx$ylim[2] + .04 * (pdfctx$ylim[2] - pdfctx$ylim[1])

        xlim  <- c(lbound[i], ubound[i])
       #if (column == "var.paleo") {
           #xlim <- c(gtzero(), pdfctx$xlim[2])
           #xlim <- c(0, 5)
       #}

        prlim   <- xlim
        range   <- xlim[2] - xlim[1]
        xlim[1] <- xlim[1] - 0.04 * range
        xlim[2] <- xlim[2] + 0.04 * range

        plot.new()
        plot.window(xlim, ylim, xaxs="i", yaxs="i")

        switch (column,
        Tcrit={
            pr <- assimctx$Tcrit_prior
            priorPlot(pr, col="light gray", lty="dashed", lwd=lwd, xlim=xlim, shade=T, negate=T)
        },
        lambda={
            pr <- assimctx$lambda_prior
        },
        var.paleo={
            pr <- invGammaPrior(assimctx$alpha, assimctx$beta)
        }, {
            pr <- uniformPrior(prlim[1], prlim[2])
        })
        priorPlot(pr, col="light gray", lty="dashed", lwd=lwd, xlim=xlim, shade=T, negate=F)

        pdfPlotWindow(pdfctx, col=col, lty=lty, lwd=lwd, xlab=xlabs[[i]], ylab="PDF", xlim=xlim, ylim=ylim, xline=2.5, yline=1, new=T)
    }

    par(mar=c(4.5, 2, 1.5, 2.50))
    plot.new()
    legend(
        "center",
        legend=c("Wide prior", NA, "Wide+expert", NA, "Wide+expert", "+paleo+obs", "+IPCC"),
        col=c(         "gray", NA, "red",         NA, "black",       NA,            NA),
        lty=c(             NA, NA, "dashed",      NA, "solid",       NA,            NA),
        lwd=c(             NA, NA, 2,             NA, 2,             NA,            NA),
        seg.len=1,
        pch=c(             15, NA, NA,            NA, NA,            NA,            NA),
        pt.cex=2.25,
        bg="white",
        bty="n"
        )
   #box()

    if (outfiles) { finDev(display=display) }
}


figMarginal <- function(assimctx=ias1, outfiles=T, filetype="pdf", display=T)
{
    figPlotMarginals("fig_marg_u", list(ias1$chain, as1$chain), assimctx=assimctx, outfiles=outfiles, filetype=filetype, display=display)
    figPlotMarginals("fig_marg_b", list(ias2$chain, as2$chain), assimctx=assimctx, outfiles=outfiles, filetype=filetype, display=display)
    figPlotMarginals("fig_marg_n", list(ias3$chain, as3$chain), assimctx=assimctx, outfiles=outfiles, filetype=filetype, display=display)
}


# Predicted AIS volume loss in 2100 with all observations
figPredict <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_predict", outfile=outfiles, height=9.7/3, filetype=filetype, mar=rep(0, 4))

    nfig <- 2
    plotLayout(matrix(1:(nfig + 1), nrow=(nfig + 1), byrow=T), heights = c(1.5, rep(10, nfig)))

    chains <- list(ipr1$prchain, ipr2$prchain, ipr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.0005, 0.9995))
    xlim   <- cictx$cis[[3]]

    par(mar=c(0, 4, 0.25, 1))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$expert_window, label="Expert assessment range", offset=0)

    col    <- c(plotGetColors(3), "black")
    lty    <- c(rep("solid", 3), "solid")
    figPdfCdf(chains=chains, col=col, lty=lty, xlim=xlim)
    legend(
        "bottomleft",
        legend=c("Expert assessment", paste(cnames, "interp+paleo+obs+IPCC")),
        bg="white",
        lty=c("solid", lty),
        col=c("black", col),
        lwd=c(1, rep(2, 3))
        )

    if (outfiles) { finDev(display=display) }
}


figCmpPredict <- function(assimctx=as1, outfiles=T, filetype="pdf", display=T)
{
    # is single=F too wide? Tony used 6 inches for his marginal PDFs
   #newDev("fig_cdf", outfile=outfiles, single=F, height=9.7/2, filetype=filetype, mar=rep(0, 4))
    newDev("fig_cdf", outfile=outfiles, width=6,  height=9.7/2, filetype=filetype, mar=rep(0, 4))

    nfig <- 6
    plotLayout(matrix(1:8, nrow=4, byrow=T), heights = c(1.5, rep(10, nfig / 2)))

    xlim <- c(0, 0.8)

    # draw arrow at top of first column
    par(mar=c(0, 4, 0.25, 1))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$expert_window, label="Expert assessment range", offset=0)

    # draw arrow at top of second column
    par(mar=c(0, 4, 0.25, 1))
    plot.new()
    plot.window(xlim, c(0, 1), xaxs="i")
    plotArrowX(xlim=assimctx$expert_window, label="Expert assessment range", offset=0)

    lty  <- c("dotted", "solid")
    col  <- plotGetColors(3)

    legends <- c("interpretation", "+paleo+obs+IPCC")

    figPdfCdf(chains=list(pr1$prchain, ipr1$prchain), col=rep(col[1], 2), lty=lty, xlim=xlim, bottom=c(2, 2))
    legend(
        "bottomleft",
        legend=c("Expert assessment", paste(cnames[1], legends)),
        bg="white",
        lty=c("solid", lty),
        col=c("black", rep(col[1], 2)),
        lwd=c(1, rep(2, 3))
        )

    figPdfCdf(chains=list(pr2$prchain, ipr2$prchain), col=rep(col[2], 2), lty=lty, xlim=xlim, bottom=c(2, 2), labels=c("c", "d"))
    legend(
        "bottomleft",
        legend=c("Expert assessment", paste(cnames[2], legends)),
        bg="white",
        lty=c("solid", lty),
        col=c("black", rep(col[2], 2)),
        lwd=c(1, rep(2, 3))
        )

    figPdfCdf(chains=list(pr3$prchain, ipr3$prchain), col=rep(col[3], 2), lty=lty, xlim=xlim, bottom=c(3, 3), labels=c("e", "f"))
    legend(
        "bottomleft",
        legend=c("Expert assessment", paste(cnames[3], legends)),
        bg="white",
        lty=c("solid", lty),
        col=c("black", rep(col[3], 2)),
        lwd=c(1, rep(2, 3))
        )

    if (outfiles) { finDev(display=display) }
}


figPlotHindcast <- function(assimctx=ias1, prctx=prdaisctx, prexp=NULL, meancol="black")
{
    present <- 2010
    xvals   <- -150000:0
    rows    <- tsGetIndicesByRange(assimctx$frc_ts, lower=xvals[1]+present, upper=present)

    cictx <- prctx$hindQuant
    ci_lo <- cictx$cis  [[1]][rows, 1]
    ci_hi <- cictx$cis  [[1]][rows, 2]
    means <- cictx$means[[1]][rows]
   #ylim  <- range(ci_lo, ci_hi)
    ylim  <- c(-18, 10)

    par(mar=c(3, 4, 1, 0.25))
    plot.new()
    plot.window(c(xvals[1], last(xvals)), ylim)
    axis(1)
    axis(2)
   #mtext(side=1, text='Year (before present)',        line=2)
   #mtext(side=2, text='Antarctic Ice Sheet SLR (m)',  line=1)
    title(xlab="Year (before present)", ylab="Antarctic Ice Sheet SLR (m)", line=2)
    box()

    # 5-95% range
    col <- "goldenrod1"
    if (T) {
        polygon(c(xvals, rev(xvals), xvals[1]), c(ci_lo, rev(ci_hi), ci_lo[1]), col=col, border=NA)
    } else {
        lines(xvals, ci_lo, col=col, lty="dotted", lwd=1)
        lines(xvals, ci_hi, col=col, lty="dotted", lwd=1)
    }

    # the zero line of SLR
    lines(c(-1e6, 1e6), c(0, 0), type='l', lty=2, col='black')

    # model mean with only expert assessment
    if (!is.null(prexp)) {
        exp_means <- prexp$hindQuant$means[[1]][rows]
        ind <- seq(from=1, to=length(xvals), length.out=1000)
        lines(xvals[ind], exp_means[ind], col="red", lty="dashed", lwd=1)
    }

    obs.years <- c(-118000, -18000, -4000, 2002) - present
    spread    <- 1000

    # paleoclimatic observations
    ts <- cbind(obs.years, assimctx$obsonly, 2*assimctx$error)
    colnames(ts) <- c("time", "SLE", "error")
    tsErrorBars(ts[1:3, ], xbeam=F, ibeam=T, shade=F, col="purple", tick=spread, lwd=1.5)

    # Kelsey likes an asterisk; cross: pch=3, cex=0.5; pch='I'?
    points(ts[1:3, "time"], ts[1:3, "SLE"], pch=8, cex=0.75, col="purple")
   #points(ts[4,   "time"], ts[4,   "SLE"], pch=8, cex=0.50, col="purple")

    # model mean with all data
    lines(xvals, means, col=meancol, lty="solid", lwd=1)

    # instrumental observation
    i <- 4
    points(obs.years[i], assimctx$obsonly[i], pch=15, cex=0.50, col="purple")

    prior <- prctx$assimctx$prior_name
    legend(
        -92500, 11,
        legend=c(paste("5-95% range,", prior, "intepreration of expert assessment+observations+paleo+IPCC"),
                 "2-sigma range, observations",
                 paste("Mean,", prior, "interpretation of expert assessment+observations+paleo+IPCC"),
                 paste("Mean,", prior, "interpretation of expert assessment only")),
        col=c(col, "purple", meancol, "red"),
        lty=c( NA, "solid",  "solid", "dotted"),
        lwd=c( NA, 1.5,          1.5, 1.5),
        pch=c( 15, NA,            NA, NA),
        seg.len=1,
        pt.cex=2.25,
        bg="white",
        bty="n"
        )
}


figHindcast <- function(assimctx=ias1, outfiles=T, filetype="pdf", display=T)
{
    newDev("fig_hindcast", outfile=outfiles, width=7, height=7, filetype=filetype, mar=rep(0, 4))

    plotLayout(matrix(1:3, ncol=1, byrow=T))

   #col <- plotGetColors(3)
    col <- rep("black", 3)

    figPlotHindcast(assimctx=assimctx, prctx=ipr1, prexp=pr1, meancol=col[1])
    labelPlot("a")

    figPlotHindcast(assimctx=assimctx, prctx=ipr2, prexp=pr2, meancol=col[2])
    labelPlot("b")

    figPlotHindcast(assimctx=assimctx, prctx=ipr3, prexp=pr3, meancol=col[3])
    labelPlot("c")

    if (outfiles) { finDev(display=display) }
}
