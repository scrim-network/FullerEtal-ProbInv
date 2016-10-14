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
#

outfiles <- T
year <- 2100
iter <- "5e+05"
filetype <- "png"
#iter <- "2e+06"

source('plot.R')
loadLibrary('RColorBrewer')


fnames <- c("uniform", "beta", "normal")
cnames <- capitalize(fnames)


if (!exists("pr1")) {
    loadChains(paste(         fnames, "_", iter, ".RData", sep=""))
    loadChains(paste("inst_", fnames, "_", iter, ".RData", sep=""), newnames=c("ias", "ipr"))
}

pr1$assimctx <- as1
pr2$assimctx <- as2
pr3$assimctx <- as3
ipr1$assimctx <- ias1
ipr2$assimctx <- ias2
ipr3$assimctx <- ias3


xlab <- paste("Projected AIS Volume Loss in", year, "[SLE m]")


plotBounds <- function(assimctx=as1, lwd=1.5)
{
    abline(v=assimctx$windows[assimctx$expert_ind, ], lty="dotted", lwd=lwd)
}


getColors <- function(n=3, alpha=255)
{
    # display.brewer.all(type="qual", colorblindFriendly=T)
   #pal <- brewer.pal(n=n, "Dark2"))
   #pal <- brewer.pal(n=n, "Paired")
    pal <- brewer.pal(n=n, "Set2")
    pal <- paste(pal, as.hexmode(alpha), sep="")
   #print(pal)

    return (pal)
}


figAisPriors <- function()
{
    newDev("fig1", outfile=outfiles, width=8.5, height=11/2, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    col <- getColors(3)
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

    caption <- paste("Figure 1. Deep uncertainty")
    mtext(caption, outer=F, line=4, side=1, font=2)

    if (outfiles) { finDev() }
}


figCmpPriors <- function()
{
    newDev("fig2", width=8.5, height=7, outfile=outfiles, filetype=filetype)

    par(mfrow=c(2, 2))
    par(mar=c(4, 3, 0, 3))

    labels   <- c("a", "b", "c")
    col      <- getColors(3)
    shadecol <- getColors(3, 48)

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
        if (assimctx$prior_name == "uniform") {
            smoothing <- 0.50
        } else {
            smoothing <- 1.25
        }
        priorPdfPlot(prctx$prchain, "2100", prior=pr, xlim=xlim, xlab=xlab, col=col[i], shadecol=shadecol[i], legends=c(fnames[i], "inversion"), smoothing=smoothing)
        labelPlot(labels[i])
    }

    # figure title
    caption <- paste("Figure 2. Probabilistic inversion works")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


figPredict <- function()
{
    newDev("fig4", outfile=outfiles, width=8.5, height=7, filetype=filetype)

    # reserve lines to use outer=T for the lower axis and label;
    # allows using 0.5 in par(fig) and getting equally sized plots;
    # one line is 0.2 inches:  par("csi") or par("cin")[2]
    #
    #par(oma=c(4, 1, 1, 1))
    # bottom, left, top, right
    par(omi=c(1.00, 0.25, 0.25, 0.25))


    # parameters common to CDF/PDF
    #

    chains <- list(ipr1$prchain, ipr2$prchain, ipr3$prchain)
    column <- as.character(2100)
   #cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    col <- getColors(3)
    lty <- rep("solid", 3)
    lwd <- 2


    # top of figure (PDFs)
    #

    par(fig=c(0, 1, 0.5, 1))
    par(mar=c(0, 3, 0, 1))
    plot.new()

    pdfctx <- pdfCalc(chains=chains, column=column, burnin=F, smoothing=c(0.50, rep(1.25, 2)))

    #    xlim=c(0, max(cictx$range)),
    plot.window(xlim=pdfctx$xlim, ylim=pdfctx$ylim, xaxs="i")

    # bottom
    axis(1, labels=F, tcl=-0.10) # bottom

    # left
    ticks <- axTicks(2)
    ticks <- c(0, last(ticks))
    axis(2, at=ticks)

    # top:  positive values for tcl put the tickmarks inside the plot
    axis(3, labels=F, tcl=-0.10)

    # right
    axis(4, at=ticks, labels=F, tcl=-0.25)

    title(ylab="Probability density", line=2)
    box()

    plotBounds()
    pdfPlot(pdfctx, col=col, lty=lty, lwd=lwd)
    legend(
        "topleft",
        legend=c(cnames, "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(rep(lwd, 3), 1.5),
        cex=0.75
        )
    labelPlot("a")



    # bottom of figure (CDFs)
    #

    par(fig=c(0, 1, 0, 0.5), new=T)
    par(mar=c(0, 3, 0, 1))

    cdfPlots(
        chains=chains,
        column=column,
        xlim=pdfctx$xlim,
        lwd=lwd, col=col, lty=lty,
        ylab="Survival (1-CDF)",
        survival=T
        )
    plotBounds()
    labelPlot("b")

    title(xlab=xlab, line=2, outer=T)


    # figure title
    caption <- paste("Figure 4. Add paleo and instrumental observations")
    mtext(caption, outer=TRUE, side=1, font=2, line=4)

    if (outfiles) { finDev() }
}


figUber <- function()
{
    newDev("uber", outfile=outfiles, width=8.5, height=8.5, filetype=filetype)

    chains <- list(as1$chain, as2$chain, as3$chain)
    lty    <- rep("solid", 3)
    lwd    <- 2
    col    <- getColors(3)


    # bottom, left, top, right
    par(mar = c(0.25, 5, 1, 0))
    layout(matrix(1:4, nrow = 2, byrow = T), widths = c(10, 3), heights = c(3, 10))


    # top PDF
    #

    plot.new()

    topctx <- pdfCalc(chains=chains, column="Tcrit", burnin=T)
    plot.window(xlim=topctx$xlim, ylim=topctx$ylim, xaxs="i")

    # left
    ticks <- axTicks(2)
    ticks <- c(0, last(ticks))
    axis(2, at=ticks)

    title(ylab="Probability density", line=2)

    pdfPlot(topctx, col=col, lty=lty, lwd=lwd)


    # legend
    #

    par(mar = c(0.25, 0.25, 0, 0))
    plot.new()
    legend(
        "center",
        legend=cnames,
        title="Prior",
        title.adj = 0.1,
        bty = "n",
        fill = col,
        border = col
        )


    # main plot
    #

    par(mar = c(4, 5, 0, 0))
    plot.new()
    #axis(3, labels = F, tck = 0.01)
    #axis(4, labels = F, tck = 0.01)
    box()


    # right PDF
    #

    par(mar = c(4, 0.25, 0, 1))
    plot.new()

    bottctx <- pdfCalc(chains=chains, column="lambda", burnin=T) # , smoothing=c(0.50, rep(1.25, 2)))

    plot.window(ylim=bottctx$xlim, xlim=bottctx$ylim, xaxs="i")

    # bottom
    ticks <- axTicks(1)
    ticks <- c(0, last(ticks))
    axis(1, at=ticks)

    title(xlab="Probability density", line=2)
    pdfPlot(bottctx, col=col, lty=lty, lwd=lwd, reverse=T)


    if (outfiles) { finDev() }
}


figCmpInst <- function()
{
    newDev("cmp_inst", outfile=outfiles, width=8.5, height=7, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain, ipr1$prchain, ipr2$prchain, ipr3$prchain)
    lty <- c(rep("dotted", 3), rep("solid", 3))
    col <- rep(getColors(3), 2)
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
        legend=c(fnames, paste("inst", fnames), "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(rep(lwd, 6), 1.5),
        cex=0.75
        )

    caption <- paste("Figure n. Add paleo and instrumental observations")
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


if (outfiles) { finDev() }
