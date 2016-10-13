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


figPlotBounds <- function(assimctx=as1, lwd=1.5)
{
    abline(v=assimctx$windows[assimctx$expert_ind, ], lty="dotted", lwd=lwd)
}


figColors <- function(n=3, alpha=255)
{
    # display.brewer.all(type="qual", colorblindFriendly=T)
   #pal <- brewer.pal(n=n, "Dark2"))
   #pal <- brewer.pal(n=n, "Paired")
    pal <- brewer.pal(n=n, "Set2")
    pal <- paste(pal, as.hexmode(alpha), sep="")
   #print(pal)

    return (pal)
}


figCmpPriors <- function()
{
    newDev("cmp_prior", width=8, height=8, outfile=outfiles, filetype=filetype)

    par(mfrow=c(2, 2))
    par(mar=c(4, 3, 0, 3))

    labels   <- c("a", "b", "c")
    col      <- figColors(3)
    shadecol <- figColors(3, 48)

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

    finDev()
}


figCmpInst <- function()
{
    newDev("cmp_inst", outfile=outfiles, width=8.5, height=11/2, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain, ipr1$prchain, ipr2$prchain, ipr3$prchain)
    lty <- c(rep("solid", 3), rep("dashed", 3))
    col <- rep(figColors(3), 2)
    lwd <- 2

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
    figPlotBounds()
    legend(
        "topright",
        legend=c(fnames, paste("inst", fnames), "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(rep(lwd, 6), 1.5),
        cex=0.75
        )

    finDev()
}


figAisPriors <- function()
{
    newDev("ais_2100_3", outfile=outfiles, width=7, height=5, filetype=filetype)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    col <- figColors(3)
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
    figPlotBounds()
    legend(
        "topright",
        legend=c(cnames, "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(lwd, 1.5),
        cex=0.75
        )

    caption <- paste("Figure 1. Deep uncertainty")
    mtext(caption, outer=F, line=4, side=1, font=2)

    finDev()
}


figPredict <- function()
{
    newDev("pred_2100", outfile=outfiles, width=8, height=8, filetype=filetype)

    par(mfrow=c(2, 1))
    #par(mar=c(4, 3, 0, 3))

    chains <- list(ipr1$prchain, ipr2$prchain, ipr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    col <- figColors(3)
    #lty <- c("solid", "dashed", "dotted")  # "dotdash"
    lty <- rep("solid", 3)
    lwd <- 2

    column=as.character(2100)

    pdfPlots(
        chains=chains,
        column=column,
        lty=lty,
        col=col,
        burnin=F,
    #    xlim=c(0, max(cictx$range)),
    #    xlim=c(-0.2, 1.1),
    #    xlab=xlab,
        lwd=lwd,
        legendloc=NULL,
        smoothing=c(0.50, rep(1.25, 2)),
        yline=2
        )
    figPlotBounds()
    legend(
        "topright",
        legend=c(cnames, "Pfeffer"),
        col=c(col, "black"),
        lty=c(lty, "dotted"),
        lwd=c(lwd, 1.5),
        cex=0.75
        )

    xlim <- par("usr")[1:2]
    cdfPlots(
        chains=chains,
        column=column,
        xlim=xlim,
        lwd=lwd, col=col, lty=lty
        )

    caption <- paste("Figure 4. Add Paleo and Instrumental Observations")
    mtext(caption, outer=F, line=4, side=1, font=2)

    finDev()
}


finDev()
