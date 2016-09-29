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
#iter <- "2e+06"

#source('plot.R')
source('pdfplot.R')
source('plotutils.R')
source('calib.R')


fnames <- c("uniform", "beta", "normal")
cnames <- capitalize(fnames)


if (!exists("as1")) {
    loadChains(paste(         fnames, "_", iter, ".RData", sep=""))
    loadChains(paste("inst_", fnames, "_", iter, ".RData", sep=""), newnames=c("ias", "ipr"))
}


xlab <- paste("Projected AIS Volume Loss in", year, "[SLE m]")


priorPlot <- function(pr, col="gray", lty="dotted", lwd=2, xlim=par("usr")[1:2], shade=F, n=1001, border=NA)
{
    x <- seq(from=xlim[1], to=xlim[2], length.out=n)
    y <- pr$dens(x, log=F)
    if (shade) {
        polygon(c(x, x[1]), c(y, y[1]), col=col, lty=lty, lwd=lwd, border=border)
    } else {
        lines(x=x, y=y, col=col, lty=lty, lwd=lwd)
    }
}


figCmpPrior <- function()
{
    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    assims <- list(as1, as2, as3)

    for (i in 1:length(chains)) {
        fname <- fnames[i]
        newDev(paste("cmp_prior_", fname, sep=""), outfile=outfiles, width=8.5, height=11/2)

        lwd <- 2
        lty <- c("solid")

        pdfPlots(
            chains=list(chains[[i]]),
            column=as.character(2100),
            lty=lty,
            legendloc=NULL,
            #legends=fname,
            #col="black",
            col="red",
            burnin=F,
            #xlim=c(0, max(cictx$range)),
            #xlim=c(-0.2, 1.1),
            xlab=xlab,
            #yline=2,
            lwd=lwd
            )

        shadecol <- rgb(255, 0, 0, alpha=32, maxColorValue=255)
        assimctx <- assims[[i]]
        pr       <- assimctx$expert_prior
        if (is.null(pr)) {
            pr   <- normPrior(assimctx$obsonly[assimctx$expert_ind], assimctx$windows[assimctx$expert_ind, 2])
        }
        priorPlot(pr, shade=T, border=NA, col=shadecol)

        legend(
            "topright",
            legend=c("inversion", "prior"),
            col=c("red", shadecol),
            lty=c(lty, NA),
            lwd=c(lwd, NA),
            pch=c(NA, 15),
            bg=c(NA, shadecol)
            )

        caption <- paste("Figure n. PDF of AIS volume loss in year", year, "for", fname, "prior")
        mtext(caption, outer=F, line=4, side=1, font=2)

        finDev()
    }
}


figCmpInst <- function()
{
    newDev("cmp_inst", outfile=outfiles, width=8.5, height=11/2)
    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain, ipr1$prchain, ipr2$prchain, ipr3$prchain)

    pdfPlots(
        chains=chains,
        column=as.character(2100),
        lty=c(rep("solid", 3), rep("dashed", 3)),
        legends=c(fnames, paste("inst", fnames)),
        col=rep(c("black", "blue", "red"), 2),
        burnin=F,
        #xlim=c(0, max(cictx$range)),
        #xlim=c(-0.2, 1.1),
        xlab=xlab,
        #yline=2,
        lwd=2
        )

    finDev()
}


figAisPriors <- function()
{
    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))

    newDev("ais_2100_3", outfile=outfiles, width=7, height=5)
    pdfPlots(
        chains=chains,
        column=as.character(2100),
        lty=c("solid", "dashed", "dotted"),  # "dotdash"),
        legends=cnames,
        col=c("black", "blue", "red"),  # , "green"),
        burnin=F,
    #    xlim=c(0, max(cictx$range)),
    #    xlim=c(-0.2, 1.1),
        xlab=xlab,
        lwd=2,
        yline=2
        )

    newDev("ais_2100_2", outfile=outfiles, width=6, height=6)
    pdfPlots(
        chains=list(pr1$prchain, pr2$prchain),
        column=as.character(2100),
        lty=c("solid", "dashed"),
        legends=c("Uniform", "Beta"),
        col=c("black", "blue"),
        burnin=F,
    #    xlim=c(0, max(cictx$range)),
        xlab=xlab,
        lwd=2,
        yline=2
        )

    finDev()
}


finDev()
