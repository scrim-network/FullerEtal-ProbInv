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

source('plot.R')


if (!exists("pr_uniform")) {
    load(paste("uniform_", iter, ".RData", sep=""))
    pr_uniform <- prdaisctx
    as_uniform <-   daisctx

    load(paste("beta_", iter, ".RData", sep=""))
    pr_beta    <- prdaisctx
    as_beta    <-   daisctx

    load(paste("normal_", iter, ".RData", sep=""))
    pr_normal  <- prdaisctx
    as_normal  <-   daisctx
}

fnames <- c("uniform", "beta", "normal")
cnames <- capitalize(fnames)


figAisPriors <- function()
{
    chains <- list(pr_uniform$prchain, pr_beta$prchain, pr_normal$prchain)
    cictx  <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))
    xlab   <- paste("Projected AIS Volume Loss in", year, "[SLE m]")

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
        chains=list(pr_uniform$prchain, pr_beta$prchain),
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
