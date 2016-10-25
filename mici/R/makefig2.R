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
# makefig2.R
#

outfiles <- F
year <- 2100


source('pdfplot.R')
source('plotutils.R')


#fnames <- c("uniform", "beta", "normal")
fnames <- "uniform"

for (i in 1:length(fnames)) {
    fname <- fnames[i]
    #load(paste(fname, ".RData", sep=""))

    bchain <- daisctx$chain[ burnedInd(daisctx$chain), ]

    newDev(paste(fname, "_pairs", sep=""), outfile=outfiles)
    plot.pairs(bchain)
    mtext(paste("Pairs for", fname, "prior"), outer=TRUE, side=1, font=2)

    pdfplot(assimctx=daisctx, caption=paste("Posterior PDFs for", fname, "prior"), outfiles=outfiles)
    file.rename("../figures/pdfs1.pdf", paste(sep="", "../figures/", fname, "_pdfs1.pdf"))
    file.rename("../figures/pdfs2.pdf", paste(sep="", "../figures/", fname, "_pdfs2.pdf"))

    newDev(paste(fname, "_margin1", sep=""), outfile=outfiles)
    plot.marginals(bchain[ , 1:9])
    mtext(paste("Figure 1. Posterior PDFs for", fname, "prior"), outer=TRUE, side=1, font=2)

    newDev(paste(fname, "_margin2", sep=""), outfile=outfiles)
    plot.marginals(bchain[ , 10:ncol(bchain)])
    mtext(paste("Figure 2. Posterior PDFs for", fname, "prior"), outer=TRUE, side=1, font=2)

    newDev(paste(fname, "_ais", sep=""), outfile=outfiles, width=6, height=6)
    pdfPlots(
        chains=list(prdaisctx$prchain),
        column=as.character(2100),
        lty="solid",
        legends=fname,
        col="black",
        burnin=T,
    #    xlim=c(0, max(cictx$range)),
        legendloc=NULL,
        xlab=paste("Projected AIS Volume Loss in", year, "[SLE m]"),
    #    yline=2,
        lwd=2
        )
}


finDev()
