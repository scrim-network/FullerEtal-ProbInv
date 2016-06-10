# Copyright 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# analysis.R

quick <- F
outfiles <- T
year <- 2100

source('pola.R')
source('predict.R')
source('pdfplot.R')
source('plotutils.R')


# need this in part because differential evolution likes to pick out
# rho=0.99 when finding the best fit (MLE)
#
set.seed(9)

polaConfigAssim()
polaRunAssim(nbatch=100000)

if (quick) {
    # a million is better, but 100,000 runs faster
    polaRunAssim(nbatch=100000)
    bchain <- polassimctx$chain[ burnedInd(polassimctx$chain), ]
} else {
    polaRunAssim(nbatch=1200000)
    bchain <- polassimctx$chain[ 200001:1200000, ]
}

if (quick) {
    polaRunPredict(year=year)
} else {
    polaRunPredict(nbatch=1000000, year=year)
}


# old-style plotting routines (graphics.off)

polaPlotPredict(outfiles=outfiles)

prchain <- cbind( prpolactx$prchain[, as.character(year) ] )
colnames(prchain) <- paste("Sea-level anomaly in year", year)
pdfplot(mcmcChain=prchain, units="mm", chartRow=1, chartCol=1,
    outfiles=outfiles, burnin=F, width=8, height=8,
    caption=paste("Probability density function of sea-level anomaly in year",
        year))
if (outfiles) {
    file.rename("../figures/pdfs1.pdf",
        paste("../figures/slr_", year, ".pdf", sep=""))
}

pdfplot(assimctx=polassimctx, outfiles=outfiles)


# new-style plotting routines (finDev)

polaPlotFit(outfiles=outfiles)

newDev("pairs", outfile=outfiles)
plot.pairs(bchain)

finDev()


# output chains

if (!quick) {
    write.table(bchain, file="../out/assim_chain.txt", row.names=F)

    #write.table(prpolactx$prchain, file="../out/pred_chain.txt", row.names=F)
    write.table(prpolactx$prchain[, as.character(year)],
        file="../out/pred_chain.txt",
        row.names=F, col.names="sea-level anomaly (mm)")
}
