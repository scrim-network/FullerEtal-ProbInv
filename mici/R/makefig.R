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


source('pdfplot.R')
source('plotutils.R')
source('roblib.R')  # burnedInd()

bchain <- daisctx$chain[ burnedInd(daisctx$chain), ]


if (0) {
load("uniform")
pr1 <- prdaisctx
load("normal")
pr2 <- prdaisctx
load("beta")
pr3 <- prdaisctx
}


if (1) {

prchain <- cbind( prdaisctx$prchain[, as.character(year) ] )
colnames(prchain) <- paste("Projected AIS Volume Loss in", year, "[SLE m]")
pdfplot(mcmcChain=prchain, units="m", chartRow=1, chartCol=1,
    outfiles=outfiles, burnin=F, width=8, height=8,
    caption=paste("Probability density function of AIS volume loss in year",
        year))

pdfplot(assimctx=daisctx, outfiles=outfiles)
    
newDev("pairs", outfile=outfiles)
plot.pairs(bchain)

#newDev("margin1", outfile=outfiles)
#plot.marginals(bchain[ , 1:9])

#newDev("margin2", outfile=outfiles)
#plot.marginals(bchain[ , 10:ncol(bchain)])

}


finDev()
