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
outfiles <- F

source('pdfplot.R')
source('plotutils.R')
source('roblib.R')  # burnedInd()

bchain <- daisctx$chain[ burnedInd(daisctx$chain), ]

pdfplot(assimctx=daisctx, outfiles=outfiles)
    
newDev("pairs", outfile=outfiles)
plot.pairs(bchain)

#newDev("margin1", outfile=outfiles)
#plot.marginals(bchain[ , 1:9])

#newDev("margin2", outfile=outfiles)
#plot.marginals(bchain[ , 10:ncol(bchain)])

finDev()
