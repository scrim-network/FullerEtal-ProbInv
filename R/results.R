# Copyright (C) 2017 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# results.R

lhs      <- F
main     <- T
analysis <- T
fast     <- T

# do the LHS run
if (lhs) {
    source('calib.R')
    daisRunLhs()
    source('figures.R')
    figLhs()
}

# load a bunch of runs
if (main | analysis) {
    source('makefig.R')
}

# load the fast dynamics run and make the figure
if (fast) {
    source('roblib.R')  # loadGlobal()

    loadFastDynamics <- function(runbase)
    {
        loaded <- T

        if (file.exists(paste(       runbase, 'n=5e6.RData', sep=""))) {
            loadGlobal(paste(        runbase, 'n=5e6.RData', sep=""))
        } else if (file.exists(paste(runbase, 'n=2e6.RData', sep=""))) {
            loadGlobal(paste(        runbase, 'n=2e6.RData', sep=""))
        } else if (file.exists(paste(runbase, 'n=5e5.RData', sep=""))) {
            loadGlobal(paste(        runbase, 'n=5e5.RData', sep=""))
        } else {
            loaded <- F
        }

        return (loaded)
    }

    loaded     <- loadFastDynamics( 'dp="u";')
    if (!loaded) {
        loaded <- loadFastDynamics('d2p="u";')
    }

    if (loaded) {
        # make the fast dynamics figure
        source('figures.R')
        figDiagFast()
    }
}

# make the other figures
if (main) {
    source('makefig.R')
    figAisPriors()
    figCmpPriors()
    figPredict()
    figInfer()
    figCmpPredict()
}

# table data
if (analysis) {
    source('calib.R')
    print(daisStats( as1))
    print(daisStats(ias1))
    print(daisStats( as2))
    print(daisStats(ias2))
    print(daisStats( as3))
    print(daisStats(ias3))
}
