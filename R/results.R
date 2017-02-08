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

    if (!exists("daisctx") || is.null(daisctx$diagChain)) {
        loaded     <- loadFastDynamics( 'dp="u";')
        if (!loaded) {
            loaded <- loadFastDynamics('d2p="u";')
        }
    } else {
        loaded     <- T
    }

    # make the fast dynamics figure
    if (loaded) {
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

    stats    <- t(sapply(list(as1, ias1, as2, ias2, as3, ias3), daisStats))
    priors   <- c(rep("uniform", 2), rep("beta", 2), rep("normal", 2))
    all_data <- c(F, T, F, T, F, T)

    tab        <- data.frame(priors=priors, all_data=all_data, stats=stats)
    names(tab) <- c("prior", "bayes", "0.05", "0.50", "0.95")
    print(tab)

    file <- "../out/quantiles.csv"
    write.table(tab, file=file, row.names=F, eol="\r\n", sep=", ")  # like a CSV with a space after the comma

   #write.csv(  tab, file=file, row.names=F)
   #write.table(tab, file=file, row.names=F)
}
