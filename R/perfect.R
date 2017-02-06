# Copyright (C) 2009, 2010, 2016, 2017 Robert W. Fuller
# email: hydrologiccycle@gmail.com
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
# perfect.R


# truth can be "mean", "max", or "mode"
truth <- function(assimctx, truth="mode", uniform=F)
{
    # create true parameters
    #
    mcmcChain <- assimctx$chain[ burnedInd(assimctx$chain), ]
    true_param <- switch(truth,
        mean=colMeans(mcmcChain),

        # TODO:  should use max of out$llik
        max=assimctx$maxLikParam,
        mode=colMode(mcmcChain)
        )

    true_mp <- true_param[  assimctx$mp_indices ]
    true_sp <- true_param[ -assimctx$mp_indices ]

    if (uniform) {

        true_mp <- (assimctx$lbound + assimctx$ubound) / 2

        # TODO:  handle log priors differently?
        #true_mp["tau"] <- 7.296413
    }

    return (list(mp=true_mp, sp=true_sp))
}


configPerfectModel <- function(assimctx, truth="mode", uniform=F, noise=T)
{
    true_p  <- truth(assimctx, truth, uniform)
    true_mp <- true_p$mp
    true_sp <- true_p$sp

    # create observations from true parameters
    # TODO:  need to set obs too?  problematic: model function truncates
    #
    assimctx$obsonly <- assimctx$modelfn(true_mp, assimctx)

    # add a noise realization to the observations;
    #
    if (noise) {
        # use c() to strip time series attribute from arima.sim()
        realiz <- c(assimctx$noise(true_sp, length(assimctx$obsonly), assimctx))
        assimctx$obsonly <- assimctx$obsonly + realiz
    }

    # start assimilation with true parameters AND save true parameters
    #
    assimctx$init_mp <- true_mp
    assimctx$init_sp <- true_sp

    print("true parameters:")
    print(c(true_mp, true_sp))

    assimctx$maxLik  <- -Inf
}


runUniformAssim <- function(nbatch=1000000, assimctx)
{
    nparms <- length(assimctx$mp_indices)

    # could use the names from init_mp here instead of lbound;
    # also, mp_indices is not strictly necessary here since statistical
    # parameters don't (currently) have bounds specified this way
    #
    names  <- names(assimctx$lbound)[assimctx$mp_indices]

    chain  <- matrix(nrow=nbatch, ncol=nparms)
    colnames(chain) <- names
    names <- paste(names, "_prior", sep="")

    for (i in 1:nparms) {

        mp_index <- assimctx$mp_indices[i]
        min=assimctx$lbound[mp_index]
        max=assimctx$ubound[mp_index]

        if (exists(names[i], envir=assimctx)) {
            prior <- get(names[i], envir=assimctx)

            # could instantiate a uniformPrior object in the alternate case
            # whereupon this would become common code;  however, that would
            # mean unnecessarily checking the range on uniform chains which
            # penalizes the common case;  this is merely a supposition,
            # until you time it:
            #
            # system.time( makeUniformChain(assimctx=allgrgisctx,nbatch=10000000) )
            #   user  system elapsed 
            # 20.040   2.638  25.870 
            #
            # system.time( makeUniformChain(assimctx=allgrgisctx,nbatch=10000000) )
            #   user  system elapsed 
            # 34.459   4.813  39.686 
            #
            chain[, i] <- prior$rand(n=nbatch)
            while (length(redo <- which(chain[, i] < min | chain[, i] > max))) {
                #print(length(redo))
                #print(chain[redo, i])
                chain[redo, i] <- prior$rand(n=length(redo))
            }
        } else {
            chain[, i] <- runif(n=nbatch, min=min, max=max)
        }
    }

    assimctx$unichain <- chain

    if (nbatch <= 20) {
        print(chain)
    }
}


runUniformPredict <- function(
    nbatch=100000, year=2200,
    assimctx,
    prctx,
    forcings,
    modelfn,
    noisefn=rep(list(noise_zeros), length(outnames)),
    outnames=c("prunichain"),
    mcmcChain=assimctx$unichain,
    xvals=assimctx$times[1]:year,
    ...
    )
{
    # yes, changing default parameters is a horror;  tried to do this
    # by computing on the language;  broke R, badly
    #
    runPredict(
        nbatch=nbatch, year=year,
        assimctx=assimctx,
        prctx=prctx,
        forcings=forcings,
        modelfn=modelfn,
        noisefn=noisefn,
        outnames=outnames,
        mcmcChain=mcmcChain,
        xvals=xvals,
        ...
        )
}
