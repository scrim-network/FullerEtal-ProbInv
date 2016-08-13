# Copyright 2009, 2010, 2016 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# fooassim.R

source("assim.R")


fooModel <- function(mp, assimctx)
{
    y <- mp["m"] * assimctx$x + mp["b"]

    return (y)
}


if (!exists("fooassimctx")) {
    fooassimctx <- env()
}


fooConfigAssim <- function(
    nobs=10,
    sd=0.05,
    mtparms=0,
    assimctx=fooassimctx
    )
{
    mp=c(m=2, b=5)

    # fooModel() uses this
    assimctx$x   <- 1:nobs

    assimctx$obsonly <- fooModel(mp, assimctx) + rnorm(length(assimctx$x), sd=sd)
    assimctx$modelfn <- fooModel
    assimctx$lbound  <- c(0,  0)
    assimctx$ubound  <- c(5, 10)

    if (mtparms) {
        lbound   <- rep(0.0, mtparms)
        ubound   <- rep(1.0, mtparms)
        newmp    <- (lbound + ubound) / 2
        names(newmp) <- paste("mtparms", 1:mtparms, sep="")

        assimctx$lbound <- c(assimctx$lbound, lbound)
        assimctx$ubound <- c(assimctx$ubound, ubound)
        mp <- c(mp, newmp)
    }

    # get initial conditions from best fit model
    configAssim(assimctx, mp, ar=0, obserr=F)
}


fooRunAssim <- function(
    nbatch=1000000,
    initial=is.null(assimctx$chain),
    assimctx=fooassimctx
    )
{
    init_mp <- assimctx$init_mp
    init_sp <- assimctx$init_sp

    nobs <- length(assimctx$obsonly)

    if (initial || ncol(assimctx$chain) != length(init_mp) + length(init_sp)) {
        print("using initial scale")

        # this isn't perfect, but works reasonably for 10-100 observations
        scale <- c(rep(2 * nobs, length(init_mp) + length(init_sp)))

        scale <- abs(c(init_mp, init_sp)) / scale
    } else {
        print("using proposal matrix")

        mult  <- 0.5
        scale <- assimProposalMatrix(assimctx$chain, mult=mult)
    }

    runAssim(assimctx, nbatch=nbatch, scale=scale)
}
