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
# lhs.R


#source("roblib.R")  # presume assim.R has already been sourced
loadLibrary("lhs")


# potential function for extrafun in assimRunLhs;
# ychain can be pre-allocated in xxxRunLhs() by caller
assimLhsSaveY <- function(i, assimctx)
{
    assimctx$lhs_ychain[i, ] <- assimctx$y
    return ()
}


assimRunLhs <- function(assimctx, nbatch=1e3, extrafun=NULL)
{
    lhsctx <- env()

    nvar  <- length(assimctx$lbound)
    chain <- randomLHS(n=nbatch, k=nvar)
    for (i in safefor(1:nvar)) {
        chain[, i] <- qunif(chain[, i], assimctx$lbound[i], assimctx$ubound[i])
    }
    colnames(chain) <- names(assimctx$lbound)

    mp_indices <- 1:length(assimctx$init_mp)
    bar <- txtProgressBar(min=1, max=nbatch, style=3)
    llik <- numeric(length=nbatch)

    # TODO:  save/restore prior values?
    assimctx$maxLik <- -Inf

    for (i in safefor(1:nbatch)) {
        params <- chain[i, ]
        mp <- params[  mp_indices ]
        sp <- params[ -mp_indices ]

        llik[i] <- logPost(mp, sp, assimctx)
        if (!is.null(extrafun)) {
            extrafun(i, assimctx)
        }
        setTxtProgressBar(bar, i)
    }

    close(bar)

    lhsctx$chain       <- chain
    lhsctx$llik        <- llik
    lhsctx$maxLik      <- assimctx$maxLik
    lhsctx$maxLikParam <- assimctx$maxLikParam

    # TODO:  save/restore prior values?
    assimctx$maxLik <- -Inf

    return (lhsctx)
}
