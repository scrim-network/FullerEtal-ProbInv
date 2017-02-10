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
# assim.R

source("roblib.R")
source("prior.R")
#source("config.R")
#source("perfect.R")
#source("fastar1.R") # ar1.sim
loadLibrary("mcmc")
loadLibrary("DEoptim")
loadLibrary("adaptMCMC")


# wrapper for the metrop function:  preserves dimensional names
# and separates model and statistical parameters
named_metrop <- function(obj, init_mp, init_sp, nbatch, blen = 1,
    nspac = 1, scale = 1, outfun, debug = FALSE, extrafun = NULL, ...)
{
    initial <- c(init_mp, init_sp)
    mp_indices <- 1:length(init_mp)

    # for recording likelihood as a vector;  see comment below for why this is nbatch+1
    llik <- numeric(length=nbatch+1)
    i    <- 0L

    # wrap obj function in order to assign names to the parameter vector,
    # separate parameters into model and statistical parameters,
    # and record likelihood
    #
    obj2 <- function(p, ...)
    {
        names(p) <- names(initial)
        mp <- p[  mp_indices ]
        sp <- p[ -mp_indices ]

        # record likelihood
        i       <<- i + 1L
        llik[i] <<- obj(mp, sp, ...)

        if (!is.null(extrafun)) {
            extrafun(i, ...)
        }

        return (llik[i])
    }

    out <- metrop(obj2, initial, nbatch, blen, nspac, scale, outfun, debug, ...)

    #print(paste("objective function called", i, "times"))

    # need to discard likelihood from first call to likelihood function;
    # metrop() uses this to get the initial likelihood;  hence, likelihood
    # function is called nbatch+1 times;  however, before we discard it,
    # determine if the first proposal was accepted or rejected
    #
    if (isTRUE(all.equal(out$batch[1, ], initial, check.names=F, check.attributes=F))) {

        # first proposal was rejected;  hence, first recorded likelihood,
        # which will be llik[2] after discarding llik[1], should be the
        # likelihood from the first call to the likelihood function,
        # rather than the likelihood from the second call
        # to the likelihood function, which was rejected
        #
        llik[2] <- llik[1]

        out$keepFirst <- T
    } else {
        out$keepFirst <- F
    }

    # discard likelihood from first call to likelihood function
    llik <- llik[ -1 ]

    # fix up recorded likelihood
    for (i in safefor(2:nbatch)) {
        if (isTRUE(all.equal(out$batch[i - 1, ], out$batch[i, ], check.names=F, check.attributes=F))) {
            llik[i] <- llik[i - 1]
        }
    }

    out$llik <- llik

    # assign names to the chain vector
    colnames(out$batch) <- names(initial)

    return (out)
}


assimFixOutput <- function(assimctx, output, adapt=assimctx$adapt)
{
    # see comments in named_metrop() that explain this craziness
    if (!adapt) {
        outtmp <- output
        if (assimctx$out$keepFirst) {
            outtmp[ 2, ] <- outtmp[ 1, ]
        }
        output <- prmatrix(nrow(outtmp) - 1, xvals=attr(outtmp, "xvals"))
        output[, 1:ncol(outtmp)] <- outtmp[-1, 1:ncol(outtmp)]
    }

    chain <- assimctx$chain
    for (i in safefor(2:nrow(chain))) {
        if (isTRUE(all.equal(chain[ i - 1, ], chain[ i, ], check.names=F, check.attributes=F))) {
            output[ i, ] <- output[ i - 1, ]
        }
    }

    return (output)
}


# wrapper for the MCMC function:  preserves dimensional names
# and separates model and statistical parameters
named_MCMC <- function(p, n, init_mp, init_sp, n.chain=ifelse(adapt, 4, 1), packages=NULL, dyn.libs=NULL,
    scale = rep(1, length(init)), adapt = !is.null(acc.rate), acc.rate = NULL, gamma = 0.5,
    n.start = 0, extrafun = NULL, ...)
{
    initial    <- c(init_mp, init_sp)
    mp_indices <- 1:length(init_mp)
    i          <- 0L

    # wrap obj function in order to assign names to the parameter vector
    # and separate parameters into model and statistical parameters
    #
    p2 <- function(params, ...)
    {
        names(params) <- names(initial)
        mp <- params[  mp_indices ]
        sp <- params[ -mp_indices ]

        llik <- p(mp, sp, ...)
        if (!is.null(extrafun)) {
            i <<- i + 1L
            extrafun(i, ...)
        }

        return (llik)
    }

    if (n.chain != 1) {
        packages <- c(packages, "mvtnorm")
        par.out <- MCMC.parallel(p2, n, initial,
                                 n.chain=n.chain, n.cpu=n.chain, packages=packages, dyn.libs=dyn.libs,
                                 scale=scale, adapt=adapt, acc.rate=acc.rate, gamma=gamma, list=T, ...)
        out     <- par.out[[1]]
        out$par <- par.out
    } else {
        out     <- MCMC         (p2, n, initial,
                                 scale=scale, adapt=adapt, acc.rate=acc.rate, gamma=gamma, list=T, n.start=n.start, ...)
    }
    out$llik   <- out$log.p
    out$batch  <- out$samples
    out$accept <- out$acceptance.rate

    #print(paste("objective function called", i, "times"))

    # assign names to the chain vector
    colnames(out$batch) <- names(initial)

    return (out)
}


# This file also contains a function for estimating a scale matrix
# for an MCMC proposal distribution based on the covariance of 
# preliminary Markov chain
#
# calculate a scale matrix to transform an iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
# J. Comp. Graph. Stat., to appear
#
assimProposalMatrix <- function(prechain,mult=1,beta=0.05)
{
    # mult = overall scale factor to adjust all step sizes
    # beta = relative influence of nugget term

    p = dim(prechain)[2]
    precov = cov(prechain)
    propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p

    mat = t(chol(propcov))
    mat = mult*mat
}


# bounded uniform prior on all parameters
assimLogPriBounds <- function(mp, sp, assimctx)
{
    inBounds <- all(mp >= assimctx$lbound, mp <= assimctx$ubound)
    lpri <- ifelse(inBounds, 0, -Inf)
}


# posterior
assimLogPost <- function(mp, sp, assimctx)
{
    lpri <- assimctx$logPri(mp, sp, assimctx)
    if (!is.finite(lpri)) {

        # save time by avoiding running the model
        return (-Inf)
    }

    lpost <- assimctx$superLogLik(mp, sp, assimctx) + lpri

    #return (ifelse(is.finite(lpost), lpost, -Inf))

    if (!is.finite(lpost)) {
        #recover()
        return (-Inf)
    } else {
        if (lpost > assimctx$maxLik) {
            assimctx$maxLik      <- lpost
            assimctx$maxLikParam <- c(mp, sp)
        }
        return (lpost)
    }
}


# an alternative to configAssim for more complex assimilations
assimInit <- function(assimctx, init_mp, init_sp, lprifn, llikfn)
{
    assimctx$init_mp     <- init_mp
    assimctx$init_sp     <- init_sp
    assimctx$logPri      <- lprifn
    assimctx$superLogLik <- llikfn

    assimctx$maxLik      <- -Inf
    assimctx$mp_indices  <- 1:length(init_mp)
}


# potential function for extrafun in runAssim;
# ychain can be pre-allocated in xxxRunAssim() by caller
assimSaveY <- function(i, assimctx)
{
    assimctx$ychain[i, ] <- assimctx$y
    return ()
}


assimRun <- function(assimctx, nbatch, nspac=1, n.chain=ifelse(adapt, 4, 1), packages=NULL, dyn.libs=NULL,
    scale=NULL, adapt=F, acc.rate = 0.234, gamma=0.5, n.start=0.01*nbatch, extrafun=NULL)
{
    # can verify changes by seeing if they produce the same chain
    #set.seed(7)

    if (is.null(scale)) {
        range_mp <- assimctx$ubound    - assimctx$lbound
        range_sp <- assimctx$ubound_sp - assimctx$lbound_sp
        scale <- abs(c(range_mp / 100, range_sp / 1000))

       #if (assimctx$ar == 2) {
            #scale["rho2"] <- scale["rho2"] / 2
       #}

        print("using initial scale:")
        print(scale)
    }

    if (adapt) {
        if (!missing(nspac)) {
            print("WARNING! adaptive MCMC ignores nspac value")
        }
        time <- system.time(
            out <- named_MCMC(p=assimLogPost, n=nbatch, init_mp=assimctx$init_mp, init_sp=assimctx$init_sp,
                n.chain=n.chain, packages=packages, dyn.libs=dyn.libs,
                scale=scale, adapt=adapt, acc.rate=acc.rate, gamma=gamma, n.start=n.start,
                extrafun=extrafun, assimctx=assimctx
                )
        )
        out$time <- time
    } else {
        if (n.chain != 1) {
            print("WARNING! non-adaptive MCMC ignores n.chain value")
        }
        out <- named_metrop(obj=assimLogPost, init_mp=assimctx$init_mp, init_sp=assimctx$init_sp,
            nbatch=nbatch, nspac=nspac, scale=scale,
            extrafun=extrafun, assimctx=assimctx
            )
    }

    print(out$time)

    # want an acceptance rate of ~25%
    print(paste("accept=", out$accept, sep=""))

    # for posterity
    assimctx$out   <- out

    # many functions outside assim.R use this
    assimctx$chain <- out$batch

    assimctx$adapt <- adapt
}


# an old name for compatibility with old code
runAssim <- assimRun


assimMinFn <- function(p, assimctx)
{
    mp <- p[  assimctx$mp_indices ]
    sp <- p[ -assimctx$mp_indices ]

    lpost <- assimLogPost(mp, sp, assimctx)
    return (-lpost)
}


# for AR(0), this should give the same results as minimizing SSE
assimMaxLikelihood <- function(assimctx, init_mp=NULL, init_sp=NULL, useDE=T, itermax=500)
{
    lbound <- c(assimctx$lbound, assimctx$lbound_sp)
    ubound <- c(assimctx$ubound, assimctx$ubound_sp)

    init_p <- c(init_mp, init_sp)
    if (is.null(init_p)) {
        init_p <- (lbound + ubound) / 2
    }

    if (useDE) {
        control <- list(CR=1.0, NP=10*length(lbound), itermax=itermax)
        fit <- named_DEoptim(
            FUN=assimMinFn,
            lower=lbound, upper=ubound,
            control=control,
            assimctx=assimctx
            )
        best_p <- fit$optim$bestmem
    } else {
        control=list()
        if (!missing(itermax)) {
            control["maxit"] <- itermax
        }
        fit <- optim(
            par=init_p,
            fn=assimMinFn,
            #method="L-BFGS-B",
            #lower=lbound, upper=ubound,
            control=control,
            assimctx=assimctx
            )
        best_p <- fit$par
    }

    return (best_p)
}


assimGelman <- function(assimctx, from=50000, by=10000, to)
{
    if (missing(to)) {
        to <- nrow(assimctx$out$par[[1]]$samples)
    }
    len    <- seq(from=from, to=to, by=by)
    stats  <- numeric(length=length(len))

    progress <- 0
    bar      <- txtProgressBar(min=len[1], max=sum(len), style=3)
    for (i in 1:length(len)) {
        chains <- list()
        for (j in 1:length(assimctx$out$par)) {
            chains[[j]] <- as.mcmc(assimctx$out$par[[j]]$samples[1:len[i], ])
        }
        mpsrf    <- gelman.diag(chains)[[2]]
        stats[i] <- mpsrf
        progress <- progress + len[i]
        setTxtProgressBar(bar, progress)
    }
    close(bar)

    assimctx$gr_len   <- len
    assimctx$gr_stats <- stats

    plot(len, stats)
}
