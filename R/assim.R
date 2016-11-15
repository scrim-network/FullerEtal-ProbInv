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

source("roblib.R")
#source("fastar1.R") # ar1.sim
loadLibrary("mcmc")
loadLibrary("mvtnorm")
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
named_MCMC <- function(p, n, init_mp, init_sp, scale = rep(1, length(init)),
    adapt = !is.null(acc.rate), acc.rate = NULL, gamma = 0.5,
    list = TRUE, n.start = 0, extrafun = NULL, ...)
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

    out <- MCMC(p2, n, initial, scale, adapt, acc.rate, gamma, list, n.start, ...)
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


# AR(1) covariance
covar_ar1 <- function(N, sigma, rho1)
{
    # Instead of the outer, you could also do
    #   toeplitz(rho1^(0:(N-1)))
    # which gives the same matrix
    #
    # sigma^2/(1-rho1^2) estimates the process variance
    #   if it doesn't converge, try estimating the innovation variance
    #   using sigma^2 in place of sigma^2/(1-rho1^2)
    #
    covar = sigma^2/(1-rho1^2) * outer(1:N, 1:N, function(i,j) rho1^abs(i-j))
    #covar = sigma^2 * outer(1:N, 1:N, function(i,j) rho1^abs(i-j))
}


# AR(2) covariance from Yule-Walker recurrence (cf. von Storch & Zwiers (2001), Section 11.1.7)
covar_ar2 <- function(N, sigma, rho1, rho2)
{
    rho = rep(NA, N)
    rho[1] = 1; rho[2] = rho1/(1-rho2)
    for (i in 3:N) {
        rho[i] = rho1*rho[i-1] + rho2*rho[i-2]    
    }

    # might want to use something other than sigma^2 as in covar_ar1?
    # would have to derive it for AR(2)
    Sigma = sigma^2 * toeplitz(rho)
}


# covariance = AR(1) + independent Gaussian observational noise
covar_ar1_obs <- function(stderror, sigma, rho1, diagerr=diag(stderror^2))
{
    N = nrow(diagerr)
    Sigma = covar_ar1(N, sigma, rho1) + diagerr
}


# covariance = AR(2) + independent Gaussian observational noise
covar_ar2_obs <- function(stderror, sigma, rho1, rho2, diagerr=diag(stderror^2))
{
    N = nrow(diagerr)
    Sigma = covar_ar2(N, sigma, rho1, rho2) + diagerr
}


# if using a gamma prior, then the chain contains variance rather than sigma
sp_sigma <- function(sp)
{
    sigma <- sp["sigma"]
    if (is.na(sigma)) {
        sigma <- sqrt(sp["var"])
    }

    return (sigma)
}


# bounded uniform prior on all parameters except sigma
lpri_sigma <- function(mp, sp, assimctx)
{
    # believe it or not, this works for runassim(ar=0),
    # when there are no rhos
    rhos = pindex("rho", sp)

    rhosum = sum(rhos)

    # check bounds
    # rwf -- changed rho checks from [0, 1) to (-1, 1)
    # rwf -- Klaus says to use -0.99 to 0.99
    # rwf -- change rho checks to reflect initial values
    #
    #rhos = append(rhos, sum(rhos))

    inBounds = all(
        rhosum >= -assimctx$rholimit, rhosum <= assimctx$rholimit,
        mp >= assimctx$lbound,    mp <= assimctx$ubound,
        sp >= assimctx$lbound_sp, sp <= assimctx$ubound_sp

        # tried adding this for DEoptim, but DEoptim simply cannot handle
        # Inf as a bound
        #, na.rm=T 
        )
    if (inBounds) {  # add priors for model parameters if non-uniform
        sigma = sp["sigma"]
        if (is.na(sigma)) {

            # inverse gamma prior for variance
            var   = sp["var"]
            alpha = assimctx$alpha
            beta  = assimctx$beta
            lpri  = (-alpha - 1)*log(var) + (-beta/var)

        } else {

            # 1/sigma^2 Jeffreys prior for variance
            lpri = -2 * log(sigma)
        }
    } else {
        lpri = -Inf  # zero prior probability
    }

    return (lpri)
}


# bounded uniform prior on all parameters
lpri_bounds <- function(mp, sp, assimctx)
{
    inBounds <- all(mp >= assimctx$lbound, mp <= assimctx$ubound)
    lpri <- ifelse(inBounds, 0, -Inf)
}


logLik <- function(mp, sp, assimctx)
{
    # dmvnorm is sensitive to the order here when used purely with residuals
    res <- assimctx$obsonly - assimctx$modelfn(mp, assimctx)

    assimctx$logLik(res, sp, assimctx)
}


# posterior
logPost <- function(mp, sp, assimctx)
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


# pure Gaussian likelihood
llik_norm <- function(res, sp, assimctx)
{
    llik <- sum(dnorm(res, sd=sp_sigma(sp), log=TRUE))  # iid normal distribution
}


# AR(n) likelihood (approximate; no initial value correction)
llik_ar <- function(res, sp, assimctx)
{
    rhos <- pindex("rho", sp)
    res  <- arwhiten(res, rhos)
    llik <- sum(dnorm(res, sd=sp_sigma(sp), log=TRUE))  # iid normal distribution
}


llik_obs <- function(res, sp, assimctx)
{
    llik <- sum(dnorm(res, sd=assimctx$error, log=T))
}


llik_ar1_obs <- function(res, sp, assimctx)
{
    Sigma <- covar_ar1_obs( , sp_sigma(sp), sp["rho1"], assimctx$diag)
    llik  <- dmvnorm(res, sigma=Sigma, log=TRUE)
}


llik_ar2_obs <- function(res, sp, assimctx)
{
    Sigma <- covar_ar2_obs( , sp_sigma(sp), sp["rho1"], sp["rho2"], assimctx$diag)
    llik  <- dmvnorm(res, sigma=Sigma, log=TRUE)
}


noise_norm <- function(sp, N, assimctx)
{
    noise <- rnorm(n=N, sd=sp_sigma(sp))
}


noise_ar <- function(sp, N, assimctx)
{
    rhos  <- pindex("rho", sp)
    noise <- arima.sim(n=N, sd=sp_sigma(sp), list( ar=rhos ))
}


noise_ar1 <- function(sp, N, assimctx)
{
    ar1.sim(N, sp["rho1"], sp_sigma(sp))
}


noise_obs <- function(sp, N, assimctx)
{
    errlen <- length(assimctx$error)
    noise <- rnorm(n=errlen, sd=assimctx$error)
    noise <- append(noise, rep(0, N - errlen))
}


noise_ar1_obs <- function(sp, N, assimctx)
{
    squares <- c(assimctx$squares, rep(0, N - length(assimctx$squares)))
    noise <- rmvnorm(1, sigma=covar_ar1_obs( , sp_sigma(sp), sp["rho1"], diag(squares)))
}


noise_ar2_obs <- function(sp, N, assimctx)
{
    squares <- c(assimctx$squares, rep(0, N - length(assimctx$squares)))
    noise <- rmvnorm(1, sigma=covar_ar2_obs( , sp_sigma(sp), sp["rho1"], sp["rho2"], diag(squares)))
}


sampleNoiseAssim <- function(assimctx, N, nbatch=100000)
{
    mcmcChain <- assimctx$chain
    nChain    <- sample(burnedInd(mcmcChain), nbatch, replace=T)

    noise <- matrix(nrow=nbatch, ncol=N)

    for (i in 1:nbatch) {
        noise[i, ] <- assimctx$noise(mcmcChain[ nChain[i], ][ -assimctx$mp_indices ], N, assimctx)
    }

    return (noise)
}


setErrorAssim <- function(assimctx, error)
{
    assimctx$error   <- error
    assimctx$squares <- error ^ 2
    assimctx$diag    <- diag(assimctx$squares)

    return (assimctx)
}


configAssim <- function(
    assimctx,
    init_mp=NULL, init_sp=NULL,
    ar=0, obserr=T, llikfn=logLik,
    itermax=500,

    # TODO:  DEoptim cannot handle Inf for a boundary
    inv_gamma_pri=F, alpha=2, beta=1, var_max=10.0,
    sigma_max=0.1,
    fixrho=F, rholimit=0.99
    )
{
    # almost all assimilations need standard deviation in the chain
    sigma <- T

    # select likelihood and noise functions based on ar and obserr parms
    if (obserr) {

        if (ar == 0) {
            assimctx$logLik <-  llik_obs
            assimctx$noise  <- noise_obs
            sigma <- F
        } else if (ar == 1) {
            assimctx$logLik <-  llik_ar1_obs
            assimctx$noise  <- noise_ar1_obs
        } else if (ar == 2) {
            assimctx$logLik <-  llik_ar2_obs
            assimctx$noise  <- noise_ar2_obs
        } else {
            stop("assimilation with observation error only implemented for AR(0:2)")
        }

    } else if (ar) {
        if (ar == 1) {
            assimctx$noise  <- noise_ar

            # TODO:  this is faster, but doesn't check for stationarity
            #assimctx$noise <- noise_ar1
        } else {
            assimctx$noise  <- noise_ar
        }
        assimctx$logLik     <-  llik_ar
    } else {
        assimctx$logLik     <-  llik_norm
        assimctx$noise      <- noise_norm
    }

    lbound_sp <- numeric()
    ubound_sp <- numeric()

    # N.B. order is relevant:  rho comes before sigma here and below
    #
    if (ar) {
        names_rhos <- paste("rho", 1:ar, sep="") 
        lbound_sp[names_rhos] <- -rholimit
        ubound_sp[names_rhos] <-  rholimit
    }

    # AR(0) uses this b/c of the economy of lpri functions
    assimctx$rholimit <- rholimit

    if (sigma) {
        if (inv_gamma_pri) {
            assimctx$alpha <- alpha
            assimctx$beta  <- beta

            lbound_sp["var"]   <- gtzero()
            ubound_sp["var"]   <- var_max
        } else {
            # sigma must be greater than zero to prevent singular matrices
            # in llik_ar1_obs() and llik_ar2_obs()
            #
            lbound_sp["sigma"] <- gtzero()
            ubound_sp["sigma"] <- sigma_max
        }
        assimctx$logPri <- lpri_sigma
    } else {
        assimctx$logPri <- lpri_bounds
    }

    assimctx$lbound_sp   <- lbound_sp
    assimctx$ubound_sp   <- ubound_sp
    assimctx$ar          <- ar
    assimctx$obserr      <- obserr
    assimctx$maxLik      <- -Inf
    assimctx$superLogLik <- llikfn

    # predict.R and assim.R use this
    assimctx$mp_indices  <- 1:length(assimctx$lbound)

    # TODO:  eradicate some other instances of setErrorAssim()
    if ("error" %in% colnames(assimctx$obs)) {
        setErrorAssim(assimctx, assimctx$obs[assimctx$obs_ind, "error"])
    } else {
        print("WARNING! not setting observation error in configAssim()")
    }


    if (!is.null(init_mp)) {

        if (is.null(init_sp)) {

            # create statistical parameters for chain as necessary
            #

            # sanity check
            if (!all(init_mp >= assimctx$lbound, init_mp <= assimctx$ubound)) {
                stop("initial model parameters out of bounds for assimilation")
            }

            init_sp <- numeric()

            # get residuals from best fit model
            res <- assimctx$obsonly - assimctx$modelfn(init_mp, assimctx)

            # get rho1 and rho2 from residuals
            if (ar) {
                pac  <- pacf(res, lag.max=ar, plot=FALSE)
                rhos <- pac$acf[ 1:ar ]
                init_sp[names_rhos] <- rhos

                # whiten residuals before calculating standard deviation
                res <- arwhiten(res, rhos)
            }

            # get sigma from standard deviation of residuals
            if (sigma) {
                if (inv_gamma_pri) {
                    init_sp["var"]   <- sd(res)^2
                } else {
                    init_sp["sigma"] <- sd(res)
                }
            }
        }

    } else {

        init_p  <- assimMaxLikelihood(assimctx, itermax=itermax)
        init_mp <- init_p[  assimctx$mp_indices ]
        init_sp <- init_p[ -assimctx$mp_indices ]
    }

    if (any(iszero(init_mp))) {
        print("WARNING! Initial model parameter has value zero (could affect scaling:)")
        print(init_mp)
    }

    if (fixrho) {
        assimctx$lbound_sp[names_rhos] <- init_sp[names_rhos]
        assimctx$ubound_sp[names_rhos] <- init_sp[names_rhos]
    }

    assimctx$init_mp <- init_mp
    assimctx$init_sp <- init_sp
}


# potential function for extrafun in runAssim;
# ychain can be pre-allocated in xxxRunAssim() by caller
assimSaveY <- function(i, assimctx)
{
    assimctx$ychain[i, ] <- assimctx$y
    return ()
}


runAssim <- function(assimctx, nbatch, nspac=1, scale=NULL,
    adapt=F, acc.rate = 0.234, gamma=0.5, n.start=0.01*nbatch, extrafun=NULL)
{
    # can verify changes by seeing if they produce the same chain
    #set.seed(7)

    if (is.null(scale)) {
        range_mp <- assimctx$ubound    - assimctx$lbound
        range_sp <- assimctx$ubound_sp - assimctx$lbound_sp
        scale <- abs(c(range_mp / 100, range_sp / 1000))

        if (assimctx$ar == 2) {
            #scale["rho2"] <- scale["rho2"] / 2
        }

        print("using initial scale:")
        print(scale)
    }

    if (adapt) {
        if (!missing(nspac)) {
            print("WARNING! adaptive MCMC ignores nspac value")
        }
        time <- system.time(
            out <- named_MCMC(p=logPost, n=nbatch,
                init_mp=assimctx$init_mp, init_sp=assimctx$init_sp, scale=scale,
                adapt=adapt, acc.rate=acc.rate, gamma=gamma, list=T, n.start=n.start,
                extrafun=extrafun, assimctx=assimctx
                )
        )
        out$time <- time
    } else {
        out <- named_metrop(
            obj=logPost,
            init_mp=assimctx$init_mp, init_sp=assimctx$init_sp,
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


uniformPrior <- function(min, max)
{
    obj <- env()

    obj$dens <- function(x, log=T) dunif(x=x, min=min, max=max, log=log)
    obj$rand <- function(n)        runif(n=n, min=min, max=max)
    obj$mode <- function()         ((min + max) / 2)
    obj$mean <- obj$mode

    return (obj)
}


betaPrior <- function(min, max, a, b)
{
    obj       <- env()
    scale     <- 1 / (max - min)
    log_scale <- log(scale)

    to_x <- function(xt) {
        x  <- (xt / scale) + min
        return (x)
    }

    to_xt <- function(x) {
        xt <- scale * (x - min)
        return (xt)
    }

    obj$dens <- function(x, log=T) {
        xt <- to_xt(x)
        if (log) {
            p  <- log_scale + dbeta(xt, shape1=a, shape2=b, log=T)
        } else {
            p  <-     scale * dbeta(xt, shape1=a, shape2=b, log=F)
        }
        return (p)
    }

    obj$rand <- function(n) {
        xt <- rbeta(n, shape1=a, shape2=b)
        return (to_x(xt))
    }

    obj$mean <- function() {
        xt <- a / (a + b)
        return (to_x(xt))
    }

    obj$mode <- function() {
        xt <- (a - 1) / (a + b - 2)
        return (to_x(xt))
    }

    obj$var <- function() {
        xt <- (a * b) / ((a + b)^2 * (a + b + 1))
        return (to_x(xt))
    }

    return (obj)
}


normPrior <- function(mean, upper)
{
    obj <- env()
    sd  <- (upper - mean) / 2  # note that this is assuming the upper is 2-sigma

    obj$rand <- function(n)        rnorm(n=n, mean=mean, sd=sd)
    obj$dens <- function(x, log=T) dnorm(x=x, mean=mean, sd=sd, log=log)
    obj$mode <- function() (mean)
    obj$mean <- function() (mean)

    return (obj)
}


gammaPrior <- function(shape, rate)
{
    obj <- env()

    obj$rand <- function(n)        rgamma(n=n, shape=shape, rate=rate)
    obj$dens <- function(x, log=T) dgamma(x=x, shape=shape, rate=rate, log=log)
    obj$mode <- function() ((shape - 1) / rate)
    obj$mean <- function()  (shape      / rate)

    return (obj)
}


lnormPrior <- function(mean, upper)
{
    obj     <- env()
    meanlog <- log(mean)
    sdlog   <- log(upper / mean) / 2  # note that this is assuming the upper is 2-sigma

    obj$rand <- function(n)        rlnorm(n=n, meanlog=meanlog, sdlog=sdlog)
    obj$dens <- function(x, log=T) dlnorm(x=x, meanlog=meanlog, sdlog=sdlog, log=log)
    obj$mode <- function() exp( meanlog - sdlog ^ 2     )
    obj$mean <- function() exp( meanlog + sdlog ^ 2 / 2 )

    return (obj)
}


llnormPrior <- function(mean, upper)
{
    obj     <- env()
    meanlog <- log(mean)
    sdlog   <- log(upper / mean) / 2  # note that this is assuming the upper is 2-sigma

    obj$rand <- function(n)        rnorm(n=n, mean=meanlog, sd=sdlog)
    obj$dens <- function(x, log=T) dnorm(x=x, mean=meanlog, sd=sdlog, log=log)
    obj$mode <- function() ( meanlog - sdlog ^ 2     )
    obj$mean <- function() ( meanlog + sdlog ^ 2 / 2 )

    return (obj)
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


runPredict <- function(
    nbatch=100000, year=2200,
    assimctx,
    prctx,
    forcings,
    modelfn,
    noisefn=list(assimctx$noise),
    outnames=c("prchain"),
    mcmcChain=assimctx$chain,
    xvals=assimctx$times[1]:year,
    replace=T,
    ...
    )
{
    prctx$assimctx <- assimctx

    nChain <- sample(burnedInd(mcmcChain), nbatch, replace=replace)

    nfrc <- length(forcings)
    samples <- matrix(nrow=nbatch, ncol=nfrc)

    # for verification against old prediction functions:
    # comment the for loop and uncomment the subsequent line
    # in order to fix the problem of the seed changing
    #
    for (i in safefor(1:nfrc)) {
        if (2 == ncol(forcings[[i]])) {
            samples[, i] <- rep(2L, nbatch)
        } else {
            samples[, i] <- sample(2:ncol(forcings[[i]]), nbatch, replace=replace)
        }
    }
    #samples[, 1] <- sample(2:ncol(forcings[[1]]), nbatch, replace=replace)

    outchains <- vector("list", length=length(outnames))
    for (i in 1:length(outnames)) {
        outchains[[i]] <- prmatrix(nbatch, xvals)
    }

    for (i in 1:nbatch) {
        p  <- mcmcChain[ nChain[i], ]
        sp <- p[ -assimctx$mp_indices ]

        # this was too expensive;  R mistakenly thought frc was being changed and kept copying forcings
        #pred <- modelfn(xvals, list(mp=p[ assimctx$mp_indices ], frc=forcings, spl=samples[i, ], ...))

        # this next block of code is MUCH faster (>= one order of magnitude)
        frc <- list()
        for (j in safefor(1:nfrc)) {
            forcing <- forcings[[j]]
            frc[[j]] <- cbind(forcing[, 1], forcing[, samples[i, j]])
            colnames(frc[[j]]) <- colnames(forcing)[ c(1, i) ]
        }
        pred <- modelfn(xvals, list(mp=p[ assimctx$mp_indices ], frc=frc, spl=rep(2, nfrc), ...))

        for (j in 1:length(outchains)) {
            noise <- noisefn[[j]](sp, nrow(pred), assimctx)
            outchains[[j]][i, ] <- pred[, j + 1 ] + noise
        }
    }

    for (i in 1:length(outnames)) {
        assign(outnames[[i]], outchains[[i]], envir=prctx)
    }

    if (nbatch <= 20) {
        print(outchains[[1]])
    }
}


noise_save <- function(sp, N, assimctx)
{
    assimctx$noise_realization <- assimctx$noise(sp, N, assimctx)
}


noise_copy <- function(sp, N, assimctx)
{
    assimctx$noise_realization
}


noise_zeros <- function(sp, N, assimctx)
{
    rep(0, N)
}


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


assimMinFn <- function(p, assimctx)
{
    mp <- p[  assimctx$mp_indices ]
    sp <- p[ -assimctx$mp_indices ]

    lpost <- logPost(mp, sp, assimctx)
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
