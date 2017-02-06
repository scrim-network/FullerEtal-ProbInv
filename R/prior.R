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
# prior.R


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
