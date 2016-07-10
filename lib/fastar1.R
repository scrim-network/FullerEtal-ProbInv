# Copyright 2010 Nathan M. Urban

source("roblib.R")
dynReload("fastar1")

# generate AR(1) time series of length n with lag-1 autocorrelation coefficient rho and innovation variance sigma^2
ar1.sim = function(n, rho, sigma)
{
    w = rnorm(n)
    
    cout = .C( "ar1sim",
                n = as.integer(n),
                rho = as.double(rho),
                sigma = as.double(sigma),
                w = as.double(w),
                r = as.double(rep(0,n))
            )
    
    return(cout$r)
}
