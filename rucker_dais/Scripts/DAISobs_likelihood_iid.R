#####################################################################################
# R Function: Daisobs_likelihood_iid.R
# -Antarctic Ice Sheet (AIS) model
#
# compute (log) likelihood for observations
# The observations are indepent and identically distributed (IID)
#
# -Author: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
#####################################################################################
# -June 10 2015
#####################################################################################

log.lik = function(p) # model.p is the dimension of model parameters 
{ 
  par=p[1:model.p]
  
  sigma.y = p[model.p+1]
  #rho.y = p[model.p+2]
  
  model.out = iceflux(par, hindcast.forcings, standards)
  y.mod = model.out
  
  llik.y  = 0
  
  #get the residuals
  r1 = median(windows[1,]) - (y.mod[120000] - mean(y.mod[SL.1961_1990]))
  r2 = median(windows[2,]) - (y.mod[220000] - mean(y.mod[SL.1961_1990]))
  r3 = median(windows[3,]) - (y.mod[234000] - mean(y.mod[SL.1961_1990]))
  r4 = median(windows[4,]) - (y.mod[240002] - mean(y.mod[SL.1961_1990]))

  resid.y = c(r1, r2, r3, r4)
  sterr.y = obs.errs #This makes the model heteroskedastic
  
  #Calculate the likelihood. The observations are not correlated. They are independent
  llik.y = sum (dnorm(resid.y, mean=rep(0,4), sd = sqrt(sigma.y + sterr.y ^2), log=TRUE))
  
  llik = llik.y # assume residuals are independent
  llik  
}
log.pri = function(p)
{
  par=p[1:model.p]
  
  sigma.y = p[model.p+1]
  #rho.y = p[model.p+2]

  in.range = all(par > bound.lower) & all(par < bound.upper)
  
  alpha_var = 2
  beta_var = 1
  var_pri = (-alpha_var - 1)*log(sigma.y) + (-beta_var/sigma.y)
  
  if(in.range) {
    lpri=0 + var_pri
  } else {
    lpri = -Inf
  }

  lpri
}

# (log) posterior distribution:  posterior ~ likelihood * prior
log.post = function(p)
{  
  lpri = log.pri(p)
  if(is.finite(lpri)) { # evaluate likelihood if nonzero prior probability
    lpost = log.lik(p) + lpri
  } else {
    lpost = -Inf  
  }
  lpost
}
 

# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
proposal.matrix = function(prechain, mult=1, beta=0.05)
{
  # mult = overall scale factor to adjust all step sizes
  # beta = relative influence of nugget term
  
  p = ncol(prechain)
  precov = cov(prechain)
  
  propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p
  propcov = mult*propcov
  
  mat = t(chol(propcov))
}

