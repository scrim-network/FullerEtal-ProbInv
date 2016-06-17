#####################################################################################
# R Function: Daisobs_likelihood_iid.R
# -Antarctic Ice Sheet (AIS) model
#
# compute (log) likelihood for observations
# The observations are indepent and identically distributed (IID)
#
# -Author: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
#####################################################################################
# -June 10 2015 R version; written in Julia May 25, 2016
#####################################################################################

function log_lik(p) # model.p is the dimension of model parameters

    par = p[1:model_p];

    sigma_y = p[model_p+1];

    model_out = iceflux(par, hindcast_forcings, standards, endinfo);
    y_mod = model_out;

    llik_y  = 0;

    #get the residuals
    r1 = median(windows[1,:]) - (y_mod[120000] - mean(y_mod[SL_1961_1990]));
    r2 = median(windows[2,:]) - (y_mod[220000] - mean(y_mod[SL_1961_1990]));
    r3 = median(windows[3,:]) - (y_mod[234000] - mean(y_mod[SL_1961_1990]));
    r4 = median(windows[4,:]) - (y_mod[240002] - mean(y_mod[SL_1961_1990]));

    resid_y = [r1, r2, r3, r4];
    sterr_y = obs_errs; #This makes the model heteroskedastic

    #Calculate the likelihood. The observations are not correlated. They are independent
    llik_y   = zeros(length(resid_y),1);
    for i in 1:length(resid_y)
      llik_y[i] = logpdf(Normal(0, sqrt(sigma_y + sterr_y[i].^2)), resid_y[i]);
    end
    llik_y = sum(llik_y);
  
    llik = llik_y; # assume residuals are independent

    return llik;
end

function log_pri(p)

    par=p[1:model_p];
  
    sigma_y = p[model_p+1];

    in_range_1 = bound_lower[1] < par[1] < bound_upper[1];
    in_range_2 = bound_lower[2] < par[2] < bound_upper[2];
    in_range_3 = bound_lower[3] < par[3] < bound_upper[3];
    in_range_4 = bound_lower[4] < par[4] < bound_upper[4];
    in_range_5 = bound_lower[5] < par[5] < bound_upper[5];
    in_range_6 = bound_lower[6] < par[6] < bound_upper[6];
    in_range_7 = bound_lower[7] < par[7] < bound_upper[7];
    in_range_8 = bound_lower[8] < par[8] < bound_upper[8];
    in_range_9 = bound_lower[9] < par[9] < bound_upper[9];
    in_range_10 = bound_lower[10] < par[10] < bound_upper[10];
    in_range_11 = bound_lower[11] < par[11] < bound_upper[11];

    in_range = [in_range_1, in_range_2, in_range_3, in_range_4, in_range_5, in_range_6, in_range_7, in_range_8, in_range_9, in_range_10, in_range_11];

    in_range = all(in_range);
  
    alpha_var = 2;
    beta_var = 1;
    var_pri = (-alpha_var - 1)*log(sigma_y) + (-beta_var/sigma_y);

    if in_range == true
        lpri=0 + var_pri;
    else
        lpri = -Inf;
    end
    return lpri;
end

# (log) posterior distribution:  posterior ~ likelihood * prior
function log_post(p)

    lpri = log_pri(p);
if isfinite(lpri); # evaluate likelihood if nonzero prior probability
        lpost = log_lik(p) + lpri;
    else
        lpost = -Inf;
    end
    return lpost;
end
 

# calculate a scale matrix to transform a iid normal distribution,
# which defines a multivariate normal proposal distribution for MCMC
# (the proposal distribution is proportional to the covariance of
# a preliminary Markov chain of the posterior distribution to sample,
# tuned to be optimally scaled if the posterior is multivariate normal,
# plus a small nugget term to ensure ergodicity)
#
# from Gareth O. Roberts and Jeffrey S. Rosenthal,
# "Examples of adaptive MCMC", unpublished
function proposal_matrix(chain; mult=1., beta=0.05)

    # mult = overall scale factor to adjust all step sizes
    # beta = relative influence of nugget term
    prechain = chain.value';
  
    p = size(prechain, 2);
    precov = cov(prechain);

    propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diagm(p)/p;
    propcov = mult*propcov;

    mat = chol(propcov)';
    return mat  #convert(Array, mat)
end

