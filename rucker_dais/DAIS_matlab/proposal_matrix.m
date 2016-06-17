%-----------------------------------------------------------------------------
% R Function: Daisobs_likelihood_iid.R (holds the proposal.atrix function)
% -Antarctic Ice Sheet (AIS) model
%
% compute (log) likelihood for observations
% The observations are indepent and identically distributed (IID)
%
% -Author: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
%-----------------------------------------------------------------------------
% -June 10 2015; R code converted to matlab March 7th
%-----------------------------------------------------------------------------

% calculate a scale matrix to transform a iid normal distribution,
% which defines a multivariate normal proposal distribution for MCMC
% (the proposal distribution is proportional to the covariance of
% a preliminary Markov chain of the posterior distribution to sample,
% tuned to be optimally scaled if the posterior is multivariate normal,
% plus a small nugget term to ensure ergodicity)
%
% from Gareth O. Roberts and Jeffrey S. Rosenthal,
% "Examples of adaptive MCMC", unpublished

function mat = proposal_matrix(prechain, mult)

beta = 0.05;  % beta = relative influence of nugget term
%mult = 1;    % mult = overall scale factor to adjust all step sizes

p = length(prechain(1,:)); %ncol(prechain);
precov = cov(prechain);

propcov = (1-beta)*2.38^2*precov/p + beta*0.1^2*diag(p)/p;
propcov = mult*propcov;

mat = transpose(chol(propcov)); % transpose and the choleski factorization


