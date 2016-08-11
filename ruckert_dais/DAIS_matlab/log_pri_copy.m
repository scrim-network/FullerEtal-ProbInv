%-----------------------------------------------------------------------------
% R Function: Daisobs_likelihood_iid.R (holds the log.pri function)
% -Antarctic Ice Sheet (AIS) model
%
% compute (log) likelihood for observations
% The observations are indepent and identically distributed (IID)
%
% -Author: Yawen Guan and Kelsey Ruckert (klr324@psu.edu)
%-----------------------------------------------------------------------------
% -June 10 2015; R code converted to matlab March 7th
%-----------------------------------------------------------------------------

function lpri = log_pri(theta,data)
par=theta(1:11);

var_y = theta(12); % the variance

% Set the upper and lower bounds 
% var_y has inverse gamma prior, so there is a lower bound at 0 but no upper bound
%parnames   = ['gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c', 'b0','slope' ,'var_y']
bound_upper = [ 4.25 ,  1     , 13.05, 0.018,  1  ,  0.06 , 1.8 ,2206.5, 142.5, 825 , 0.00075,   Inf];
bound_lower = [ 0.5  ,  0     , 4.35 , 0.006,0.175,  0.02 , 0.6 , 735.5,  47.5, 725 , 0.00045 ,    0];

in_range = all(theta > bound_lower) & all(theta < bound_upper);

% Set the inverse gamme prior for the variance
alpha_var = 2; 
beta_var = 1;
var_pri = 0;

if in_range;
    var_pri = (-alpha_var - 1)*log(var_y) + (-beta_var/var_y);
    lpri=0 + var_pri;
else;
    lpri = -Inf;
end

lpri =lpri;
