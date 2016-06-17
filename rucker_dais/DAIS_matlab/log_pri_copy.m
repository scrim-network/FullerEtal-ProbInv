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

sigma_y = theta(12); %sigma_y is the variance

% Best Case (Case #4) from Shaffer (2014)
IP = [2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006];

%Set the upper and lower bounds
bound_lower = IP - (IP*0.5)    ; bound_upper = IP + (IP*0.5);
bound_lower(1:2) = [1/2, 0]   ; bound_upper(1:2) = [17/4, 1]; %Set bounds for gamma and alpha
bound_lower(10:11) = [725, 0.00045]   ; bound_upper(10:11) = [825, 0.00075]; %Set bounds for bo and s
%bound_lower(12) = 0 ; bound_upper(12) = 1; % Prior uniform range for sigma (the variance)

in_range = all(par > bound_lower) & all(par < bound_upper);

alpha_var = 2; beta_var = 1;
var_pri = (-alpha_var - 1)*log(sigma_y) + (-beta_var/sigma_y);

if in_range;
lpri=0 + var_pri;
else;
lpri = -Inf;
end

lpri =lpri;
