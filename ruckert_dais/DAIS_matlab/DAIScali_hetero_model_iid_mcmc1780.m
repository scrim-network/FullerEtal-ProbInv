%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  - file = "DAIScali_hetero_model_iid_mcmcRversion.R"
%  - Code written: August 2015, updated March 2016
%  - Author: Kelsey Ruckert (klr324@psu.edu)
%
%  -This program runs a Markov Chain Monte Carlo analysis of the DAIS model
%       assuming heteroskedastic errors and IID residuals as described in Ruckert et al. (2016).
%       For further description and references, please read the paper.
%
% THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
% NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
% BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
% APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
% AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Set the seed.
%rand('seed', 1234);  % For reproducibility
rand('seed', 1780);  % For reproducibility

% Read in hindcast/ forcing data, read in standard values, and read in AIS dates both specific and ranges so there is no use of magic numbers.
run DAIS_data.m

%-------------------- List of physcial model parameters --------------------%
% We will set the initial parameters to specifications from Shaffer [2014]
%Define the parameters:
% [1] Gamma = 2 			  %sensitivity of ice flow to sea level
% [2] alpha = 0.35 			  %sensitivity of ice flow to ocean subsurface temperature
% [3] mu = 8.7    			  %Profile parameter related to ice stress [m^(1/2)]
% [4] nu = 0.012   		      %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
% [5] P0 = 0.35    			  %Precipitation at 0C [m of ice/yr]
% [6] kappa = 0.04 			  %Relates precipitation to temperature [K^-1]
% [7] f0 = 1.2                %Constant of proportionality for ice speed [m/yr]
% [8] h0 = 1471               %Initial value for runoff line calculation [m]
% [9] c = 95                  %Second value for runoff line calculation [m]
% [10] b0 = 775               %Height of bed at the center of the continent [m]
% [11] slope = 0.0006         %Slope of the bed

%Create matrices for projections and hindcasts. Create a matrix that has 4 columns and 240300 rows
project_forcings(1:240300,1)=Ta; project_forcings(1:240300,2)=Toc;
project_forcings(1:240300,3)=GSL; project_forcings(1:240300,4)=SL;

% 240010x4 matrix
hindcast_forcings(1:240010,1)=Ta(1:240010); hindcast_forcings(1:240010,2)=Toc(1:240010);
hindcast_forcings(1:240010,3)=GSL(1:240010); hindcast_forcings(1:240010,4)=SL(1:240010);

% Set initial parameters to the best case (Case #4) from Shaffer (2014)
IP = [2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006];

%Source the function with the standards and the initial parameters (IP) to
%get the best estimated AIS volume loss with respect to the present day in sea level equivalence (SLE):
standards = [Tf, rho_w, rho_i, rho_m, Toc_0, Rad0, Volo];

% Load in the physical model "DAIS_ICEFLUX_model" and run with Case #4 parameters:
AIS_melt = DAIS_IceFlux_model(IP, hindcast_forcings, standards);

% Estimate the AIS volume loss projection for Case #4 with respect to the present day in sea level equivalence (SLE):
Project_melt = DAIS_IceFlux_model(IP, project_forcings, standards);
%%
%--------------------  Setup observational constraint info --------------------%
% These windows are presented in Shaffer (2014) and calculated from Shepherd et al. (2012)
% # Calculate the uncertainty with the +/- 2 standard error

% Create a vector with each observation year
%            120kyr,  20Kyr,   6kyr,   2002
obs_years = [120000, 220000, 234000, 240002];

% Accummulate the sea-level equivalent in meters from 1992 to the year 2002
% using the 1992 to 2011 trend from Shepherd et al. 2012; -71 +/- 53 Gt per yr.
% Conversion: 360 Gt = 1 mm SLE
estimate_SLE_rate = abs(-71/360)/1000;
time_years = 2002-1992;
mid_cum_SLE_2002 = estimate_SLE_rate*time_years;

estimate_SLE_error = sqrt(time_years)*abs(-53/360)/1000; % 1-sigma error
% (*sqrt(10) because 10 years of potentially accumulated error:
%  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
%                = 10*year X error^2)
SE2_2002 = estimate_SLE_error*2; %2-sigma error

% Add and subtract the 2 standard error to the mean value
positive_2SE = mid_cum_SLE_2002 + SE2_2002; 
negative_2SE = mid_cum_SLE_2002 - SE2_2002; 

% Create observational constraint windows.
upper_wind = [6.0, -6.9, -1.25, positive_2SE];
lower_wind = [1.8, -15.8, -4.0, negative_2SE];
windows(1:4,1) = lower_wind;
windows(1:4,2) = upper_wind;

% Determine observational error from windows: half-width of window = uncertainty; assume all windows are 2*stdErr (last one actually is)
obs_errs = (windows(:,2)-windows(:,1))*.5;       
%%
%-------------------- CALCULATE RESIDUALS and INITIAL SIGMA VALUE --------------------%
resid(1:length(obs_years)) = NaN;     %Create a vector of the residuals

% Estimate the residuals: modification from equation (1)
for i = 1:length(obs_years)
resid(i) = (median(windows(i,:))-(AIS_melt(obs_years(i))-mean(AIS_melt(SL_1961_1990))));
end

% Calculate the variance, sigma^2
sigma = std(resid)^2;

%%
%--------------------  SETUP MCMC --------------------%
%Set up the priors; the upper and lower bounds 
% var_y has inverse gamma prior, so there is a lower bound at 0 but no upper bound
%parnames   = ['gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c', 'b0','slope' ,'var_y']
bound_upper = [ 4.25 ,  1     , 13.05, 0.018,  1  ,  0.06 , 1.8 ,2206.5, 142.5, 825 , 0.00075,   Inf];
bound_lower = [ 0.5  ,  0     , 4.35 , 0.006,0.175,  0.02 , 0.6 , 735.5,  47.5, 725 , 0.00045 ,    0];

% Set up the constraint structure for likelihood function
constraint.LIG = median(windows(1,:));
constraint.LGM = median(windows(2,:));
constraint.MH = median(windows(3,:));
constraint.IP = median(windows(4,:));
constraint

% SL_1961_1990 = 239961:239990;

%print out obs_errs
obs_errs

% Setup hindcast forcings to be called during MCMC calibration
global data; 
data = hindcast_forcings;

% Call physical and statistical model
logmodelprior=@log_pri_copy; % prior.
loglike=@log_lik_calibration_copy; % log likelihood.

% Use Shaffer [2014] Case #4 parameters as the initial starting values for
% parameters.
minit = [IP, sigma];

% Set the step size, number of iterations, and amount of thinning.
nsimu = 1000;  % number of simulations
%step = minit/150;
%step = [0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15]/5;
step = [0.05, 0.01, 0.15, 0.035, 0.1, 0.01, 0.1, 50, 10, 30, 0.0005, 0.1]/5;
thin=5;

%%
%-------------- Run the MCMC chain --------------%
%rand('seed', 1234);  % For reproducibility
rand('seed', 1780);  % For reproducibility

mmc=mcmc(minit,loglike,logmodelprior,step,nsimu,thin);
mmc1 = mmc;

%Estimate a more appropriate step size:
% step2 = proposal_matrix(mmc1,0.5);    % step2(1:12) = diag(scale);

% New starting value:
minit2 = mmc1(length(mmc1),:);

rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

nsimu2 = 1.2e6;
mmc = mcmc(minit2, loglike, logmodelprior, step, nsimu2, thin);
mmc2 = mmc;

% Save chain output:
%save('DAIS_MCMCchain_1234', 'mmc1','mmc2')
save('DAIS_MCMCchain_1780', 'mmc1','mmc2')

%--------------------  Analysis of the MCMC chain produced --------------------%
% m(1:100,:)=[]; %crop drift
% plotmatrix(m);
%
%% Plot some figures with the chain:
%figure(1); clf
%mcmcplot(mmc2,[],[],'chainpanel')

%figure(2); clf
%mcmcplot(mmc2,[],[],'dens')

%figure(3); clf
%mcmcplot(mmc2,[],[],'hist')

%save('DAIS_Matlab_MCMCcalibration','-v7.3') %Save the workspace into a R readable file

% The rest of the analysis will be run in R
%----------------------------------------------------------------------%

