%-----------------------------------------------------------------------------
%------- Now let's calibrate the model parameters
%------- using iid model
%%%-----HETEROSKEDASTIC
%-----------------------------------------------------------------------------

clear all
close all

rand('seed', 1234);  % For reproducibility

% step 1 define the boundary for parameters
run DAIS_data.m

%-------------------- SET INITIAL PARAMETERS --------------------%
% We will set the initial parameters to specifications from Shaffer [2014]
%Define the parameters:
% [1] gamma = 2 			  %sensitivity of ice flow to sea level
% [2] alpha = 0.35 			  %sensitivity of ice flow to ocean subsurface temperature
% [3] mu = 8.7    			  %Profile parameter related to ice stress [m^(1/2)]
% [4] eta = 0.012   		  %Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
% [5] Po = 0.35    			  %Precipitation at 0C [m of ice/yr]
% [6] kappa = 0.04 			  %Relates precipitation to temperature [K^-1]
% [7] fo = 1.2                %Constant of proportionality for ice speed [m/yr]
% [8] ho = 1471               %Initial value for runoff line calculation [m]
% [9] co = 95                 %Second value for runoff line calculation [m]
% [10] bo = 775               %Height of bed at the center of the continent [m]
% [11] s = 0.0006             %Slope of the bed

%Create a matrix that has 4 columns and 240300 rows
project_forcings(1:240300,1)=TA; project_forcings(1:240300,2)=TO;
project_forcings(1:240300,3)=GSL; project_forcings(1:240300,4)=SL;

% 240010x4 matrix
hindcast_forcings(1:240010,1)=TA(1:240010); hindcast_forcings(1:240010,2)=TO(1:240010);
hindcast_forcings(1:240010,3)=GSL(1:240010); hindcast_forcings(1:240010,4)=SL(1:240010);


% Best Case (Case #4) from Shaffer (2014)
IP = [2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006];

%Source the function with the standards and the initial parameters (IP) to
%get the best estimated AIS volume loss with respect to the present day in sea level equivalence (SLE):
standards = [Tice, eps1, del, eps2, TOo, Volo, Roa, R];
endinfo = [endsat, enddate];

%source("Scripts/DAIS_IceFlux_model.R")
AIS_melt = DAIS_IceFlux_model(IP, hindcast_forcings, standards, endinfo);

%set the end dates to the year 2300 to get future projections
endsat = 240298;
enddate = 240300;
endinfo = [endsat, enddate];
Project_melt = DAIS_IceFlux_model(IP, project_forcings, standards, endinfo);

%Set the end dates back to the hindcast period:
endsat = 240000;
enddate = 240010;
endinfo = [endsat, enddate];

%--------------------  CALCULATE RESIDUALS (PRIOR SIGMA) --------------------%
% For this model the residuals are based off of the windowing approach and are weighted
% These windows are presented in Shaffer (2014) and calculated from Shepherd et al. (2012)
% # Calculate the uncertainty with the +/- 2 standard error

% Create a vector with each observation year
% 120kyr, 20Kyr, 6kyr, 2002
obs_years = [120000, 220000, 234000, 240002];

% 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
% We want the cumulative sea-level equivalent in meters for the year 2002
% 360Gt = 1mm SLE
estimate_SLE_rate = abs(-71/360)/1000;
time_years = 2002-1992;
mid_cum_SLE_2002 = estimate_SLE_rate*time_years;

estimate_SLE_error = abs(-53/360)/1000; %1- sigma error
SE2_2002 = estimate_SLE_error*2; %2-sigma error

positive_2SE = mid_cum_SLE_2002 + SE2_2002; % Add the 2 standard error to the mean value
negative_2SE = mid_cum_SLE_2002 - SE2_2002; % Subtract the 2 standard error to the mean value

upper_wind = [6.0, -6.9, -1.25, positive_2SE];
lower_wind = [1.8, -15.8, -4.0, negative_2SE];
windows(1:4,1) = lower_wind;
windows(1:4,2) = upper_wind;

obs_errs = [abs(median(windows(1,:))-windows(1,1)); abs(median(windows(2,:))-windows(2,1));
            abs(median(windows(3,:))-windows(3,1)); SE2_2002];

%Set up equation to find the residuals and then the prior sigma
resid(1:length(obs_years)) = NaN;     %Create a vector of the residuals

for i = 1:length(obs_years)
resid(i) = (median(windows(i,:))-(AIS_melt(obs_years(i))-mean(AIS_melt(SL_1961_1990))));
end

sigma = std(resid); %calculate the standard deviation (sigma);

%--------------------  RUN MCMC --------------------%
%Set up the priors and ranges for the MCMC
%Set the upper and lower bounds
bound_lower = IP - (IP*0.5)    ; bound_upper = IP + (IP*0.5);
bound_lower(1:2) = [1/2, 0]   ; bound_upper(1:2) = [17/4, 1]; %Set bounds for gamma and alpha
bound_lower(10:11) = [725, 0.00045]   ; bound_upper(10:11) = [825, 0.00075]; %Set bounds for bo and s

bound_lower(12) = 0 ; bound_upper(12) = 1; % Prior uniform range for sigma (the variance)

% step 2 define number of model parameters
model_p = 11;
%parnames=['gamma', 'alpha', 'mu', 'eta', 'po', 'kappa', 'fo', 'ho', 'co', 'bo', 's', 'sigma_y'];
% step 3 source the physical model and statistical model
%source("Scripts/DAIS_IceFlux_model.R")
%source("Scripts/DAISobs_likelihood_iid.R")

%Shaffer [2014] best guess parameters
p = [IP, sigma];
p0 = [2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.6];

%----------------- Try to find good initial parameters ------------%
% NOTE: This does not work yet so skip for now

%options.MaxIterations = 10;
%p0=fminsearch(ssfun,p0,options,data);
%round(p0,4)
%clear options
%install.packages("mcmc")
%library(mcmc)

%%
%----------------- set up mcmc info ------------%
%step = [0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15]/100;  % Old step size: step = [0.001, 0.0001, 0.001, 0.00001, 0.0001, 0.00001, 0.001, 0.5, 0.1, 0.5, 0.000001, 0.001]
%NI = 10;  %test number of iterations NI = 1E6; %number of iterations
%burnin = linspace(1, 0.01*NI, 0.01*NI);
%------------------------------------------------
% Calibration with mcmcrun (matlab) instead of metrop (R)
%
%MCMCRUN Metropolis-Hastings MCMC simulation for nonlinear Gaussian models
%function [results,chain,s2chain,sschain, hchain]=mcmcrun(model,data,params,options,res)
% properties:
%  multiple y-columns, sigma2-sampling, adaptation,
%  Gaussian prior, parameter limits, delayed rejection, dram
%
% [RESULTS,CHAIN,S2CHAIN,SSCHAIN] = MCMCRUN(MODEL,DATA,PARAMS,OPTIONS)
% MODEL   model options structure
%    model.ssfun    -2*log(likelihood) function
%    model.priorfun -2*log(pior) prior function
%    model.sigma2   initial error variance
%    model.N        total number of observations
%    model.S20      prior for sigma2
%    model.N0       prior accuracy for S20
%    model.nbatch   number of datasets
%
% DATA the data, passed directly to ssfun. The structure of DATA is given
%      by the user. A possible 'time' variable must be given in the first column of
%      xdata. Note that only data.xdata is needed for model simulations.
%      In addition, DATA may include any user defined structure needed by
%      |modelfun| or |ssfun|
%
% PARAMS  theta structure
%   {  {'par1',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      {'par2',initial, min, max, pri_mu, pri_sig, targetflag, localflag}
%      ... }
%
% OPTIONS mcmc run options
%    options.nsimu            number of simulations
%    options.qcov             proposal covariance
%    options.method           'dram','am','dr' or 'mh'
%    options.adaptint         interval for adaptation, if 'dram' or 'am' used
%                             DEFAULT adaptint = 100
%    options.drscale          scaling for proposal stages of dr
%                             DEFAULT 3 stages, drscale = [5 4 3]
%    options.updatesigma      update error variance. Sigma2 sampled with updatesigma=1
%                             DEFAULT updatesigma=0
%    options.verbosity        level of information printed
%    options.waitbar          use graphical waitbar?
%    options.burnintime       burn in before adaptation starts
%
% Output:
%  RESULTS   structure that contains results and information about
%            the simulations
%  CHAIN, S2CHAIN, SSCHAIN
%           parameter, sigma2 and sum-of-squares chains
%
%------------------------------------------------
%%
% For the MCMC run we need the sum of squares function. For the
% plots we shall also need a function that returns the model.
% Both the model and the sum of squares functions are
% easy to write as one line anonymous functions using the @
% construct.

% Set up the structure for likelihood function
constraint.LIG = median(windows(1,:));
constraint.LGM = median(windows(2,:));
constraint.MH = median(windows(3,:));
constraint.IP = median(windows(4,:));
constraint

% SL_1961_1990 = 239961:239990;

%print out obs_errs
obs_errs

data = hindcast_forcings;

modelfun = @DAIS_calibration_model;
ssfun    = @log_lik_calibration;


nsimu = 1e5;  % number of simulations
npar  = 12;     % dimension of the unknown


params = {
    {'gamma', p0(1), bound_lower(1), bound_upper(1)}
    {'alpha', p0(2), bound_lower(2), bound_upper(2)}
    {'mu', p0(3), bound_lower(3), bound_upper(3)}
    {'eta', p0(4), bound_lower(4), bound_upper(4)}
    {'Po', p0(5), bound_lower(5), bound_upper(5)}
    {'kappa', p0(6), bound_lower(6), bound_upper(6)}
    {'fo', p0(7), bound_lower(7), bound_upper(7)}
    {'ho', p0(8), bound_lower(8), bound_upper(8)}
    {'co', p0(9), bound_lower(9), bound_upper(9)}
    {'bo', p0(10), bound_lower(10), bound_upper(10)}
    {'s', p0(11), bound_lower(11), bound_upper(11)}
    {'sigma_y', p0(12), bound_lower(12), bound_upper(12)}
};

%%
%-------------- Run the MCMC chain --------------%
% The |model| structure holds information about the model. Minimally
% we need to set |ssfun|.

model.ssfun  = ssfun;

% The |options| structure has settings for the MCMC run. We need at
% least the number of simulations in |nsimu|.

options.nsimu = nsimu;
options.burnintime = 0.01*nsimu;

% Generate the MCMC chain. The actual MCMC simulation run is done using the function
% |mcmcrun|.
[results,chain] = mcmcrun(model,data,params,options);

save('DAIS_MCMCchain', 'chain')

% R syntax below
% #Run the MCMC chain
% dais.out.heter = metrop(log.post, p0, nbatch=NI, scale=step)
% dh.chain = dais.out.heter$batch
% dais.out.heter$accept
% #Calculate the parameter acceptance rate
% acceptrate = dais.out.heter$accept * 100
% #Echo the acceptance rate. Should be ~ 25%
% cat("Accept rate =", acceptrate, "%\n")

%--------------------  Analysis of the MCMC chain produced --------------------%
%% Plot some figures with the chain:
%figure(1); clf
%mcmcplot(chain,[],results.names,'chainpanel')

%figure(2); clf
%mcmcplot(chain,[],results.names,'dens')

%figure(3); clf
%mcmcplot(chain,[],results.names,'hist')

%figure(4); clf
%plot(chain(options.burnintime:nsimu,1))

% Test for MCMC chain convergence in R:
% library(coda)
% conv = results[burnin:NI,]
% heidel.diag(results, eps=0.1, pvalue=0.05)
% geweke.diag(conv, frac1=0.1, frac2=0.5)
% raftery.diag(conv, q=0.025, r=0.005, s=0.95, converge.eps=0.001)

%# To check if subset is sufficient:
%They should be roughly similiar
%par(mfrow=c(3,2))
%plot(density(results[(NI/2):NI,1]), main="gamma", xlab="")
%lines(density(sschain[,1]), col="red")
%plot(density(results[(NI/2):NI,2]), main="alpha", xlab="")
%lines(density(sschain[,2]), col="red")
%plot(density(results[(NI/2):NI,3]), main="mu", xlab="")
%lines(density(sschain[,3]), col="red")
%plot(density(results[(NI/2):NI,4]), main="eta", xlab="")
%lines(density(sschain[,4]), col="red")
%plot(density(results[(NI/2):NI,5]), main="Po", xlab="")
%lines(density(sschain[,5]), col="red")
%plot(density(results[(NI/2):NI,6]), main="kappa", xlab="")
%lines(density(sschain[,6]), col="red")
%plot(density(results[(NI/2):NI,7]), main="fo.y", xlab="")
%lines(density(sschain[,7]), col="red")
%plot(density(results[(NI/2):NI,8]), main="ho", xlab="")
%lines(density(sschain[,8]), col="red")

%Calculate the new best estimates from the mcmc chain
% NOTE: results may be in different format aka parameters in rows instead of columns
mean_dais_par = [mean(chain(:,1)); mean(chain(:,2)); mean(chain(:,3)); mean(chain(:,4));
                 mean(chain(:,5)); mean(chain(:,6)); mean(chain(:,7)); mean(chain(:,8));
                 mean(chain(:,9)); mean(chain(:,10)); mean(chain(:,11))]

%Create New mean estimate hindcast
new_dais.mcmc_est =  DAIS_IceFlux_model(mean_dais_par, hindcast_forcings, standards, endinfo);

%Create New best projection
endsat = 240298;
enddate = 240300;
endinfo = [endsat, enddate];
new_dais_mcmc_proj =  DAIS_IceFlux_model(mean_dais_par, project_forcings, standards, endinfo);

%save.image(file = "Workspace/DAIS_MCMC_calibration.RData")
% Estimating with all 10 million runs is not neccasary if the chains have
% converged. A subset of every 500th number should be sufficient.
subset_N = (nsimu-(nsimu*0.01))/2000; % #calculate the every nth number to get a subset of 2000
%subset_N = nsimu/2000;
R_subset = round(subset_N,0)-1;
sschain = chain(options.burnintime:R_subset:nsimu,:);
%sschain = chain(1:R_subset:nsimu,:);
subset_length = length(sschain(:,1));

%--------------------  PROJECT AIS SEA-LEVEL equivalence --------------------%
%Calculate all possible hindcasts from the subset parameter estimates.
par_mcmc(:,1)=sschain(:,1);
par_mcmc(:,2)=sschain(:,2);
par_mcmc(:,3)=sschain(:,3);
par_mcmc(:,4)=sschain(:,4);
par_mcmc(:,5)=sschain(:,5);
par_mcmc(:,6)=sschain(:,6);
par_mcmc(:,7)=sschain(:,7);
par_mcmc(:,8)=sschain(:,8);
par_mcmc(:,9)=sschain(:,9);
par_mcmc(:,10)=sschain(:,10);
par_mcmc(:,11)=sschain(:,11);

%## Superimpose the bias onto the model
%## True world = model + bias + error
bias_mcmc = sschain(:,12);

proj_dais_fits(1:subset_length, 1:enddate) = NaN;
for i=1:subset_length
proj_dais_fits(i,:) = DAIS_IceFlux_model(par_mcmc(i,:), project_forcings, standards, endinfo);
end

%## Superimpose the bias onto the model
proj_mcmc_w_bias(1:subset_length, 1:enddate) = NaN;
for i=1:subset_length
proj_mcmc_w_bias(i,:) = proj_dais_fits(i,:) + bias_mcmc(i)^2;
end

proj_mcmc_1961_1990(1:subset_length, 1:enddate) = NaN;
for i=1:subset_length
proj_mcmc_1961_1990(i,:) = proj_mcmc_w_bias(i,:)-mean(proj_mcmc_w_bias(i,SL_1961_1990));
end

%------------------ Save the workspace --------------------------------%
%save('DAIS_MCMCchain', 'chain')
save('DAIS_Matlab_MCMCcalibration','-v7.3') %Save the workspace into a R readable file

% The rest of the analysis will be run in R
%----------------------------------------------------------------------%

