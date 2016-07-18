#-----------------------------------------------------------------------------
#------- Now let's calibrate the model parameters
#------- using iid model
###-----HETEROSKEDASTIC
#-----------------------------------------------------------------------------

srand(1234);  # For reproducibility
#srand(1780);  # For reproducibility

# Install packages:
# Pkg.add("DataFrames")
# Pkg.add("Distributions")
# Pkg.add("DataFramesMeta")
# Pkg.add("Optim")
# Pkg.add("Lora")
using DataFrames
using Distributions
using DataFramesMeta
using Optim
using Lora

# step 1 define the boundary for parameters
include("Data/DAIS_data.jl")

#-------------------- SET INITIAL PARAMETERS --------------------#
# We will set the initial parameters to specifications from Shaffer [2014]
#Define the parameters:
# gamma = 2 			  #sensitivity of ice flow to sea level
# alpha = 0.35 			  #sensitivity of ice flow to ocean subsurface temperature
# mu = 8.7    			  #Profile parameter related to ice stress [m^(1/2)]
# eta = 0.012   		  #Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
# Po = 0.35    			  #Precipitation at 0C [m of ice/yr]
# kappa = 0.04 			  #Relates precipitation to temperature [K^-1]
# fo = 1.2                #Constant of proportionality for ice speed [m/yr]
# ho = 1471               #Initial value for runoff line calculation [m]
# co = 95                 #Second value for runoff line calculation [m]
# bo = 775               #Height of bed at the center of the continent [m]
# s = 0.0006             #Slope of the bed

#Create a matrix that has 4 columns and 240300 rows
project_forcings = [TA TO GSL SL];

# 240010x4 matrix
hindcast_forcings = [TA[1:240010] TO[1:240010] GSL[1:240010] SL[1:240010]];

# Best Case (Case #4) from Shaffer (2014)
IP = [2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006];

#Source the function with the standards and the initial parameters (IP) to
#get the best estimated AIS volume loss with respect to the present day in sea level equivalence (SLE):
standards = [Tice, eps1, del, eps2, TOo, Volo, Roa, R];
endinfo = [endsat, enddate];

include("Scripts/DAIS_IceFlux_model.jl")
AIS_melt = iceflux(IP, hindcast_forcings, standards, endinfo);

#set the end dates to the year 2300 to get future projections
endsat = 240298;
enddate = 240300;
endinfo = [endsat, enddate];
Project_melt = iceflux(IP, project_forcings, standards, endinfo);

#Set the end dates back to the hindcast period:
endsat = 240000;
enddate = 240010;
endinfo = [endsat, enddate];

#--------------------  CALCULATE RESIDUALS (PRIOR SIGMA) --------------------#
# These windows are presented in Shaffer (2014) and calculated from Shepherd et al. (2012)
# # Calculate the uncertainty with the +/- 2 standard error

# Create a vector with each observation year
# 120kyr, 20Kyr, 6kyr, 2002
obs_years = [120000, 220000, 234000, 240002];

# 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
# We want the cumulative sea-level equivalent in meters for the year 2002
# 360Gt = 1mm SLE
estimate_SLE_rate = abs(-71/360)/1000;
time_years = 2002-1992;
mid_cum_SLE_2002 = estimate_SLE_rate*time_years;

estimate_SLE_error = abs(-53/360)/1000; #1- sigma error
SE2_2002 = estimate_SLE_error*2; #2-sigma error

positive_2SE = mid_cum_SLE_2002 + SE2_2002; # Add the 2 standard error to the mean value
negative_2SE = mid_cum_SLE_2002 - SE2_2002; # Subtract the 2 standard error to the mean value

upper_wind = [6.0, -6.9, -1.25, positive_2SE];
lower_wind = [1.8, -15.8, -4.0, negative_2SE];
windows = [lower_wind upper_wind];

obs_errs = [abs(median(windows[1,:])-windows[1,1]); abs(median(windows[2,:])-windows[2,1]);
            abs(median(windows[3,:])-windows[3,1]); SE2_2002];

#Set up equation to find the residuals and then the prior sigma
resid = zeros(length(obs_years),1);   #Create a vector of the residuals

for i in 1:length(obs_years)
  resid[i] = (median(windows[i,:])-(AIS_melt[obs_years[i]]-mean(AIS_melt[SL_1961_1990])));
end

sigma = std(resid)^2; #calculate the variance (sigma^2);

##
#--------------------  SETUP MCMC --------------------#
#Set up the priors; the upper and lower bounds
bound_lower = IP - (IP*0.5)    ; bound_upper = IP + (IP*0.5);
bound_lower[1:2] = [1/2, 0]   ; bound_upper[1:2] = [17/4, 1]; #Set bounds for gamma and alpha
bound_lower[10:11] = [725, 0.00045]   ; bound_upper[10:11] = [825, 0.00075]; #Set bounds for bo and s

#bound_lower(12) = 0 ; bound_upper(12) = 1; % Prior uniform range for sigma (the variance)

#Shaffer [2014] best guess parameters
p_1 = [IP; sigma];
p0 = [2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.8];

# step 2 define number of model parameters
model_p=11
parnames= ["gamma", "alpha", "mu", "eta", "po", "kappa", "fo", "ho", "co", "bo", "s", "sigma_y"]
# step 3 source the physical model and statistical model

include("Scripts/DAIS_IceFlux_model.jl")
include("Scripts/DAISobs_likelihood_iid.jl")

#function f(p)
#  -log_post(p)
#end

#p0 = optimize(f, p);
#p0 = Optim.minimizer(p0);

#R"source("Scripts/DAIS_IceFlux_model.R")"
#R"source("Scripts/DAISobs_likelihood_iid.R")"
#R"optim($p0, function(p) - log.post(p))$par"

###########################
### Define the parameter via BasicContMuvParameter (it is a continuous multivariate variable)
### The input arguments for BasicContMuvParameter are:
### 1) the variable key,
### 2) the log-target

plogtarget(p::Vector{Float64}) = log_post(p)
p = BasicContMuvParameter(:p, logtarget=plogtarget)
model = likelihood_model(p, false)

step = [0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15]/100
nsteps = 100;
burnin = (0.01*nsteps)
sampler = MH(step)
mcrange = BasicMCRange(nsteps = 100, burnin = 1, thinning = 5)

v0 = Dict(:p=>p_1)

job = BasicMCJob(model, sampler, mcrange, v0, tuner=VanillaMCTuner(verbose=true))
run(job)

### Get simulated values
chain = output(job)
results = chain.value;
mean(chain)
###########################

plogtarget(z::Vector{Float64}) = -dot(z, z)
model = likelihood_model(BasicContMuvParameter(:p, logtarget=log_post), false)

state = BasicContMuvParameterState(Int64[2.1,0.29,8.0,0.015,0.4,0.04,1.0,1450.0,90,770.0,0.0005,0.8], [:accept], Float32, [false])



job = BasicMCJob(model, MH(step), BasicMCRange(nsteps, burnin), Dict(:p=>p0))


job = GibbsJob(model, MH(step), BasicMCRange(nsteps, burnin), Dict(:p=>p0))
mcmc_results = run(job)
results = output(job)

#step = c(0.001, 0.0001, 0.001, 0.00001, 0.0001, 0.00001, 0.001, 0.5, 0.1, 0.5, 0.000001, 0.001)
# NI = 900
NI = 1.2E6 #number of iterations
burnin = seq(1, 0.01*NI, 1)

dais.out.heter = metrop(log.post, p0, nbatch=NI, scale=step)
results = dais.out.heter$batch



% Set up the structure for likelihood function
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

%----------------- Estimate good initial parameters ------------%
%Initial parameters from optim in R
load OtimizedInitialParameters.txt;

minit(1,:) = OtimizedInitialParameters(:,2);
minit(12)=[1.16];
% minit=[2.2413, 0.3682, 8.0021, 0.0156, 0.4647, 0.0486, 1.2850, 1558.7492, 89.9974, 805.6928, 0.0006, 0.5782];
nsimu = 1000;  % number of simulations
step = minit/150;
thin=5;

%%
%-------------- Run the MCMC chain --------------%
rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

mmc=mcmc(minit,loglike,logmodelprior,step,nsimu,thin);
mmc1 = mmc;

%Estimate a more appropriate step size:
scale=proposal_matrix(mmc1,0.5);
step2(1:12) = diag(scale);

% New starting value:
minit2 = mmc1(length(mmc1),:);

rand('seed', 1234);  % For reproducibility
%rand('seed', 1780);  % For reproducibility

nsimu2 = 1.2e6;
mmc = mcmc(minit2, loglike, logmodelprior, step2, nsimu2, thin);
mmc2 = mmc;

save('DAIS_MCMCchain_1234', 'mmc1','mmc2')
%save('DAIS_MCMCchain_1780', 'mmc1','mmc2')

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

