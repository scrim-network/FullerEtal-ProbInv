#################################################################################
#
#  - file = "DAIScali_hetero_model_iid_mcmcRversion.R"
#  - Code written: August 2015, updated March 2016
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program runs a Markov Chain Monte Carlo analysis of the DAIS model
#       assuming heteroskedastic errors and IID residuals as described in Ruckert et al. (2016).
#       For further description and references, please read the paper.
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOT IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL,
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
###################################################################################
rm(list =ls()) # Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)
library(mcmc)

#Set the seed.
set.seed(1)

## Run multiple times with different seeds to check for convergence &
## robustness of the results:
# set.seed(1234)
# set.seed(1)

# Read in hindcast/ forcing data, read in standard values, and read in AIS dates both specific and ranges so there is no use of magic numbers.
source("Data/DAIS_data_C.R")

############################## List of physcial model parameters ##############################
# We will set the initial parameters to specifications from Shaffer [2014]
#Define the parameters:
# [1] gamma = 2 				#sensitivity of ice flow to sea level
# [2] alpha = 0.35 			#sensitivity of ice flow to ocean subsurface temperature
# [3] mu = 8.7    			#Profile parameter related to ice stress [m^(1/2)]
# [4] nu = 0.012   			#Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
# [5] P0 = 0.35    			#Precipitation at 0C [m of ice/yr]
# [6] kappa = 0.04 			#Relates precipitation to temperature [K^-1]
# [7] f0 = 1.2          #Constant of proportionality for ice speed [m/yr]
# [8] h0 = 1471         #Initial value for runoff line calculation [m]
# [9] c = 95            #Second value for runoff line calculation [m]
# [10] b0 = 775         #Height of bed at the center of the continent [m]
# [11] slope = 0.0006   #Slope of the bed

# Create matrices for projections and hindcasts.
project.forcings = matrix(c(Ta, Toc, GSL, SL), ncol=4, nrow=240300)
hindcast.forcings = matrix(c(Ta[1:240010], Toc[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

# SEt initial parameters to the best case (Case #4) from Shaffer (2014)
IP = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Source the function with the standards and the initial parameters (IP) to
#get the best estimated AIS volume loss 

# Load in the physical model. Calls C model by default and set vector of standard values:
source("models.R") 
standards = c(Tf, rho_w, rho_i, rho_m, Toc_0, Rad0)

# Estimate the AIS volume loss hindcast for Case #4 with respect to the present day in sea level equivalence (SLE):
AIS_melt = iceflux(IP, hindcast.forcings, standards)

# Estimate the AIS volume loss projection for Case #4 with respect to the present day in sea level equivalence (SLE):
Project_melt = iceflux(IP, project.forcings, standards)

############################## Setup observational constraint info ##############################
# For this model the residuals are based off of the windowing approach.
# These windows are presented in Shaffer (2014) and calculated from figure 5 in Shepherd et al. (2012).

# Accummulate the sea-level equivalent in meters from 1992 to the year 2002
# using the 1992 to 2011 trend from Shepherd et al. 2012; -71 +/- 53 Gt per yr.
# Conversion: 360 Gt = 1 mm SLE
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002 - 1992
mid.cum.SLE_2002 = estimate.SLE.rate * time.years

# Accumulate the 1 sigma error for the year 2002 and estimate the 2-sigma error:
estimate.SLE.error = abs(-53/360)/1000 #1- sigma error
estimate.SLE.error = sqrt(time.years)*abs(-53/360)/1000 # 1-sigma error
# (*sqrt(10) because 10 years of potentially accumulated error:
#  total error^2 = year 1 error^2 + year 2 error^2 + ... year 10 error^2
#                = 10*year X error^2)
SE2_2002 = estimate.SLE.error*2 # 2-sigma error

# Add and subtract the 2 standard error to the mean value
positive_2SE = mid.cum.SLE_2002 + SE2_2002 
negative_2SE = mid.cum.SLE_2002 - SE2_2002

# Create observational constraint windows.
upper.wind = c(6.0, -6.9, -1.25, positive_2SE )
lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
windows = matrix(c(lower.wind, upper.wind), nrow = 4, ncol=2)

# Determine observational error from windows: half-width of window = uncertainty; assume all windows are 2*stdErr (last one actually is)
obs.errs = (windows[,2]-windows[,1])*.5
 
# Create a vector with each observation year.
#            120 kyr, 20 Kyr,  6 kyr,  2002
obs.years = c(120000, 220000, 234000, 240002)

############################## CALCULATE RESIDUALS and INITIAL SIGMA VALUE ##############################
resid <- rep(NA,length(obs.years))

# Estimate the residuals: modification from equation (1)
for(i in 1:length(obs.years)){
    resid[i] <- (median(windows[i,]) - (AIS_melt[obs.years[i]] - mean(AIS_melt[SL.1961_1990])))
}

# Calculate the variance, sigma^2
variance = sd(resid)^2

############################## RUN MCMC CALIBRATION #######################################
# Set up priors.
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
print(bound.lower)
print(bound.upper)

# var.y has inverse gamma prior, so there is a lower bound at 0 but no upper bound
parnames    = c('gamma','alpha','mu'  ,'nu'  ,'P0' ,'kappa','f0' ,'h0'  ,'c'  , 'b0','slope' ,'var.y')
bound.upper = c( 4.25 ,  1     , 13.05, 0.018,0.525,  0.06 , 1.8 ,2206.5, 142.5, 825 , 0.00075,   Inf)
bound.lower = c( 0.5  ,  0     , 4.35 , 0.006,0.175,  0.02 , 0.6 , 735.5,  47.5, 725 , 0.00045 ,    0)

# Specify the number of model parameters.
# Variance is a statistical parameter and is not counted in the number.
model.p=11

# Load the likelihood model assuming heteroskedastic observation errors and non-correlated residuals.
source("Scripts/DAISobs_likelihood_iid.R")

# Optimize the likelihood function to estimate initial starting values.
p = c(IP, variance) # Shaffer [2014] Case #4 parameters
p0 = c(2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.6) # Random guesses
p0 = optim(p0, function(p) - log.post(p))$par
print(round(p0,4))

# Set the step size and number of iterations.
step = c(0.1, 0.015, 0.2, 0.035, 0.1, 0.01, 0.1, 50, 10, 25, 0.0005, 0.1)/5
#step = c(0.05, 0.01, 0.15, 0.035, 0.1, 0.01, 0.1, 50, 10, 30, 0.0005, 0.1)/5
#step = c(0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15)/5
# NI = 5E4
NI = 4E6

# Run MCMC calibration.
DAIS_mcmc_output1780 = metrop(log.post, p0, nbatch=NI, scale=step)
DAIS_chains1780 = DAIS_mcmc_output1780$batch

# Print the acceptance rate as a percent. Should be ~ 25%
acceptrate = DAIS_mcmc_output1780$accept * 100
cat("Accept rate =", acceptrate, "%\n")

## Save workspace image - you do not want to re-simulate all those!
save.image(file = "Scratch/Workspace/DAIS_calib_MCMC_C1.RData")

############################## RUN ADAPTIVE MCMC #######################################
# # NOT FULLY CODED YET
# # load robust adaptive Metropolis
# library(adaptMCMC)
#
# # Set optimal acceptance rate as # parameters->infinity (Gelman et al, 1996; Roberts et al, 1997)
# accept.mcmc = 0.234
#
# # Set number of iterations, burnin, and rate of adaptation (between 0.5 and 1, lower is faster adaptation)
# niter.mcmc = 5E4											# number of iterations for MCMC
# gamma.mcmc = 0.5											#
# burnin = 0.01*NI
#
# # Specify when to stop adapting (niter*1 => don't stop) and step size
# # stopadapt.mcmc = round(niter.mcmc*1.0)
# # step.mcmc = (bound.upper-bound.lower)*.05
#
# ## Actually run the calibration.
# amcmc.out1 = MCMC(log.post, niter.mcmc, parameters0, scale=step.mcmc, adapt=TRUE, acc.rate=accept.mcmc,gamma=gamma.mcmc, list=TRUE,
#                   n.start=round(0.01*niter.mcmc),parnames.in=parnames, bound.lower.in=bound.lower, bound.upper.in=bound.upper,
#                                 obs.in=obs.targets  , obs.err.in=obs.err        , obs.step.in=obs.years     ,
#                                 alpha.var=alpha.var , beta.var=beta.var         ,
#                                 slope.Ta2Tg.in=slope.Ta2Tg , intercept.Ta2Tg.in=intercept.Ta2Tg ,
#                                 Tg.in=Tg.recon      , SL.in=SL )
#
# DAIS_chains = amcmc.out1$samples
#
########################## Analysis of the mean of MCMC chains produced  #################
# Identify the burn-in period and subtract it from the chains.
#burnin = seq(1, 0.01*NI, 1) # 1% burnin
#DAIS_chains_burnin <- DAIS_chains[-burnin,]

# Find mean of the estimated parameters.
#mean.parameters = c(mean(DAIS_chains_burnin[,1]), mean(DAIS_chains_burnin[,2]), mean(DAIS_chains_burnin[,3]), mean(DAIS_chains_burnin[,4]),
#                  mean(DAIS_chains_burnin[,5]), mean(DAIS_chains_burnin[,6]), mean(DAIS_chains_burnin[,7]), mean(DAIS_chains_burnin[,8]),
#                  mean(DAIS_chains_burnin[,9]), mean(DAIS_chains_burnin[,10]), mean(DAIS_chains_burnin[,11]), mean(DAIS_chains_burnin[,12]))
#print(mean.parameters)

# Estimate mean bias:
#bias.mean = sqrt(mean.parameters[12])

# Estimate model hindcast from parameter means.
#mean.hind.simulation = iceflux(mean.parameters[1:11], hindcast.forcings, standards)
#mean.hind.anomaly = mean.hind.simulation - mean(mean.hind.simulation[SL.1961_1990])
#mean.hindcast = mean.hind.anomaly + rnorm(length(mean.hind.simulation), mean=0, sd=bias.mean)

# Estimate model projection from parameter means.
#mean.proj.simulation = iceflux(mean.parameters[1:11], project.forcings, standards)
#mean.proj.anomaly = mean.proj.simulation - mean(mean.proj.simulation[SL.1961_1990])
#mean.projection = mean.proj.anomaly + rnorm(length(mean.proj.anomaly), mean=0, sd=bias.mean)

############################## ANALYZE CHAINS #######################################
# Load function to calculate pdfs of each DAIS parameter 'parameter.pdfs'
#source('Scripts/plot_PdfCdfSf.R')

# Find the probability density function for each of the estimated parameters
#dais_parameter_PDFs = parameter.pdfs(DAIS_chains_burnin)

# Thin the chain to a subset; ~2,500 is sufficient.
#subset_N = NI/5000
#R_subset = round(subset_N,0)
#sub_chain = DAIS_chains_burnin[seq(1, length(DAIS_chains_burnin[,1]), R_subset),]
#subset_length = length(sub_chain[,1])

# Check for simularities between full chain and the subset.
#par(mfrow=c(4,3))
#for(i in 1:12){
#  plot(density(DAIS_chains_burnin[ ,i]), type="l",
#       xlab=paste('Parameter =',' ', parnames[i], sep=''), ylab="PDF", main="")
#  lines(density(sub_chain[ ,i]), col="red")
#}

# Set up parameter matrix.
#par.mcmc = mat.or.vec(subset_length, 11)
#for(i in 1:subset_length) {
#  par.mcmc[i,1] = sub_chain[i,1]
#  par.mcmc[i,2] = sub_chain[i,2]
#  par.mcmc[i,3] = sub_chain[i,3]
#  par.mcmc[i,4] = sub_chain[i,4]
#  par.mcmc[i,5] = sub_chain[i,5]
#  par.mcmc[i,6] = sub_chain[i,6]
#  par.mcmc[i,7] = sub_chain[i,7]
#  par.mcmc[i,8] = sub_chain[i,8]
#  par.mcmc[i,9] = sub_chain[i,9]
#  par.mcmc[i,10] = sub_chain[i,10]
#  par.mcmc[i,11] = sub_chain[i,11]
#}

################### HINDCAST AIS CONTRIBUTIONS ##########################
# enddate = 240010

# # Loop over the physical model to generate a distribution of model simulations.
# hindcast_sim = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length) {
#   hindcast_sims[i,] = iceflux(par.mcmc[i,], hindcast.forcings, standards)
# }
#
# # Estimate model simulation anomalies to the mean 1961-1990 period.
# dais_hindcast_anomaly = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#     dais_hindcast_anomaly[i,] = hindcast_sims[i,] - mean(hindcast_sims[i,SL.1961_1990])
# }
# 
# # Estimate bias from the variance
# bias.mcmc = sqrt(sub_chain[,12]) # standard deviation
#
# # Estimate the hindcasts: add the residuals (bias) onto the model simulations. Equations 1 & 2
# # True world = model + bias + error
# hind.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#   hind.mcmc.1961_1990[i,] = dais_hindcast_anomaly[i,] + rnorm(enddate, mean = 0, sd = bias.mcmc[i])
# }
################### PROJECT AIS CONTRIBUTIONS ##########################
#enddate = 240300

# Loop over the physical model to generate a distribution of model simulations.
#projection_sims = mat.or.vec(subset_length, enddate)
#for(i in 1:subset_length) {
#  projection_sims[i,] = iceflux(par.mcmc[i,], project.forcings, standards)
#}

# # Estimate model simulation anomalies to the mean 1961-1990 period.
#proj.mcmc.anomaly = mat.or.vec(subset_length, enddate)
#for(i in 1:subset_length){
#  proj.mcmc.anomaly[i,] = projection_sims[i,] - mean(projection_sims[i,SL.1961_1990])
#}

# Estimate bias from the variance
#bias.mcmc = sqrt(sub_chain[,12]) # standard deviation

# # Estimate the projections: add the residuals (bias) onto the model simulations. Equations 1 & 2
# # True world = model + bias + error
#proj.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
#for(i in 1:subset_length){
#  proj.mcmc.1961_1990[i,] = proj.mcmc.anomaly[i,] + rnorm(enddate, mean = 0, sd = bias.mcmc[i])
#}

#--------------------- Estimate PDFs, CDFs, and SFs in certain years --------------------------
# Function to find SLE values in certain years 'fn.prob.proj'
#year.pcs = c(120000, 220000, 234000, 240002, 240050, 240100, 240300)

#mcmc.prob_proj <- fn.prob.proj(proj.mcmc.1961_1990, year.pcs, subset_length, un.constr=T)

# Calculate the pdf, cdf, and sf of AIS melt estimates in:
#LIG.sf.mcmc <- plot.sf(mcmc.prob_proj[,1], make.plot=F) # 120,000 BP (Last interglacial)
#LGM.sf.mcmc <- plot.sf(mcmc.prob_proj[,2], make.plot=F) # 20,000 BP (Last glacial maximum)
#MH.sf.mcmc <- plot.sf(mcmc.prob_proj[,3], make.plot=F) # 6,000 BP (Mid-holocene)
#sf.2002.mcmc <- plot.sf(mcmc.prob_proj[,4], make.plot=F) # 2002 (Observed trend from 1993-2011)
#sf.2050.mcmc <- plot.sf(mcmc.prob_proj[,5], make.plot=F) # 2050
#sf.2100.mcmc <- plot.sf(mcmc.prob_proj[,6], make.plot=F) # 2100
#sf.2300.mcmc <- plot.sf(mcmc.prob_proj[,7], make.plot=F) # 2300

#Find out parameter relationships; set up a matrix
#d.pos_parameters = sub_chain
#colnames(d.pos_parameters, do.NULL = FALSE)
#colnames(d.pos_parameters) = c("gamma", "alpha", "mu", "nu", "p0", "kappa", "f0", "h0", "c","b0", "slope", "sigma")

#save.image(file = "Scratch/Workspace/DAIS_MCMC_R_C_calibration1780.RData")

######################## For Further Analysis Plot graphes ######################
#source("MCMC_plots.R")

