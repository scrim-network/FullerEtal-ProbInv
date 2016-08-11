#########################################################
#------- Now let's calibrate the model parameters
#------- using AR model
###-----HETEROSKEDASTIC
#########################################################
rm(list =ls()) #Clear global environment
library(compiler)
library(R.matlab)
enableJIT(3)
enableJIT(3)

#Set the seed
set.seed(1234)

# step 1 define the boundary for parameters
source("Data/DAIS_data_C.R")

############################## SET INITIAL PARAMETERS ##############################
# We will set the initial parameters to specifications from Shaffer [2014]
#Define the parameters:
# [1] gamma = 2 				#sensitivity of ice flow to sea level
# [2] alpha = 0.35 			#sensitivity of ice flow to ocean subsurface temperature
# [3] mu = 8.7    				#Profile parameter related to ice stress [m^(1/2)]
# [4] eta = 0.012   			#Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
# [5] Po = 0.35    			#Precipitation at 0C [m of ice/yr] 
# [6] kappa = 0.04 			#Relates precipitation to temperature [K^-1]
# [7] fo = 1.2                #Constant of proportionality for ice speed [m/yr]
# [8] ho = 1471               #Initial value for runoff line calculation [m]
# [9] co = 95                 #Second value for runoff line calculation [m]
# [10] bo = 775                #Height of bed at the center of the continent [m]
# [11] s = 0.0006              #Slope of the bed

# Create matrices for projections and hindcasts.
project.forcings = matrix(c(Ta, Toc, GSL, SL), ncol=4, nrow=240300)
hindcast.forcings = matrix(c(Ta[1:240010], Toc[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

# SEt initial parameters to the best case (Case #4) from Shaffer (2014)
IP = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Source the function with the standards and the initial parameters (IP) to
#get the best estimated AIS volume loss

# Load in the physical model. Calls C model by default and set vector of standard values:
source("models.R")
standards = c(Tf, rho_w, rho_i, rho_m, Toc_0, Rad0) # Volo is specified in the model

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

############################## SETUP MCMC #######################################
############################## RUN MCMC CALIBRATION #######################################
# Set up priors.
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
print(bound.lower)             ; print(bound.upper)

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

############################## RUN MCMC #######################################
# MCMC calibration is run in matlab for speed. Read in the chains produced by the
# matlab code.
mat_chains = readMat("DAIS_matlab/DAIS_MCMCchain_1234.mat")
results = mat_chains$mmc2

NI = length(results[,1]) #number of iterations
burnin = (1.2e6*0.01)+1
results = results[burnin:NI,]
NI = length(results[,1])

# results = dh.chain
mean.parameters = c(mean(results[,1]), mean(results[,2]), mean(results[,3]), mean(results[,4]),
                  mean(results[,5]), mean(results[,6]), mean(results[,7]), mean(results[,8]),
                  mean(results[,9]), mean(results[,10]), mean(results[,11]), mean(results[,12]))
print(mean.parameters)

############################## ANALYZE CHAINS #######################################
# Load function to calculate pdfs of each DAIS parameter 'parameter.pdfs'
source('Scripts/plot_PdfCdfSf.R')

# Find the probability density function for each of the estimated parameters
dais_parameter_PDFs = parameter.pdfs(results)

# Estimate mean bias:
bias.mean = sqrt(mean.parameters[12])

# Estimate model hindcast from parameter means.
mean.hind.simulation = iceflux(mean.parameters[1:11], hindcast.forcings, standards)
mean.hind.anomaly = mean.hind.simulation - mean(mean.hind.simulation[SL.1961_1990])
mean.hindcast = mean.hind.anomaly + rnorm(length(mean.hind.simulation), mean=0, sd=bias.mean)

# Estimate model projection from parameter means.
mean.proj.simulation = iceflux(mean.parameters[1:11], project.forcings, standards)
mean.proj.anomaly = mean.proj.simulation - mean(mean.proj.simulation[SL.1961_1990])
mean.projection = mean.proj.anomaly + rnorm(length(mean.proj.anomaly), mean=0, sd=bias.mean)

# Estimating with all the runs is not neccasary if the chains have
# converged. A subset of 2000 iterations should be sufficient.
#subset_N = NI/2500 #thinthe chain by every nth number to get a subset of 2000
#R_subset = round(subset_N,0)
#sschain = results[seq(1, length(results[,1]), R_subset),]
#subset_length = length(sschain[,1])

subset_length = 3500
sub_chain = mat.or.vec(subset_length, 12)
for(i in 1:12){
    sub_chain[,i] = sample(results[,i], subset_length)
}

# Check for simularities between full chain and the subset.
pdf(file="simplecheck_matlab.pdf")
par(mfrow=c(4,3))
for(i in 1:12){
    plot(density(results[ ,i]), type="l",
    xlab=paste('Parameter =',' ', parnames[i], sep=''), ylab="PDF", main="")
    lines(density(sub_chain[ ,i]), col="red")
}
dev.off()

# Set up parameter matrix.
par.mcmc = mat.or.vec(subset_length, 11)
for(i in 1:subset_length) {
    par.mcmc[i,1] = sub_chain[i,1]
    par.mcmc[i,2] = sub_chain[i,2]
    par.mcmc[i,3] = sub_chain[i,3]
    par.mcmc[i,4] = sub_chain[i,4]
    par.mcmc[i,5] = sub_chain[i,5]
    par.mcmc[i,6] = sub_chain[i,6]
    par.mcmc[i,7] = sub_chain[i,7]
    par.mcmc[i,8] = sub_chain[i,8]
    par.mcmc[i,9] = sub_chain[i,9]
    par.mcmc[i,10] = sub_chain[i,10]
    par.mcmc[i,11] = sub_chain[i,11]
}

################### HINDCAST AIS CONTRIBUTIONS ##########################
#Set the end dates back to the hindcast period:
# end = 240000
# enddate = 240010
# 
# hetero.dais.fit=mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length) {
#   hetero.dais.fit[i,] = iceflux(par.mcmc[i,], hindcast.forcings, standards)
# }

### Superimpose the bias onto the model
### True world = model + bias + error
# bias.mcmc = sqrt(sschain[,12]) #standard deviation
# dais.mcmc = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#     dais.mcmc[i,] = hetero.dais.fit[i,] - mean(hetero.dais.fit[i,SL.1961_1990])
# }
# 
# dais.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#   dais.mcmc.1961_1990[i,] = dais.mcmc[i,] + rnorm(enddate,mean=0,sd=bias.mcmc[i])
# }
################### PROJECT AIS CONTRIBUTIONS ##########################
enddate = 240300

# Loop over the physical model to generate a distribution of model simulations.
projection_sims = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length) {
    projection_sims[i,] = iceflux(par.mcmc[i,], project.forcings, standards)
}

# # Estimate model simulation anomalies to the mean 1961-1990 period.
proj.mcmc.anomaly = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length){
    proj.mcmc.anomaly[i,] = projection_sims[i,] - mean(projection_sims[i,SL.1961_1990])
}

# Estimate bias from the variance
bias.mcmc = sqrt(sub_chain[,12]) # standard deviation

# # Estimate the projections: add the residuals (bias) onto the model simulations. Equations 1 & 2
# # True world = model + bias + error
proj.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length){
    proj.mcmc.1961_1990[i,] = proj.mcmc.anomaly[i,] + rnorm(enddate, mean = 0, sd = bias.mcmc[i])
}

#--------------------- Estimate PDFs, CDFs, and SFs in certain years --------------------------
# Function to find SLE values in certain years 'fn.prob.proj'
year.pcs = c(120000, 220000, 234000, 240002, 240050, 240100, 240300)

mcmc.prob_proj <- fn.prob.proj(proj.mcmc.1961_1990, year.pcs, subset_length, un.constr=T)

# Calculate the pdf, cdf, and sf of AIS melt estimates in:
LIG.sf.mcmc <- plot.sf(mcmc.prob_proj[,1], make.plot=F) # 120,000 BP (Last interglacial)
LGM.sf.mcmc <- plot.sf(mcmc.prob_proj[,2], make.plot=F) # 20,000 BP (Last glacial maximum)
MH.sf.mcmc <- plot.sf(mcmc.prob_proj[,3], make.plot=F) # 6,000 BP (Mid-holocene)
sf.2002.mcmc <- plot.sf(mcmc.prob_proj[,4], make.plot=F) # 2002 (Observed trend from 1993-2011)
sf.2050.mcmc <- plot.sf(mcmc.prob_proj[,5], make.plot=F) # 2050
sf.2100.mcmc <- plot.sf(mcmc.prob_proj[,6], make.plot=F) # 2100
sf.2300.mcmc <- plot.sf(mcmc.prob_proj[,7], make.plot=F) # 2300

#Find out parameter relationships; set up a matrix
d.pos_parameters = sub_chain
colnames(d.pos_parameters, do.NULL = FALSE)
colnames(d.pos_parameters) = c("gamma", "alpha", "mu", "nu", "p0", "kappa", "f0", "h0", "c","b0", "slope", "variance")

save.image(file = "Scratch/Workspace/DAIS_MCMC_Matlabcalibration_1234_C.RData")

######################## For Further Analysis Plot graphes ######################
#source("heteroskedastic_dais_plots.R")

