#########################################################
#------- Now let's calibrate the model parameters
#------- using AR model
###-----HETEROSKEDASTIC
#########################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)

#Set the seed
set.seed(1234)

# step 1 define the boundary for parameters
source("Data/DAIS_data.R")

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

project.forcings = matrix(c(TA,TO,GSL,SL), ncol=4, nrow=240300)
hindcast.forcings = matrix(c(TA[1:240010], TO[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

# Best Case (Case #4) from Shaffer (2014)
IP = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Source the function with the standards and the initial parameters (IP) to
#get the best estimated AIS volume loss with respect to the present day in sea level equivalence (SLE):
standards = c(Tice,eps1, del, eps2, TOo, Volo, Roa, R)
source("Scripts/DAIS_IceFlux_model.R")
AIS_melt = iceflux(IP, hindcast.forcings, standards)

#set the end dates to the year 2300 to get future projections
end = 240298
enddate = 240300
Project_melt = iceflux(IP, project.forcings, standards)

#Set the end dates back to the hindcast period:
end = 240000
enddate = 240010

############################## CALCULATE RESIDUALS (PRIOR SIGMA) ##############################
# For this model the residuals are based off of the windowing approach and are weighted
#These windows are presented in Shaffer (2014) and calculated from figure 5 in Shepherd et al. (2012)
# ice_sat_data = read.csv("shepherd_AIS_observations.csv", skip=1)
# # Calculate the uncertainty with the +/- 2 standard error
# positive_2SE = ice_sat_data[1,12] + ice_sat_data[1,14] # Add the 2 standard error to the mean value
# negative_2SE = ice_sat_data[1,12] - ice_sat_data[1,14] # Subtract the 2 standard error to the mean value

# upper.wind = c(6, -8, -2, positive_2SE )
# lower.wind = c(2.5, -17, -4, negative_2SE)
# windows = matrix(c(lower.wind, upper.wind), nrow = 4, ncol=2)

# obs.errs = c(abs(median(windows[1,])-windows[1,1]), abs(median(windows[2,])-windows[2,1]),
#              abs(median(windows[3,])-windows[3,1]), ice_sat_data[1,14])
# 
# # Create a vector with each observation year
# obs.years = c(120000, 220000, 234000, 240002)

#1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
# We want the cumulative sea-level equivalent in meters for the year 2002
# 360Gt = 1mm SLE
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002-1992
mid.cum.SLE_2002 = estimate.SLE.rate*time.years

estimate.SLE.error = abs(-53/360)/1000 #1- sigma error
SE2_2002 = estimate.SLE.error*2 #2-sigma error

positive_2SE = mid.cum.SLE_2002 + SE2_2002 # Add the 2 standard error to the mean value
negative_2SE = mid.cum.SLE_2002 - SE2_2002 # Subtract the 2 standard error to the mean value

upper.wind = c(6.0, -6.9, -1.25, positive_2SE )
lower.wind = c(1.8, -15.8, -4.0, negative_2SE)
windows = matrix(c(lower.wind, upper.wind), nrow = 4, ncol=2)

obs.errs = c(abs(median(windows[1,])-windows[1,1]), abs(median(windows[2,])-windows[2,1]),
             abs(median(windows[3,])-windows[3,1]), SE2_2002)

# Create a vector with each observation year
#120kyr, 20Kyr, 6kyr, 2002
obs.years = c(120000, 220000, 234000, 240002)

#Set up equation to find the residuals and then the prior sigma
resid <- rep(NA,length(obs.years)) #Create a vector of the residuals
for(i in 1:length(obs.years)){
    resid[i] <- (median(windows[i,])-(AIS_melt[obs.years[i]]-mean(AIS_melt[SL.1961_1990]))) #/sd(windows[i,])
	}
sigma = sd(resid) #calculate the standard deviation (sigma)

############################## RUN MCMC #######################################
#Set up the priors and ranges for the MCMC
#Set the upper and lower bounds
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
bound.lower[1:2] = c(1/2, 0)   ; bound.upper[1:2] = c(17/4, 1) #Set bounds for gamma and alpha
bound.lower[10:11] = c(725, 0.00045)   ; bound.upper[10:11] = c(825, 0.00075) #Set bounds for bo and s

bound.lower[12] = 0 ; bound.upper[12] = 1 # Prior uniform range for sigma (the variance)

# step 2 define number of model parameters
model.p=11
parnames=c("gamma", "alpha", "mu", "eta", "po", "kappa", "fo", "ho", "co", "bo", "s", "sigma.y")
# step 3 source the physical model and statistical model
source("Scripts/DAIS_IceFlux_model.R")
source("Scripts/DAISobs_likelihood_iid.R")

#Shaffer [2014] best guess parameters
p = c(IP, sigma)
p0 = c(2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, 0.6)
p0 = optim(p0, function(p) - log.post(p))$par
print(round(p0,4))
#install.packages("mcmc")
library(mcmc)

step = c(0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15)/100
#step = c(0.001, 0.0001, 0.001, 0.00001, 0.0001, 0.00001, 0.001, 0.5, 0.1, 0.5, 0.000001, 0.001)
# NI = 900
NI = 2.2E6 #number of iterations
n.cores = 64 #number of cores #needs to be same as specified in shell #use only half of the cores in the node
burnin = seq(1, 0.01*NI, 1)
library(doParallel)
registerDoParallel(cores = n.cores)

results = foreach(i=1:n.cores, .combine=rbind) %dopar%{
  dais.out.heter=metrop(log.post, p0, nbatch=NI/n.cores, scale=step)
  dais.out.heter$batch
}

# #Run the MCMC chain
# dais.out.heter = metrop(log.post, p0, nbatch=NI, scale=step)
# dh.chain = dais.out.heter$batch
# dais.out.heter$accept
# #Calculate the parameter acceptance rate
# acceptrate = dais.out.heter$accept * 100
# #Echo the acceptance rate. Should be ~ 25%
# cat("Accept rate =", acceptrate, "%\n")

#If you don't want to use a burnin then uncomment out the following lines and these codes will do the same:
#mcmc.out2 = metrop(mcmc.out1, nbatch=1e5, scale=proposal.matrix(prechain1,mult=0.5))
#prechain2 = mcmc.out2$batch
#mcmc.out2$accept
#mcmc.out= metrop(mcmc.out2, nbatch=1e6, scale=proposal.matrix(prechain2,mult=0.5))
#chain = mcmc.out$batch
#mcmc.out$accept

########################## Analysis of the MCMC chain produced  #################
# Test for MCMC chain convergence:
# library(coda)
# conv = results[burnin:NI,]
# heidel.diag(results, eps=0.1, pvalue=0.05)
# geweke.diag(conv, frac1=0.1, frac2=0.5)
# raftery.diag(conv, q=0.025, r=0.005, s=0.95, converge.eps=0.001)

## To check if subset is sufficient: 
#They should be roughly similiar
#par(mfrow=c(3,2))
#plot(density(results[(NI/2):NI,1]), main="gamma", xlab="")
#lines(density(sschain[,1]), col="red")
#plot(density(results[(NI/2):NI,2]), main="alpha", xlab="")
#lines(density(sschain[,2]), col="red")
#plot(density(results[(NI/2):NI,3]), main="mu", xlab="")
#lines(density(sschain[,3]), col="red")
#plot(density(results[(NI/2):NI,4]), main="eta", xlab="")
#lines(density(sschain[,4]), col="red")
#plot(density(results[(NI/2):NI,5]), main="Po", xlab="")
#lines(density(sschain[,5]), col="red")
#plot(density(results[(NI/2):NI,6]), main="kappa", xlab="")
#lines(density(sschain[,6]), col="red")
#plot(density(results[(NI/2):NI,7]), main="fo.y", xlab="")
#lines(density(sschain[,7]), col="red")
#plot(density(results[(NI/2):NI,8]), main="ho", xlab="")
#lines(density(sschain[,8]), col="red")

#Calculate the new best estimates from the mcmc chain
# results = dh.chain
mean.dais.par = c(mean(results[,1]), mean(results[,2]), mean(results[,3]), mean(results[,4]),
                  mean(results[,5]), mean(results[,6]), mean(results[,7]), mean(results[,8]),
                  mean(results[,9]), mean(results[,10]), mean(results[,11]))
print(mean.dais.par)

##Calculate the probability density of each parameter estimated from MCMC:
#save.image(file = "DAIS_convergence_v4.RData")

# Function to calculate the pdfs of each parameter 'parameter.pdfs'
source('Scripts/plot_PdfCdfSf.R')
ais.mcmc.parPDF = parameter.pdfs(results)
# ais.pdfgamma <- density(results[,1]) ; ais.pdfalpha <- density(results[,2])
# ais.pdfmu <- density(results[,3]) ; ais.pdfeta <- density(results[,4])
# ais.pdfpo <- density(results[,5]) ; ais.pdfkappa <- density(results[,6])
# ais.pdffo <- density(results[,7]) ; ais.pdfho <- density(results[,8])
# ais.pdfco <- density(results[,9]) ; ais.pdfsigma <- density(results[,10])

#Create New mean estimate hindcast
standards = c(Tice,eps1, del, eps2, TOo, Volo, Roa, R)
new.dais.mcmc.est = iceflux(mean.dais.par, hindcast.forcings, standards)

#Create New best projection
end = 240298
enddate = 240300
new.dais.mcmc.proj = iceflux(mean.dais.par, project.forcings, standards)

save.image(file = "Workspace/DAIS_MCMC_calibration.RData")
# Estimating with all 10 million runs is not neccasary if the chains have
# converged. A subset of every 500th number should be sufficient.
subset_N = (NI-(NI*0.01))/2000 #calculate the every nth number to get a subset of 2000
sschain = results[seq(length(burnin), length(results[,1]), subset_N),]
subset_length = length(sschain[,1])

########################## HINDCAST AIS SEA-LEVEL equivalence ################################
#Lets calculate all possible hindcasts from the subset parameter estimates.
par.mcmc=mat.or.vec(subset_length, 11)
for(i in 1:subset_length) {
  par.mcmc[i,1]=sschain[i,1]
  par.mcmc[i,2]=sschain[i,2]
  par.mcmc[i,3]=sschain[i,3]
  par.mcmc[i,4]=sschain[i,4]
  par.mcmc[i,5]=sschain[i,5]
  par.mcmc[i,6]=sschain[i,6]
  par.mcmc[i,7]=sschain[i,7]
  par.mcmc[i,8]=sschain[i,8]
  par.mcmc[i,9]=sschain[i,9]
  par.mcmc[i,10]=sschain[i,10]
  par.mcmc[i,11]=sschain[i,11]
}

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
bias.mcmc = sschain[,12]
# dais.mcmc = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#     dais.mcmc[i,] = hetero.dais.fit[i,] + rnorm(1,mean=0,sd=bias[i])
# }
# 
# dais.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#   dais.mcmc.1961_1990[i,] = dais.mcmc[i,]-mean(dais.mcmc[i,SL.1961_1990])
# }
################### PROJECT AIS CONTRIBUTIONS ##########################
#Set the end dates to the projection period:
end = 240298
enddate = 240300

proj.dais.fits = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length) {
    proj.dais.fits[i,] = iceflux(par.mcmc[i,], project.forcings, standards)
}

### Superimpose the bias onto the model
proj.mcmc.w.bias = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length){
    proj.mcmc.w.bias[i,] = proj.dais.fits[i,] + rnorm(1,mean=0,sd=bias.mcmc[i])
}

proj.mcmc.1961_1990 = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length){
  proj.mcmc.1961_1990[i,] = proj.mcmc.w.bias[i,]-mean(proj.mcmc.w.bias[i,SL.1961_1990])
}
#------------------ Save the workspace --------------------------------#
#save.image(file = "Workspace/DAIS_MCMC_calibration.RData")
#--------------------- Estimate PDFs, CDFin 2100 & 2050 --------------------------
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

# #Plot the probability distribution of sea-level estimates in 2300
# ais.prob_proj=mat.or.vec(k,1)
# ais.prob_proj=proj.mcmc.1961_1990[,enddate]
# ais.pdf <- density(ais.prob_proj)
# ais.cdf = ecdf(ais.prob_proj)
# ais.cmcdf=ecdf(ais.prob_proj*100)
# ais.survival=1-ais.cmcdf(1:k) 
# 
# #Plot the probability distribution of sea-level estimates in 2100
# ais.prob_proj2100=mat.or.vec(k,1)
# ais.prob_proj2100=proj.mcmc.1961_1990[,240100] #The year 2100 is the 240100 number in the sequence
# ais.pdf2100 <- density(ais.prob_proj2100)
# ais.cdf2100 = ecdf(ais.prob_proj2100)
# ais.cmcdf2100=ecdf(ais.prob_proj2100*100)
# ais.survival2100=1-ais.cmcdf2100(1:k) 
# 
# #Plot the probability distribution of sea-level estimates in 2050
# ais.half_prob_proj=mat.or.vec(k,1)
# ais.half_prob_proj=proj.mcmc.1961_1990[,240050] #The year 2050 is the 240050 number in the sequence
# ais.half_pdf <- density(ais.half_prob_proj)
# ais.half.cdf = ecdf(ais.half_prob_proj)
# ais.cmhalf.cdf=ecdf(ais.half_prob_proj*100)
# ais.survival.half=1-ais.cmhalf.cdf(1:k) 

#Lets find out how the parameter relationships
#But lets also add in the initial value to the mix
#set up a matrix
d.pos_parameters = sschain
colnames(d.pos_parameters, do.NULL = FALSE)
colnames(d.pos_parameters) = c("gamma", "alpha", "mu", "eta", "po", "kappa", "fo", "ho", "co","bo", "s", "sigma.y")

save.image(file = "Workspace/DAIS_MCMC_calibration.RData")

###Calculate all possible residuals and find the sum of the squared error (SSE)
#d.all.residuals = mat.or.vec(k,nyears.obs)
#windowa = c(2.5,6)
#windowb = c(-8,-17)
#windowc = c(-2,-4)
#windowd = c(0.002,0.006)
#LI = (4.25-(AIS_melt$SLE[120000]-mean(AIS_melt$SLE[SL.1961_1990])))/sd(windowa)
#LG = (-12.5-(AIS_melt$SLE[220000]-mean(AIS_melt$SLE[SL.1961_1990])))/sd(windowb)
#H = (-3-(AIS_melt$SLE[234000]-mean(AIS_melt$SLE[SL.1961_1990])))/sd(windowc)
#pres = (0.004-(AIS_melt$SLE[240009]-mean(AIS_melt$SLE[SL.1961_1990])))/sd(windowd)
#resid = c(LI, LG, H, pres)
#d.sse = rep(NA,k)
#for(i in 1:k){
#  d.all.residuals[i,] = slr - homo.mcmc.fit[i,]
#  h.sse[i] = -sum(h.all.residuals[i,]^2)
#}
######################## For Further Analysis Plot graphes ######################
#source("heteroskedastic_dais_plots.R")

