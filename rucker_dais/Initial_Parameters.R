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
#These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
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
sigma = sd(resid)^2 #calculate the variance (sigma^2)

############################## RUN OPTIMIZATION #######################################
#Set up the priors; the upper and lower bounds
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
bound.lower[1:2] = c(1/2, 0)   ; bound.upper[1:2] = c(17/4, 1) #Set bounds for gamma and alpha
bound.lower[10:11] = c(725, 0.00045)   ; bound.upper[10:11] = c(825, 0.00075) #Set bounds for bo and s

# bound.lower[12] = 0 ; bound.upper[12] = 1 # Prior uniform range for sigma (the variance)

# step 2 define number of model parameters
model.p=11
parnames=c("gamma", "alpha", "mu", "eta", "po", "kappa", "fo", "ho", "co", "bo", "s", "sigma.y")
# step 3 source the physical model and statistical model
source("Scripts/DAIS_IceFlux_model.R")
source("Scripts/DAISobs_likelihood_iid.R")

#Shaffer [2014] best guess parameters
p = c(IP, sigma)
p0 = c(2.1, 0.29, 8, 0.015, 0.4, 0.04, 1.0, 1450, 90, 770, 0.0005, sigma-0.1)
p0 = optim(p0, function(p) - log.post(p))$par
print(round(p0,4))

write.csv(p0, file = "DAIS_matlab/OtimizedInitialParameters.csv")

############################## END #######################################
