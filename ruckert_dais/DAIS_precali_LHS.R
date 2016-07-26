###################################
# file: DAIS_precali_LHS.R
###################################
# Author and copyright: Kelsey Ruckert
# Pennsylvania State University
# klr324@psu.edu
###################################
# Latin Hypercube Sampling of parameters
# generates data frame with PDFs
# This LHS precalibrates the parameters used in the DAIS Model
###################################

rm(list =ls()) #Clear global environment

# Setup -------------------------------------------------------------------
#install.packages('lhs')
#install.packages('compiler')
#install.packages('pscl')
require(lhs)
library(pscl) # install inverse gamma distribution
library(compiler)
enableJIT(3)
enableJIT(3)

#set the seed
set.seed(1234)

# step 1 define the boundary for parameters
source("Data/DAIS_data.R")

#Define the parameters:
# [1] gamma #sensitivity of ice flow to sea level
# [2] alpha #sensitivity of ice flow to ocean subsurface temperature
# [3] mu    #Profile parameter related to ice stress [m^(1/2)]
# [4] eta   #Relates balance gradient to precipitation [m^(-1/2) yr^(-1/2)]
# [5] Po    #Precipitation at 0C [m of ice/yr] 
# [6] kappa #Relates precipitation to temperature [K^-1]
# [7] fo = 1.2                #Constant of proportionality for ice speed [m/yr]
# [8] ho = 1471               #Initial value for runoff line calculation [m]
# [9] co = 95                  #Second value for runoff line calculation [m]
# [10] bo = 775                #Height of bed at the center of the continent [m]
# [11] s = 0.0006              #Slope of the bed

project.forcings = matrix(c(TA,TO,GSL,SL), ncol=4, nrow=240300)
hindcast.forcings = matrix(c(TA[1:240010], TO[1:240010], GSL[1:240010], SL[1:240010]), ncol=4, nrow=240010)

# Best Case (Case #4) from Shaffer (2014)
IP = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Set the upper and lower bounds
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
bound.lower[1:2] = c(1/2, 0)   ; bound.upper[1:2] = c(17/4, 1) #Set bounds for gamma and alpha
bound.lower[10:11] = c(725, 0.00045)   ; bound.upper[10:11] = c(825, 0.00075) #Set bounds for bo and s

#bound.lower[12] = 0 ; bound.upper[12] = 1 # Prior uniform range for sigma (the variance)

#Source the function with the standards and the new parameters to
#get the best estimated AIS melt SLE:
standards = c(Tice,eps1, del, eps2, TOo, Volo, Roa, R)
source("Scripts/DAIS_IceFlux_model.R")
AIS_melt = iceflux(IP, hindcast.forcings, standards)

#set the end dates to the year 2300 to get future projections
end = 240298
enddate = 240300
Project_melt = iceflux(IP, project.forcings, standards)

####################### Set up windows and optimize parameters #################
#These windows are presented in Shaffer (2014) and Shepherd et al. (2012)
#1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
# We want the cumulative sea-level equivalent in meters for the year 2002
# 360Gt = 1mm SLE
estimate.SLE.rate = abs(-71/360)/1000
time.years = 2002-1992
mid.cum.SLE_2002 = estimate.SLE.rate*time.years

estimate.SLE.error = sqrt(time.years)*abs(-53/360)/1000 #1- sigma error
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

# Optimize the parameters to find the best hindcast and projection (slow ~50 minutes on single CPU)
end = 240000
enddate = 240010
lower = bound.lower[1:11] 
upper = bound.upper[1:11]
iter=50  # specify how many times the command should be run
source("Scripts/optimize_DAIS_model.R")
tresult = optim(IP, modelfn, gr=NULL, lower = lower, upper = upper,
                method = "L-BFGS-B", control = list(maxit = iter, 
                                                    trace=TRUE), 
                hindcast.forcings, standards, windows, obs.years)
print(tresult$par) #print out the estimates
#------------------ Save the workspace --------------------------------#
save.image(file = "Workspace/DAIS_precalibration_LHS.RData")
#----------------------------------------------------------------------#
params = c(tresult$par[1],tresult$par[2],tresult$par[3],tresult$par[4],tresult$par[5],tresult$par[6],
           tresult$par[7],tresult$par[8],tresult$par[9],tresult$par[10],tresult$par[11])

#Project the best fit from the optimized parameters:
end = 240298
enddate = 240300
best.project = iceflux(params, project.forcings, standards)


# Set up the function for Latin Hypercube Sampling  -----------------------
n_samples = 500 # For 12 parameters 500 samples should be sufficient
n_parameters = 12
# LHS Function
fn.LHS <- function(n_samples, x, y, bound.upper, bound.lower) {
  x <- maximinLHS(n_samples, n_parameters)
  y <- x
  alpha.var = 2
  beta.var = 1
  y[,1] <- qunif(x[,1], bound.lower[1], bound.upper[1])
  y[,2] <- qunif(x[,2], bound.lower[2], bound.upper[2]) 
  y[,3] <- qunif(x[,3], bound.lower[3], bound.upper[3])
  y[,4] <- qunif(x[,4], bound.lower[4], bound.upper[4])
  y[,5] <- qunif(x[,5], bound.lower[5], bound.upper[5]) 
  y[,6] <- qunif(x[,6], bound.lower[6], bound.upper[6])
  y[,7] <- qunif(x[,7], bound.lower[7], bound.upper[7])
  y[,8] <- qunif(x[,8], bound.lower[8], bound.upper[8]) 
  y[,9] <- qunif(x[,9], bound.lower[9], bound.upper[9])
  y[,10] <- qunif(x[,10], bound.lower[10], bound.upper[10])
  y[,11] <- qunif(x[,11], bound.lower[11], bound.upper[11]) 
  y[,12] <- qigamma(x[,12], alpha.var, beta.var) #  y[,12] <- qunif(x[,12], bound.lower[12], bound.upper[12])
  return(as.data.frame(y))
}

# Output ------------------------------------------------------------------

Parameters <- fn.LHS(n_samples, x, y, bound.upper, bound.lower)
colnames(Parameters, do.NULL = FALSE)
colnames(Parameters) <- c("gamma", "alpha", "mu", "eta", "po", "kappa","fo", "ho", "co","bo","s","sigma")

sample_length = n_samples #500
par=mat.or.vec(sample_length, 11)

for(i in 1:sample_length) {
  par[i,1] = Parameters[i,1]
  par[i,2] = Parameters[i,2]
  par[i,3] = Parameters[i,3]
  par[i,4] = Parameters[i,4]
  par[i,5] = Parameters[i,5]
  par[i,6] = Parameters[i,6]
  par[i,7] = Parameters[i,7]
  par[i,8] = Parameters[i,8]
  par[i,9] = Parameters[i,9]
  par[i,10] = Parameters[i,10]
  par[i,11] = Parameters[i,11]
}

end = 240298
enddate = 240300
lhs.dais.models = mat.or.vec(sample_length, enddate)
source("Scripts/DAIS_IceFlux_model.R")
standards = c(Tice, eps1, del, eps2, TOo, Volo, Roa, R)

for(i in 1:sample_length) {
  lhs.dais.models[i,] = iceflux(par[i,], project.forcings, standards)
}

# ### Superimpose the bias onto the model
# ### True world = model + bias
# bias = Parameters[,12]
# lhs.dais.mpbias = mat.or.vec(sample_length, enddate)
# for(i in 1:sample_length){
#   lhs.dais.mpbias[i,] = lhs.dais.models[i,] + bias[i]^2
# }

### Superimpose the bias onto the model
### True world = model + bias
bias = sqrt(Parameters[,12])
lhs.dais.mpanom = mat.or.vec(sample_length, enddate)
for(i in 1:sample_length){
  lhs.dais.mpanom[i,] = lhs.dais.models[i,] - mean(lhs.dais.models[i,SL.1961_1990])
}

dais.pre.cali = mat.or.vec(sample_length, enddate)
for(i in 1:sample_length){
  dais.pre.cali[i,] = lhs.dais.mpanom[i,] + rnorm(enddate,mean=0,sd=bias[i])
}
#------------------ Save the workspace --------------------------------#
save.image(file = "Workspace/DAIS_precalibration_LHS.RData")
# Write csv of SLE values for the targeted years -------------------------------------------------------------------
surv.targ = matrix(c(dais.pre.cali[1:sample_length,120000], dais.pre.cali[1:sample_length,220000],dais.pre.cali[1:sample_length,234000], 
                     dais.pre.cali[1:sample_length,240002]), nrow=sample_length, ncol=4)
colnames(surv.targ, do.NULL = FALSE)
colnames(surv.targ) <- c("Last Interglacial", "Last Glacial Max", "Holocene", "93-2011 Trend")
write.csv(surv.targ, file="Random_out/surviving_targets_LHS.csv")

surv.targ.nonoise = matrix(c(lhs.dais.mpanom[1:sample_length,120000], lhs.dais.mpanom[1:sample_length,220000], 
                             lhs.dais.mpanom[1:sample_length,234000], lhs.dais.mpanom[1:sample_length,240002]), 
                           nrow=sample_length, ncol=4)

#------------------------ Find the runs that pass through each constraint ---------------------------
source("Scripts/surviveTargetfunc.R")
surLIG = surviveTarget(windows[1,], surv.targ[,1])
surLGM = surviveTarget(windows[2,], surv.targ[,2])
surMH = surviveTarget(windows[3,], surv.targ[,3])
sur9311trend = surviveTarget(windows[4,], surv.targ[,4])

# If sur9311trend is true then uncomment the lines with sur9311trend or sur.all

#Find the runs that pass through the instrumental constraint
# total.sur = c(surLIG,surLGM,surMH,sur9311trend)
# template <- table(as.vector(total.sur))
# all = names(template)[template == max(template)]
# sur.all = as.numeric(all)

#Find the runs that pass through the LIG and the instrumental constraint (Do they better represent expert Judgement)
# LIG_trend.sur = c(surLIG,sur9311trend)
# LIG_trend.template <- table(as.vector(LIG_trend.sur))
# all.LIG.trend = names(LIG_trend.template)[LIG_trend.template == max(LIG_trend.template)]
# sur.LIG.trend = as.numeric(all.LIG.trend)

percent.include = c((sample_length/sample_length)*100, (length(surLIG)/sample_length)*100, 
                    (length(surLGM)/sample_length)*100, (length(surMH)/sample_length)*100, 
                    (length(sur9311trend)/sample_length)*100)
constraints = c("No constraints","Last integlacial","Last glacial maximum","Mid-Holocene","Instrumental period")
table.parts = matrix(c(constraints, percent.include), nrow=5, ncol=2)
write.csv(table.parts, file="Random_out/constraint_trend_percent.csv")

# No noise ----------------------------------------------------------------
surLIG_NN = surviveTarget(windows[1,], surv.targ.nonoise[,1])
surLGM_NN = surviveTarget(windows[2,], surv.targ.nonoise[,2])
surMH_NN = surviveTarget(windows[3,], surv.targ.nonoise[,3])
sur9311trend_NN = surviveTarget(windows[4,], surv.targ.nonoise[,4])

# If sur9311trend is true then uncomment the lines with sur9311trend or sur.all

#Find the runs that pass through all of the constraints
total.sur_NN = c(surLIG_NN, surLGM_NN, surMH_NN, sur9311trend_NN)
template_NN <- table(as.vector(total.sur_NN))
all_NN = names(template_NN)[template_NN == max(template_NN)]
sur.all_NN = as.numeric(all_NN)

#Find the runs that pass through the LIG and the instrumental constraint (Do they better represent expert Judgement)
LIG_trend.sur_NN = c(surLIG_NN, sur9311trend_NN)
LIG_trend.template_NN <- table(as.vector(LIG_trend.sur_NN))
all.LIG.trend_NN = names(LIG_trend.template_NN)[LIG_trend.template_NN == max(LIG_trend.template_NN)]
sur.LIG.trend_NN = as.numeric(all.LIG.trend_NN)

percent.include = c((sample_length/sample_length)*100, (length(surLIG_NN)/sample_length)*100, 
                    (length(surLGM_NN)/sample_length)*100, (length(surMH_NN)/sample_length)*100, 
                    (length(sur9311trend_NN)/sample_length)*100)
table.parts.nonoise = matrix(c(constraints, percent.include), nrow=5, ncol=2)
write.csv(table.parts.nonoise, file="Random_out/NoNoise_constraint_trend_percent.csv")

#--------------------- Parameter Pairsplot ------------------------------------
# Find the parameters that fit each of the constraints
LIG.parameters = Parameters[surLIG[1:length(surLIG)],]
LGM.parameters = Parameters[surLGM[1:length(surLGM)],]
MH.parameters = Parameters[surMH[1:length(surMH)],]
# instr.parameters = Parameters[sur9311trend[1:length(sur9311trend)],]
# all.const.parameters = Parameters[sur.all[1:length(sur.all)],]

length.parameters = length(MH.parameters[,1]) + length(LIG.parameters[,1]) + # Calculate the vector length
  length(LGM.parameters[,1])# + length(instr.parameters[,1]) + length(all.const.parameters[,1])
constr.parameters = mat.or.vec(length.parameters, 12) # set up an empty matrix 
for(i in 1:12){ #Combine all the surviving parameters. The parameters will stay in the order in which they were added
  constr.parameters[,i] = c(MH.parameters[,i], LIG.parameters[,i], LGM.parameters[,i])#, instr.parameters[,i], all.const.parameters[,i])  
}
# Create a vector to classify each line in the matrix to identify which parameter survived which target
color.ident = c(rep(1,length(MH.parameters[,1])), rep(2,length(LIG.parameters[,1])), 
                rep(3,length(LGM.parameters[,1])))#, rep(4,length(instr.parameters[,1])), rep(5,length(all.const.parameters[,1])))

colnames(constr.parameters, do.NULL = FALSE) # Name the columns in the matrix
colnames(constr.parameters) = c("Gamma","Alpha","Mu","Eta","Po","Kappa","fo", "ho", "co","bo", "s","sigma")

#--------------------- Estimate Parameter PDFs ----------------------------------
# Function to calculate hte pdfs of each parameter 'parameter.pdfs'
source('Scripts/plot_PdfCdfSf.R')

uncon.parameter.pdf = parameter.pdfs(Parameters)
LIG.parameter.pdf = parameter.pdfs(LIG.parameters)
LGM.parameter.pdf = parameter.pdfs(LGM.parameters)
MH.parameter.pdf = parameter.pdfs(MH.parameters)
# instr.parameter.pdf = parameter.pdfs(instr.parameters)
# all.parameter.pdf = parameter.pdfs(all.const.parameters)

#--------------------- Estimate PDFs, CDFin 2100 & 2050 --------------------------
# Function to find SLE values in certain years 'fn.prob.proj'
year.pcs = c(120000, 220000, 234000, 240002, 240050, 240100, 240300)

unprob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, n_samples, un.constr=T)
LIGprob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, surLIG)
LGMprob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, surLGM)
MHprob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, surMH)
# present.prob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, sur9311trend)
# all.prob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, sur.all)
# LIG.trend.prob_proj <- fn.prob.proj(dais.pre.cali, year.pcs, sur.LIG.trend)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 124,000 BP (Last interglacial)
un.sflig <- plot.sf(unprob_proj[,1], make.plot=F); LIG.sflig <- plot.sf(LIGprob_proj[,1], make.plot=F)
LGM.sflig <- plot.sf(LGMprob_proj[,1], make.plot=F); MH.sflig <- plot.sf(MHprob_proj[,1], make.plot=F)
# present.sflig <- plot.sf(present.prob_proj[,1], make.plot=F); all.sflig <- plot.sf(all.prob_proj[,1], make.plot=F)
# LIGtrend.sflig <- plot.sf(LIG.trend.prob_proj[,1], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 20,000 BP (Last glacial maximum)
un.sflgm <- plot.sf(unprob_proj[,2], make.plot=F); LIG.sflgm <- plot.sf(LIGprob_proj[,2], make.plot=F)
LGM.sflgm <- plot.sf(LGMprob_proj[,2], make.plot=F); MH.sflgm <- plot.sf(MHprob_proj[,2], make.plot=F)
# present.sflgm <- plot.sf(present.prob_proj[,2], make.plot=F); all.sflgm <- plot.sf(all.prob_proj[,2], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 6,000 BP (Mid-holocene)
un.sfmh <- plot.sf(unprob_proj[,3], make.plot=F); LIG.sfmh <- plot.sf(LIGprob_proj[,3], make.plot=F)
LGM.sfmh <- plot.sf(LGMprob_proj[,3], make.plot=F); MH.sfmh <- plot.sf(MHprob_proj[,3], make.plot=F)
# present.sfmh <- plot.sf(present.prob_proj[,3], make.plot=F); all.sfmh <- plot.sf(all.prob_proj[,3], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 2002 (Observed trend from 1993-2011)
un.sf2002 <- plot.sf(unprob_proj[,4], make.plot=F); LIG.sf2002 <- plot.sf(LIGprob_proj[,4], make.plot=F)
LGM.sf2002 <- plot.sf(LGMprob_proj[,4], make.plot=F); MH.sf2002 <- plot.sf(MHprob_proj[,4], make.plot=F)
# present.sf2002 <- plot.sf(present.prob_proj[,4], make.plot=F); all.sf2002 <- plot.sf(all.prob_proj[,4], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 2050
un.sf2050 <- plot.sf(unprob_proj[,5], make.plot=F); LIG.sf2050 <- plot.sf(LIGprob_proj[,5], make.plot=F)
LGM.sf2050 <- plot.sf(LGMprob_proj[,5], make.plot=F); MH.sf2050 <- plot.sf(MHprob_proj[,5], make.plot=F)
# present.sf2050 <- plot.sf(present.prob_proj[,5], make.plot=F); all.sf2050 <- plot.sf(all.prob_proj[,5], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 2100
un.sf2100 <- plot.sf(unprob_proj[,6], make.plot=F); LIG.sf2100 <- plot.sf(LIGprob_proj[,6], make.plot=F)
LGM.sf2100 <- plot.sf(LGMprob_proj[,6], make.plot=F); MH.sf2100 <- plot.sf(MHprob_proj[,6], make.plot=F)
# present.sf2100 <- plot.sf(present.prob_proj[,6], make.plot=F); all.sf2100 <- plot.sf(all.prob_proj[,6], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

# Calculate the pdf, cdf, and sf of AIS melt estimates in 2300
un.sf2300 <- plot.sf(unprob_proj[,7], make.plot=F); LIG.sf2300 <- plot.sf(LIGprob_proj[,7], make.plot=F)
LGM.sf2300 <- plot.sf(LGMprob_proj[,7], make.plot=F); MH.sf2300 <- plot.sf(MHprob_proj[,7], make.plot=F)
# present.sf2300 <- plot.sf(present.prob_proj[,7], make.plot=F); all.sf2300 <- plot.sf(all.prob_proj[,7], make.plot=F)
# LIGtrend.sflgm <- plot.sf(LIG.trend.prob_proj[,2], make.plot=F)

#------------------ Save the workspace --------------------------------#
save.image(file = "Workspace/DAIS_precalibration_LHS.RData")
#----------------------------------------------- END ------------------------------------------------------------#


