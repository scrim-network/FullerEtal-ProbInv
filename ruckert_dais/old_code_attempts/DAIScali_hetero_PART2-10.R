#########################################################
#------- Now let's calibrate the model parameters
#------- using iid model
###-----HETEROSKEDASTIC
#########################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)

#Set the seed
set.seed(1234)

# step 1 load in previous mcmc chain
load(file = "Workspace/DAIS_MCMC_part1.RData")
#load(file = "Workspace/DAIS_MCMC_part2.RData")
#load(file = "Workspace/DAIS_MCMC_part3.RData")
#load(file = "Workspace/DAIS_MCMC_part4.RData")
#load(file = "Workspace/DAIS_MCMC_part5.RData")
#load(file = "Workspace/DAIS_MCMC_part6.RData")
#load(file = "Workspace/DAIS_MCMC_part7.RData")
#load(file = "Workspace/DAIS_MCMC_part8.RData")
#load(file = "Workspace/DAIS_MCMC_part9.RData")

############################## SET INITIAL PARAMETERS ##############################
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

standards = c(Tice,eps1, del, eps2, TOo, Volo, Roa, R)

#Set the end dates back to the hindcast period:
end = 240000
enddate = 240010

############################## Setup MCMC #######################################
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

# We will set the initial parameters to the last set in the MCMC chain
p1 = dh.chain[NI,]
#p2 = dh.chain1[NI,]
#p3 = dh.chain2[NI,]
#p4 = dh.chain3[NI,]
#p5 = dh.chain4[NI,]
#p6 = dh.chain5[NI,]
#p7 = dh.chain6[NI,]
#p8 = dh.chain7[NI,]
#p9 = dh.chain8[NI,]

library(mcmc)

step = c(0.1, 0.015, 0.2, 0.025, 0.1, 0.01, 0.1, 50, 10, 20, 0.0005, 0.15)/100
#step = c(0.001, 0.0001, 0.001, 0.00001, 0.0001, 0.00001, 0.001, 0.5, 0.1, 0.5, 0.000001, 0.001)
# NI = 300
NI = 1E5 #number of iterations
#n.cores = 64 #number of cores #needs to be same as specified in shell #use only half of the cores in the node
burnin = seq(1, 0.01*NI, 1)

#Run the MCMC chain using p1
dais.out.heter1 = metrop(log.post, p1, nbatch=NI, scale=step)
dh.chain1 = dais.out.heter1$batch
dais.out.heter1$accept
#Echo the acceptance rate. Should be ~ 25%

#Run the MCMC chain using p2
#dais.out.heter2 = metrop(log.post, p2, nbatch=NI, scale=step)
#dh.chain2 = dais.out.heter2$batch
#dais.out.heter2$accept

#Run the MCMC chain using p3
#dais.out.heter3 = metrop(log.post, p3, nbatch=NI, scale=step)
#dh.chain3 = dais.out.heter3$batch
#dais.out.heter3$accept

#Run the MCMC chain using p4
#dais.out.heter4 = metrop(log.post, p4, nbatch=NI, scale=step)
#dh.chain4 = dais.out.heter4$batch
#dais.out.heter4$accept

#Run the MCMC chain using p5
#dais.out.heter5 = metrop(log.post, p5, nbatch=NI, scale=step)
#dh.chain5 = dais.out.heter5$batch
#dais.out.heter5$accept

#Run the MCMC chain using p6
#dais.out.heter6 = metrop(log.post, p6, nbatch=NI, scale=step)
#dh.chain6 = dais.out.heter6$batch
#dais.out.heter6$accept

#Run the MCMC chain using p7
#dais.out.heter7 = metrop(log.post, p7, nbatch=NI, scale=step)
#dh.chain7 = dais.out.heter7$batch
#dais.out.heter7$accept

#Run the MCMC chain using p8
#dais.out.heter8 = metrop(log.post, p8, nbatch=NI, scale=step)
#dh.chain8 = dais.out.heter8$batch
#dais.out.heter8$accept

#Run the MCMC chain using p9
#dais.out.heter9 = metrop(log.post, p9, nbatch=NI, scale=step)
#dh.chain9 = dais.out.heter9$batch
#dais.out.heter9$accept

# Save workspace
save.image(file = "Workspace/DAIS_MCMC_part2.RData")
#save.image(file = "Workspace/DAIS_MCMC_part3.RData")
#save.image(file = "Workspace/DAIS_MCMC_part4.RData")
#save.image(file = "Workspace/DAIS_MCMC_part5.RData")
#save.image(file = "Workspace/DAIS_MCMC_part6.RData")
#save.image(file = "Workspace/DAIS_MCMC_part7.RData")
#save.image(file = "Workspace/DAIS_MCMC_part8.RData")
#save.image(file = "Workspace/DAIS_MCMC_part9.RData")
#save.image(file = "Workspace/DAIS_MCMC_part10.RData")

############################## END #######################################
