#########################################################
#------- Calibrate the model parameters
#------- using AR model
###-----HETEROSKEDASTIC
#########################################################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)

#Set the seed
set.seed(1234)


library(deSolve)
source("roblib.R")
dynReload("dais", srcname=c("dais.c", "r.c"), extrasrc="r.h")


# step 1 define the boundary for parameters
source("Data/DAIS_data.R")
source("Scripts/put_fig_letter.r")

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
case1 = c(1, 0, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)
case2 = c(2, 0, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)
case3 = c(1, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)
case4 = c(2, 0.35, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)
case5 = c(3.5, 0.45, 8.7, 0.012, 0.35, 0.04, 1.2, 1471, 95, 775, 0.0006)

#Source the function with the standards and the initial parameters (IP) to
# return (Vol):
standards = c(Tice,eps1, del, eps2, TOo, Volo, Roa, R)
source("Scripts/DAIS_IceFlux_model.R")
ptm <- proc.time()
AIS_case1 = iceflux(case1, hindcast.forcings, standards)
proc.time() - ptm
AIS_case2 = iceflux(case2, hindcast.forcings, standards)
AIS_case3 = iceflux(case3, hindcast.forcings, standards)
AIS_case4 = iceflux(case4, hindcast.forcings, standards)
AIS_case5 = iceflux(case5, hindcast.forcings, standards)

############################## CALCULATE RESIDUALS (PRIOR SIGMA) ##############################
#These windows are presented in Shaffer (2014) and the 1992 to 2011 trend from Shepherd et al. 2012 is -71 +/- 53 Gt per yr
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
    resid[i] <- (median(windows[i,])-(AIS_case4[obs.years[i]]-mean(AIS_case4[SL.1961_1990]))) #/sd(windows[i,])
	}
sigma = sd(resid)^2 #calculate the standard deviation (sigma)

############################## SETUP MCMC #######################################
#Set up the priors: the upper and lower bounds
bound.lower = case4 - (case4*0.5)    ; bound.upper = case4 + (case4*0.5)
bound.lower[1:2] = c(1/2, 0)   ; bound.upper[1:2] = c(17/4, 1) #Set bounds for gamma and alpha
bound.lower[10:11] = c(725, 0.00045)   ; bound.upper[10:11] = c(825, 0.00075) #Set bounds for bo and s

#bound.lower[12] = 0 ; bound.upper[12] = 1 # Prior uniform range for sigma (the variance)

Ta <<- hindcast.forcings[,1]     
SL <<- hindcast.forcings[,4]      
Toc <<- hindcast.forcings[,2] 

source("dais.R")

ptm <- proc.time()

print(system.time(for (i in 1:100)
brick_case1 = dais(
  b0    = case1[10],
  slope = case1[11],
  mu    = case1[3],
  h0    = case1[8],
  c     = case1[9],
  P0    = case1[5],
  kappa = case1[6],
  nu    = case1[4],
  f0    = case1[7],
  gamma = case1[1],
  alpha = case1[2],
  Tf    = Tice,
  ro_w  = Dsw*1000,
  ro_i  = Dice*1000,
  ro_m  = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa,
  tstep = 1)))

mp <- c(
  b0    = case1[10],
  slope = case1[11],
  mu    = case1[3],
  h0    = case1[8],
  c     = case1[9],
  P0    = case1[5],
  kappa = case1[6],
  nu    = case1[4],
  f0    = case1[7],
  Gamma = case1[1],
  alpha = case1[2],
  Tf    = Tice,
  rho_w = Dsw*1000,
  rho_i = Dice*1000,
  rho_m = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa
)

np     <- length(Ta)
Rad    <- numeric(length=np)               # Radius of ice sheet
Vais   <- numeric(length=np)               # Ice volume
SLE    <- numeric(length=np)               # Sea-level equivalent [m]

print(system.time(for (i in 1:100) .Call("daisOdeC", list(mp=mp, frc=hindcast.forcings, out=list(Rad, Vais, SLE)), PACKAGE = "dais")))

if (1) {
    vol1 <- Vais
    .Call("daisOdeC", list(mp=mp, frc=hindcast.forcings, out=list(Rad, Vais, SLE)), PACKAGE = "dais")
    vol2 <- Vais
    print(any(vol1 != vol2))
}


# brick_case1_sle = 57*(1-brick_case1/Volo)
proc.time() - ptm

brick_case2 = dais(
  b0    = case2[10],
  slope = case2[11],
  mu    = case2[3],
  h0    = case2[8],
  c     = case2[9],
  P0    = case2[5],
  kappa = case2[6],
  nu    = case2[4],
  f0    = case2[7],
  gamma = case2[1],
  alpha = case2[2],
  Tf    = Tice,
  ro_w  = Dsw*1000,
  ro_i  = Dice*1000,
  ro_m  = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa,
  tstep = 1) 

# brick_case2_sle = 57*(1-brick_case2/Volo)

brick_case3 = dais(
  b0    = case3[10],
  slope = case3[11],
  mu    = case3[3],
  h0    = case3[8],
  c     = case3[9],
  P0    = case3[5],
  kappa = case3[6],
  nu    = case3[4],
  f0    = case3[7],
  gamma = case3[1],
  alpha = case3[2],
  Tf    = Tice,
  ro_w  = Dsw*1000,
  ro_i  = Dice*1000,
  ro_m  = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa,
  tstep = 1) 

# brick_case3_sle = 57*(1-brick_case3/Volo)

brick_case4 = dais(
  b0    = case4[10],
  slope = case4[11],
  mu    = case4[3],
  h0    = case4[8],
  c     = case4[9],
  P0    = case4[5],
  kappa = case4[6],
  nu    = case4[4],
  f0    = case4[7],
  gamma = case4[1],
  alpha = case4[2],
  Tf    = Tice,
  ro_w  = Dsw*1000,
  ro_i  = Dice*1000,
  ro_m  = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa,
  tstep = 1) 

# brick_case4_sle = 57*(1-brick_case4/Volo)

brick_case5 = dais(
  b0    = case5[10],
  slope = case5[11],
  mu    = case5[3],
  h0    = case5[8],
  c     = case5[9],
  P0    = case5[5],
  kappa = case5[6],
  nu    = case5[4],
  f0    = case5[7],
  gamma = case5[1],
  alpha = case5[2],
  Tf    = Tice,
  ro_w  = Dsw*1000,
  ro_i  = Dice*1000,
  ro_m  = Drock*1000,
  Toc_0 = TOo,
  Rad0  = Roa,
  tstep = 1) 

large = matrix(c(AIS_case1,AIS_case2,AIS_case3,AIS_case4,AIS_case5,brick_case1, 
                 brick_case2, brick_case3, brick_case4, brick_case5), nrow=length(Ta), ncol=10)
write.csv(large, file="OriginalVSbrick.csv")


ptm <- proc.time()
brick_case1_sle = 57*(1-brick_case1/Volo)
proc.time() - ptm

############################## PLOT VOLUME #######################################
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlGnBu")

boxpalette <- brewer.pal(11,"Spectral")

plot(date[1:240010], AIS_case4, type="l", col="blue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume [m]", xaxt="n")
lines(date[1:240010], brick_case4, col="lightblue", lwd=2)
lines(date[1:240010], AIS_case3, col="green", lwd=2)
lines(date[1:240010], brick_case3, col="lightgreen", lwd=2)
lines(date[1:240010], AIS_case2, col="magenta", lwd=2)
lines(date[1:240010], brick_case2, col="pink", lwd=2)
lines(date[1:240010], AIS_case1, col="darkred", lwd=2)
lines(date[1:240010], brick_case1, col="red", lwd=2)
lines(date[1:240010], AIS_case5, col="darkturquoise", lwd=2)
lines(date[1:240010], brick_case5, col="cyan", lwd=2)

ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))


