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
    resid[i] <- (median(windows[i,])-(AIS_melt[obs.years[i]]-mean(AIS_melt[SL.1961_1990]))) #/sd(windows[i,])
	}
sigma = sd(resid)^2 #calculate the standard deviation (sigma)

############################## SETUP MCMC #######################################
#Set up the priors: the upper and lower bounds
bound.lower = IP - (IP*0.5)    ; bound.upper = IP + (IP*0.5)
bound.lower[1:2] = c(1/2, 0)   ; bound.upper[1:2] = c(17/4, 1) #Set bounds for gamma and alpha
bound.lower[10:11] = c(725, 0.00045)   ; bound.upper[10:11] = c(825, 0.00075) #Set bounds for bo and s

#bound.lower[12] = 0 ; bound.upper[12] = 1 # Prior uniform range for sigma (the variance)

############################## PLOT #######################################
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlGnBu")

boxpalette <- brewer.pal(11,"Spectral")

pdf(file="Ruckertetal_dais_Fig1.pdf", family="Helvetica", width=6.7, height=8.1, pointsize=12)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig1_dais_mcmcLHS.pdf", family="Helvetica",height=2.7, width=6.7,pointsize=11)
#par(mfrow=c(1,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
# Last interglacial 240 kyr Bp - 2010 AD
plot(date, Project_melt-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-18,windows[1,2]), xaxt="n")

lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2],median(windows[1,]), col="black", pch=8, cex=0.75)

ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

axis(1, at = date[kyrbp_25]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

#pdf(file="NEWnFig25kyr.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[kyrbp_25:240010], Project_melt[kyrbp_25:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue",
lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]",ylim=c(windows[2,1],0), xaxt="n")

lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=2)
points(last.glacialmax[2],median(windows[2,]), col="black", pch=8, cex=0.75)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="b.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_6]:date[240010], labels = FALSE, col=boxpalette[2], lwd=2)

# Mid-Holocene6 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig6kyr.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[kyrbp_6:240010], Project_melt[kyrbp_6:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(windows[3,1],0))

lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=2)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=2)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=2)
points(holocene[2],median(windows[3,]), col="black", pch=8, cex=0.75)
put.fig.letter(label="c.", location="topleft", font=2)

box(col = boxpalette[2])
axis(1, at = date[AD_1880]:date[240010], labels = FALSE, col=boxpalette[3], lwd=2)

# Present day 1880 - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

#pdf(file="NEWnFig1880.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
plot(date[AD_1880:240010], Project_melt[AD_1880:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(-0.10,0.05), xaxt="n")

err_positive = windows[4,2]
err_negative = windows[4,1]
present = 2
arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
points(present,median(windows[4,]), col="black", pch=8, cex=0.75)
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="d.", location="topleft", font=2)

box(col = boxpalette[3])
axis(1, at = date[240000]:date[240020], labels = FALSE, col=boxpalette[4], lwd=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
#par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(date[240000:enddate], Project_melt[240000:enddate]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", xaxt="n")

err_positive = windows[4,2]
err_negative = windows[4,1]
arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[4])

dev.off()

############################## END #######################################
