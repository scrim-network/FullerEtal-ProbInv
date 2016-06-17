##################################################################################################
#
#  -file = "Paper2_plots_v2.R"   Code written September 2014
#  - Author: Kelsey Ruckert (klr324@psu.edu)
#
#  -This program loads in the LHS and mcmc workspace
#
#   -NOTE: The graphs will be saved as pdf files in the current working directory.
#       Large files will also be save as a jpeg file.
#
###################################################################################################

load("Workspace/DAIS_precalibration_LHS.RData") # Load in the saved workspace from LHS precalibration

load("Workspace/DAIS_MCMC_calibration.RData") # Load in the saved workspace from MCMC calibration

load("Workspace/DAIS_MCMC_Matlabcalibration.RData")

#install.packages('ash')
library(ash)
#install.packages('fields')
library(fields)
#install.packages('RColorBrewer')
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlGnBu")

boxpalette <- brewer.pal(11,"Spectral")

source("Scripts/put_fig_letter.r")
source("Scripts/plot_rangefn.R")

#------------------------------------- Create mean estimate projection and output R, H, & F ----------------------
source("Scripts/iceflux.mult_func_outRHF.R")
#Project the best fit from the optimized parameters:
end = 240298
enddate = 240300
RHF_outputs = iceflux_RHF(mean.dais.par, project.forcings, standards)

#------------------------------------- Create mean estimate projection -------------------------------------------
#pdf(file="Figures/SuppFigures/Ruckertetal_daisRHF_Sfig2a.pdf", family="Helvetica", height=5.4, width=6.7,pointsize=11)
png(file="Figures/SuppFigures/Ruckertetal_daisRHF_Sfig2a.tif", family="Helvetica", width=6.7, 
    height=5.4, units="in",pointsize=12, res=300)
par(mfrow=c(2,2)) # set figure dimensions
# Last interglacial 240 kyr Bp - 2010 AD
plot(date[90000:140000], RHF_outputs$SLE[90000:140000]-mean(RHF_outputs$SLE[SL.1961_1990]),
     typ="l", col="goldenrod2", lwd=2, xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", 
     ylim=c(-1,6), xaxt="n")
abline(h=0, lty=2, col="black")

lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2],median(windows[1,]), col="black", pch=8, cex=0.75)
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="a.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$Rad[90000:140000], col="blue", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Radius [m]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="b.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$WatDepth[90000:140000], col="black", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Water depth at the grounding line (H) [m]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="c.", location="topleft", font=2)

plot(date[90000:140000], RHF_outputs$Flow[90000:140000], col="red", lwd=2,typ="l",
     xlab="Date [Kyr BP]", ylab="Ice flux at the grounding line (F) [m/yr]", xaxt="n")
ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
put.fig.letter(label="d.", location="topleft", font=2)
dev.off()

# pdf(file="WaterdepthminusSL.pdf", family="Helvetica", height=2.7, width=3.5,pointsize=7)
# par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 3, 3, 3.5))
# plot(date[90000:140000], SL[90000:140000], col="purple", lwd=2,typ="l",
#      xlab="Date [Kyr BP]", ylab="Global sea-level anomaly [m]", xaxt="n")
# abline(h=0, lty=2, col="black")
# ticks=c(-150000,-140000,-130000,-120000,-110000,-100000)
# axis(side=1, at=ticks, labels=expression(-150, -140, -130,-120,-110,-100))
# text(-140000, -100, "Sea-level", col="purple")
# text(-145000, -20, "H - SL", col="black")
# 
# par(new=TRUE)
# plot(date[90000:140000], RHF_outputs$WatDepth[90000:140000] - SL[90000:140000], col="black", lwd=2,typ="l",
#      xlab="", ylab="",yaxt="n", xaxt="n")
# 
# axis(side=4)
# mtext(4, text="Water depth at the grounding line - SL [m]", line=2)
# dev.off()

#------------------------------------- Transparent Color Function -------------------------------------------
makeTransparent<- function(somecolor, alpha=100){
  someColor = someColor
  newColor<-col2rgb(someColor)
  apply(newColor,2 ,
        function(curcoldata)
        {rgb(red=curcoldata[1],
             green=curcoldata[2],
             blue=curcoldata[3], alpha=alpha,
             maxColorValue=255)})
}

#------------------------------------- Calculate 90% and 99% credible intervals -------------------------------------------
mcmc_5 <-
  mcmc_95 <-
  mcmc_2point5 <-
  mcmc_975 <-
  mcmc_point5 <-
  mcmc_std <-
  mcmc_mean<-
  mcmc_995 <-rep(NA,enddate) # 240300 years
for(i in 1:enddate){
  mcmc_5[i] = quantile(proj.mcmc.1961_1990[,i],0.05) # 90%
  mcmc_95[i] = quantile(proj.mcmc.1961_1990[,i],0.95)
 
  mcmc_2point5[i] = quantile(proj.mcmc.1961_1990[,i],0.025) # 95%
  mcmc_975[i] = quantile(proj.mcmc.1961_1990[,i],0.975)
  
  mcmc_point5[i] = quantile(proj.mcmc.1961_1990[,i],0.005) # 99%
  mcmc_995[i] = quantile(proj.mcmc.1961_1990[,i],0.995)
  
  mcmc_std[i] = sd(proj.mcmc.1961_1990[,i]) # 1 standard deviation
  mcmc_mean[i] = mean(proj.mcmc.1961_1990[,i]) #mean
}

x_mcmc_90=c(mcmc_5, rev(mcmc_95)); y_mcmc_90=c(date, rev(date))

mcmc_plusSTD = mcmc_mean + mcmc_std
mcmc_minusSTD = mcmc_mean - mcmc_std
plus_minusSTD = c(mcmc_minusSTD[1:240100], rev(mcmc_plusSTD[1:240100])); y_STD=c(date[1:240100], rev(date[1:240100]))

#------------------------------------- Figure 4 -------------------------------------------
someColor = c("red", "goldenrod2", "blue")
studies.color <- makeTransparent(someColor, 100)

png(file="Figures/Fig_4.tif", family="Helvetica", units="in", width=3.4, height=3, pointsize=11, res=300)
par(mfrow=c(1,1), mgp=c(1.5,.5,0),mar=c(4, 4, 2, 1))
plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Year", ylab="AIS Volume loss [SLE m]", 
     xlim=c(-10,106),ylim=c(-0.1,1.35), xaxt="n")

polygon(y_STD, plus_minusSTD, col=studies.color[2], border=NA)
lines(date[1:240100], mcmc_mean[1:240100], col="goldenrod2")
abline(v=100, lty=2)
text(80, 0.2, "0.05 ± 0.05", col="goldenrod2", cex=0.75)

# DeConto and Pollard 2016
highPliocene = c((1.05-0.3), 1.05, (1.05+0.3))
lowPliocene = c((0.64-0.49), 0.64, (0.64+0.49))

place.where = c(104, 108)
width = 2.5  # Width of bars

polygon(x = c((place.where[1]), (place.where[1]),  (place.where[1] - width), (place.where[1]-width)),
          y = c(highPliocene[1], highPliocene[3], highPliocene[3], highPliocene[1]), 
          border=NA, col = studies.color[1])
points((place.where[1]-1), highPliocene[2], col="red", pch="_")
text(80, 0.9, "1.05 ± 0.30", col="red", cex=0.75)

polygon(x = c((place.where[2]), (place.where[2]),  (place.where[2] - width), (place.where[2]-width)),
        y = c(lowPliocene[1], lowPliocene[3], lowPliocene[3], lowPliocene[1]), 
        border=NA, col = studies.color[3])
points((place.where[2]-1), lowPliocene[2], col="blue", pch="_")
text(80, 0.6, "0.64 ± 0.49", col="blue", cex=0.75)

ticks=c(0,20,40,60,80,100)
axis(side=1, at=ticks, labels=expression(2000,2020,2040,2060,2080,2100))

legend.names = c(expression(paste("±1", sigma, " This study")), expression(paste("±1", sigma, " DeConto & Pollard 2016")), 
                 expression(paste("±1", sigma, " DeConto & Pollard 2016")))
legend("topleft", legend=legend.names, pch=15, col=c(studies.color[2], studies.color[1], studies.color[3]), bty="n", cex=0.85)
dev.off()
#------------------------------------- Figure 1 -------------------------------------------
# MCMC hindcasts & projection
png(file="Figures/Ruckertetal_dais90MCMC_Fig1a.tif", family="Helvetica", width=6.7, 
    height=8.1, units="in",pointsize=12, res=300)
# pdf(file="Figures/Ruckertetal_dais90MCMC_Fig1.pdf", family="Helvetica", width=6.7, height=8.1, pointsize=12)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

#pdf(file="NEWnFig1_dais_mcmcLHS.pdf", family="Helvetica",height=2.7, width=6.7,pointsize=11)
#par(mfrow=c(1,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1], proj.mcmc.1961_1990[1,1], type="l", col="powderblue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[1], date[240010]), ylim=c(-21,7), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2", border=NA)
lines(date, mcmc_2point5, lty=3, col=boxpalette[1])
lines(date, mcmc_975, lty=3, col=boxpalette[1])
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)
abline(h=0, lty=3, col="black")

lines(date,proj.mcmc.1961_1990[990,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[1225,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[2005,] , col="dimgray", lwd=0.75)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray", lwd=0.75)
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
plot(date[kyrbp_25], proj.mcmc.1961_1990[1,kyrbp_25], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]",
     xlim=c(date[kyrbp_25], date[240010]), ylim=c(-26,0.5), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2", border=NA)
lines(date, mcmc_2point5, lty=3, col=boxpalette[1])
lines(date, mcmc_975, lty=3, col=boxpalette[1])
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date,proj.mcmc.1961_1990[990,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[1225,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[2005,] , col="dimgray", lwd=0.75)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
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
plot(date[kyrbp_6], proj.mcmc.1961_1990[1,kyrbp_6], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[kyrbp_6], date[240010]),ylim=c(-5,0.5))

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2", border=NA)
lines(date, mcmc_2point5, lty=3, col=boxpalette[1])
lines(date, mcmc_975, lty=3, col=boxpalette[1])
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date,proj.mcmc.1961_1990[990,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[1225,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[2005,] , col="dimgray", lwd=0.75)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
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
plot(date[AD_1880], proj.mcmc.1961_1990[1,AD_1880], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(date[AD_1880], date[240010]),ylim=c(-0.04,0.02), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2", border=NA)
lines(date, mcmc_2point5, lty=3, col=boxpalette[1])
lines(date, mcmc_975, lty=3, col=boxpalette[1])
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date,proj.mcmc.1961_1990[990,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[1225,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[2005,] , col="dimgray", lwd=0.75)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
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
plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", 
     xlim=c(-10,290),ylim=c(-0.1,3), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2", border=NA)
lines(date, mcmc_2point5, lty=3, col=boxpalette[1])
lines(date, mcmc_975, lty=3, col=boxpalette[1])
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date,proj.mcmc.1961_1990[990,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[1225,] , col="dimgray", lwd=0.75)
lines(date,proj.mcmc.1961_1990[2005,] , col="dimgray", lwd=0.75)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
err_positive = windows[4,2]
err_negative = windows[4,1]
arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[4])

plot.new()

#pdf(file="NEWnFiglegend2.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
#par(mfrow=c(1,1),mgp=c(1.5,.5,0)) # set figure dimensions
legend("center", c("Data constraint", "Mean hindcast & projection", "90% credible interval","99% credible interval"),
       lty=c(1,1,NA,2), pch=c(NA,NA,15,NA), col=c("black", "dimgray", "goldenrod2", "black"))
dev.off()

#------------------------------------- Figure 2 -------------------------------------------

#Set up previous study ranges for the year 2100 using RCP 8.5
#the previous ranges are in cm so divide by 100 to get meters
pfeffer = c(12.8, 14.6, 61.9)/100 #Low2, low1, and high estimate in Pfeffer et al. 2008
Bamber = c(-2, 14, 83)/100 # 5%, Median, 95% estimates in Bamber and Aspinall 2013
IPCC_AR5 = c(-15, 4, 23)/100 # 5%, Median, 95% estimates in IPCC AR5
Little = c(-8, 2.4, 13.3)/100 # 5%, Median, 95% estimates in Little et al. 2013
Kopp = c(-11, 4, 33)/100 # 5%, Median, 95% estimates in Kopp et al. 2014
present = 2

#------------------------------------- Figure 3 -------------------------------------------

#pdf(file="Figures/Ruckertetal_daisLIG2100pdf_Fig2aa.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
png(file="Figures/Ruckertetal_daisMcmc2100pdf_Fig2aa.tif", family="Helvetica", width=6.7, 
    height=5.4, units="in",pointsize=12, res=300)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions

plot(LIG.sf.mcmc$pdf, main="",lwd=3, col="goldenrod2", xlab="Projected AIS Volume loss during", sub="Last interglacial [SLE m]", 
     ylab="Probability Density",xlim=c(min(LIG.sf.mcmc$pdf$x),14), ylim=c(-0.85,2), yaxt="n")

lines(LIG.sf.mcmc$pdf, col="goldenrod2", lwd=2)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,1], probs, width.size = 0.75, where.at = -0.5, tick.length = 0.05, line.width = 2, color = "goldenrod2")
#add.hor.box(all.prob_proj[,1], probs, width.size = 0.35, where.at = -0.5, tick.length = 0.02, line.width = 2, color = mypalette[9])

abline(h=1, lty=2)
text(10, 1+0.25, cex=0.75, "Last Interglacial period
Data Constraint")
text(10, 1-0.25, cex=0.75, "Model Inversion
This study")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=1.5, color="black")
put.fig.letter(label="a.", location="topleft", font=2)

par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(sf.2100.mcmc$pdf, main="", col="goldenrod2",lwd=2,
     xlab="Projected AIS Volume loss in 2100 [SLE m]", 
     ylab="Probability Density", 
     xlim=c(IPCC_AR5[1],Bamber[3]), ylim=c(-7,53), yaxt="n") 
lines(sf.2100.mcmc$pdf, col="goldenrod2", lwd=2)

probs = c(0.05, 0.95)
add.hor.box(mcmc.prob_proj[,6], probs, width.size = 8, where.at = -3, tick.length = 0.35, line.width = 2, color = "goldenrod2")
#add.hor.box(all.prob_proj[,6], probs, width.size = 3.5, where.at = -6, tick.length = 0.25, line.width = 2, color = mypalette[9])

abline(h=30, lty=2)
text(0.6, 30+3, cex=0.75, "Expert Assessments")
text(0.6, 30-5, cex=0.75, "Model Inversion
This study")
plotrange(Little[1], Little[2], Little[3], year=F, height=34, color="pink")
plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=38, color="purple")
plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=42, color="orange")
plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=46, color="gray")
plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=50, color="red")

put.fig.letter(label="b.", location="topleft", font=2)

plot.new()
legend("left", c("MCMC", "90% C.I. (Little et al. 2013)","90% C.I. (IPCC AR5)", 
                 "90% C.I. (Kopp et al. 2014)", "90% C.I. (Pfeffer et al. 2008)","90% C.I. (Bamber & Aspinall 2013)"),
       lty=c(NA,1,1,1,1,1), pch=c(15,8,8,8,8,8), lwd=2, bty="n", 
       col=c("goldenrod2","pink","purple","orange","gray","red"))
dev.off()

###################################### SUPPLEMENTARY FIGURES ############################################
#------------------------------------- Supplementary Figure 1 -------------------------------------------

# # LHS hindcasts & projection
# width=1920, height=1080
png(file="Figures/SuppFigures/SuppFig3_dais.tif", family="Helvetica", width=6.7, height=8.1, units="in",pointsize=12, res=300)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
# jpeg(file="nSuppFig3_dais_mcmcLHS.jpeg", family="Helvetica", width=1590, height=1920, units="px", pointsize=40)
# par(mfrow=c(3,2), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 2))

# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1:240010], AIS_melt-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-30,10), xaxt="n")
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
abline(h=0, lty=2, col="black")
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=1.5)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=1.5)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=1.5)
ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

axis(1, at = date[kyrbp_25]:date[240010], labels = FALSE, col=boxpalette[1], lwd=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[kyrbp_25:240010], AIS_melt[kyrbp_25:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", 
     lwd=1,xlab="Date [Kyr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-30,7), xaxt="n")
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
# abline(h=0, lty=2, col="gray")
lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=1.5)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=1.5)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=1.5)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="b.", location="topleft", font=2)

box(col = boxpalette[1])
axis(1, at = date[kyrbp_6]:date[240010], labels = FALSE, col=boxpalette[2], lwd=2)

# Mid-Holocene6 kyr Bp - 2010 AD 
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(date[kyrbp_6:240010], AIS_melt[kyrbp_6:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-8,4))
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=1.5)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=1.5)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=1.5)
put.fig.letter(label="c.", location="topleft", font=2)

box(col = boxpalette[2])
axis(1, at = date[AD_1880]:date[240010], labels = FALSE, col=boxpalette[3], lwd=2)

# Present day 1880 - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))

plot(date[AD_1880:240010], Project_melt[AD_1880:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", ylim=c(-0.10,0.05), xaxt="n")
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
err_positive = windows[4,2]
err_negative = windows[4,1]
present = 2
arrows(present, err_positive, present, err_negative, length=0, lwd=1.5, col="black")
points(present,median(windows[4,]), col="black", pch=8, cex=0.75)
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="d.", location="topleft", font=2)

box(col = boxpalette[3])
axis(1, at = date[240000]:date[240020], labels = FALSE, col=boxpalette[4], lwd=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

plot(date[240000:enddate], Project_melt[240000:enddate]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=1,
     xlab="Date [yr BP]", ylab="AIS Volume loss [SLE m]", xlim=c(-10,290),ylim=c(-0.5,5), xaxt="n")
for(i in 1:sample_length){
  lines(date, dais.pre.cali[i,], col="gray91", lwd=1)
}
for(i in 1:length(surMH)){
  lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=1)
}
for(i in 1:length(surLIG)){
  lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=1)
}
for(i in 1:length(surLGM)){
  lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=1)
}
for(i in 1:length(sur9311trend)){
  lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=1)
}
for(i in 1:length(sur.all)){
  lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=1)
}
lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=1)
err_positive = windows[4,2]
err_negative = windows[4,1]
arrows(present, err_positive, present, err_negative, length=0, lwd=1.5, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="e.", location="topleft", font=2)

box(col = boxpalette[4])

plot.new()
legend("center", c("Data constraints", "Optim best-fit", "Last interglacial constraint", "Last glacial maximum constraint","Mid-Holocene constraint", 
                   "Instrumental period (1992-2011) constraint", "All constraints"),lty=c(1,1,NA,NA,NA,NA,NA,NA), pch=c(NA,NA,15,15,15,15,15,15),
       col=c("black", "black", mypalette[3], mypalette[5], mypalette[1], mypalette[7], mypalette[9]))
dev.off()

#------------------------------------- Supplementary Figure 2 -------------------------------------------
# someColor = c(mypalette[3], mypalette[5], mypalette[7])
# trans_con_colors = makeTransparent(someColor,100)
# 
# pdf(file="nSuppFig2_dais_mcmcLHS.pdf", family="Helvetica", pointsize=12)
# pairs(constr.parameters, main="", pch=20, 
#       col = c(mypalette[1], trans_con_colors[1], trans_con_colors[2], trans_con_colors[3], "black")[color.ident], upper.panel=NULL)
# 
# par(xpd=TRUE)
# legend(0.5, 1, c("Last interglacial constraint", "Last glacial maximum constraint","Mid-Holocene constraint", "Instrumental period (1993-2011)
# constraint", "All constraints"), bty="n",fill=c(trans_con_colors[1], trans_con_colors[2], mypalette[1], trans_con_colors[3], "black"))
# dev.off()

#------------------------------------- Supplementary Figure 3 -------------------------------------------
pairs.image <- function(x) {
  pairs(x, panel=function(x,y) {
    foo <- bin2(cbind(x,y),nbin=c(50,50))
    foo <- ash2(foo,m=c(8,8))
    image(foo,add=T,xlab="",ylab="",col=topo.colors(1000))
  }, upper.panel = NULL)
}

pdf(file="Figures/SuppFig5_dais_mcmc.pdf", family="Helvetica", pointsize=12)
#png(file="Figures/SuppFig5_dais_mcmc.tif", family="Helvetica",pointsize=12, res=300)
pairs.image(d.pos_parameters)
par(fig = c(.5, 1, .5, 1), new=TRUE) 

image.plot(zlim=c(0,5), legend.only=TRUE,col=topo.colors(1000), legend.shrink = 0.75, legend.lab=NULL,
           legend.width = 0.8, horizontal = TRUE, 
           axis.args=list(at=c(0,5),labels=c("Less likely", "More likely")))
dev.off()
#------------------------------------- Supplementary Figure 3 -------------------------------------------

#pdf(file="Figures/SuppFigures/suppFig4_dais_LHS.pdf", family="Helvetica",height=5.4, width=6.7,pointsize=11)
png(file="Figures/SuppFigures/suppFig4_dais_LHS.tif", family="Helvetica", width=6.7, 
    height=5.4, units="in",pointsize=12, res=300)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions

plot(un.sflig$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS Volume loss during", sub="Last interglacial [SLE m]", 
     ylab="Probability Density",xlim=c(min(un.sflig$pdf$x),20), ylim=c(-0.92,3), yaxt="n")

lines(LIG.sflig$pdf, col=mypalette[3], lwd=2) 
lines(LGM.sflig$pdf, col=mypalette[5], lwd=2) 
lines(MH.sflig$pdf, col=mypalette[1], lwd=2) 
lines(present.sflig$pdf, col=mypalette[7], lwd=2)
lines(all.sflig$pdf, col=mypalette[9], lwd=2)

probs = c(0.05, 0.95)
add.hor.box(all.prob_proj[,1], probs, width.size = 0.2, where.at = -0.1, tick.length = 0.01, line.width = 2, color = mypalette[9])
add.hor.box(present.prob_proj[,1], probs, width.size = 0.2, where.at = -0.28, tick.length = 0.01, line.width = 2, color = mypalette[7])
add.hor.box(MHprob_proj[,1], probs, width.size = 0.2, where.at = -0.44, tick.length = 0.01, line.width = 2, color = mypalette[1])
add.hor.box(LGMprob_proj[,1], probs, width.size = 0.2, where.at = -0.6, tick.length = 0.01, line.width = 2, color = mypalette[5])
add.hor.box(LIGprob_proj[,1], probs, width.size = 0.2, where.at = -0.76, tick.length = 0.01, line.width = 2, color = mypalette[3])
add.hor.box(unprob_proj[,1], probs, width.size = 0.2, where.at = -0.90, tick.length = 0.01, line.width = 2, color = "gray91")


abline(h=2, lty=2)
text(14, 2+0.3, cex=0.75, "Last Interglacial period
Data Constraint")
text(14, 2-0.3, cex=0.75, "Model Inversion
This study")

plotrange(windows[1,2], (windows[1,1] + windows[1,2])/2, windows[1,1], year=F, height=2.4, color="black")

put.fig.letter(label="a.", location="topleft", font=2)

par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(un.sf2100$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS Volume loss in 2100 [SLE m]", 
     ylab="Probability Density", 
     xlim=c(IPCC_AR5[1],Bamber[3]), ylim=c(-18,78), yaxt="n") 

lines(LIG.sf2100$pdf, col=mypalette[3], lwd=2) 
lines(LGM.sf2100$pdf, col=mypalette[5], lwd=2) 
lines(MH.sf2100$pdf, col=mypalette[1], lwd=2) 
lines(present.sf2100$pdf, col=mypalette[7], lwd=2)
lines(all.sf2100$pdf, col=mypalette[9], lwd=2)

probs = c(0.05, 0.95)
add.hor.box(all.prob_proj[,6], probs, width.size = 3.5, where.at = -3, tick.length = 0.25, line.width = 2, color = mypalette[9])
add.hor.box(present.prob_proj[,6], probs, width.size = 3.5, where.at = -6, tick.length = 0.25, line.width = 2, color = mypalette[7])
add.hor.box(MHprob_proj[,6], probs, width.size = 3.5, where.at = -9, tick.length = 0.25, line.width = 2, color = mypalette[1])
add.hor.box(LGMprob_proj[,6], probs, width.size = 3.5, where.at = -12, tick.length = 0.25, line.width = 2, color = mypalette[5])
add.hor.box(LIGprob_proj[,6], probs, width.size = 3.5, where.at = -15, tick.length = 0.25, line.width = 2, color = mypalette[3])
add.hor.box(unprob_proj[,6], probs, width.size = 3.5, where.at = -18, tick.length = 0.25, line.width = 2, color = "gray91")

abline(h=55, lty=2)
text(0.6, 55+4, cex=0.75, "Expert Assessments")
text(0.6, 55-6, cex=0.75, "Model Inversion
This study")
plotrange(Little[1], Little[2], Little[3], year=F, height=59, color="pink")
plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=63, color="purple")
plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=67, color="orange")
plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=71, color="gray")
plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=75, color="red")

put.fig.letter(label="b.", location="topleft", font=2)

plot.new()
legend("left", c("Unconstrained L.H. fits (500)","LIG constraint fits (285)", "LGM constraint fits (177)",
                 "MH constraint fits (299)","Instrumental constraint fits (18)","All L.H. constrained fits (2)",
                 "90% C.I. (Little et al. 2013)","90% C.I. (IPCC AR5)","90% C.I. (Kopp et al. 2014)", 
                 "90% C.I. (Pfeffer et al. 2008)","90% C.I. (Bamber & Aspinall 2013)"), cex=0.75, 
       lty=c(NA,NA,NA,NA,NA,NA,1,1,1,1,1), pch=c(15,15,15,15,15,15,8,8,8,8,8), lwd=1.5, bty="n", 
       col=c("gray91", mypalette[3], mypalette[5], mypalette[1], mypalette[7],mypalette[9],"pink",
             "purple","orange","gray","red"))
dev.off()

#------------------------------------- Supplementary Figure 3 -------------------------------------------
cummelt.mean.mcmcProj <- mean.dais.mcmcProj - mean(mean.dais.mcmcProj[SL.1961_1990])
meanMelt.rate <-rep(NA, enddate)
for(i in 1:enddate){
  meanMelt.rate[i] = cummelt.mean.mcmcProj[i+1] - cummelt.mean.mcmcProj[i]
}
 # x = AIS air temp; y = AIS ocean temp
degree = 3
fit.temps=lm(project.forcings[,2] ~ poly(project.forcings[,1], degree, raw=TRUE))
smoothoceantemp = predict(fit.temps)

fit.melt=lm(project.forcings[1:(enddate-1),2] ~ poly(meanMelt.rate[1:(enddate-1)], degree, raw=TRUE))
smoothoceanmelt = predict(fit.melt)

jpeg(file="nSuppFig6_dais_mcmcLHS.jpeg", family="Helvetica", width=1590, height=1200, units="px", pointsize=20)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
plot(project.forcings[,1], smoothoceantemp, xlab="Antarctica surface air temp. [C]",
     typ="l", ylab="High-latitude subocean temp. [C]")
points(project.forcings[,1], project.forcings[,2], col="skyblue", pch=20)
put.fig.letter(label="a.", location="topleft", font=2)

par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(meanMelt.rate[1:(enddate-1)], smoothoceanmelt, xlab="AIS melt rate [SLE m]", typ="l",
     ylab="High-latitude subocean temp. [C]")
points(meanMelt.rate, project.forcings[,2], col="skyblue", pch=20)
put.fig.letter(label="b.", location="topleft", font=2)

# 
# plot.new()
# legend("center", c("Maximum LIG year of AIS contribution ","2100", "Last interglacial period","Median MCMC fit"), 
#        pch=c(20,20,NA,NA), lty=c(NA,NA,1,1), bty="n", col=c("red","blue","red", "black"))
dev.off()

###################################### OTHER FIGURES ############################################
#------------------------------------- Other Figure 1 -------------------------------------------
#iceFlux_at_groundL = fo*((1-alpha)+alpha*((TO[i]-Tice)/(TOo-Tice))^2)/((s*Roa-bo)^(gamma-1))
iceFlux_at_groundL = median.dais.par[7]*((1-median.dais.par[2])+median.dais.par[2]*
                                           ((project.forcings[,2]-Tice)/(TOo-Tice))^2)/((s*Roa-bo)^(median.dais.par[1]-1))

# The max estimate in the 90% credible interval during the LIG
ais_during_LIG <- mcmc_95[90000:140000]
max.LIG <- max(ais_during_LIG)
max.year <- which.max(ais_during_LIG)
max.year <- 90000 + max.year

pdf(file="OtherFig1_dais_mcmcLHS.pdf", family="Helvetica", height=2.7, width=6.7,pointsize=9)
par(mfrow=c(1,2),mgp=c(3, 1, 0), mar=c(5,4,1,1))
plot(project.forcings[,2], iceFlux_at_groundL, xlab="High-latitude ocean subsurface 
     temperature anomaly [C]", typ="l", ylab="Ice flux across the grounding line")
points(project.forcings[240100,2], iceFlux_at_groundL[240100],pch=20, col="blue")
points(project.forcings[max.year,2], iceFlux_at_groundL[max.year],pch=20, col="red")
abline(v=1, h=0.0016, lty=2)
text(3, 0.0016+0.0002, cex=0.75, "Missing Marine Ice 
     Sheet Instability (MISI)")
#text(2.5, 0.0016+0.0001, cex=0.75, "Missing MISI")
#put.fig.letter(label="c.", location="topleft", font=2)

plot.new()
legend("center", c("Maximum LIG year of AIS contribution ","2100", "Median MCMC fit"), 
       pch=c(20,20,NA), cex=0.75, lty=c(NA,NA,1), bty="n", col=c("red","blue","black"))
dev.off()

#------------------------------------- Other Figure 2 -------------------------------------------
pdf(file="OtherFig2_dais_mcmcLHS.pdf", family="Helvetica", height=5.4, width=6.7,pointsize=11)
par(mfrow=c(2,2),mgp=c(1.5,.5,0), mar=c(3.5,4,1,1)) # set figure dimensions
plot(date, iceFlux_at_groundL, xlab="Date [Kyr BP]", typ="l", ylab="Ice flux across the grounding line", xaxt="n")
abline(h=0.0016, lty=2)
lines(c(-130000,-110000), c(0.0016, 0.0016), lwd=2,col="red")
ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)
points(date[max.year], iceFlux_at_groundL[max.year],pch=20, col="red")
# text(-185000, 0.0016+0.0002, cex=0.75, "Missing Marine Ice 
# Sheet Instability (MISI)")
text(-185000, 0.0016+0.0001, cex=0.75, "Missing MISI")

par(mgp=c(1.5,.5,0), mar=c(3.5, 3, 1, 2))
plot(date[240000:240300], iceFlux_at_groundL[240000:240300], xlab="Date [yr BP]", typ="l", xaxt="n",
     ylab="Ice flux across the grounding line", ylim=c(0.0010,0.0035))
abline(h=0.0016, lty=2)
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="b.", location="topleft", font=2)
points(date[240100], iceFlux_at_groundL[240100],pch=20, col="blue")
text(175, 0.0016+0.0001, cex=0.75, "Missing MISI")

plot.new()
legend("center", c("Maximum LIG year of AIS contribution ","2100", "Last interglacial period","Median MCMC fit"), 
       pch=c(20,20,NA,NA), lty=c(NA,NA,1,1), bty="n", col=c("red","blue","red", "black"))
dev.off()

###################################### END ############################################
