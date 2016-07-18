##################################################################################################
#
#  -file = "Paper2_plots.R"   Code written September 2014
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

#install.packages('RColorBrewer')
library(RColorBrewer)
mypalette <- brewer.pal(9,"YlGnBu")

source("Scripts/put_fig_letter.r")
source("Scripts/plot_rangefn.R")

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
  mcmc_point5 <-
  mcmc_995 <-rep(NA,enddate) # 240300 years
for(i in 1:enddate){
  mcmc_5[i] = quantile(proj.mcmc.1961_1990[,i],0.05) #heteroskedastic SLR values
  mcmc_95[i] = quantile(proj.mcmc.1961_1990[,i],0.95)
  
  mcmc_point5[i] = quantile(proj.mcmc.1961_1990[,i],0.005) # Bootstrap SLR values
  mcmc_995[i] = quantile(proj.mcmc.1961_1990[,i],0.995)
}

x_mcmc_90=c(mcmc_5, rev(mcmc_95)); y_mcmc_90=c(date, rev(date))

#------------------------------------- Figure 1 -------------------------------------------
# MCMC hindcasts & projection

pdf(file="Fig1_dais_mcmcLHS.pdf", family="Helvetica", width=6.7, height=8.1, pointsize=12)
par(mfrow=c(3,2), mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))

# Last interglacial 240 kyr Bp - 2010 AD
plot(date[1], proj.mcmc.1961_1990[1,1], type="l", col="powderblue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]", 
     xlim=c(date[1], date[240010]), ylim=c(-16,7), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)
abline(h=0, lty=3, col="black")

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray", lwd=0.75)
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
ticks=c(-200000,-150000,-100000,-50000,0)
axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
put.fig.letter(label="a.", location="topleft", font=2)

# Last glacial maximum 25 kyr Bp - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[kyrbp_25], proj.mcmc.1961_1990[1,kyrbp_25], type="l", col="powderblue", 
     lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]",
     xlim=c(date[kyrbp_25], date[240010]), ylim=c(windows[2,1],0.5), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=2)
lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=2)
ticks=c(-25000,-20000,-15000,-10000,-5000,0)
axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
put.fig.letter(label="b.", location="topleft", font=2)

# Mid-Holocene6 kyr Bp - 2010 AD 
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
plot(date[kyrbp_6], proj.mcmc.1961_1990[1,kyrbp_6], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", 
     xlim=c(date[kyrbp_6], date[240010]),ylim=c(-4,0.025))

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=2)
lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=2)
lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=2)
put.fig.letter(label="c.", location="topleft", font=2)

# Present day 1880 - 2010 AD
par(mgp=c(1.5,.5,0), mar=c(4, 3, 2, 2))
plot(date[AD_1880], proj.mcmc.1961_1990[1,AD_1880], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", 
     xlim=c(date[AD_1880], date[240010]),ylim=c(-0.01,0.01), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
err_positive = windows[4,2]
err_negative = windows[4,1]
present = 2
arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(-120,-100,-50,0)
axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
put.fig.letter(label="d.", location="topleft", font=2)

# Future 2010 - 2300 AD
par(mgp=c(1.5,.5,0), mar=c(4, 4, 2, 1))
plot(date[240000], proj.mcmc.1961_1990[1,240000], type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", 
     xlim=c(-10,290),ylim=c(-0.01,1.8), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
err_positive = windows[4,2]
err_negative = windows[4,1]
arrows(present, err_positive, present, err_negative, length=0, lwd=2, col="black")
ticks=c(0,50,100,150,200,250,300)
axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
put.fig.letter(label="e.", location="topleft", font=2)

plot.new()
legend("center", c("Reconstruction/observation", "Median MCMC fit", "90% credible interval","99% credible interval"),
       lty=c(1,1,NA,2), pch=c(NA,NA,15,NA), col=c("black", "dimgray", "goldenrod2", "black"))
dev.off()
#------------------------------------- Figure 2 -------------------------------------------
# projection/pdf in 2100

#Set up previous study ranges for the year 2100 using RCP 8.5
#the previous ranges are in cm so divide by 100 to get meters
pfeffer = c(12.8, 14.6, 61.9)/100 #Low2, low1, and high estimate in Pfeffer et al. 2008
Bamber = c(-2, 14, 83)/100 # 5%, Median, 95% estimates in Bamber and Aspinall 2013
IPCC_AR5 = c(-15, 4, 23)/100 # 5%, Median, 95% estimates in IPCC AR5
Little = c(-8, 2.4, 13.3)/100 # 5%, Median, 95% estimates in Little et al. 2013
Kopp = c(-11, 4, 33)/100 # 5%, Median, 95% estimates in Kopp et al. 2014
present = 2

pdf(file="Fig2a_dais_mcmcLHS.pdf", family="Helvetica", width=6.7, height=2.7, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1)) # set figure dimensions

plot(all.sf2100$pdf, main="", col=mypalette[9],lwd=2,
     xlab="Projected AIS Volume loss in SLE in 2100 [m]", 
     ylab="Probability Density [Dimensionless]", 
     xlim=c(IPCC_AR5[1],Bamber[3]), ylim=c(0,78)) 
lines(sf.2100.mcmc$pdf, col="goldenrod2", lwd=2)

probs = c(0.05, 0.95)
add.hor.box(all.prob_proj[,6], probs, width.size = 3.5, where.at = 55, tick.length = 0.25, line.width = 2, color = mypalette[9])
add.hor.box(mcmc.prob_proj[,6], probs, width.size = 3.5, where.at = 58, tick.length = 0.25, line.width = 2, color = "goldenrod2")

abline(h=60, lty=2)
text(0.6, 60+2, cex=0.75, "Expert Assessment")
plotrange(Little[1], Little[2], Little[3], year=F, height=62, color="pink")
plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=65, color="purple")
plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=68, color="orange")
plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=71, color="gray")
plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=74, color="red")

plot.new()
par(mgp=c(1.5,.5,0), mar=c(4, 3, 3, 2)) # set figure dimensions
legend("left", c("L.H. all constraints (5)", "MCMC","90% C.I. Little et al. 2013","90% C.I. IPCC AR5", 
                 "90% C.I. Kopp et al. 2014", "90% C.I. Pfeffer et al. 2008","90% C.I. Bamber & Aspinall 2013"),
       lty=c(NA,NA,1,1,1,1,1), pch=c(15,15,NA,NA,NA,NA,NA), lwd=2, bty="n", col=c(mypalette[9],"goldenrod2","pink","purple","orange","gray","red"))
dev.off()

#------------------------------------- Figure 2b -------------------------------------------

pdf(file="Fig2b_dais_mcmcLHS.pdf", family="Helvetica", width=6.7, height=2.7, pointsize=9)
par(mfrow=c(1,2), mgp=c(1.5,.5,0), mar=c(4, 4, 3, 1)) # set figure dimensions
plot(un.sf2100$pdf, main="",lwd=3, col="gray91", xlab="Projected AIS Volume loss in SLE in 2100 [m]", 
     ylab="Probability Density [Dimensionless]", 
     xlim=c(IPCC_AR5[1],Bamber[3]), ylim=c(0,40)) 

lines(LIG.sf2100$pdf, col=mypalette[3], lwd=2) 
lines(LGM.sf2100$pdf, col=mypalette[5], lwd=2) 
lines(MH.sf2100$pdf, col=mypalette[1], lwd=2) 
lines(present.sf2100$pdf, col=mypalette[7], lwd=2)

probs = c(0.05, 0.95)
add.hor.box(LIGprob_proj[,6], probs, width.size = 3.5, where.at = 21, tick.length = 0.25, line.width = 2, color = mypalette[3])
add.hor.box(LGMprob_proj[,6], probs, width.size = 3.5, where.at = 19, tick.length = 0.25, line.width = 2, color = mypalette[5])
add.hor.box(MHprob_proj[,6], probs, width.size = 3.5, where.at = 16, tick.length = 0.25, line.width = 2, color = mypalette[1])
add.hor.box(present.prob_proj[,6], probs, width.size = 3.5, where.at = 13, tick.length = 0.25, line.width = 2, color = mypalette[7])

abline(h=23, lty=2)
text(0.6, 23+2, cex=0.75, "Expert Assessment")
plotrange(Little[1], Little[2], Little[3], year=F, height=25, color="pink")
plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=F, height=28, color="purple")
plotrange(Kopp[1], Kopp[2], Kopp[3], year=F, height=31, color="orange")
plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=F, height=34, color="gray")
plotrange(Bamber[1], Bamber[2], Bamber[3], year=F, height=37, color="red")

plot.new()
par(mgp=c(1.5,.5,0), mar=c(4, 3, 3, 2)) # set figure dimensions
legend("left", c("Unconstrained L.H. fits (500)","LIG best fits (239)", "LGM best fits (223)","MH best fits (230)",
                 "Instrumental best fits (102)","90% C.I. Little et al. 2013","90% C.I. IPCC AR5", 
                 "90% C.I. Kopp et al. 2014", "90% C.I. Pfeffer et al. 2008","90% C.I. Bamber & Aspinall 2013"),
       lty=c(NA,NA,NA,NA,NA,1,1,1,1,1), pch=c(15,15,15,15,15,NA,NA,NA,NA,NA), lwd=2, bty="n", col=c("gray91", mypalette[3], mypalette[5], 
                                                                                                    mypalette[1], mypalette[7],"pink","purple","orange","gray","red"))
dev.off()

#---------------------------------Calculate max estimates in LIG and 2100------------------------------------------------
# The max estimate in the 90% credible interval during the LIG
ais_during_LIG <- mcmc_95[90000:140000]
max.LIG <- max(ais_during_LIG)
max.year <- which.max(ais_during_LIG)
max.year <- 90000 + max.year

# Check max estimate (they should both be the same number)
print(max.LIG)
print(mcmc_95[max.year])

# The max estimate in the 90% credible interval during the LIG
max.2100 <- mcmc_95[240100]

someColor = "red"
trans = makeTransparent(someColor,100)
#------------------------------------- Figure 3 -------------------------------------------

pdf(file="Fig3_dais_mcmcLHS.pdf", family="Helvetica", width=10, height=2.7, pointsize=12)
par(mfrow=c(1,3), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 0))
plot(date[110000], AIS_melt[90000]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
     xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]", 
     xlim=c(date[110000], date[140000]),ylim=c(-4.5,8), xaxt="n")

polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)

lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")
lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=2)
lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=2)
lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=2)
points(last.interglacial[2], (windows[1,1] + windows[1,2])/2, pch=19)

# add in disrepancy box
LIG_x=c(rep(max.LIG, length(last.interglacial)), rep(windows[1,2], length(last.interglacial)))
LIG_y=c(last.interglacial, rev(last.interglacial))
polygon(LIG_y, LIG_x, col=trans)

ticks=c(-125000,-120000,-115000,-110000,-105000)
axis(side=1, at=ticks, labels=expression(-125,-120,-115,-110,-105))

par(mgp=c(1.5,.5,0), mar=c(4, 0, 3, 4)) # set figure dimensions
plot(date[240090], Project_melt[240090]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
     xlab="Date [yr BP]", xlim=c(date[240090], date[240110]),
     ylim=c(IPCC_AR5[1],Bamber[3]), xaxt="n", yaxt="n") # ylab="AIS Volume loss in SLE [m]",ylim=c(-0.5,1), xaxt="n")
polygon(y_mcmc_90, x_mcmc_90, col="goldenrod2")
lines(date, mcmc_point5, lty=2)
lines(date, mcmc_995, lty=2)
lines(date, new.dais.mcmc.proj-mean(new.dais.mcmc.proj[SL.1961_1990]), col="dimgray")

ticks=c(90,95,100,105,110)
axis(side=1, at=ticks, labels=expression(2090,2095,2100,2105,2110))
axis(4)

# add in disrepancy box
x_2100=c(rep(max.2100, length(98:103)), rep(Bamber[3], length(98:103)))
y_2100=c(98:103, rev(98:103))
polygon(y_2100, x_2100, col=trans)

#add expert judgement
plotrange(Little[1], Little[2], Little[3], year=98, height=F, color="pink")
plotrange(IPCC_AR5[1], IPCC_AR5[2], IPCC_AR5[3], year=99, height=F, color="purple")
plotrange(Kopp[1], Kopp[2], Kopp[3], year=101, height=F, color="orange")
plotrange(pfeffer[1], pfeffer[2], pfeffer[3], year=102, height=F, color="gray")
plotrange(Bamber[1], Bamber[2], Bamber[3], year=103, height=F, color="red")

plot.new()
legend("left", c("LIG reconstruction", "Median MCMC fit", "99% C.I.","90% C.I.",
                    "Disrepency","90% C.I. Little et al. 2013","90% C.I. IPCC AR5", 
                 "90% C.I. Kopp et al. 2014", "90% C.I. Pfeffer et al. 2008", 
                 "90% C.I. Bamber & Aspinall 2013"),
       lty=c(1,1,2,NA,NA,1,1,1,1,1), pch=c(NA,NA,NA,15,15,NA,NA,NA,NA,NA), col=c("black", "dimgray", "black", "goldenrod2",trans,"pink",
                                                                           "purple","orange","gray","red"))

dev.off()
###################################### SUPPLEMENTARY FIGURES ############################################
#------------------------------------- Supplementary Figure 1 -------------------------------------------

# # LHS hindcasts & projection
# 
# jpeg(file="Fig1_dais_mcmcLHS.jpeg", family="Helvetica", width=700, height=600, units="px", pointsize=20)
# par(mfrow=c(3,2), mgp=c(1.5,.5,0),mar=c(4, 4, 3, 2))
# 
# # Last interglacial 240 kyr Bp - 2010 AD
# plot(date[1:240010], AIS_melt-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
#      xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(-30,10), xaxt="n")
# for(i in 1:sample_length){
#   lines(date, dais.pre.cali[i,], col="gray91", lwd=2)
# }
# for(i in 1:length(surMH)){
#   lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=2)
# }
# for(i in 1:length(surLIG)){
#   lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=2)
# }
# for(i in 1:length(surLGM)){
#   lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=2)
# }
# for(i in 1:length(sur9311trend)){
#   lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=2)
# }
# for(i in 1:length(sur.all)){
#   lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=2)
# }
# abline(h=0, lty=2, col="black")
# lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=2)
# lines(last.interglacial, c(windows[1,1],windows[1,1],windows[1,1]), col="black", lwd=4)
# lines(last.interglacial, c(windows[1,2],windows[1,2],windows[1,2]), col="black", lwd=4)
# lines(c(last.interglacial[2],last.interglacial[2]), c(windows[1,1],windows[1,2]), col="black", lwd=4)
# ticks=c(-200000,-150000,-100000,-50000,0)
# axis(side=1, at=ticks, labels=expression(-200,-150,-100,-50,0))
# put.fig.letter(label="a.", location="topleft", font=2, offset=c(0.05, -0.1))
# 
# # Last glacial maximum 25 kyr Bp - 2010 AD
# plot(date[kyrbp_25:240010], AIS_melt[kyrbp_25:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", 
#      lwd=2,xlab="Date [Kyr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(-30,7), xaxt="n")
# for(i in 1:sample_length){
#   lines(date, dais.pre.cali[i,], col="gray91", lwd=2)
# }
# for(i in 1:length(surMH)){
#   lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=2)
# }
# for(i in 1:length(surLIG)){
#   lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=2)
# }
# for(i in 1:length(surLGM)){
#   lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=2)
# }
# for(i in 1:length(sur9311trend)){
#   lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=2)
# }
# for(i in 1:length(sur.all)){
#   lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=2)
# }
# lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=2)
# # abline(h=0, lty=2, col="gray")
# lines(c(last.glacialmax[2],last.glacialmax[2]), c(windows[2,1],windows[2,2]), col="black", lwd=4)
# lines(last.glacialmax, c(windows[2,1],windows[2,1],windows[2,1]), col="black", lwd=4)
# lines(last.glacialmax, c(windows[2,2],windows[2,2],windows[2,2]), col="black", lwd=4)
# ticks=c(-25000,-20000,-15000,-10000,-5000,0)
# axis(side=1, at=ticks, labels=expression(-25,-20,-15,-10,-5,0))
# put.fig.letter(label="d.", location="topleft", font=2, offset=c(0.05, -0.1))
# 
# # Mid-Holocene6 kyr Bp - 2010 AD 
# plot(date[kyrbp_6:240010], AIS_melt[kyrbp_6:240010]-mean(AIS_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
#      xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(-8,4))
# for(i in 1:sample_length){
#   lines(date, dais.pre.cali[i,], col="gray91", lwd=2)
# }
# for(i in 1:length(surLIG)){
#   lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=2)
# }
# for(i in 1:length(surLGM)){
#   lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=2)
# }
# for(i in 1:length(sur9311trend)){
#   lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=2)
# }
# for(i in 1:length(surMH)){
#   lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=2)
# }
# for(i in 1:length(sur.all)){
#   lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=2)
# }
# lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=2)
# lines(c(holocene[2],holocene[2]), c(windows[3,1],windows[3,2]), col="black", lwd=4)
# lines(holocene, c(windows[3,1],windows[3,1],windows[3,1]), col="black", lwd=4)
# lines(holocene, c(windows[3,2],windows[3,2],windows[3,2]), col="black", lwd=4)
# put.fig.letter(label="b.", location="topleft", font=2, offset=c(0.05, -0.1))
# 
# # Present day 1880 - 2010 AD
# plot(date[AD_1880:240010], Project_melt[AD_1880:240010]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
#      xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", ylim=c(-0.10,0.05), xaxt="n")
# for(i in 1:sample_length){
#   lines(date, dais.pre.cali[i,], col="gray91", lwd=2)
# }
# for(i in 1:length(surMH)){
#   lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=2)
# }
# for(i in 1:length(surLIG)){
#   lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=2)
# }
# for(i in 1:length(surLGM)){
#   lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=2)
# }
# for(i in 1:length(sur9311trend)){
#   lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=2)
# }
# for(i in 1:length(sur.all)){
#   lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=2)
# }
# lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=2)
# err_positive = windows[4,2]
# err_negative = windows[4,1]
# present = 2
# arrows(present, err_positive, present, err_negative, length=0, lwd=4, col="black")
# ticks=c(-120,-100,-50,0)
# axis(side=1, at=ticks, labels=expression(1880,1900,1950,2000))
# put.fig.letter(label="e.", location="topleft", font=2, offset=c(0.05, -0.1))
# 
# # Future 2010 - 2300 AD
# plot(date[240000:enddate], Project_melt[240000:enddate]-mean(Project_melt[SL.1961_1990]), type="l", col="powderblue", lwd=2,
#      xlab="Date [yr BP]", ylab="AIS Volume loss in SLE [m]", xlim=c(-10,290),ylim=c(-0.5,5), xaxt="n")
# for(i in 1:sample_length){
#   lines(date, dais.pre.cali[i,], col="gray91", lwd=2)
# }
# for(i in 1:length(surMH)){
#   lines(date, dais.pre.cali[surMH[i],], col=mypalette[1], lwd=2)
# }
# for(i in 1:length(surLIG)){
#   lines(date, dais.pre.cali[surLIG[i],], col=mypalette[3], lwd=2)
# }
# for(i in 1:length(surLGM)){
#   lines(date, dais.pre.cali[surLGM[i],], col=mypalette[5], lwd=2)
# }
# for(i in 1:length(sur9311trend)){
#   lines(date, dais.pre.cali[sur9311trend[i],], col=mypalette[7], lwd=2)
# }
# for(i in 1:length(sur.all)){
#   lines(date, dais.pre.cali[sur.all[i],], col=mypalette[9], lwd=2)
# }
# lines(date, best.project-mean(best.project[SL.1961_1990]), col="coral3", lwd=2)
# err_positive = windows[4,2]
# err_negative = windows[4,1]
# arrows(present, err_positive, present, err_negative, length=0, lwd=4, col="black")
# ticks=c(0,50,100,150,200,250,300)
# axis(side=1, at=ticks, labels=expression(2000,2050,2100,2150,2200,2250,2300))
# put.fig.letter(label="c.", location="topleft", font=2, offset=c(0.05, -0.1))
# 
# plot.new()
# legend("center", c("Reconstruction/observation", "Optim best-fit", "Last interglacial constraint", "Last glacial maximum constraint","Mid-Holocene constraint", 
#                    "Instrumental period (1993-2011) constraint", "All constraints"),lty=c(1,1,NA,NA,NA,NA,NA,NA), pch=c(NA,NA,15,15,15,15,15,15),
#        col=c("black", "black", mypalette[3], mypalette[5], mypalette[1], mypalette[7], mypalette[9]))
# dev.off()

#------------------------------------- Supplementary Figure 2 -------------------------------------------

#------------------------------------- Supplementary Figure 3 -------------------------------------------

#------------------------------------- Supplementary Figure 4 -------------------------------------------

#------------------------------------- Supplementary Figure 5 -------------------------------------------

###################################### END ############################################