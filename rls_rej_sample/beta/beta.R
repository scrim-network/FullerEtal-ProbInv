############################################
##  file: main.R
##   - Application to sea level rise (SLR) projections 
#
######################################################
##  Authors: Klaus Keller (klaus@psu.edu) and Ryan Sriver 
##  copyright by the authors
######################################################
##  distributed under the GNU general public license
##  no warranty
##  This is PRELIMINARY RESEARCH CODE, USE AT YOUR OWN RISK AND ONLY
##  AFTER YOU CONVINCED YOURSELF THAT THIS IS OK TO USE FOR YOUR APPLICATION 
##################################################
# version: 3
# last changes:
# August 22 2011
#/Users/klaus/Desktop/wip/sea-level-rise/windowing_analysis/codes_Ryan_Aug_22/beta_ann_negc
# Mon Aug 22 12:05:47 EDT 2011
##################################################
# how to run:
# - save the file in a directory
# - go to the directory with this file
# - open R
# - install packages "timeDate" and "timeSeries" if needed  
# - type 'source('main.R')' inside the R prompt
# - look at the ./output directory
# - open the new pdf file(s) to analyze the results
# - the file array_raw.txt is a table of the resulting model parameters for the prior 
# - the file array_cauch.txt is a table of the resulting model parameters for the posterior 
#	- after rejection resampling to fit the chain to theoretical slr distribution
# - the entries are (index, a, b, c, tstar, cstar, projected slr in 2100 in mm)
############################################

rm(list = ls(all = TRUE))
set.seed(1) # set the seed to be reproducible

library(timeSeries)

# Define the number of bootstrap/mc samples
N=55000  # Used throughout analysis

# read in the global data from Jevrejava
# see /Users/klaus/Desktop/wip_local/sea-level/data
raw.global<- scan(file="../data/jevrejava_yearly.txt",
        what=numeric(), quiet=T)
slr.global.data <- matrix(raw.global, ncol=2, byrow=T)

# fit a simple polynominal model to the data and extract the year 2000 prediion
years.global=slr.global.data[,1]
years.global=years.global-(2011)  # normalize time series to current year 2011
  nyears.global=length(years.global)
slr.global=slr.global.data[,2]

# 2nd order polynomial ~ a + bx + cx^2
fit.global=lm(slr.global ~ years.global + I(years.global^2))

# scale the observations such that the slr anomalie is zero in the year 2000
predict.global.2000=predict(fit.global,newdata=data.frame(years.global=-11))
slr.global=slr.global-predict.global.2000		# data
predict.hind=predict(fit.global)-predict.global.2000	# poly fit

### Calculate Residuals during observed time series (data - polynomial fit)  ###
res=slr.global-predict.hind
res_stdev=sd(res[length(res)-50:length(res)]) # stdev for the final 50 years of data
print(res_stdev)

# write-out the residuals
residuals.annual=mat.or.vec(2,nyears.global) #(nr,nc)
residuals.annual[1,]=years.global+2011
residuals.annual[2,]=res
write.table(residuals.annual,
            file="./output/residuals.annual.ascii",
            row.names=FALSE,
            col.names=FALSE
            )
           

### Bootstrap the residuals ###
### confine to post-1900 portion of time series ###
white.boot = mat.or.vec(N, nyears.global) # create matrix (nr,nc)
white.boot_sd = rep(NA,N)

  for(i in 1:N) {
    white.boot[i,1:nyears.global] = sample(res[94:nyears.global],size=nyears.global,replace=TRUE)
    white.boot_sd[i] = sd(white.boot[i,])
  }

### create new residuals from bootstraps with original AR structure
pac <- pacf(res[94:nyears.global], lag.max=5, plot=FALSE)  # apply partial auto-correlation to determine correlation coefficients
  print(pac$acf)

res.boot=mat.or.vec(N, nyears.global) #(nr,nc)

  for(n in 1:N) {
    for(i in 2:nyears.global) {
        res.boot[n,i] = pac$acf[1]*res.boot[n,i-1] + rnorm(1,mean=0,sd=white.boot_sd[n]) 
    }
  }

# write out boot strap samples of auto-correlated residuals
#res.boot.t=t(res.boot)
#res.boot.array=mat.or.vec(nyears.global,N+1) #(nr,nc)
# res.boot.array[,1]=years.global+2011
# res.boot.array[,2:1001]=res.boot.t
#write.table(res.boot.array,file="./output/residual_sample.txt",quote=FALSE,row.names=FALSE,col.names=FALSE)

# write out just a simple time-series
#res.boot.simple=mat.or.vec(nyears.global,N) #(nr,nc)
#res.boot.simple[1:N]=res.boot
#write.csv(res.boot.simple,file="./output/just_residuals_sample.csv")               


### apply new AUTOCORRELATED residuals to polynomial fit
slr.boot=mat.or.vec(N, nyears.global) #(nr,nc)

  for(i in 1:N) {
    slr.boot[i,]=predict.hind+res.boot[i,]
  }

### plot observations and polynomial fit with superimposed autocorrelated bootstraps ###
pdf(file="./output/hindcast.pdf")  # write to pdf, define a pdf file to write to
plot(years.global+2011,slr.global,type="p",col="green",pch=1,
     xlab=("Year"),
     ylab=("mean sea-level anomaly (mm) with respect to the year 2000"),
     xlim=c(1800,2000),
     ylim=c(-400,50))
title("")

# add bootstrap samples
nexamp=3
nrand=sample(N,nexamp)
  for(i in 1:nexamp) {
    lines(years.global+2011, slr.boot[nrand[i],1:nyears.global], col="red", lwd=1)
  }
points(years.global+2011,slr.global,type="p",col="green",pch=1)
lines(years.global+2011,predict.hind,col="black",lty="solid")

legend(x="topleft",
       legend=c("Observations (Jevrejava et al (2006)", "polynomial best fit", "best fit + bootstrap residuals"),
       col=c("green","black","red"),
       lty=c("solid","solid","solid")
       )

dev.off()

# define the extrapolation, predict, and center on 2000
years.extra=seq(-211,89, length=300) # using jevrejava annual data 1800-2100 (t0=2011)
n.years.extra=length(years.extra)
predict.global=predict(fit.global,newdata=data.frame(years.global=seq(-211,89,length=300)))-predict.global.2000


### calculate polynomial coefficients for projections ###
boot.fit_coef=mat.or.vec(N, 3)
boot.fit_predict=mat.or.vec(N, n.years.extra)
for(i in 1:N) {
 fit.new=lm(slr.boot[i,] ~ years.global + I(years.global^2))
  boot.fit_coef[i,1]=fit.new$coefficients[1]
  boot.fit_coef[i,2]=fit.new$coefficients[2]
  boot.fit_coef[i,3]=fit.new$coefficients[3]
  boot.fit_predict[i,]=predict(fit.new,newdata=data.frame(years.global=seq(-211,89,length=300))) - 
	predict(fit.new,newdata=data.frame(years.global=-11))
  }


# define the scenario models
slrmodel <- function(times, slr.normal, c, t.star)
  # this function adds the scenarios for unresolved processes using a
  # time t.star at which the change happen and an additional sea level
  # rate c after t.star
  # example call: foo=slrmodel(years.extra,predict.global,1,2050)
{
  n.points=length(times)
  slr.scenario=slr.normal
   for(i in 1:n.points)
     {
      if (times[i] > t.star)
        {
        slr.scenario[i]=slr.normal[i]+(times[i]-t.star)*c
        }
     }
    return (slr.scenario)
}

# produce an example scenario projections
scenario.global=slrmodel(years.extra,predict.global,17,4)


# produce many slr scenarios drawing from the prior
###################################################
c.sample=mat.or.vec(N,1) #(nr,nc)
t.star.sample=mat.or.vec(N,1)
scenario.mc=mat.or.vec(N,n.years.extra)
for(n in 1:N)
  {
  # draw random samples from a uniform distribution
  c.sample[n]<-runif(1, min=-15, max=35)
  t.star.sample[n]<-runif(1, min=0, max=89)
  # produce the scenario
  scenario.mc[n,]=slrmodel(years.extra,boot.fit_predict[n,],c.sample[n],t.star.sample[n])
  }

scenario.poly.only=scenario.mc   # define the polynomial projections (used as input for pdfs)


### calculate projected residuals ###
res.boot_proj=mat.or.vec(N, n.years.extra) #(nr,nc)

 for(n in 1:N) {
    for(i in 2:n.years.extra) {
        res.boot_proj[n,i] = pac$acf[1]*res.boot_proj[n,i-1]  + rnorm(1,mean=0,sd=white.boot_sd[n])
    }
  }



### superimpose residuals on polynomial fit projections ###
  for(i in 1:N) {
    scenario.mc[i,]=scenario.mc[i,]+res.boot_proj[i,]
  }



# cut the scenarios within the window
ub.2100=2508
lb.2100=255
# see notes Aug. 16th.
#(pfeffer + thermosteric + reg. uncert.) )

# define arrays for all parameters and slr time series
array=mat.or.vec(N,6)
colnames(array) = c("a","b","c","t.star","c.star","slr")
scenario.filter=mat.or.vec(N,n.years.extra)

for(n in 1:N)
  {
    if (scenario.poly.only[n,n.years.extra] < ub.2100 && scenario.poly.only[n,n.years.extra] > lb.2100)
      {
       array[n,1:3]=boot.fit_coef[n,]
       array[n,4]=t.star.sample[n] + 2011
       array[n,5]=c.sample[n]
       array[n,6]=scenario.poly.only[n,n.years.extra]

       scenario.filter[n,]=scenario.mc[n,]
      } 
    else 
     {
      array[n,]=NA	# set out of bound entries to missing values
      scenario.filter[n,]=NA
     }
  }  


array=removeNA(array)
scenario.filter=removeNA(scenario.filter)  # remove missing values
  print(array)


# projections figure
pdf(file="./output/projections.pdf")  
plot(years.global+2011,slr.global,type="p",col="green",pch=1,
     xlab=("Year"),
     ylab=("mean sea-level anomaly (mm) with respect to the year 2000"),
     xlim=c(1800,2110),
     ylim=c(-300,2600))
title("")

# add projected bootstrap samples
nexamp=30
nrand=sample(N/2,nexamp)  # decrease sampling due to cutting off outliers
  for(i in 1:nexamp) {
    lines(years.extra+2011, scenario.filter[nrand[i],], col="grey", lwd=0.5)
  }
points(years.global+2011,slr.global,type="p",col="green",pch=1)
lines(years.extra+2011,predict.global,col="black",lty="solid")

# add CA scenarios
#arrows(2100,780,2100,1768,col="gray",lty="solid",code=3);
# add sriver et al version of the Pfeffer et al scenarios
arrows(2100,lb.2100,2100,ub.2100,col="blue",lty="solid",code=3);

# add the zero line top guide the eye
abline(h=0,lty="dotted")

legend(x="topleft",
       legend=c("Observations (Jevrejava et al (2006)",
	'polynomial best fit projection',
         'model scenarios',
         'Pfeffer et al (2008), expanded'),
       col=c("green","black","grey","blue"),
       lty=c("solid","solid","solid","solid")
       )
dev.off()

write.table(array,file="./output/array_raw.txt",quote=FALSE,col.names=FALSE)

###############################################
### Thin the model output to fit theoretical
### distribution of 2100 slr
###   - using rejection sampling
################################################


###  beta distriution  ###

a=lb.2100 # define slr bounds
b=ub.2100
n.bin=50  # number of bins for rejection sampling

s1 = 2    # define shape parameters (yields good fit to expert priors -- pfeffer, sriver, etc.)
s2 = 3

xs <- seq(0, 1, length=n.bin)	# define standard beta function on interval [0,1]
xt <- seq(a, b, length=n.bin)   # transform distribution to slr bounds

beta_pdf<-dbeta(xs, s1, s2)
  print(beta_pdf)

### check that pdf integrates to 1 ###
beta_wt_test=beta_pdf*(xs[2]-xs[1])  # standard distribution
  print(sum(beta_wt_test))
#beta_wt=beta_pdf*(xt[2]-xt[1])/(b-a) # transformed distribution
beta_wt=beta_pdf/sum(beta_pdf)
  print(sum(beta_wt))

# get the uniform pdf
x.sample.CA=seq(0,2500,length=1000)
lb.CA=10
ub.CA=2060
uniform.fit.CA=dunif(x.sample.CA, min=lb.CA, max=ub.CA, log = FALSE)

################################################

### read in model output ###
array = read.table("./output/array_raw.txt")
x.bin=n.bin-1

### read through slr and bin the values accordingly ###
nrows=length(array[,1])
bin.array=array(NA, dim=c(x.bin,nrows,7))

for (i in 1:x.bin)
{
 for(n in 1:nrows)                                         
  {
    if (array[n,7] > xt[i]  && array[n,7] < xt[i+1])
      {
       bin.array[i,n,1]=array[n,1]
       bin.array[i,n,2]=array[n,2]
       bin.array[i,n,3]=array[n,3]
       bin.array[i,n,4]=array[n,4]
       bin.array[i,n,5]=array[n,5]
       bin.array[i,n,6]=array[n,6]
       bin.array[i,n,7]=array[n,7]
      }
  }
}


print(bin.array[1,,])

### thin the binned data to fit the beta distribution ###

### test on first sample ###
bin.sub=bin.array[1,,]
bin.sub=removeNA(bin.sub)
  print(bin.sub)

sample.size=beta_wt[1]*nrows
  print(sample.size)
sample.ind=sample(bin.sub[,1],sample.size+1,replace=TRUE)  # add one sample to ensure no empty bins 
  print(sample.ind)
beta.dist=array(NA, dim=c(length(sample.ind),6))
colnames(beta.dist) = c("a","b","c","t.star","c.star","slr")
  print(beta.dist)

for(i in 1:length(sample.ind))
  {
#   beta.dist[i,1]=array[sample.ind[i],1]
   beta.dist[i,1]=array[sample.ind[i],2]
   beta.dist[i,2]=array[sample.ind[i],3]
   beta.dist[i,3]=array[sample.ind[i],4]
   beta.dist[i,4]=array[sample.ind[i],5]
   beta.dist[i,5]=array[sample.ind[i],6]
   beta.dist[i,6]=array[sample.ind[i],7]
  }
  print(beta.dist)

### repeat for rest of chain ###
for (n in 2:x.bin)
  {
   bin.sub_tmp=bin.array[n,,]
   bin.sub_tmp=removeNA(bin.sub_tmp)
   sample.size_tmp=beta_wt[n]*nrows
   sample.ind_tmp=sample(bin.sub_tmp[,1],sample.size_tmp+1,replace=TRUE)
   beta.dist_tmp=array(NA, dim=c(length(sample.ind_tmp),6))

  for(i in 1:length(sample.ind_tmp))
    {
#     beta.dist_tmp[i,1]=array[sample.ind_tmp[i],1]
     beta.dist_tmp[i,1]=array[sample.ind_tmp[i],2]
     beta.dist_tmp[i,2]=array[sample.ind_tmp[i],3]
     beta.dist_tmp[i,3]=array[sample.ind_tmp[i],4]
     beta.dist_tmp[i,4]=array[sample.ind_tmp[i],5]
     beta.dist_tmp[i,5]=array[sample.ind_tmp[i],6]
     beta.dist_tmp[i,6]=array[sample.ind_tmp[i],7]
    }

   beta.dist=rbind(beta.dist,beta.dist_tmp)  #combines bins into single array
  }

print(beta.dist)

### plot pdf ###
d.array=hist(array[,7], breaks=xt,plot=FALSE)
d.beta.dist=hist(beta.dist[,6],breaks=xt,plot=FALSE) # use histograms for easier weighting

pdf(file="./output/fitted_beta.pdf",height=10, width=6) 
par(mfcol=c(2,1))

# panel A: plot the beta distribution fit to Pfeffer

plot(xt,beta_wt,type="l",col="black",
     lwd=4,xlim=c(-50,3000),
     ylim=c(0,0.05),
     xlab="Projected Sea-Level Rise in 2100 [mm]",
     ylab="Probability density function", main="(a) Extended Scenario of Pfeffer et al (2008)")
# add the zero line
 lines(xt,beta_pdf*0,lwd=1,lty=2,col="gray")
# add an arrow for the range of the distribution
 y.loc.arrow.a=9.0e-4;
#arrows(lb.2100,y.loc.arrow.a,ub.2100,y.loc.arrow.a,col="black",lty="solid",code=3);
 abline(v=lb.2100,lwd=1,col="blue",lty=2)
 abline(v=ub.2100,lwd=1,col="blue",lty=2)

#lines(d.array$mids,d.array$counts/length(array[,7]),lwd=3,lty=2,col="blue")

# add normalized fitted beta distribution
  lines(d.beta.dist$mids,d.beta.dist$counts/length(beta.dist[,6]),lwd=3,lty=2,col="red")

# check that normalized distributions sum to one
  print(sum(d.array$counts)/length(array[,7]))
  print(sum(d.beta.dist$counts)/length(beta.dist[,6]))
  print(sum(beta_wt))



legend(1600,0.045, inset=.2,
c("bounds","Pdf fit","Model Approximation"),
   lwd=2, col=c("blue","black","red"), cex=.7)

# panel B: plot the uniform dist fit to the CA scenarios
plot(x.sample.CA,uniform.fit.CA,type="l",col="black", lwd=4,
     xlim=c(-50,3000),
     ylim=c(0,6e-4),
     xlab="Projected Sea-Level Rise in 2100 [mm]",
     ylab="Probability density function", main="(b) Extended Scenario of Co-CAT (2010)")
# add the zero line
lines(x.sample.CA,uniform.fit.CA*0,lwd=1,lty=2,col="gray")
#arrows(lb.2100,y.loc.arrow.a,ub.2100,y.loc.arrow.a,col="black",lty="solid",code=3);
abline(v=lb.CA,lwd=1,col="blue",lty=2)
abline(v=ub.CA,lwd=1,col="blue",lty=2)
legend(2100,5.5e-4, inset=.2,
c("bounds","Pdf fit"),
   lwd=2, col=c("blue","black"), cex=.7)


dev.off()   #write file




# Nathan's resources for marginals and pair plots
source("../plotutils.R")

figure(pdf(file="./output/marginals.pdf"))
plot.marginals(beta.dist) 
dev.off()

figure(pdf(file="./output/pairs.pdf"))
plot.pairs(beta.dist)
dev.off()

write.table(beta.dist,file="./output/array_beta.txt",quote=FALSE,col.names=FALSE)


# zoom into the joint distribution of t.star and c.star 
figure(pdf(file="./output/c*vst*.pdf"))
plot(beta.dist[,4],beta.dist[,5],pch=16,cex=0.2,col="black",
     xlab="t* estimates [Years after 2011]",
     ylab="c* estimate [mm/a]", main=" ")
dev.off()
