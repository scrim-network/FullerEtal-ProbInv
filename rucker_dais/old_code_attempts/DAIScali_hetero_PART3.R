########################## Analysis of the MCMC chain produced  #################
rm(list =ls()) #Clear global environment
library(compiler)
enableJIT(3)
enableJIT(3)

#Set the seed
set.seed(1234)

load("Workspace/DAIS_MCMC_calibrationOrig.RData") # Load in the saved workspace from MCMC calibration

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
# new.dais.mcmc.est = iceflux(mean.dais.par, hindcast.forcings, standards)

#Create New best projection
end = 240298
enddate = 240300
new.dais.mcmc.proj = iceflux(mean.dais.par, project.forcings, standards)

#save.image(file = "Workspace/DAIS_MCMC_calibration.RData")
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

################### PROJECT AIS CONTRIBUTIONS ##########################
#Set the end dates to the projection period:
end = 240298
enddate = 240300

proj.dais.fits = mat.or.vec(subset_length, enddate)
for(i in 1:subset_length) {
  proj.dais.fits[i,] = iceflux(par.mcmc[i,], project.forcings, standards)
}
save.image(file = "Workspace/DAIS_MCMC_caliPART3.RData")

# ### Superimpose the bias onto the model
### True world = model + bias + error
bias.mcmc = sschain[,12]

# proj.mcmc.w.bias = mat.or.vec(subset_length, enddate)
# for(i in 1:subset_length){
#   proj.mcmc.w.bias[i,] = proj.dais.fits[i,] + bias.mcmc[i]^2
# }

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
#save.image(file = "Workspace/DAIS_MCMC_caliPART3.RData")
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

#Lets find out how the parameter relationships
#But lets also add in the initial value to the mix
#set up a matrix
d.pos_parameters = sschain
colnames(d.pos_parameters, do.NULL = FALSE)
colnames(d.pos_parameters) = c("gamma", "alpha", "mu", "eta", "po", "kappa", "fo", "ho", "co","bo", "s", "sigma.y")

save.image(file = "Workspace/DAIS_MCMC_caliPART3.RData")

######################## For Further Analysis Plot graphes ######################
#source("heteroskedastic_dais_plots.R")
