set.seed(9);
source("daisassim.R")
daisConfigAssim()
daisRunFit()
daisRunAssim(nbatch=1e4)
#save.image("DAIS_MCMC_Rversioncalibration.RData")
