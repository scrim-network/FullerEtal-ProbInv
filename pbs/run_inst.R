cmd    <- paste(commandArgs(), collapse=" ")
prior  <-            sub(".*inst_(.*)_.*",       "\\1", cmd)
nbatch <- as.numeric(sub(".*inst_.*_(.*)\\.R.*", "\\1", cmd))

source("calib.R")
daisConfigAssim(prior=prior, instrumental=T, paleo=T)
daisRunFit()
daisRunAssim(nbatch=nbatch)
npredict <- ifelse(nbatch > 5e5, 1e5, 1e4)
npredict <- nbatch / 5
daisRunPredict(nbatch=npredict)
save.image(paste("inst_", prior, "_", nbatch, ".RData", sep=""))
