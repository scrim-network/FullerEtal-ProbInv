set.seed(9)
source("calib.R")
daisConfigAssim(heaviside=T)  # all_predict=T)
daisRunFit()
daisRunAssim(nbatch=1e4)
daisRunPredict()
save.image("test.RData")
