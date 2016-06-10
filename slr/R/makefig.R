# Copyright 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# makefig.R


# RAND spaghetti

load("../runs/paper/pred_tau900_spagh")
source('newfig.R')
figRandSpag(outfiles=T)

rm("alleyctx", "allgrgisctx", "grctx", "prallgrgisctx", "savedir", "tau")


# Uber figure and RAND plots

load("../runs/paper/pred_tau300")
pr1 <- prallgrgisctx
load("../runs/paper/pred_tau900")
pr2 <- prallgrgisctx
load("../runs/paper/pred_tau2700")
pr3 <- prallgrgisctx

#source('newfig.R')
source('klausfig.R')

figUber(outfiles=T)
figRand(outfiles=T)
figPredict(outfiles=T)
figHindcast(outfiles=T)

rm("alleyctx", "allgrgisctx", "grctx", "pr1", "pr2", "pr3", "prallgrgisctx", "savedir", "tau")


# identifiability problem in Tau and predictions

chainload("../runs/paper/ar1/prperf", oldnames=c("grinassimctx", "prgrinctx"), newnames=c("assim", "pr"))

source('newfig.R')

figIdent(outfiles=T)

rm( "assim1", "assim2", "assim3", "assim4", "assim5",
    "pr1", "pr2", "pr3", "pr4", "pr5",
    "grctx", "perf",
    "savedir", "subdir")


# influence of mass balance and paleo priors on uncertainty and predictions

load("../runs/paper/pred_tau900_no_gis")
rename("prallgrgisctx", "prbase")

load("../runs/paper/pred_tau900_paleo_only")
rename("prallgrgisctx", "prpaleo")

load("../runs/paper/pred_tau900_no_paleo")
rename("prallgrgisctx", "prrignot")

load("../runs/paper/pred_tau900")
rename("prallgrgisctx", "prall")

source('newfig.R')

newDev("gis_2100", outfile=T, width=6, height=6)

chains <- list(prall$prchain, prbase$prchain, prpaleo$prchain, prrignot$prchain)
cictx <- ciCalc(chains=chains, xvals=2100, probs=c(0.005, 0.995))
pdfPlots(
    chains=chains,
    column=as.character(2100),
    lty=c("solid", "dashed", "dotted", "dotdash"),
    legends=c("All data", "Tide gage", "Tide gage & paleo", "Tide gage & TMB"),
    col=c("black", "blue", "red", "green"),
    burnin=F,
    xlim=c(0, max(cictx$range)),
    lwd=2
    )

rm("alleyctx", "allgrgisctx", "grctx", "prall", "prbase", "prpaleo", "prrignot", "savedir", "tau", "chains", "cictx")


# fit to IPCC GCMs

source('ipcc.R')

mri <- env()
grinConfigAssim(assimctx=mri, ipccFile="mri_cgcm2_3_2a.txt", ipcc=T, ar=0)

miub <- env()
grinConfigAssim(assimctx=miub, ipccFile="miub_echo_g.txt", ipcc=T, ar=0)

inm <- env()
grinConfigAssim(assimctx=inm, ipccFile="inmcm3_0.txt", ipcc=T, ar=0)

ccc <- env()
grinConfigAssim(assimctx=ccc, ipccFile="cccma_cgcm3_1.txt", ipcc=T, ar=0)

aom <- env()
grinConfigAssim(assimctx=aom, ipccFile="giss_aom.txt", ipcc=T, ar=0)

ipccPlotFits(assimlist=list(mri, miub, inm, ccc, aom), outfiles=T)

rm("aom", "ccc", "grctx", "grinassimctx", "inm", "miub", "mri", "prgrinctx")


# pairs plot

load("../runs/paper/pred_tau900")
source("plotutils.R")
newDev("pairs", outfile=T, height=8.5, width=8.5)
plot.pairs(allgrgisctx$chain)

rm("alleyctx", "allgrgisctx", "grctx", "prallgrgisctx", "savedir", "tau")


# fit to paleo-data

load("../runs/paper/pred_tau900")
source("allgrgis.R")
allgrgisRunSeqPredict(startTemp=-12)

newDev("seqpred", outfile=T, width=6, height=6)
alleyPlotSeqGis(outfiles=F, newdev=F)

rm("alleyctx", "allgrgisctx", "grctx", "prallgrgisctx", "savedir", "tau")


# latest Klausian figure

load("../runs/paper/pred_tau900")
source("klausfig.R")
figProbForce(outfiles=T)
rm("alleyctx", "allgrgisctx", "grctx", "prallgrgisctx", "savedir", "tau")


finDev()
