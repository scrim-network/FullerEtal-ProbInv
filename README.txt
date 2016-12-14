HOWTO:  Running the code for "Probabilistic Inversion of an Antarctic Ice Sheet Model"

Copyright 2016, Robert William Fuller <hydrologiccycle@gmail.com>.  The intent is to release this code under the GPLv3 after the corresponding article is accepted by a journal.  Until the author has publicly released this code, please treat this code as confidential.  Also included is code that is Copyright 2016 by Kelsey Ruckert, Yawen Guan, and Tony Wong.  They are greatfully acknowledge for their contributions.  For details see the file "calib.R" in the source code.

This code has been tested on MacOS El Capitan, Fedora 24 Workstation, and Redhat Enterprise Linux 6.  The following versions of R have been used with the code:  3.3.0, 3.2.1, and 3.1.3.  You will need to install the following R packages:  adaptMCMC, DEoptim, fields, KernSmooth, lhs, mcmc, mvtnorm, and RColorBrewer.  Detailed instructions follow.

It will be much simpler to run the assimilation on a host with the Public Batch System (PBS) installed, as well as ample memory to accomodate multiple simultaneous jobs.  These instructions are written with these assumptions in mind.  If you are a member of SCRiM, napa is the recommended host.

INSTALL the source:

cmd:  ssh napa
cmd:  cd ~
cmd:  git clone /home/scrim/rwf136/pi ~/pi
cmd:  git clone git@bitbucket.org:aub48/robustslr.git ~/bitbucket
cmd:  cd ~/bitbucket
cmd:  git checkout FastDynamics
cmd:  cd ~/pi/R
cmd:  R
from R:  install.packages(c('adaptMCMC', 'DEoptim', 'fields', 'KernSmooth', 'lhs', 'mcmc', 'mvtnorm', 'RColorBrewer'))
from R:  source('calib.R')
If there are any errors, contact the author before proceeding.
Exit R using ctrl-D.

RUN the assimilation:

cmd:  cd ~/pi/pbs
cmd:  ./run
cmd:  ./run_diagnose2
Monitor the PBS jobs with the command:  qstat
Console output from the jobs is saved in ~/pi/out.  These may be used to monitor the progress of MCMC.
Wait for all 7 jobs to complete.  There will be 7 ".RData" save files in ~/pi/R.
cmd:  cd ~/pi/R
cmd:  R
from R:  source('makefig.R'); figAisPriors(); figCmpPriors(); figPredict(); figUber(); figCmpPredict()
from R:  source('calib.R'); daisRunLhs(); source('figures.R'); figLhs()
from R:  load('d2p="u";n=5e6.RData'); source('figures.R'); figDiagFast()
The figures are in ~/pi/figures.

Please direct questions to Robert W. Fuller, hydrologiccycle@gmail.com.