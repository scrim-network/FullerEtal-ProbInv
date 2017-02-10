HOWTO: Running the code for the paper with the working title "Combining Expert Assessments with Paleo- and Instrumental Observations to Produce Probabilistic Projections about the Antarctic Ice Sheet"

Copyright (C) 2016, 2017 Robert William Fuller <hydrologiccycle@gmail.com>. The intent is to release this code under the GPLv3 after the corresponding article is accepted by a journal. Until the author has publicly released this code, please treat this code as confidential. Also included is code that is Copyright (C) 2016 by Kelsey Ruckert, Yawen Guan, and Tony Wong. They are greatfully acknowledge for their contributions. For details see the file "calib.R" in the source code.

This code has been tested on MacOS El Capitan, Fedora 24 Workstation, and Redhat Enterprise Linux 6. The following versions of R have been used with the code: 3.3.2, 3.3.0, 3.2.1, and 3.1.3. You will need to install the following R packages: adaptMCMC, DEoptim, fields, KernSmooth, lhs, mcmc, and RColorBrewer. Detailed instructions follow.

It will be much simpler to run the assimilation on a host with the Public Batch System (PBS) installed, as well as ample memory to accomodate multiple simultaneous jobs. These instructions are written with these assumptions in mind. If you are a member of SCRiM, napa is the recommended host. This code has also been run on PSU's ACI-b batch processing cluster.


Install and build the source:

cmd:  git clone git@github.com:scrim-network/FullerEtal-ProbInv.git ~/pi
cmd:  cd ~/pi/R
cmd:  R
from R:  install.packages(c('adaptMCMC', 'DEoptim', 'fields', 'KernSmooth', 'lhs', 'mcmc', 'RColorBrewer'))
from R:  source('calib.R')

If there are any errors, contact the author before proceeding.  Exit R using ctrl-D.


Run the assimilation:

cmd:  cd ~/pi/pbs
cmd:  ./run
cmd:  ./run_diagnose
cmd:  qstat (to confirm that the jobs are running)

Console output from the jobs is saved in ~/pi/out. These may be used to monitor the progress of MCMC. Wait for all 7 jobs to complete. There will be 7 ".RData" save files in ~/pi/R.


Make the figures and table:

cmd:  cd ~/pi/R
cmd:  Rscript results.R

The figures are in ~/pi/figures. The table is in ~/pi/out/quantiles.csv.


Please direct questions to Robert W. Fuller <hydrologiccycle@gmail.com>.
