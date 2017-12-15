HOWTO: Running the code for the paper titled "Probabilistic inversion of expert assessments to inform projections about Antarctic ice sheet responses"

Copyright (C) 2016, 2017 Robert William Fuller <hydrologiccycle@gmail.com>. Also included is code that is Copyright (C) 2016 by Kelsey Ruckert, Yawen Guan, and Tony Wong. They are gratefully acknowledged for their contributions. For details see the file "calib.R" in the source code.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program.  If not, see <http://www.gnu.org/licenses/>.

This code has been tested on MacOS El Capitan, Fedora 24 Workstation, and Redhat Enterprise Linux 6. The following versions of R have been used with the code: 3.3.2, 3.3.0, 3.2.1, and 3.1.3. You will need to install the following R packages: adaptMCMC, DEoptim, entropy, fields, KernSmooth, lhs, mcmc, MCMCpack, and RColorBrewer. Detailed instructions follow.

It will be much simpler to run the assimilation on a host with the Public Batch System (PBS) installed, as well as ample memory to accomodate multiple simultaneous jobs. These instructions are written with these assumptions in mind. If you are a member of SCRiM, napa is the recommended host. This code has also been run on PSU's ACI-b batch processing cluster.


Install and build the source:

cmd:  git clone git@github.com:scrim-network/FullerEtal-ProbInv.git ~/pi
cmd:  cd ~/pi/R
cmd:  R
from R:  install.packages(c('adaptMCMC', 'DEoptim', 'entropy', 'fields', 'KernSmooth', 'lhs', 'mcmc', 'MCMCpack', 'RColorBrewer'))
from R:  source('calib.R')

If there are any errors, contact the author before proceeding.  Exit R using ctrl-D.


Run the assimilation:

cmd:  cd ~/pi/pbs
cmd:  ./run
cmd:  ./run_diagnose
cmd:  qstat (to confirm that the jobs are running)

Console output from the jobs is saved in ~/pi/out. These may be used to monitor the progress of MCMC. Wait for all 8 jobs to complete. There will be 8 ".RData" save files in ~/pi/R.


Make the figures and table:

cmd:  cd ~/pi/R
cmd:  Rscript results.R

The figures are in ~/pi/figures. The table data is in the files ~/pi/out/quantiles.csv and ~/pi/out/quantiles_tony.csv.


Please direct questions to Robert W. Fuller <hydrologiccycle@gmail.com>.
