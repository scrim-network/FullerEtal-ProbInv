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

source("grinassim.R")

# set seed to ensure all 5 perfect model experiments have same truth values
#set.seed(42)
set.seed(12)

grinConfigAssim()
grinRunAssim(nbatch=100000)
grinRunAssim()

perf <- as.numeric(sub(".*perf([[:digit:]]+).*", "\\1", paste(commandArgs(), collapse=" ")))
savedir <- "../runs/paper/"
subdir  <- "ar1/"
dir.create(paste(savedir, subdir, sep=""))

# need different noise realization for each perfect model experiment
set.seed(perf)

grinConfigPerfectModel(uniform=T, useUrban=T)
grinRunAssim(nbatch=100000, initial=T)
grinRunAssim()
#save.image(paste(savedir, subdir, "asperf", perf, sep=""))
grinRunPredict()
save.image(paste(savedir, subdir, "prperf", perf, sep=""))
