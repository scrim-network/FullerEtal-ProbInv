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

savedir <- "../runs/paper/"
dir.create(savedir)

source('allgrgis.R')
allgrgisConfigAssim(fp=c(tau=900, gis_scl=1.7, gis_s0=0), rignot_ts=F)
allgrgisRunAssim(nbatch=100000)
allgrgisRunAssim(nbatch=100000)
allgrgisRunAssim()
save.image(paste(savedir, "assim_tau900_paleo_only", sep=""))
allgrgisRunPredict()
save.image(paste(savedir, "pred_tau900_paleo_only", sep=""))
