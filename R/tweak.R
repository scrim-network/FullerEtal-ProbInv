# Copyright (C) 2017 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# tweak.R
#

iter   <- "5e6"
fnames <- c("uniform", "beta", "normal")
files  <- c(paste('ep="', substr(fnames,    1, 1), '";n=', iter, ".RData", sep=""),
            paste('ip="', substr(fnames,    1, 1), '";n=', iter, ".RData", sep=""),
            paste('dp="', substr(fnames[1], 1, 1), '";n=', iter, ".RData", sep=""))

for (file in files) {
    load(file)
    print(paste("processing", file))

source('calib.R')
daisRunHindcast()

    save.image(file)
    rm(  daisctx,
       prdaisctx)
}
