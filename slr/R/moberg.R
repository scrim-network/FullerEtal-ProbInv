# Copyright 2009, 2010 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# moberg.R

source('data.R')
source('ts.R')

ar4bias <- T

moberg <- loadMoberg( ar4bias=ar4bias)
crut2  <- loadHadcrut(ar4bias=ar4bias, dataset="nh")$annual
crut3  <- loadHadcrut(ar4bias=ar4bias, file="../data/brohan/data/www.cru.uea.ac.uk/cru/data/temperature/hadcrut3nh.txt")$annual

moberg  <- tsTrim(moberg,   1856, 1979)
crut2   <- tsTrim(crut2,    1856, 1979)
crut3   <- tsTrim(crut3,    1856, 1979)

print(colMeans(moberg))
print(colMeans(crut2))

# hadcrut 3 is not compatible
print(colMeans(crut3))
