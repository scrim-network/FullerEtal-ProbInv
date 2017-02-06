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
# results.R

# do the LHS run
source('calib.R')
daisRunLhs()
source('figures.R')
figLhs()

# load a bunch of runs
source('makefig.R')
load('dp="u";n=5e6.RData')

# make the fast dynamics figure
source('figures.R')
figDiagFast()

# make the other figures
source('makefig.R')
figAisPriors()
figCmpPriors()
figPredict()
figInfer()
figCmpPredict()
