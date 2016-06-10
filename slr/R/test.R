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
# test.R

source('assim.R')

# truth can be "mean", "max", or "mode"
getEqTemp <- function(seq_lvl, assimctx=grinassimctx, truth="mode")
{
    mp <- truth(assimctx, truth, F)$mp
    temp <- (seq_lvl - mp["b"]) / mp["a"]

    print(temp)
}
