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
# makeascii.R


save_year <- function(assimctx=prallgrgisctx, chain=assimctx$prchain, year="2100", basename=NULL)
{
    #xvals <- attr(prallgrgisctx$prchain, "xvals")
    preds <- chain[, year]
    filename <- paste(sep="", "../runs/paper/", basename, "_", year, ".txt")
    write.table(preds, file=filename, row.names=F, col.names=F)
}


save_grin <- function(tau, year="2100")
{
    basename <- paste(sep="", "tau", tau)
    load(paste(sep="", "../runs/paper/pred_", basename))
    save_year(assimctx=prallgrgisctx, basename=basename, year=year)
}


save_grin(300)
save_grin(900)
save_grin(2700)
