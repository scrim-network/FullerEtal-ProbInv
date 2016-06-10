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
# marlos.R
# written by Robert W. Fuller on 090403

source("data.R")

refYear <- 1870

jev    <-    loadJev(units_mm=T, lower=refYear, upper=refYear)$annual
church <- loadChurch(units_mm=T, lower=refYear, upper=refYear)$annual

eol="\r\n"
sep=", "
col.names=c("year", "gmsl (mm)")

write.table(jev[, 1:2],    file="jev.txt",    row.names=FALSE, col.names=col.names, eol=eol, sep=sep)
write.table(church[, 1:2], file="church.txt", row.names=FALSE, col.names=col.names, eol=eol, sep=sep)
