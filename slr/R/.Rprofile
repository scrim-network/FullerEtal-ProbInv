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


formLibPath <- function(path)
{
    # archdir is of the form x86_64-pc-linux-gnu-library
    archdir <- paste(R.version$platform, "-library", sep="")
    fqp     <- paste(path, "/", archdir, sep="")

    return (fqp)
}


addLibPath <- function(...)
{
    args <- list(...)
    for (path in args) {
        fqp <- formLibPath(path)
        if (FALSE == file.exists(fqp)) {
            dir.create(fqp, recursive=TRUE)
        }
        .libPaths(c(fqp, .libPaths()))
    }
}


.First <- function()
{
    .robSeed <<- as.integer(Sys.time())
    set.seed(.robSeed)

    options(help_type="html")
}


#source("~/.Rprofile")
.Second <- .First

.First <- function()
{
    .Second()

    # this doesn't work too well -- breaks help
    #addLibPath(".")

    # but this is fine
    addLibPath(getwd())

    .robPath <<- .libPaths()[1]
}
