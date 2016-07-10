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
# restd.R

source("roblib.R")


ci_to_q <- function(ci)
{
    q <- (1.0 - ci) / 2

    return (c(q, 1.0 - q))
}


get_sd <- function(range, ci)
{
    # p is equivalent to ci_to_q(ci)[2]
    p <- 0.5 + ci / 2
    sd <- (range[2] - range[1]) / 2 / qnorm(p)

    return (sd)
}


standardize <- function(oldRange, oldci, newci)
{
    sd <- get_sd(oldRange, oldci)
    mean <- mean(oldRange)

    q <- ci_to_q(newci)
    newRange <- qnorm(q, mean=mean, sd=sd)

    return (newRange)
}
