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
# ts.R
# written by Robert W. Fuller on 090622


tsGetIndices <- function(ts, lower, upper, factor)
{
    if (missing(lower)) {
        lower <- ts[ 1, "time" ]
    }
    if (missing(upper)) {
        upper <- ts[ nrow(ts), "time" ]
    }
    if (missing(factor)) {
        factor <- 1
    }

    # close the upper boundary of the interval to simplify computations
    upper <- upper + 1

    indices <- which( (ts[, "time"] >= lower) & (ts[, "time"] < upper) )

    len <- length(indices)
    if (len != ((upper - lower) * factor)) {
        print(len)
        stop("did not find correct number of values in interval")
    }

    return (indices)
}


# ar4 reference period is 1980-1999, i.e. [1980, 2000)
tsBias <- function(ts, factor, lower=1980, upper=1999, cols=2:ncol(ts))
{
    indices <- tsGetIndices(ts, factor, lower, upper)
    len <- length(indices)

    for (col in cols) {
        bias <- sum(ts[indices, col]) / len
        ts[, col] <- ts[, col] - bias
    }

    return (ts) 
}


# could use this for tsBias(), but we would lose locality of reference
# for matrices that don't fit in the processor cache such as that
# returned by loadUrban()
#
tsGetBias <- function(ts, factor, lower=1980, upper=1999, cols=2:ncol(ts))
{
    indices <- tsGetIndices(ts, factor, lower, upper)
    len <- length(indices)

    #bias <- matrix(nrow=1, ncol=max(cols))
    bias <- numeric(length=max(cols))
    for (col in cols) {
        #bias[1, col] <- sum(ts[indices, col]) / len
        bias[col] <- sum(ts[indices, col]) / len
    }

    return (bias)
}


tsFindByDate <- function(ts, date, col=2, exact=T)
{
    index <- findInterval(date, ts[, "time"], all.inside=TRUE)
    if (length(index) != 1) {
        stop("found more than one matching date")
    }

    if (exact) {
        assert(ts[index, "time"] == date)
    } else {
        # operator "<=" ensures 0.5 rounds up
        if (abs(ts[index + 1, "time"] - date) <= abs(ts[index, "time"] - date)) {
            index <- index + 1
        }
    }

    return (ts[index, col])
}


tsTrim <- function(ts, startYear, endYear, col="time")
{
    # close upper boundary
    endYear <- endYear + 1

    return (ts[ (ts[, col] >= startYear & ts[, col] < endYear), ])
}


tsTrimForcing <- function(forcing, ts)
{
    indices <- ts[, "time"] - forcing[1, "time"] + 1
    return (forcing[indices, ])
}


# method can be "lowess", "linear", "simple"
tsDriftCorrect <- function(ts, control, col=2, truncate=F, method="lowess", f=0.25)
{
    switch (method,
        lowess={
            smoothed <- lowess(x=control[, "time"], y=control[, col], f=f)
            control[, col] <- smoothed$y
        },
        linear={
            x <- control[, "time"]
            y <- control[, col]
            fit <- lm(y ~ x)
            #plot(control, type="l", col="blue")
            control[, col] <- predict(fit)
            #lines(control, col="red")
            #stop("foo")
        },
        simple={
            # simple subtraction
        }, {
          stop("unknown method in tsDriftCorrect()")  
        })

    common <- intersect(ts[, "time"], control[, "time"])
    if (length(common) != nrow(ts)) {
        if (truncate) {
            ind <- ts[, "time"] %in% control[, "time"]
            ts  <- ts[ind, ]
        } else {
            stop("control does not cover range of time series")
        }
    }

    ind <- control[, "time"] %in% ts[, "time"]

    #dev.new()
    #plot(control[ind,])
    #print(control[ind,])

    ts[, col] <- ts[, col] - control[ind, col]

    return (ts)
}
