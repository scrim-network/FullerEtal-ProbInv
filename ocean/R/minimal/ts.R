# ts.R
# written by Robert W. Fuller on 090622


tsTrim <- function(ts, startYear, endYear, col="time")
{
    # close upper boundary
    endYear <- endYear + 1

    return (ts[ (ts[, col] >= startYear & ts[, col] < endYear), ])
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
