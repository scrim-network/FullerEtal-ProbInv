# data.R
# written by Robert W. Fuller on 090321

source("months.R")
source("ts.R")


subLoadAnnual <- function(gm, ar4bias, ...)
{
    # get gm in average form and fix up month fraction in monthly form

    gmAvg <- matrix(nrow=(nrow(gm) / 12), ncol=ncol(gm))
    colnames(gmAvg) <- colnames(gm) <- c("time", "temp")

    row <- 1
    for (i in 1:nrow(gmAvg)) {

        year <- floor(gm[row, "time"])
        endRow <- row + 11
        if (year != floor(gm[endRow, "time"])) {
            stop("error in calculating average")
        }

        temp <- sum(gm[ row:endRow, "temp" ]) / 12

        # TODO:  this is very much like subLoadSeaLevel();
        # differences:  error column, column names, mm unit conversion;
        # these could be parameterized thereby replacing both functions
        # with one function;  mm unit conversion could be parameterized
        # in the same manner as convFn in subLoadTemp()
        #
        #error  <- sum(gm[ row:endRow, "error"  ]) / 12
        gmAvg[i, ] <- c(year, temp)

        # fix up month fraction to be consistent with others
        gm[ row:endRow, "time" ] <- year + months

        row <- row + 12
    }

    if (ar4bias) {
        # specify cols="temp" to avoid applying bias to observation error
        gm    <- tsBias(gm,    12, ..., cols="temp")
        gmAvg <- tsBias(gmAvg,  1, ..., cols="temp")
    }

    return (list( annual=gmAvg, monthly=gm ))
}


loadIpccTemp <- function(file="../data/ipcc/gmst/mri_cgcm2_3_2a.txt", ar4bias=T, ...)
{
    error <- F
    gmstRaw <- tryCatch(scan(file, what=numeric(), quiet=T), error=function(e) { error <<- T; return (e) })
    if (error) {
        print(gmstRaw)
        stop(paste("could not scan file", file))
    }

    gmst <- matrix(gmstRaw, ncol=2, byrow=T)

    return (subLoadAnnual(gm=gmst, ar4bias=ar4bias, ...))
}


loadIpccTempDriftCorrectNathan <- function(file)
{
    #temp  <- loadIpccTemp(paste("../data/20c3m/run1/",   file, sep=""), ar4bias=F)$annual

    temp  <- loadIpccTemp(paste("../data/sresa1b/run1/", file, sep=""), ar4bias=F)$annual
    ctl   <- loadIpccTemp(paste("../data/picntrl/run1/", file, sep=""), ar4bias=F)$annual

    #print(ctl)
    ctl   <- tsTrim(ctl, 2001, 2099)

    drift <- tsDriftCorrect(temp, ctl, truncate=T, method="linear", f=0.75)

    #frc   <- tsBias(drift, 1)
    #print(drift)

    return (drift)
}
