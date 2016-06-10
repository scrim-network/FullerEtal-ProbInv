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
# data.R
# written by Robert W. Fuller on 090321

source("months.R")
source("ts.R")


subLoadTemp <- function(gmstRaw, convFn, ar4bias, ..., calcAvg=(ncol(gmstRaw) < 14))
{
    # get gmst data in same form as gmsl

    rows <- nrow(gmstRaw)

    gmst    <- matrix(nrow=rows * 12, ncol=2)
    gmstAvg <- matrix(nrow=rows,      ncol=ncol(gmst))

    row <- 1
    for (i in 1:rows) {
        year  <- gmstRaw[ i, 1 ]
        temps <- convFn(gmstRaw[ i, 2:13 ])

        gmst[ row:(row+11), ] <- cbind(year + months, temps)

        if (calcAvg) {
            gmstAvg[ i, ] <- c(year, sum(temps) / 12)
        } else {
            gmstAvg[ i, ] <- c(year, gmstRaw[ i, 14 ])
        }

        row <- row + 12
    }

    colnames(gmstAvg) <- colnames(gmst) <- c("time", "temp")

    if (ar4bias) {
        gmst    <- tsBias(gmst,    12, ...)
        gmstAvg <- tsBias(gmstAvg,  1, ...)
    }

    return (list( annual=gmstAvg, monthly=gmst ))
}


loadHadcrut <- function(ar4bias=T, dataset="variance", file=NULL, ...)
{
    if (is.null(file)) {
        file <- switch(dataset,
            global="../data/brohan/data/www.cru.uea.ac.uk/cru/data/temperature/hadcrut3gl.txt",
            variance="../data/brohan/data/www.cru.uea.ac.uk/cru/data/temperature/hadcrut3vgl.txt",
            nh="../data/jones/tavenh2v.dat"
            )
    }
    gmstRaw <- read.table(file, fill=TRUE)
    gmst <- as.matrix(gmstRaw)
    gmst <- gmst[ seq(1, nrow(gmstRaw), by=2), ]

    return (subLoadTemp(gmst=gmst, convFn=identity, ar4bias=ar4bias, ...))
}


loadVinther <- function(ar4bias=T, ...)
{
    gmstRaw <- scan("../data/vinther/www.cru.uea.ac.uk/cru/data/greenland/swgreenlandave.dat",
                    skip=14, what=numeric(), quiet=T)
    gmst <- matrix(gmstRaw, ncol=13, byrow=T)

    return (subLoadTemp(gmst=gmst, convFn=function(x) { x / 10 },
                        ar4bias=ar4bias, ...))
}


loadRignot <- function(pure=T)
{
    file <- ifelse(pure, "../data/rignot/puredata.txt", "../data/rignot/data.txt")
    gisSlr <- scan(file, what=numeric(), quiet=T, sep=",")
    gisSlr <- matrix(gisSlr, ncol=3, byrow=T)

    # 360 Gton per mm SLR, aka SLE (Sea Level Equivalent)
    gisSlr <- cbind(gisSlr[, 1], gisSlr[, 2:3] / -360 / 1000, gisSlr[, 2:3])
    colnames(gisSlr) <- c("time", "SLE", "error", "mass", "mass error")
    gisSlr[, "error"] <- abs(gisSlr[, "error"])

    return (gisSlr)
}


subLoadPrior <- function(basename, temp_bias)
{
    filename <- paste("../data/", basename, "/data.txt", sep="")
    fig <- scan(filename, what=numeric(), quiet=T, sep=",", skip=2)
    fig <- matrix(fig, ncol=6, byrow=T)
    colnames(fig) <- c("temp", "SLE", "temp_lbound", "temp_ubound", "SLE_lbound", "SLE_ubound")

    cols <- c("temp", "temp_lbound", "temp_ubound")
    fig[, cols] <- fig[, cols] + temp_bias

    return (fig)
}


loadAlley <- function(temp_bias=0)
{
    return (subLoadPrior("alley", temp_bias))
}


loadArcher <- function(temp_bias=0)
{
    return (subLoadPrior("archer", temp_bias))
}


subLoadSeaLevel <- function(gmsl, units_mm, ar4bias, ...)
{
    # get gmsl in average form and fix up month fraction in monthly form

    gmslAvg <- matrix(nrow=(nrow(gmsl) / 12), ncol=ncol(gmsl))
    colnames(gmslAvg) <- colnames(gmsl) <- c("time", "sealvl", "error")

    row <- 1
    for (i in 1:nrow(gmslAvg)) {

        year <- floor(gmsl[row, "time"])
        endRow <- row + 11
        if (year != floor(gmsl[endRow, "time"])) {
            stop("error in calculating average sea level")
        }

        sealvl <- sum(gmsl[ row:endRow, "sealvl" ]) / 12
        error  <- sum(gmsl[ row:endRow, "error"  ]) / 12
        gmslAvg[i, ] <- c(year, sealvl, error)

        # fix up month fraction to be consistent with others
        gmsl[ row:endRow, "time" ] <- year + months

        row <- row + 12
    }

    # convert to meters
    if (!units_mm) {
        cols <- c("sealvl", "error")
        gmsl   [, cols] <- gmsl   [, cols] / 1000
        gmslAvg[, cols] <- gmslAvg[, cols] / 1000
    }

    if (ar4bias) {
        # specify cols="sealvl" to avoid applying bias to observation error
        gmsl    <- tsBias(gmsl,    12, ..., cols="sealvl")
        gmslAvg <- tsBias(gmslAvg,  1, ..., cols="sealvl")
    }

    return (list( annual=gmslAvg, monthly=gmsl ))
}


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


loadJev <- function(startYear=1808, endYear=2001, units_mm=F, ar4bias=T, ...)
{
    gmslRaw <- scan("../data/jevrejeva/www.pol.ac.uk/psmsl/author_archive/jevrejeva_etal_gsl/gslJC2006.txt",
                    skip=17, what=numeric(), quiet=T)
    gmsl <- matrix(gmslRaw, ncol=5, byrow=T)

    # discard extra columns
    gmsl <- gmsl[, c(1, 4, 5)]

    # trim partial years
    gmsl <- tsTrim(gmsl, startYear, endYear, col=1)

    return (subLoadSeaLevel(gmsl=gmsl, units_mm=units_mm, ar4bias=ar4bias, ...))
}


loadJevAmsterdam <- function(units_mm=F, ar4bias=T, ...)
{
    gmslRaw <- scan("../data/jevrejeva/www.pol.ac.uk/psmsl/author_archive/jevrejeva_etal_1700/gslGRL2008.txt",
                    skip=14, what=numeric(), quiet=T)
    gmsl <- matrix(gmslRaw, ncol=3, byrow=T)
    colnames(gmsl) <- c("time", "sealvl", "error")

    if (!units_mm) {
        cols <- c("sealvl", "error")
        gmsl   [, cols] <- gmsl   [, cols] / 1000
    }

    if (ar4bias) {
        # specify cols="sealvl" to avoid applying bias to observation error
        gmsl <- tsBias(gmsl,  1, ..., cols="sealvl")
    }

    return (gmsl)
}


loadMoberg <- function(ar4bias=T, instrumental=F, ...)
{
    gmstRaw <- scan("../data/moberg/ftp.ncdc.noaa.gov/pub/data/paleo/contributions_by_author/moberg2005/nhtemp-moberg2005.txt", skip=93, what=numeric(), quiet=T)
    gmstRaw <- matrix(gmstRaw, ncol=9, byrow=T)
    gmst <- cbind(gmstRaw[, 1], gmstRaw[, 2])
    colnames(gmst) <- c("time", "temp")

    # supplement paleo data with instrumental record to which paleo data
    # was calibrated;  See Moberg et al. and Jones and Moberg for details
    #
    inst <- loadHadcrut(ar4bias=F, dataset="nh")$annual

    # period of overlap is 1856-1979
    if (instrumental) {
        gmst <- tsTrim(gmst, -Inf, 1855)
    } else {
        inst <- tsTrim(inst, 1980, Inf)
    }

    gmst <- rbind(gmst, inst)

    if (ar4bias) {
        gmst <- tsBias(gmst,  1, ...)
    }

    return (gmst)
}


loadChurch <- function(units_mm=F, ar4bias=T, ...)
{
    gmslRaw <- scan("../data/church and white/data/church_white_grl_gmsl.txt",
                    what=numeric(), quiet=T)
    gmsl <- matrix(gmslRaw, ncol=3, byrow=T)

    return (subLoadSeaLevel(gmsl=gmsl, units_mm=units_mm, ar4bias=ar4bias, ...))
}


subLoadIpccSlrAnnual <- function(gmsl, ar4bias, ...)
{
    colnames(gmsl) <- c("time", "sealvl", "error")
    gmsl[, "time"] <- floor(gmsl[, "time"])
    if (ar4bias) {
        gmsl <- tsBias(gmsl, 1, ...)
    }

    return (list( annual=gmsl, monthly=NULL ))
}


loadIpccSeaLevel <- function(file="../data/ipcc/hindcast/giss_aom.txt", ar4bias=T, ...)
{
    gmslRaw <- scan(file, what=numeric(), quiet=T)
    gmsl <- matrix(gmslRaw, ncol=2, byrow=T)

    # TODO:  fix this hack (error=0)
    gmsl <- cbind(gmsl, rep(0, nrow(gmsl)))

    if (gmsl[1, 1] - floor(gmsl[1, 1]) > 0.5) {
        return (subLoadIpccSlrAnnual(gmsl, ar4bias=ar4bias, ...))
    }

    return (subLoadSeaLevel(gmsl=gmsl, units_mm=T, ar4bias=ar4bias, ...))
}


loadDomingues <- function(file="../data/domingues/tsl_three_yearly.txt", units_mm=F, ar4bias=T, ...)
{
    gmslRaw <- scan(file, what=numeric(), quiet=T)
    gmsl <- matrix(gmslRaw, ncol=3, byrow=T)

    colnames(gmsl) <- c("time", "sealvl", "error")
    gmsl[, "time"] <- floor(gmsl[, "time"])

    if (!units_mm) {
        gmsl[, "sealvl"] <- gmsl[, "sealvl"] / 1000
        gmsl[, "error"]  <- gmsl[, "error"]  / 1000
    }

    if (ar4bias) {
        gmsl <- tsBias(gmsl, 1, ...)
    }

    return (gmsl)
}


loadUrban <- function(startYear=1850, endYear=2300, ar4bias=T, ...)
{
    # reading the file takes about 2 seconds
    raw <- scan("../data/urban/temps.dat", what=numeric(), quiet=T)
    cols <- endYear - startYear + 1
    dim(raw) <- c(cols, length(raw) / cols)

    gmst <- cbind(startYear:endYear, raw)
    colnames(gmst) <- c("time", rep("", ncol(raw)))
    #colnames(gmst)[1] <- "time"

    # eliminate columns with bad predictions
    bad <- which(F == is.finite(gmst), arr.ind=T)[, "col"]
    if (length(bad)) {
        gmst <- gmst[, -bad]
    }

    if (ar4bias) {
        gmst <- tsBias(gmst, 1, ...)
    }

    return (gmst)
}


loadUrbanAvg <- function(...)
{
    urban <- loadUrban(...)
    gmst <- cbind(urban[, "time"], rowMeans(urban[, 2:ncol(urban)]))
    colnames(gmst) <- c("time", "temp")

    return (gmst)
}


loadNichols <- function()
{
    raw <- scan("../data/nichols/data.txt", comment.char="#", quiet=T)
    mat <- matrix(raw, ncol=3, byrow=T)

    # calculate rate of sea level rise in m/year
    rise  <- (mat[2:nrow(mat), 2] - mat[1, 2]) * 7.0/20.05
    years <- (mat[2:nrow(mat), 1] - mat[1, 1])
    rate  <- rise / years

    # subtract control year
    rate  <- rate[2:length(rate)] - rate[1]

    # tabulate it
    tab <- cbind(mat[3:nrow(mat), 3], rate)
    colnames(tab) <- c("scenario", "rate")

    return (tab)
}


loadIpccSeaLevelDriftCorrect <- function(file, truncate=F)
{
    hind    <- loadIpccSeaLevel(paste("../data/ipcc/hindcast/", file, sep=""), ar4bias=F)$annual
    fore    <- loadIpccSeaLevel(paste("../data/ipcc/predict/",  file, sep=""), ar4bias=F)$annual
    control <- loadIpccSeaLevel(paste("../data/ipcc/control/",  file, sep=""), ar4bias=F)$annual
    both    <- rbind(hind, fore)
    correct <- tsDriftCorrect(both, control, truncate=truncate, f=0.25)
    obs     <- tsBias(correct, 1)

    return (list(obs=obs, endYear=last(hind[, "time"])))
}


loadIpccTempDriftCorrect <- function(file)
{
    temp  <- loadIpccTemp(paste("../data/ipcc/gmst/",         file, sep=""), ar4bias=F)$annual
    ctl   <- loadIpccTemp(paste("../data/ipcc/control_gmst/", file, sep=""), ar4bias=F)$annual
    drift <- tsDriftCorrect(temp, ctl, truncate=T, f=0.75)
    frc   <- tsBias(drift, 1)

    return (frc)
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
