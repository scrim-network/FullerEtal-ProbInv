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
# grinsted.R
# written by Robert W. Fuller on 090322
# ported from Robert W. Fuller's original matlab code
# see corresponding matlab code for function descriptions for now

source("roblib.R")
source("data.R")
loadLibrary("deSolve")
loadLibrary("DEoptim")

loadSlrModel()


grinstedOde <- function(t, y, parms)
{
    mp      <- parms$mp
    gmst    <- parms$frc[[1]]
    gmstCol <- parms$spl

    s <- y[1]
    temp <- tsFindByDate(gmst, t, gmstCol)
    ds <- (mp["a"] * temp + mp["b"] - s) / mp["tau"]

    return (list(ds))
}


# Grinsted's model fit:  s0=-0.249, a=3.10, b=3.68, tau=1193
grinsted <- function(times, parms)
{
    if (length(parms$ep)) {
        parms$mp <- replace(parms$mp, names(parms$ep), parms$ep)
        parms$ep <- NULL
    }

    y <- parms$mp["s0"]
    names(y) <- "sealvl"

    if (parms$sw["log_tau"]) {
        parms$mp["tau"] <- exp(parms$mp["tau"])
    }

    # ode( y=y, times=times, func=grinstedOde,    parms=parms, method="euler" )
    euler( y=y, times=times, func="grinstedOdeC", parms=parms, dllname="slrmodel", initfunc="grinstedOdeInit" )
}


grinstedMinFn <- function(mp, fitctx)
{
    model <- grinsted(fitctx$times, list(mp=mp, frc=list(fitctx$frc), spl=2, sw=fitctx$sw, ep=fitctx$ep))

    # searching for optimal parameters for model function
    err <- sse(fitctx$obsonly, model[fitctx$mod_ind, "sealvl"])
    if (is.na(err)) {
        return (Inf)
    }

    return (err)
}


configTimes <- function(startTime, endTime, frc, obs, envir=env(), reconstruct=F, obscol=2)
{
    # must have a forcing to run the model
    times  <- frc[, "time"]

    # default: run model where there are both forcings and observations
    common <- intersect(times, obs[, "time"])

    if (missing(startTime)) {
        startTime <- ifelse(reconstruct, times[1], common[1])
        #print(paste("default startTime to", startTime))
    } else if (!(startTime %in% times)) {
        stop("no forcing for startTime")
    }

    if (missing(endTime)) {
        endTime <- ifelse(reconstruct, last(times), last(common))
        #print(paste("default endTime to", endTime))
    } else if (!(endTime %in% times)) {
        stop("no forcing for endTime")
    }

    assert(startTime < endTime)

    # times where the model is run (mod_times would be a better name)
    envir$times   <- times[ times >= startTime & times <= endTime ]

    # index into observations which intersect with where model is run
    envir$obs_ind <- which(obs[, "time"] %in% envir$times)

    # just the observations where the model is run
    envir$obsonly <- obs[envir$obs_ind, obscol]

    # times that go with obsonly
    envir$obstime <- obs[envir$obs_ind, "time"]

    # index into model output corresponding to observations
    envir$mod_ind <- which(envir$times   %in% envir$obstime)

    # index into forcings corresponding to observations (where model is run)
    envir$frc_ind <- which(frc[, "time"] %in% envir$obstime)

    envir$obs     <- obs
    envir$frc     <- frc

    return (envir)
}


removeFixedParms <- function(x, fp)
{
    if (!is.null(fp)) {
        ind <- which(names(x) %in% names(fp))
        if (length(ind)) {
            return (x[ -ind ])
        }
    }

    return (x)
}


grinstedConfig <- function(fitctx, startTime, endTime, log_tau=F, paleo=F, historical=F, ipcc=F, fp=NULL, ipccFile="miub_echo_g.txt")
{
    if (ipcc) {

        sealvl <- loadIpccSeaLevelDriftCorrect(ipccFile, truncate=T)
        obs    <- sealvl$obs
        frc    <- loadIpccTempDriftCorrect(ipccFile)

        #plot(obs)

        startTime <- obs[1, "time"]
        endTime   <- sealvl$endYear
    } else {
        if (paleo) {
            obs <- loadJevAmsterdam()
            frc <- loadMoberg()
        } else {
            obs <- loadJev()$annual
            frc <- loadHadcrut()$annual
        }
    }

    configTimes(startTime, endTime, frc, obs, envir=fitctx, reconstruct=historical)


    # a priori constraints from Grinsted's manuscript

    # Grinsted was using 1871 starting sea level for 1850 run:
    # if Grinsted started with Church & White, this is likely vestigial
    #
    #s0center <- -0.21
    s0center <- tsFindByDate(fitctx$obs, fitctx$times[1], exact=!historical)

    # how does Grinsted derive his four standard errors? (+- 0.25)
    # from earlier than 1850 based on the value he uses:
    # 1850 should have the value 0.166
    #
    stderror <- 4 * tsFindByDate(fitctx$obs, fitctx$times[1], "error", exact=!historical)

    # a manuscript that Grinsted cites states that sea level did not vary
    # more than +- 0.25 meters "from 2000 to 100 years before present"
    #

    # for reconstructions
    #if (fitctx$times[1] < fitctx$obs[1, "time"]) {
    #    stderror <- stderror + 0.25
    #}
    if (historical) {
        s0center <- 0
        stderror <- 1
    }

    # this is roughly the variability from the ordinary assimilation (!historical & !paleo)
    if (ipcc) {
        stderror <- 0.16
    }

    # N.B. this line makes dead code of preceding lines
    #
    #stderror <- 0.25

    # a minimum of 1 year for tau avoids numerical instability problems
    # with a timestep of 1 year
    #                                   s0     a       b     tau
    #fitctx$lbound <- c(s0center - stderror,  0.5, gtzero(),   50)

    # TODO:  are these the priors we want?
    fitctx$lbound <- c(s0center - stderror,  gtzero(), gtzero(),   1)
    fitctx$ubound <- c(s0center + stderror, 30.0,    5.0,   3000)

    #fitctx$ubound <- c(s0center + stderror, 10.0,    5.0,   10000)
    #fitctx$ubound <- c(s0center + stderror, 30.0,    24.0,   10000)
    names(fitctx$lbound) <- names(fitctx$ubound) <- c("s0", "a", "b", "tau")

    # plotting functions use this
    fitctx$units <- c("m", "m/C", "m", "Years", "", "", "")

    fitctx$ep <- numeric()
    names(fitctx$ep) <- character()

    fitctx$sw <- logical()
    fitctx$sw["log_tau"] <- log_tau
    if (log_tau) {
        fitctx$lbound["tau"] <- log(fitctx$lbound["tau"])
        fitctx$ubound["tau"] <- log(fitctx$ubound["tau"])
    }

    configFixedParms(fitctx, fp)
}


grinstedRunFit <- function(fitctx=grctx, startYear, endYear, log_tau=F, paleo=F, historical=F, ipcc=F, fp=NULL, ipccFile="miub_echo_g.txt")
{
    grinstedConfig(fitctx=fitctx,
        startTime=startYear, endTime=endYear,
        log_tau=log_tau, paleo=paleo, historical=historical, ipcc=ipcc, fp=fp, ipccFile=ipccFile
        )

    # advice:  increase NP and decrease F together to improve convergence
    # advice:  CR and F are generally in the range of [0.5, 1.0] for most problems
    #
    # defaults give fit of 0.0448
    # NP=500 (default  50), gave a slightly better fit (0.0447), but ran for a very long time
    #  F=0.5 (default 0.8), gave a slightly better fit (0.0447), and ran in about the same time
    # CR=0   (default 0.5), terrible
    # CR=0.9 (default 0.5, recommended 0.9), got to 0.0448 much faster, and best fit so far! (0.0440)
    #
    # CR gave the biggest difference
    #
    # combinations:
    # F=0.5, CR=0.9 performed similarly to CR=0.9
    # F=0.5, NP=250 is doing slightly better than changing F only, BUT SLOW (0.0443)
    # CR=0.9, NP=250 is converging a few steps ahead of CR=0.9 alone, BUT SLOW
    # CR=1.0, F=0.25, NP=500, itermax=500 is worse than CR=1.0 alone (at iteration 290:  0.0441)
    #
    # CR=0.75 did not converge as fast
    # CR=0.9 converged around 370 iterations (-0.2302 0.8052 0.4674 159.6999)
    # CR=1.0 converged around 300 iterations (-0.2302 0.8052 0.4674 159.6999)
    #
    # try F=0.1 performs like crap...
    #
    # conclusion:  CR=1.0 with itermax=500
    #
    control <- list(CR=1.0, itermax=500)

    fitctx$fit <- named_DEoptim(
        FUN=grinstedMinFn,
        lower=fitctx$lbound, upper=fitctx$ubound,
        control=control,
        fitctx=fitctx
        )
}


if (!exists("grctx")) {
    grctx <- env()
    #grinstedRunFit()
}
