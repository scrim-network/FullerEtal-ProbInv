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
# predict.R
# written by Robert W. Fuller on 090725

source("pdfplot.R")
source("data.R") # archerPlotPredict() and rignotPlotPredict()
source("roblib.R") # thin_chain()


prplots <- function(years=c(2050, 2100, 2150, 2200), prctx=prallgrgisctx, ...)
{
    chains <- cbind( prctx$prchain[, as.character(years) ] )
    colnames(chains) <- gmslLab(years)

    #units <- rep("m", length(years))
    units <- rep("", length(years))
    pdfplot(mcmcChain=chains, units=units, chartRows=2, chartCols=2, burnin=F, width=8, height=8, ...)
}


prplot <- function(year=2100, prctx=prallgrgisctx, ...)
{
    chain <- cbind( prctx$prchain[, as.character(year) ] )
    colnames(chain) <- gmslLab(year)

    pdfplot(mcmcChain=chain, units="", chartRow=1, chartCol=1, burnin=F, width=8, height=8,
        caption=paste("Probability density function of sea-level anomaly in year",
            year), ...)
}


plotInterval <- function(pred, mean=rep(0, length(noise)), ci=.68, ...)
{
    lower <- (1 - ci) / 2
    probs <- c(lower, 1 - lower)

    xvals <- as.numeric(colnames(pred))
    interval <- colQuantile(pred, probs)
    interval[, 1] <- interval[, 1] + mean
    interval[, 2] <- interval[, 2] + mean
    interval <- cbind(xvals, interval)
    tsColLines(interval, ...)
}


prnoisesub <- function(
    ..., legends, col, lty,
    startYear=assimctx$times[1], endYear=last(assimctx$times),
    assimctx,
    outfiles=F
    )
{
    newDev("pnoise", outfiles)

    xvals <- startYear:endYear

    times <- assimctx$obstime
    if (is.null(times)) {
        times <- assimctx$times
    }
    gmsl  <- cbind(times, assimctx$obsonly, assimctx$error)
    colnames(gmsl) <- c("time", "sealvl", "error")
    gmsl  <- tsTrim(gmsl, startYear, endYear)
    slr   <- gmsl[, "sealvl"]
    err   <- gmsl[, "error"]

    xlim <- c(startYear, endYear)
    ylim <- c(min(slr) - 0.05, max(slr) + 0.05)
    emptyPlot(xlim, ylim, "Year", gmslLab())

    col <- c("black", "black", col)
    lty <- c(NA,      "solid", lty)

    # observations
    points(xvals, slr, pch=3, cex=0.5, col=col[1])
 
    # standard error
    lines(xvals, slr + err, col=col[2], lty=lty[2])
    lines(xvals, slr - err, col=col[2], lty=lty[2])

    args <- list(...)
    for (i in 1:length(args)) {
        noise <- args[[ i ]]
        noise <- noise[, as.character(startYear:endYear)]
        plotInterval(noise, slr, col=col[ i + 2 ], lty=lty[ i + 2 ])
    }

    legend(
        "bottomright",
        legend=c("Observations", "Published error", legends),
        lty=lty,
        col=col,
        lwd=c( NA, 2,  2,  2),
        pch=c( 3, NA, NA, NA)
        )

    mtext("Figure 1. Noise intervals for global mean sea level", outer=TRUE, side=1)

    if (outfiles) { finDev() }
}


# fancy <- grinassimctx$prnoise
prnoiseplot <- function(fancyNoise=fancy, simpleNoise=simple, assimctx=grinassimctx, ...)
{
    prnoisesub(
        fancyNoise, simpleNoise,
        legends=c("AR + obs (68%)", "AR only (68%)"),
        col=    c("blue",           "red"),
        lty=    c("dashed",         "dotted"),
        assimctx=assimctx,
        ...
        )
}


prqplot <- function(prctx=prallgrgisctx, ...)
{
    prqplotsub(
        prctx$prchain,
        legends=c("Mean", "95% credible interval"),
        col=c("blue", "red"),
        lty=c("solid", "dashed"),
        prctx=prctx,
        ...
        )
}


prqplotsub <- function(
    ..., legends, col, lty,
    prctx=NULL,
    args=list(...),
    N=ifelse(plotci, 2 * length(args), length(args)),
    lwd   =rep( 2, N),
    pch   =rep(NA, N),
    pt.bg =rep(NA, N),
    pt.lwd=rep(NA, N),
    pt.cex=rep( 1, N),
    xmin=as.numeric(colnames(args[[1]])[1]),
    xmax=as.numeric(last(colnames(args[[1]]))),
    obs=T, shade=T, invert=F, obserr=F,
    caption="Figure 1. Posterior predictive for global mean sea level",
    outfiles=F,
    xvals= as.numeric(colnames(args[[1]]) [
           as.numeric(colnames(args[[1]])) >= xmin
         & as.numeric(colnames(args[[1]])) <= xmax ]),
    xcaption="Year",
    ycaption=gmslLab(),
    plotci=T,
    ymin=Inf,
    ymax=-Inf,
    height=6,
    width=6,
    newdev=T
    )
{
    if (newdev) {
        newDev("ppred", outfiles, height=height, width=width)
    }

    missing_ymin <- missing(ymin)
    missing_ymax <- missing(ymax)

    prmean <- list()
    ci     <- list()

    for (i in 1:length(args)) {

        chain <- args[[ i ]]
        if (!missing(xmin) || !missing(xmax)) {
            chain <- chain[, as.character(xvals)]
        }

        prmean[[i]] <- cbind(xvals, colMeans(chain))
        #prmean[[i]] <- cbind(xvals, colMode(chain))
        ci[[i]]     <- cbind(xvals, colQuantile(chain))
        if (missing_ymin) {
            ymin <- c(min(ymin, ci[[i]][, "0.025"]))
        }
        if (missing_ymax) {
            ymax <- c(max(ymax, ci[[i]][, "0.975"]))
        }
    }

    xlim <- c(xvals[1], last(xvals))
    ylim <- c(ymin, ymax)
    emptyPlot(xlim, ylim, xcaption, ycaption)

    if (obs) {
        legends <- append(legends, "Observations")
        lty     <- append(lty,     NA)
        col     <- append(col,     "purple")
        lwd     <- append(lwd,     NA)
        pch     <- append(pch,      3)
        pt.bg   <- append(pt.bg,   NA)
        pt.lwd  <- append(pt.lwd,  NA)
        pt.cex  <- append(pt.cex, 1.0)
    }

    # assimilation range as gray shaded area;
    # draw this first so that it is in the back in terms of Z-order
    #
    if (shade) {
        shadecol <- rgb(0.8, 0.8, 0.8)

        if (invert) {
            if (xmin < prctx$assimctx$times[1]) {
                xshade(xmin, prctx$assimctx$times[1], col=shadecol)
            }
            if (xmax > last(prctx$assimctx$times)) {
                xshade(last(prctx$assimctx$times), xmax, col=shadecol)
            }
            legends <- append(legends, "Prediction range")
        } else {
            xshade(prctx$assimctx$times[1], last(prctx$assimctx$times), col=shadecol)
            legends <- append(legends, "Assimilation range")
        }

        # a sneaky way of doing boxes is to use pch=22 with a border color (col="black"),
        # a fill color (pt.bg), a larger character size (pt.cex=2.0), and a smaller border (pt.lwd=0.5)
        lty     <- append(lty,     NA)
        col     <- append(col,     "black")
        lwd     <- append(lwd,     NA)
        pch     <- append(pch,     22)
        pt.bg   <- append(pt.bg,   shadecol)
        pt.lwd  <- append(pt.lwd, 0.5)
        pt.cex  <- append(pt.cex, 2.0)
    }

    # draw borders AFTER xshade() so that shading doesn't overwrite borders
    box()

    # observations as plus signs
    if (obs) {
        times <- prctx$assimctx$obstime
        if (is.null(times)) {
            times <- prctx$assimctx$times
        }

        if (obserr) {
            ts <- cbind(times, prctx$assimctx$obsonly, prctx$assimctx$error)
            colnames(ts) <- c("time", "y", "error")
            tsErrorBars(ts, col="purple", shade=F)
        }

        cex <- ifelse(any(par("mfrow") > 1, length(xvals) > 200), 0.5, 1)

        # a temporary hack for showing all observations when using prqplot()
        # for GCM assimilation;  not needed anymore with ipccPlotFits()
        #
        #points(prctx$assimctx$obs, pch=3, cex=cex, col="purple")
        points(times, prctx$assimctx$obsonly, pch=3, cex=cex, col="purple")
    }

    j <- 1
    for (i in 1:length(args)) {

        # mean of runs as a line
        lines(prmean[[i]], col=col[j], lty=lty[j], lwd=lwd[j])
        j <- j + 1

        if (plotci) {
            # credible interval as dashed lines
            tsColLines(ci[[i]], col=col[j], lty=lty[j], lwd=lwd[j])
            j <- j + 1
        }
    }

    legend(
        "topleft",
        legend=legends,
        lty=lty,
        lwd=lwd,
        col=col,

        # this forces the Z-order
        bg="white",

        pch=pch,
        pt.bg=pt.bg,
        pt.lwd=pt.lwd,
        pt.cex=pt.cex
        )

    if (newdev) {
        mtext(caption, outer=TRUE, side=1)
    }

    if (newdev && outfiles) { finDev() }
}


# prcompare <- prgrinctx
prqplot2 <- function(prctx=prallgrgisctx, ...)
{
    prqplotsub(
        prctx$prchain,
        prcompare$prchain,
        legends=c("Mean AR(1) + obs error", "95% credible interval AR(1)", "Mean AR(2) + obs error", "95% credible interval AR(2)"),
        col=c("blue", "blue", "red", "red"),
        lty=c("solid", "dashed", "solid", "dotted"),
        prctx=prctx,
        ...
        )
}


archerPlotPredict <- function(
    prctx=prallgrgisctx,
    total=F,
    outfiles=F, newdev=T,
    ...
    )
{
    if (newdev) {
        newDev("seqpred", outfiles)
    }

    if (total) {
        name <- "All components"
        prchain <- prctx$seq_otherchain + prctx$seq_gischain
    } else {
        name <- "Non-GIS component"
        prchain <- prctx$seq_otherchain
    }

    prqplotsub(
        prchain,
        legends=c(paste(c("Mean", "95% CI"), name), "Archer prior"),
        col=c(rep("mediumblue", 2), "seagreen"),
        lty=c("solid", "dashed", "solid"),
        prctx=prctx,
        obs=F,
        shade=F,
        xvals=as.numeric(colnames(prchain)),
        xcaption="Temperature anomaly (C)",
        ycaption="Equilibrium sea-level anomaly (m)",
        caption="Figure 1. Posterior predictive for equilibrium sea level",
        newdev=F,
        outfiles=outfiles,
        ymin=-150,
        ymax=100,
        ...
        )

    # TODO:  merge with alleyPlotPrior()?
    #alleyPlotPrior()

    mcmcChain <- prctx$assimctx$chain
    bias <- mean(mcmcChain[ burnedInd(mcmcChain), "gis_temp"])
    fig613 <- loadArcher(temp_bias=-bias)

    # pch=16, 19, 20 (tried 3 and 10 as well);  tried 2-3 for cex as well
    points(
        x=fig613[, "temp"], y=fig613[, "SLE"],
        col="seagreen", lwd=2, pch=16, cex=1.5
        )
    rect(
        fig613[, "temp_lbound"], fig613[, "SLE_lbound"],
        fig613[, "temp_ubound"], fig613[, "SLE_ubound"],
        border="seagreen", lwd=2
        )

    abline(h=0, lwd=1, lty="dotted")
    abline(v=0, lwd=1, lty="dotted")

    if (newdev && outfiles) { finDev() }
}


rignotPlotPredict <- function(prctx=prallgrgisctx, outfiles=F, newdev=T, obserr=F, plotci=T, ...)
{
    if (newdev) {
        newDev("gispred", outfiles)
    }

    prqplotsub(
        prctx$ds_gis,
        legends=c("Mean", "95% CI", "Observ."),
        col=c("blue", "red", "purple"),
        lty=c("solid", "dashed", NA),
        pch=c(NA, NA, 3),
        prctx=prctx,
        obs=F,
        plotci=plotci,
        shade=F,
        newdev=F,
        outfiles=F,
        xmin=1957,
        xmax=2007,
        ycaption=slrGreenlandLab(),
        ...
        )

    massBal <- loadRignot(pure=T)

    if (obserr) {
        # we do 2*error for the 95% CI (published errors are 1 sigma, which is a 68% CI)
        ts <- cbind(prctx$assimctx$gis_times[prctx$assimctx$gis_ind], prctx$assimctx$gis_obs, ifelse(plotci, 2, 1) * prctx$assimctx$gis_err$error)
        colnames(ts) <- c("time", "SLE", "error")
        tsErrorBars(ts, shade=F, col="purple", xbeam=T, ibeam=T)
    } else {
        points(massBal[, "time"], massBal[, "SLE"], pch=3, cex=1, col="purple")
    }

    if (newdev && outfiles) { finDev() }
}


prpanel <- function(
    prctx=prallgrgisctx,
    caption="Figure 1. Posterior predictive for equilibrium sea level",
    outfiles=F
    )
{
    newDev("gispred", outfiles)

    # bottom, left, top, right (margins in inches)
    #par(mfrow=c(3, 2), omi=c(0.25, 0.25, 0.25, 0.25))
    par(mfrow=c(3, 2))

    # PDF for 2100
    prplot(prctx=prctx, ylim=c(0, 3.5), xlim=c(0, 2), newdev=F, outfiles=F)

    # forecast by component
    prqplotsub(
        prctx$gischain,
        prctx$otherchain,
        legends=c("Mean GIS", "95% CI", "Mean other", "95% CI"),
        col=c("blue", "blue", "red", "red"),
        lty=c("solid", "dashed", "solid", "dotted"),
        prctx=prctx,
        obs=F,
        shade=F,
        xmin=last(prctx$assimctx$times) + 1,
        xmax=2100,
        ymin=-.01,
        ymax=0.80,
        newdev=F,
        outfiles=F
        )

    # hindcast
    prqplot(
        xmin=prctx$assimctx$times[1],
        xmax=last(prctx$assimctx$times),
        prctx=prctx,
        shade=F,
        outfiles=F,
        ymin=-.25,
        ymax=0.05,
        newdev=F
        )

    # equilibrium sea level
    prqplotsub(
        prctx$seq_gischain,
        prctx$seq_otherchain,
        legends=c("Mean GIS", "95% CI", "Mean other", "95% CI"),
        col=c("blue", "blue", "red", "red"),
        lty=c("solid", "dashed", "solid", "dotted"),
        prctx=prctx,
        obs=F,
        shade=F,
        #xvals=as.numeric(colnames(prctx$seq_gischain)),
        xcaption="Temperature anomaly (Celsius)",
        ycaption="Equilibrium sea-level anomaly (m)",
        ymin=-2,
        ymax=15,
        newdev=F,
        outfiles=F
        )

    # fit to Rignot's mass loss data
    rignotPlotPredict(prctx=prctx, ymin=-0.5e-4, ymax=1.1e-3, outfiles=F, newdev=F)

    # fit to Alley's prior
    alleyPlotSeqGis(prctx=prctx, ymin=-1, ymax=7.5, outfiles=F, newdev=F)

    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


bool <- function(b)
{
    return (substring(as.character(b), 1, 1))
}


prnathan <- function(fname, prctx=prallgrgisctx, outfiles=F)
{
    assimctx <- prctx$assimctx

    allgrgisRunSeqPredict()

    subcaption <- paste(sep="",
        " tide_ts=",      bool(assimctx$sw["tide_ts"]),
        ", rignot_ts=",   bool(assimctx$sw["rignot_ts"]),
        ", alley_prior=", bool(assimctx$sw["alley_prior"]),
        ", tau2_prior=",  bool(assimctx$sw["tau2_prior"])
        )

    subfname <- paste(sep="",
        bool(assimctx$sw["tide_ts"]),
        bool(assimctx$sw["rignot_ts"]),
        bool(assimctx$sw["alley_prior"]),
        bool(assimctx$sw["tau2_prior"])
        )

    pdfplot(assimctx=assimctx, caption=paste("Posterior PDFs:", subcaption), outfiles=outfiles, height=7.5)
    prpanel(prctx=prctx, caption=paste("Posterior predictive:", subcaption), outfiles=outfiles)

    file.rename("../figures/gispred.pdf", paste(sep="", "../figures/pred-", subfname, ".pdf"))
    file.rename("../figures/pdfs1.pdf",   paste(sep="", "../figures/pdf-",  subfname, ".pdf"))
}


prrangecompare <- function(newdev=T, outfiles=F)
{
    if (newdev) {
        newDev("box", outfiles)
    }

    ipcc     <- c(0.18, 0.76)
    rahms    <- c(0.50, 1.40)
    grinsted <- c(0.30, 1.59)
    siddall  <- c(0.07, 0.82)
    #us       <- c(0.52, 1.37)
    us <- prallgrgisctx$prchain[, "2100"]

    #par(las=1)

    boxplot(us, siddall, grinsted, rahms, ipcc,
        width=      rep(1, 5),
        boxwex=     c(.1, rep(.0001, 4)),
        staplewex=  c(1, rep(1000, 4)),
        horizontal=T,
        names=c("This study", "Siddall", "Grinsted", "Rahmst.", "IPCC4"),
        xaxp=c(.2, 1.6, 7),
        outline=F,
        xlab=gmslLab(2100)
        )

    if (newdev && outfiles) { finDev() }
}


wais_nichols <- function(years=2000:2100, scenario=2130)
{
    startYear <- 2030
    totalRise <- 5

    ts <- matrix(nrow=length(years), ncol=2)
    colnames(ts) <- c("time", "sealvl")
    ts[, "time"] <- years

    if (!is.finite(scenario)) {
        ts[, "sealvl" ] <- rep(0, nrow(ts))
    } else {
        rate <- totalRise / (scenario - startYear)

        norise  <- which(ts[, "time"] < startYear)
        linear  <- which(ts[, "time"] >= startYear & ts[, "time"] <= scenario)
        maxrise <- which(ts[, "time"] > scenario)
        ts[ norise,  "sealvl" ] <- 0
        ts[ linear,  "sealvl" ] <- (ts[linear, "time"] - startYear) * rate
        ts[ maxrise, "sealvl" ] <- totalRise
    }
    
    return (ts)
}


rate_wais <- function(scenario=2130)
{
    startYear <- 2030
    totalRise <- 5

    if (is.finite(scenario)) {
        rate <- totalRise / (scenario - startYear) * 100
    } else {
        rate <- 0
    }

    return (rate)
}


add_wais <- function(chain=prallgrgisctx$prchain, scenario=2130)
{
    addme <- wais_nichols(years=as.numeric(colnames(chain)), scenario=scenario)[, 2]

    #chain[1:nrow(chain), ] <- chain[1:nrow(chain), ] + addme

    # don't ask why...
    for (i in 1:nrow(chain)) {
        chain[i, ] <- chain[i, ] + addme
    }

    return (chain)
}


prChainRate <- function(chain=prallgrgisctx$prchain)
{
    rows <- nrow(chain)
    cols <- ncol(chain)

    xvals <- attr(chain, "xvals")
    xvals <- xvals[ -cols ]
    rates <- prmatrix(rows, xvals)

    for (i in 1:rows) {
        rates[i, ] <- chain[i, 2:cols] - chain[i, 1:(cols-1)]
    }

    return (rates)
}


prChainThin <- function(
    ctx=prallgrgisctx,
    names=c("prchain", "otherchain", "gischain", "ds_gis", "seq_gischain", "seq_otherchain", "ds_total"),
    nthin=10000
    )
{
    for (name in names) {
        chain <- get(name, envir=ctx)

        xvals <- attr(chain, "xvals")
        chain <- thin_chain(chain, nthin)
        attr(chain, "xvals") <- xvals

        assign(name, chain, envir=ctx)
    }
}


prob_exceed <- function(chain=prallgrgisctx$prchain, threshold=(48*2.54/100))
{
    cols <- ncol(chain)
    p_exceed <- numeric(cols)
    
    for (i in 1:cols) {
        cdf <- ecdf(chain[, i])
        p_exceed[i] <- 1 - cdf(threshold)
    }

    names(p_exceed) <- colnames(chain)

    return (p_exceed)
}
