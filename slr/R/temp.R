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
# temp.R
# written by Robert W. Fuller 0907xx

source('data.R')
source('plot.R')


# 4.b. on list of what Nathan wants
tempPlotPacf <- function(outfiles=F, startYear=1852, endYear=2005)
{
    newDev("pacfGreen", outfiles)

    gmst <- tsTrim(loadHadcrut()$annual, startYear, endYear)
    grst <- tsTrim(loadVinther()$annual, startYear, endYear)

    grst <- grst[, "temp"]
    gmst <- gmst[, "temp"]

    fit <- lm(grst ~ gmst)
    print(summary(fit))

    print(paste("correlation =", cor(grst, gmst)))

    pacf(
        residuals(fit),
        main="",
        ylab="Partial ACF",
        xlab="Years"
        )

    mtext("Figure 1. Southwest Greenland vs. global mean surface temp (Celsius)", outer=TRUE, side=1)

    if (outfiles) { finDev() }

    #fit2 <- lm(grst ~ poly(gmst, degree=2))
    #print(summary(fit2))
}


# 4.a. on list of what Nathan wants
tempPlotScatter <- function(outfiles=F, startYear=1852, endYear=2005, newDev=T, ar4bias=T)
{
    if (newDev) {
        newDev("grstvsgmst", outfiles)
    }

    gmst <- tsTrim(loadHadcrut(ar4bias=ar4bias)$annual, startYear, endYear)
    grst <- tsTrim(loadVinther(ar4bias=ar4bias)$annual, startYear, endYear)

    lmPlot(gmst[,"temp"], grst[,"temp"], xlab="Global mean surface temperature (C)", ylab="S&W Greenland temp (C)")

    legend("bottomright",
       legend=c("Linear fit", "Observations"),
       lty=c("solid", NA),
       lwd=c(2, NA),
       pch=c(NA, 1)
       )

    if (newDev && outfiles) { finDev() }
}


# 1 on list of what Nathan wants
tempPlotMass <- function(outfiles=F, startRow=2)
{
    newDev("massvstemp", outfiles)

    mass <- loadRignot(pure=F)
    gmst <- loadHadcrut()$annual
    forcing <- tsTrimForcing(gmst, mass)

    y <- mass   [startRow:nrow(mass),    "mass"]
    x <- forcing[startRow:nrow(forcing), "temp"]

    fit1 <- lm(y ~ x)
    fit2 <- lm(y ~ poly(x, degree=2))

    lty <- c("dashed", "solid", NA)
    lwd <- c(1, 1, NA)

    plot(x, y, xlab="Global mean surface temperature (Celsius)", ylab="GIS mass balance (Gton/a)")
    abline(fit1, lty=lty[1], lwd=lwd[1])
    curve(predict(fit2, newdata=data.frame(x=x)), add=T, lty=lty[2], lwd=lwd[2])

    legend("bottomleft",
       legend=c("Linear fit", "Quadratic fit", "Observations"),
       lty=lty,
       lwd=lwd,
       pch=c(NA, NA, 1)
       )

    if (outfiles) { finDev() }
}


tempPlotNathan <- function(outfiles=F)
{
    newDev("ntemp", outfiles, filetype="png", height=7, width=7)
    par(mfcol=c(2, 2), omi=c(0.5, 0, 0, 0))

    grst <- loadVinther(ar4bias=F)$annual
    grst <- tsTrim(grst, startYear=1852, endYear=2005)
    plot(grst, type="l", xaxp=c(1850, 2010, (2010-1850)/10), xlab="Year", ylab="S&W Greenland temp (C)",
         sub="Figure 1. South and west Greenland temperature")

    gmst <- loadHadcrut(ar4bias=F)$annual
    plot(gmst, type="l", xaxp=c(1850, 2010, (2010-1850)/10), xlab="Year", ylab="Global mean temp (C)",
         sub="Figure 2. Global mean surface temperature")

    tempPlotScatter(newDev=F, ar4bias=F)
    title(sub="Figure 3. Greenland vs. global temperature")

    #mtext("Figure 1. South and west Greenland temperature", outer=TRUE, side=1)

    if (outfiles) { finDev() }
}
