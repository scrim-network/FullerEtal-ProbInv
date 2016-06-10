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
# plots.R
# written by Robert W. Fuller on 090323

source("grinsted.R")
source("plot.R")
source("roblib.R") # arwhiten


runplots <- function(fitctx=grctx, outfiles=FALSE, whitenonly=FALSE)
{
    if (is.null(fitctx$fit$optim$bestmem)) {
        fitctx <- env()
        grinstedRunFit(fitctx=fitctx)
    }

    parms <- paste("Fitted parameters:",
        paste(names(grctx$fit$optim$bestmem),
            "=", format(grctx$fit$optim$bestmem, digits=4, trim=T),
            sep="", collapse=", "))


    # calculate residuals
    #

    slr <- grinsted(fitctx$times, list(mp=fitctx$fit$optim$bestmem, frc=list(fitctx$frc), spl=2, sw=fitctx$sw))

    # to be consistent: observations - model
    res <- fitctx$obsonly - slr[fitctx$mod_ind, "sealvl"]


    newDev("wres", outfiles)
    par(mfcol=c(2, 2), omi=c(0.5, 0, 0, 0))


    # plot whitened residuals
    #

    #newDev("wresid", outfiles)
    wres <- acfwhiten(res)
    lmPlot(
        x=fitctx$obstime[ -(1:2) ], xlab="Year",
        y=wres, ylab="Whitened sea-level residuals (m)"
        )
    mtext(parms, outer=TRUE, side=1)


    # plot whitened residuals as a PDF
    #

    #newDev("wrespdf", outfiles)
    normHist(
        wres,
        main="",
        xlab="Whitened sea-level residuals (m)",
        probability=T
        )
    mtext(parms, outer=TRUE, side=1)


    # plot partial auto-correlation of whitened residuals
    #

    #newDev("wrespacf", outfiles)
    pacf(
        wres,
        main="",
        ylab="Whitened partial ACF",
        xlab="Years"
        )
    mtext(parms, outer=TRUE, side=1)


    # plot whitened residuals at t+1 vs. t
    #

    #newDev("wserial", outfiles)
    wresAtT1 <- wres[ -1 ]
    wresAtT0 <- wres[ -(length(wres)) ]
    lmPlot(
        x=wresAtT0, xlab="Whitened sea-level residuals ( t )",
        y=wresAtT1, ylab="Whitened sea-level residuals ( t + 1 )"
        )
    mtext(parms, outer=TRUE, side=1)


    if (!whitenonly) {

        newDev("fitres", outfiles)
        par(mfcol=c(2, 2), omi=c(0.5, 0, 0, 0))


        # plot model fit
        #

        #newDev("fit", outfiles)

        # plot observations
        emptyPlot(
            xlim=range(fitctx$obstime), ylim=range(fitctx$obsonly),
            xlab="Year", ylab=gmslLab()
            )
        tsErrorBars(fitctx$obs[fitctx$obs_ind, ], shade=F)
        points(fitctx$obs[fitctx$obs_ind, ])

        # plot model
        lines(slr, lwd=2)
        legend("topleft",
               c("Model", "Data"),
               lty=c( 1, -1),
               lwd=c( 2, -1),
               pch=c(-1,  1)
               )
        mtext(parms, outer=TRUE, side=1)


        # plot residuals as a function of driving force (temperature)
        #
        #newDev("restemp", outfiles)
        lmPlot(
            x=fitctx$frc[fitctx$frc_ind, "temp"], xlab="GMST anomaly (Celsius)",
            y=res, ylab="Sea-level residuals"
            )
        mtext(parms, outer=TRUE, side=1)

if (1) {
        # check for normal distribution of whitened residuals
        #

        #newDev("resqq", outfiles)

        #qqnorm(res, main="", sub=parms)
        #mtext(
        #    "Sea-level residuals normal Q-Q plot",
        #    outer=TRUE,
        #    side=1) 

        qqnorm(wres, main="", sub="Sea-level residuals normal Q-Q plot")
        qqline(wres)
} else {
        # display the Moberg reconstruction
        #

        emptyPlot(
            xlim=range(fitctx$times),
#            ylim=range(fitctx$obsonly, slr[, "sealvl"]),
            ylim=range(slr[, "sealvl"]),
            xlab="Year", ylab=gmslLab())

        # plot observations
        #points(fitctx$obs[fitctx$obs_ind, ])

        # plot model
        lines(slr, lwd=1)
        legend("topleft",
               #c("Model", "Data"),
               c("Model"),
               lty=c( 1, -1),
               lwd=c( 2, -1),
               pch=c(-1,  1)
               )
}
        mtext(parms, outer=TRUE, side=1)


        # plot auto-correlation of residuals
        #

        #newDev("resacf", outfiles)
        acf(
            res,
            main="",
            xlab="Years"
            )
        mtext(parms, outer=TRUE, side=1)


        newDev("resdist", outfiles)
        par(mfcol=c(2, 2), omi=c(0.5, 0, 0, 0))


        # plot residuals
        #

        #newDev("resid", outfiles)
        lmPlot(
            x=fitctx$obstime, xlab="Year",
            y=res, ylab="Sea-level residuals (m)"
            )
        mtext(parms, outer=TRUE, side=1)


        # plot residuals as a PDF
        #

        #newDev("respdf", outfiles)
        normHist(
            res,
            main="",
            xlab="Sea-level residuals (m)",
            probability=T
            )
        mtext(parms, outer=TRUE, side=1)


        # plot partial auto-correlation of residuals
        #

        #newDev("respacf", outfiles)
        pacf(
            res,
            main="",
            xlab="Years"
            )
        mtext(parms, outer=TRUE, side=1)


        # plot residuals at t+1 vs. t
        #

        #newDev("serial", outfiles)

        resAtT1 <- res[ -1 ]
        resAtT0 <- res[ -(length(res)) ]
        lmPlot(
            x=resAtT0, xlab="Sea-level residuals ( t )",
            y=resAtT1, ylab="Sea-level residuals ( t + 1 )"
            )
        mtext(parms, outer=TRUE, side=1)
    }


    if (outfiles) { finDev() }
}

runplots()
