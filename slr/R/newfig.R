# Copyright 2010, 2016 Robert W. Fuller <hydrologiccycle@gmail.com>
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
# newfig.R

source("predict.R")
source("roblib.R")


# from run_tau900.R:
#prallgrgisctx$prrates <- rateCalc(prallgrgisctx$prchain,
#    xvals=(last(prallgrgisctx$assimctx$times) + 1):2100,
#    rows=sample(1:nrow(prallgrgisctx$prchain), 10000, replace=F))

# fit a slope to the prior 30 years of predictions for each year
rateCalc <- function(prchain,
    xvals=attr(prchain, "xvals")[ -(1:(years-1)) ], years=30,
    rows=sample(1:nrow(prchain), 1000, replace=F)
    )
{
    rates <- prmatrix(length(rows), xvals)
    cols  <- which(attr(prchain, "xvals") %in% xvals)
    x     <- 1:years

    # do this by column in order to only calculate idx in the outer loop;
    # this may be sub-optimal if length(rows) is large
    #
    for (j in 1:length(cols)) {
        col <- cols[j]
        idx <- (col - years + 1):col

        for (i in 1:length(rows)) {
            row <- rows[i]

            ts  <- prchain[ row, idx ]
            fit <- lm(ts ~ x)
            rates[i, j] <- coef(fit)[2]
        }
    }

    return (rates)
}


figRandSpag <- function(
    prctx=prallgrgisctx,
    xvals=(last(prctx$assimctx$times) + 1):2100,
    realiz=10,
    rate="average",
    outfiles=F
    )
{
    newDev("rand4", outfiles)


    # reserve lines to use outer=T for the lower axis and label;
    # allows using 0.5 in par(fig) and getting equally sized plots;
    # one line is 0.2 inches:  par("csi") or par("cin")[2]
    #
    #par(oma=c(4, 1, 1, 1))
    par(omi=c(0.75, 0.25, 0.25, 0.25))


    # top of figure

    par(fig=c(0, 1, 0.5, 1))
    par(mar=c(0, 3, 0, 1))
    plot.new()

    cictx <- ciCalc(prctx$prchain, xvals=xvals)
    plot.window(xlim=range(xvals), ylim=cictx$range, yaxs="i")

    axis(2) # left
    axis(1, labels=F, tcl=-0.10) # bottom
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right
    title(ylab=gmslLab(), line=2)

    samples <- sample(1:nrow(prctx$prchain), realiz + 1, replace=F)
    for (row in samples) {
        lines(xvals, prctx$prchain[row, cictx$cols], lty="solid", col="lightgray", lwd=1)
    }
    ciPlot(cictx)
    lines(xvals, prctx$prchain[realiz + 1, cictx$cols], lty="solid", col="black", lwd=1)
    box()

    legend(
        "topleft",
        legend=c("Mean of simple model", "95% credible interval", paste(realiz, "sample runs"), "Example run"),
        lty=c("solid", "dashed", "solid", "solid"),
        lwd=c(2, 2, 1, 1),
        col=c("red", "red", "lightgray", "black")

        # this forces the Z-order
        #bg="white",
        )
    labelPlot("a")

    # bottom of figure

    par(fig=c(0, 1, 0, 0.5), new=T)
    par(mar=c(0, 3, 0, 1))
    plot.new()

    switch(rate,
        average={
            if (is.null(prctx$prrates)) {
                prctx$prrates <- rateCalc(prctx$prchain, xvals)
            }
            rates <- prctx$prrates
        },

        # these are not the same thing, even if ds_total includes noise;
        # subtracting one year from the prior, as prChainRate() does, subtracts
        # some of the auto-correlated noise;  hence, prChainRate() produces
        # tighter credible intervals
        #

        nonoise={
            rates <- prctx$ds_total
        },
        noise={
            rates <- prChainRate(prctx$prchain)
        }, {
            stop("unknown rate in prrandspag()")
        })

    rates <- rates * 1000
    cictx <- ciCalc(rates, xvals=xvals)
    plot.window(xlim=range(xvals), ylim=cictx$range, yaxs="i")

    axis(1, outer=T) # bottom
    axis(2) # left
    axis(3, labels=F, tcl=-0.10) # top
    axis(4, labels=F, tcl=-0.25) # right

    title(ylab="30a mean rate of global mean sea-level rise (mm/a)", line=2)
    title(xlab="Year", line=2, outer=T)

    for (row in samples) {
        if (rate == "average") {
            ts <- rateCalc(prctx$prchain, xvals, rows=row)
            ts <- ts * 1000
        } else {
            ts <- rates[row, cictx$cols]
        }
        lines(xvals, ts, lty="solid", col="lightgray", lwd=1)
    }

    ciPlot(cictx)
    lines(xvals, ts, lty="solid", col="black", lwd=1)
    box()

    labelPlot("b")

    if (outfiles) { finDev() }
}


figUber <- function(
    prctx=prallgrgisctx,
    #caption="Figure 1. Hindcasts and posterior predictive PDFs for global mean sea level",
    caption=NULL,
    outfiles=F
    )
{
    newDev("randpred", width=8, height=8, outfiles)

    par(mfcol=c(2, 2))
    par(mar=c(4, 3, 0, 3))

    taus <- c(pr1$assimctx$ep["tau"], pr2$assimctx$ep["tau"], pr3$assimctx$ep["tau"])
    legends <- expression()
    for (tau in taus) {
        legends <- append(legends, bquote(paste(tau[1]==.(tau), " a")))
    }

    col     <- c("black", "blue", "red")
    lty     <- c("solid", "dotdash", "longdash")

    column <- as.character(2100)

    pdflim <- c(0, 1.1)


    # hindcast
    #

    times <- prctx$assimctx$times
    cictx <- ciCalc(pr1$prchain, pr2$prchain, pr3$prchain, xvals=times)
    emptyPlot(xlim=range(times), ylim=cictx$range, xlab="Year", ylab=gmslLab())
    ts <- cbind(times, prctx$assimctx$obsonly, prctx$assimctx$error)
    colnames(ts) <- c("time", "y", "error")
    tsErrorBars(ts, col="purple", shade=F)
    ciPlot(cictx, col=col, lty=lty, plotci=F)
    legend(
        "topleft",
        legend=c(legends, "Jevrejeva et al. 2006"),
        lty=c(lty, NA),
        lwd=c(rep(2, 3), NA),
        col=c(col, "purple"),
        pch=c(rep(NA, 3), 3)
        )
    labelPlot("a")


    # Rignot fit
    #

    massBal <- loadRignot(pure=T)
    times <- massBal[, "time"]
    times <- (times[1] - 1):(last(times) + 1)
    cictx <- ciCalc(pr1$ds_gis, pr2$ds_gis, pr3$ds_gis, xvals=times)
    emptyPlot(xlim=range(times), ylim=cictx$range, xlab="Year", ylab=slrGreenlandLab())
    legend(
        "topleft",
        legend=c(legends, "Rignot et al. 2008"),
        lty=c(lty, NA),
        lwd=c(rep(2, 3), NA),
        col=c(col, "purple"),
        pch=c(rep(NA, 3), 3)
        )

    # could use xbeam=T, ibeam=T below, but would need to match to size of cex?
    points(massBal[, "time"], massBal[, "SLE"], pch=3, cex=0.5, col="purple")

    ts <- cbind(prctx$assimctx$gis_times[prctx$assimctx$gis_ind], prctx$assimctx$gis_obs, prctx$assimctx$gis_err$error)
    colnames(ts) <- c("time", "SLE", "error")
    tsErrorBars(ts, shade=F, col="purple")
    ciPlot(cictx, col=col, lty=lty, plotci=F)

    labelPlot("b")


    pdfPlots(
        pr1$prchain, pr2$prchain, pr3$prchain,
        column=column,
        legends=legends,
        burnin=F,
        xlim=pdflim,
        col=col, lty=lty,
        yline=2
        )
    labelPlot("c")


    cdfPlots(
        pr1$prchain, pr2$prchain, pr3$prchain,
        column=column,
        xlim=pdflim,
        col=col, lty=lty
        )
    labelPlot("d")


    # figure title
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


figRand <- function(year=2100, outfiles=F)
{
    column=as.character(year)

    pdflim <- c(0, 5)

    low_wais  <- "Never"
    med_wais  <- 2330
    high_wais <- 2130

    # order here should match order in calls to add_wais()
    taus <- c(
        pr3$assimctx$ep["tau"],
        pr1$assimctx$ep["tau"],
        pr2$assimctx$ep["tau"],
        pr3$assimctx$ep["tau"],
        pr1$assimctx$ep["tau"]
        )
    rates <- format(c(
        rate_wais(low_wais),
        rate_wais(low_wais),
        rate_wais(med_wais),
        rate_wais(high_wais),
        rate_wais(high_wais)
        ), digits=2)

    legends <- expression()
    for (i in 1:length(taus)) {
        legends <- append(legends, bquote(paste(WAIS==.(rates[i]), " m/c., ", tau[1]==.(taus[i]), " a")))
    }

    chains <- list(
        add_wais(pr3$prchain, low_wais),
        add_wais(pr1$prchain, low_wais),
        add_wais(pr2$prchain, med_wais),
        add_wais(pr3$prchain, high_wais),
        add_wais(pr1$prchain, high_wais)
        )

    # my scenarios
    col <- c("red", "black", "blue", "red", "black")

    # wais scenarios
    lty <- c("solid", "solid", "dotdash", "longdash", "longdash")

    #lty <- rep("solid", 5)
    lwd <- 2


    newDev("rand2", outfiles)
    par(mfrow=c(2, 1))

    pdfPlots(
        chains=chains,
        column=column,
        legends=legends,
        burnin=F,
        col=col, lty=lty, lwd=lwd,
        xlim=pdflim,
        legendloc="top",
        yline=2
        )

    cdfPlots(
        chains=chains,
        column=column,
        xlim=pdflim,
        lwd=lwd,
        col=col, lty=lty
        )


    newDev("rand3", outfiles)
    par(mfrow=c(2, 1))


    # probability of exceeding 48 inches
    xlim <- c(last(pr1$assimctx$times) + 1, 2100)
    emptyPlot(
        xlim=xlim, ylim=c(0, 1),
        xlab="Year", ylab="Probability of exceeding 48 inches",
        xpad=F
        )
    years <- xlim[1]:xlim[2]
    for (i in 1:length(chains)) {
        exceed <- prProbExceed(chains[[i]])
        lines(years, exceed[ as.character(years) ], lty=lty[i], col=col[i], lwd=lwd)
    }
    legend(
        "topleft",
        legend=legends,
        lty=lty,
        col=col,
        lwd=rep(lwd, length(lty))
        )


    # calculate slr rate
    rate_chains <- list()
    for (i in 1:length(chains)) {
        rate_chains[[i]] <- prChainRate(chains[[i]]) * 1000
    }

    column <- as.character(2050)
    pdfPlots(
        chains=rate_chains,
        column=column,
        legends=legends,
        burnin=F,
        col=col, lty=lty, lwd=lwd,
        legendloc=NULL,
        xlab=paste("Rate of global mean sea-level rise in year", column, "(mm/a)"),
        yline=2
        )
    rm(rate_chains)

    if (outfiles) { finDev() }
}


figIdent <- function(
    #caption="Figure 2. Identifiability problem in posterior PDFs for perfect model, AR(1)",
    caption=NULL,
    plotmeans=F,
    column=as.character(2100),
    outfiles=F
    )
{
    newDev("ident", width=8.5, height=4, outfiles)

    par(mfcol=c(1, 2))

    # bott, left, top, right
    par(mar=c(4, 2, 0, 2))


    pdfPlots(
        assim1$chain, assim2$chain, assim3$chain, assim4$chain, assim5$chain,
        column="tau",
        lty=rep("solid", 5), legends=paste("Run no.", 1:5),
        col=c("black", "red", "blue", "green", "purple", "yellow", "gray"),

        # note that Tau[1] is an option also
        xlab=expression(paste(tau[1], " (years)")),

        burnin=T, truth=T,
        trueval=assim1$init_mp["tau"], xlim=c(assim1$lbound["tau"], assim1$ubound["tau"]),
        plotmeans=plotmeans,
        #legendloc="topleft",
        legendloc=NULL
        )
    labelPlot("a", line=2)

    chains <- list(pr1$prchain, pr2$prchain, pr3$prchain, pr4$prchain, pr5$prchain)
    cictx <- ciCalc(chains=chains, xvals=as.numeric(column), probs=c(0.0005, 0.9995))
    pdfPlots(
        chains=chains,
        lty=rep("solid", 5), legends=paste("Run no.", 1:5),
        col=c("black", "red", "blue", "green", "purple", "yellow", "gray"),
        burnin=F, truth=T,
        trueval=pr1$assimctx$true_predict,
        column=column,
        xlim=cictx$range,
        plotmeans=plotmeans,
        legendloc=NULL
        )
    labelPlot("b", line=2)


    # figure title:  TODO:  left justify (adj=0)?
    mtext(caption, outer=TRUE, side=1, font=2)

    if (outfiles) { finDev() }
}


# per Klaus: change Our Study to This Study
figCompareGrin <- function(..., legends=c("Grinsted", "This study"))
{
    pdfPlots(
        prgrinctx$prchain, prallgrgisctx$prchain,
        legends=legends, col=c("red", "blue"), lty=c("dashed", "solid"),
        burnin=F,
        ...)
}
