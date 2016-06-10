# ipcc.R
# written by Robert W. Fuller in 2010

source('data.R')
source('plot.R')
source('roblib.R')


ipccPlotPredict <- function(xaxis=c(2001, 2100), outfiles=F)
{
    ipccPlotOceanTemp(path="../data/sresa1b/run1/", caption="Figure 1. Comparison of IPCC predictions (SRESA1B)", xaxis=xaxis, legendloc="topleft", outfiles=outfiles, graphFile="ipcc_pred")
}


ipccPlotControl <- function(xaxis=c(2001, 2100), outfiles=F)
{
    ipccPlotOceanTemp(path="../data/picntrl/run1", caption="Figure 1. Comparison of IPCC pre-industrial control runs", xaxis=xaxis, outfiles=outfiles, graphFile="ipcc_ctl")
}


ipccPlotDriftCorrect <- function(xaxis=c(2001, 2100), outfiles=F)
{
    ipccPlotOceanTemp(path="../data/picntrl/run1", caption="Figure 1. Comparison of drift-corrected predictions", xaxis=xaxis, legendloc="topleft", correct=T, outfiles=outfiles, graphFile="ipcc_drift")
}


ipccPlotOceanTemp <- function(
    path="../data/20c3m/run1",
    caption="Figure 1. Comparison of IPCC hindcasts",
    xaxis=NULL, yaxis=NULL,
    legendloc="bottomright",
    correct=F,
    newdev=T,
    outfiles=F,
    graphFile="ipcc_hind"
    )
{
    if (newdev) {
        newDev(graphFile, outfiles)
    }

    preds <- list()
    xrange <- numeric()
    yrange <- numeric()

    filenames <- notDir(list.files(path, full.names=T))
    for (i in 1:length(filenames)) {
        filename <- filenames[i]

        if (correct) {

            filename <- basename(filenames[i])
            ts <- loadIpccTempDriftCorrectNathan(filename)

        } else {
            ts <- loadIpccTemp(filename, ar4bias=F)$annual
            if (F) {
                ts[, "temp"] <- ts[, "temp"] - ts[1, "temp"]
            }
        }

        preds[[i]] <- ts

        xrange <- range(xrange, ts[, "time"])
        yrange <- range(yrange, ts[, "temp"])
    }

    color <- rainbow(i)

    if (is.null(xaxis)) {
        xaxis <- xrange
    }
    if (is.null(yaxis)) {
        yaxis <- yrange
    }

    emptyPlot(xlim=xaxis, ylim=yaxis,
        xlab="Year", ylab="Ocean temperature (C)")

    for (i in 1:length(preds)) {
        ts <- preds[[i]]
        lines(ts[, "time"], ts[, "temp"], col=color[i])
    }

    modelnames <- sub("\\..*", "", basename(filenames))
    legend(
        legendloc,
        legend=modelnames,
        lty="solid",
        col=color,
        lwd=rep(1, length(preds))
        )

    mtext(caption, outer=TRUE, side=1)

    if (newdev && outfiles) { graphics.off() }
}
