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
# plot.R
# written by Robert W. Fuller on 090809


#graphics.off()

newDev <- function(fname, outfile, height=11, width=8.5, filetype="pdf", horiz=F)
{
    fname <- paste("../figures/", fname, sep="")
    if (outfile) {
        switch(filetype, 
            png={
                fname <- paste(sep="", fname, ".png")
                png(fname, height=height, width=width, units="in", res=300)
            },
            pdf={
                fname <- paste(sep="", fname, ".pdf")
                # a4 and a4r are options for Europe
                pdf(fname, paper=ifelse(horiz, "USr", "letter"), onefile=F, height=height, width=width)
            },
            eps={
                fname <- paste(sep="", fname, ".eps")
                postscript(fname, horizontal=horiz, onefile=F, height=height, width=width)
            }, {
                stop("unimplemented filetype in newDev()")
            })
    } else {
        dev.new()
    }

    # there are problems copying to postscript and pdf devices;  hence, copy from
    # postscript and pdf devices;  this has to be enabled
    #
    if (!dev.interactive() && interactive()) {
        dev.control("enable")
    }

    # bottom, left, top, right (margins in inches)
    par(omi=c(0.25, 0.25, 0.25, 0.25))
}


# finalize devices
finDev <- function()
{
    devs <- dev.list()
    for (src in devs) {
        dev.set(which=src)

        # if printing to the device, close it
        if (!dev.interactive()) {

            # if there is a user, display what is printed
            if (interactive()) {
                dev.copy(device=dev.new)
            }

            dev.off(which=src)
        }
    }
}


lmPlot <- function(x, y, ...)
{
    plot(x, y, ...)
    abline(coef(lm(y ~ x)), lwd=2)
}


normHist <- function(x, ...)
{
    out <- hist(x, ...)
    box()

    # this is the standard, but ugly solution
    #mu <- mean(x)
    #sigma <- sd(x)
    #curve(dnorm(x, mean=mu, sd=sigma), add=T)

    xnorm <- seq(min(x), max(x), length.out=1000)
    lines(xnorm, dnorm(xnorm, mean=mean(x), sd=sd(x)))

    return (out)
}


tsColLines <- function(x, ...)
{
    for (col in safefor(2:ncol(x))) {
        lines(x[, 1], x[, col], ...)
    }
}


xshade <- function(xleft, xright, col=rgb(0.8, 0.8, 0.8))
{
    plotxy  <- par()$usr
    ybottom <- plotxy[3]
    ytop    <- plotxy[4]

    if (missing(xleft)) {
        xleft <- plotxy[1]
    }
    if (missing(xright)) {
        xright <- plotxy[2]
    }

    rect(xleft, ybottom, xright, ytop, col=col, border=F)
}


# avoid the problem of having to switch between low and high level
# plotting functions as objects are added and removed from a plot
#
emptyPlot <- function(xlim, ylim, xlab=NULL, ylab=NULL, rhs=T, box=T, top=F, xpad=T, line=2)
{
    plot.new()
    plot.window(xlim, ylim, xaxs=ifelse(xpad, "r", "i"))
    axis(1)
    axis(2)
    if (rhs) {
        axis(4, labels=F, tcl=-0.25)
    }
    if (top) {
        # positive values for tcl put the tickmarks inside the plot
        axis(3, labels=F, tcl=-0.10)
    }
    title(xlab=xlab, ylab=ylab, line=line)
    if (box) {
        box()
    }
}


labelPlot <- function(letter, line=3, where="topleft")
{
    switch(where,
        topleft={
            at <- par("usr")[4]
        },
        left={
            at <- NA
        }, {
            stop("unknown location in labelPlot()")
        })

    # left:  side=2
    # horizontal:  las=1
    # above axis label:  line=3
    # bold:  font=2
    # locate at y-axis:  at=par("usr")[4]
    #
    # font=2
    mtext(letter, side=2, las=1, line=line, at=at, cex=1.25, font=2)
}


tsErrorBars <- function(ts, shade=T, lines=!shade, xbeam=F, ibeam=!xbeam, obscol=2, lwd=1, col="gray", tick=0.5)
{
    x     <- ts[, "time"]
    upper <- ts[, obscol] + ts[, "error"]
    lower <- ts[, obscol] - ts[, "error"]

    if (lines) {
        segments(x, lower, x, upper, col=col, lwd=lwd)
        if (ibeam) {
            segments(x - tick, c(lower, upper), x + tick, c(lower, upper), col=col, lwd=lwd)
            #segments(x - tick, upper, x + tick, upper, col=col, lwd=lwd)
        }
        if (xbeam) {
            segments(x - tick, ts[, obscol], x + tick, ts[, obscol], col=col, lwd=lwd)
        }
    }

    if (shade) {
        polygon(c(x, rev(x), x[1]), c(lower, rev(upper), lower[1]), col=col, border=NA)
    }
}
