# plot.R
# written by Robert W. Fuller on 090809


graphics.off()

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

    # bottom, left, top, right (margins in inches)
    par(omi=c(0.25, 0.25, 0.25, 0.25))
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
