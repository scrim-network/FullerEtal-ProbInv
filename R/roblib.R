# Copyright (C) 2009, 2010, 2016, 2017 Robert W. Fuller
# email: hydrologiccycle@gmail.com
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
# roblib.R


lsj <- function(name=".GlobalEnv")
{
    vars  <- ls(name)
    types <- sapply(sapply(vars, get), typeof)
    types <- types[ types != "closure" ]
    vars  <- names(types)

    return (vars)
}


mpmatch <- function(x, table)
{
    return (which(x == substring(table, 1, nchar(x))))
}


pindex <- function(x, table, namefn=names)
{
    return (table[ mpmatch(x, namefn(table)) ])
}


named_DEoptim <- function(FUN, lower, upper, control = list(), ...)
{
    FUN2 <- function(p, ...)
    {
        names(p) <- names(lower)
        FUN(p, ...)
    }

    out <- DEoptim(FUN2, lower, upper, control, ...)

    return (out)
}


named_c <- function(..., pos=parent.frame())
{
    varnames <- as.character(match.call(expand.dots=FALSE)$...)
    vec <- c(...)
    names(vec) <- varnames

    return (vec)
}


# AR(n) whitening (approximate; no initial value correction)
arwhiten <- function(res, rhos)
{
    nres  <- length(res)
    nrhos <- length(rhos)
    w <- res[ (nrhos + 1): nres ]
    for (i in 1:nrhos) {
        w = w - rhos[i] * res [ (nrhos + 1 - i):(nres - i) ]
    }

    return (w)
}


acfwhiten <- function(res, order=2)
{
    pac <- pacf(res, plot=FALSE)
    return (arwhiten(res, pac$acf[1:order]))
}


sse <- function(observed, model)
{
    sum( (observed - model) ^ 2 )
}


last <- function(x)
{
    return (x[ length(x) ])
}


env <- function(..., hash=T, parent=emptyenv())
{
    as_env(list(...), hash=hash, parent=parent)
}


as_env <- function(srcList, hash=T, parent=emptyenv())
{
    newEnv <- new.env(hash=hash, parent=parent)

    vars <- names(srcList)
    N    <- length(vars)
    for (i in safefor(1:N)) {
        if (!length(vars[i])) {
            stop("variables must be named")
        }
        assign(vars[i], srcList[[ i ]], envir=newEnv)
    }

    return (newEnv)
}


copy_env <- function(srcEnv, hash=T)
{
    newEnv <- new.env(hash=hash, parent=parent.env(srcEnv))

    vars <- ls(srcEnv)
    for (var in vars) {
        assign(var, get(var, envir=srcEnv), newEnv)
    }
    
    return (newEnv)
}


extraParms <- function(fn, ..., hash=T, envir=parent.frame())
{
    newEnv <- as_env(list(...), hash=hash, parent=environment(fn))
    myexpr <- substitute({
        environment(fn) <- newEnv
    })
    eval(myexpr, envir=envir)
}


dllName <- function(basename)
{
    return (paste(basename, .Platform$dynlib.ext, sep=""))
}


dynLoad <- function(basename, ..., srcname=paste(sep="", basename, ".c"), extrasrc=NULL, makevars=NULL)
{
    libName <- dllName(basename)
    error <- F

    # are any of the source files newer than the library?  if so, rebuild
    if (file.exists(libName)) {
        if (any(file.info(c(srcname, extrasrc))$mtime > file.info(libName)$mtime)) {
            error <- T
        }
    }

    # try to load the library;  catch any errors
    if (!error) {
        rc <- tryCatch(dyn.load(libName, ...), error=function(e) { error <<- T; return(e) })
    }

    # rebuild the library if it did not load, or if the source files are out of date
    if (error) {
        cmd <- paste(makevars, "R CMD SHLIB --preclean -o", libName, paste(srcname, collapse=" "))
        rc <- system(cmd, intern=F)
        if (rc != 0) {
            stop(paste("could not build library", libName))
        }
        rc <- dyn.load(libName, ...)
    }

    return (rc)
}


dynUnload <- function(basename)
{
    tryCatch(dyn.unload(dllName(basename)), error=identity)
}


dynReload <- function(basename, ...)
{
    dynUnload(basename)
    dynLoad(basename, ...)
}


loadSlrModel <- function()
{
    dynReload("slrmodel", srcname=c("slrmodel.c", "r.c"), extrasrc="r.h")
}


safefor <- function(seq)
{
    if (last(seq) < seq[1]) {
        seq <- NULL
    }

    return (seq)
}


prmatrix <- function(nbatch, xvals)
{
    prchain <- matrix(nrow=nbatch, ncol=length(xvals))
    colnames(prchain) <- as.character(xvals)
    attr(prchain, "xvals") <- xvals

    return (prchain)
}


# global data pollutes lsj()
gtzero <- function()
{
    # gtzero is used where the constraint is > 0
    return (1e-16)
}


iszero <- function(x)
{
    return (x >= 0 & x <= gtzero())
}


rmif <- function(..., list=character(), envir=parent.frame(), inherits=FALSE)
{
    varnames <- as.character(match.call(expand.dots=FALSE)$...)
    for (name in c(list, varnames)) {
        if (exists( name, envir=envir, inherits=inherits)) {
            rm(list=name, envir=envir, inherits=inherits)
        }
    }
}


burnlen <- function(chain)
{
    return (min(250000, nrow(chain) / 4))
}


burnedInd <- function(chain)
{
    start <- burnlen(chain) + 1

    return (start:nrow(chain))
}


burninInd <- function(chain)
{
    return (1:burnlen(chain))
}


rowXxx <- function(x, f, ...)
{
    rows <- nrow(x)
    vals <- numeric(length=rows)
    names(vals) <- rownames(x)
    for (row in safefor(1:rows)) {
        vals[row] <- f(x[row, ], ...)
    }

    return (vals)
}


colXxx <- function(x, f, ...)
{
    cols <- ncol(x)
    vals <- numeric(length=cols)
    names(vals) <- colnames(x)
    for (col in safefor(1:cols)) {
        vals[col] <- f(x[, col], ...)
    }

    return (vals)
}


#rowMean <- function(x, ...) { rowXxx(x, mean, ...) }
#colMean <- function(x, ...) { colXxx(x, mean, ...) }


colMode <- function(x, ...) { colXxx(x, fMode, ...) }

fMode <- function(col, ...)
{
    dens <- density(col, ...)
    ind  <- which.max(dens$y)
    return (dens$x[ind])
}


assert <- function(assertion, text="expression is FALSE")
{
    if (!assertion) {
        #text <- deparse(substitute(assertion))
        stop(text)
    }
}


rename <- function(oldname, newname, envir=parent.frame(), inherits=FALSE)
{
    old <- get(oldname,  envir=envir, inherits=inherits)
    assign(newname, old, envir=envir, inherits=inherits)
    rm(list=oldname,     envir=envir, inherits=inherits)
}


loadChains <- function(filenames, oldnames=c("daisctx", "prdaisctx"), newnames=c("as", "pr"),
                       envir=as.environment(".GlobalEnv"))
{
    for (i in 1:length(filenames)) {
        load(filenames[i], envir=envir)
        for (j in safefor(1:length(oldnames))) {
            rename(oldnames[j], paste(newnames[j], i, sep=""), envir=envir)
        }
    }
}


# chainload("../runs/paper/ar1/prperf", oldnames=c("grinassimctx", "prgrinctx"), newnames=c("gr", "pr"))
chainload <- function(basename, oldnames=NULL, newnames=NULL, envir=as.environment(".GlobalEnv"))
{
    n <- 1
    while(T) {

        filename <- paste(basename, n, sep="")
        if (!file.exists(filename)) {
            break;
        }
        load(filename, envir=envir)
        for (i in safefor(1:length(oldnames))) {
            rename(oldnames[i], paste(newnames[i], n, sep=""), envir=envir)
        }

        n <- n + 1
    }
}


thinChain <- function(chain, nthin=10000)
{
    rows <- nrow(chain)

    # would random be better?  probably not:  predictions are already
    # randomly drawn from the assimilation chain;  the assimilation chain
    # should be sampled uniformly to avoid missing excursions
    # in parameter space
    #
    chain <- chain[ seq(1, rows, len=min(nthin, rows)), ]

    return (chain)
}


sampleChain <- function(chain, nbatch, replace=T)
{
    return (chain[ sample(nrow(chain), nbatch, replace=replace), ])
}


acceptRate <- function(chain)
{
    rows   <- nrow(chain)
    accept <- 0
    for (i in safefor(2:rows)) {
        if (!isTRUE(all.equal(chain[ (i), ], chain[ (i - 1), ], check.names=F, check.attributes=F))) {
            accept <- accept + 1
        }
    }

    return (accept / rows)
}


notDir <- function(filenames)
{
    fileinfo  <- file.info(filenames)
    filenames <- rownames(fileinfo[fileinfo[, "isdir"] == F, ])

    return (filenames)
}


gmslLab <- function(year=NULL)
{
    if (is.null(year)) {
        return ("Global mean sea-level anomaly (m)")
    } else {
        return (paste("Global mean sea-level anomaly in year", year, "(m)"))
    }
}


slrGreenlandLab <- function()
{
    return ("Sea-level rise from Greenland ice (m)")
}


# DEPRECATED:  this is re-defined below
# this depends on formLibPath() and/or .robPath existing in .Rprofile;
loadLibrary <- function(package)
{
    if (FALSE == file.exists(.robPath)) {
        fqp <- formLibPath(getwd())
        print(paste("updating library path from", .robPath, "to", fqp))
        .robPath <<- fqp
    }
    path <- paste(.robPath, "/", package, sep="")
    if (length(notDir(path))) {

        # need to build the package
        install.packages(package, lib=.robPath, dependencies=TRUE, repos=paste("file:", getwd(), sep=""), type="source")
    }

    return (require(package, character.only=T))
}


loadLibrary <- function(package)
{
    library(package, character.only=T)
}


configFixedParms <- function(assimctx, fp)
{
    if (!is.null(fp)) {
        ind <- which(names(assimctx$lbound) %in% names(fp))
        if (length(ind)) {
            assimctx$lbound <- assimctx$lbound[ -ind ]
            assimctx$ubound <- assimctx$ubound[ -ind ]
            assimctx$units  <- assimctx$units [ -ind ]
        }

        # allow overriding things like max_sle
        #assimctx$ep <- append(assimctx$ep, fp)
        assimctx$ep <- replace(assimctx$ep, names(fp), fp)
    }
}


bool <- function(b)
{
    return (substring(as.character(b), 1, 1))
}


prThinChains <- function(
    prctx=prallgrgisctx,
    names=c("prchain", "otherchain", "gischain", "ds_gis", "seq_gischain", "seq_otherchain", "ds_total"),
    nthin=10000
    )
{
    for (name in names) {
        chain <- get(name, envir=prctx)
        xvals <- attr(chain, "xvals")

        chain <- thinChain(chain, nthin)

        attr(chain, "xvals") <- xvals
        assign(name, chain, envir=prctx)
    }
}


prTrimChain <- function(prchain=prdaisctx$prchain, lower, upper)
{
    xvals <- attr(prchain, "xvals")
    if (missing(lower)) {
        lower <- xvals[1]
    }
    if (missing(upper)) {
        upper <- last(xvals)
    }

    indices <- which( (xvals >= lower) & (xvals <= upper) )
    xvals   <- xvals[indices]
    prchain <- prchain[ , indices]
    attr(prchain, "xvals") <- xvals

    return (prchain)
}


prTrimChains <- function(prctx=prdaisctx, names="prchain", lower, upper)
{
    for (name in names) {
        chain <- get(name,  envir=prctx)
        chain <- prTrimChain(chain, lower, upper)
        assign(name, chain, envir=prctx)
    }
}


capitalize <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep="")
}


# TODO:  separate some parts of roblib.R into stats.R and chains.R?
rejectSample <- function(x, tgt_dense_fn)
{
    # need to cut the artifacts at the edges.  less bandwidth does it.
    # need a bigger gridsize for less smoothing.  canonical did not help.
    #
    kernel <- "box"
    bw     <- 0.25 * dpik(x, kernel=kernel)
    d      <- bkde(x, kernel=kernel, bandwidth=bw, gridsize=1001L)
    d_fn   <- approxfun(d)

    # instrumental density
    h_x   <- d_fn(x)

    # target density
    f_x   <- tgt_dense_fn(x, log=F)

    # make sure h(x) envelopes f(x)
    c     <- max(f_x) / min(h_x)

    # uniform random draw
    u     <- runif(length(x), 0, 1)

    # rejection sample formula
    index <- (u < f_x / (h_x * c))

    return (index)
}


condMeans <- function(x, y, xlim=range(x), nbins=401L, na.rm=T)
{
    nBreaks  <- nbins + 1  # +1 for the end points
    breaks_x <- seq(xlim[1], xlim[2], length.out=nBreaks)
    midpts_x <- (breaks_x[ 1:(nBreaks - 1) ] + breaks_x[ 2:nBreaks ]) / 2
    binned_x <- .bincode(x, breaks_x, include.lowest=T)
    counts_y <- tabulate(binned_x, nbins)

    sums_y <- numeric(length=nbins)
    for (bin in 1:nbins) {
        sums_y[ bin ] <- sum( y[ which(binned_x == bin) ] )
    }

    means_y <- sums_y / counts_y

    if (na.rm) {
        finite   <- which(is.finite(means_y))
        midpts_x <- midpts_x[ finite ]
        means_y  <- means_y [ finite ]
    }

    return (list(x=midpts_x, y=means_y))
}


constrain <- function(x, xlim)
{
    lo <- which(x < xlim[1])
    x[ lo ] <- xlim[1]
    hi <- which(x > xlim[2])
    x[ hi ] <- xlim[2]

    return (x)
}
