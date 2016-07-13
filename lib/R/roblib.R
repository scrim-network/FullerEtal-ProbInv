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
# roblib.R


lsj <- function(name=".GlobalEnv")
{
    vars  <- ls(name)
    types <- sapply(sapply(vars, get), typeof)
    types <- types[ types != "closure" ]
    vars  <- names(types)

    return (vars)
}


colQuantile <- function(x, probs=c(0.025, 0.975), ...)
{
    cols <- ncol(x)
    q    <- prmatrix(cols, probs)
    for (col in safefor(1:cols)) {
        q[col, ] <- quantile(x[, col], probs=probs, ...)
    }

    return (q)
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


dynLoad <- function(basename, ..., srcname=paste(sep="", basename, ".c"), extrasrc=NULL)
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
        cmd <- paste("R CMD SHLIB --preclean -o", libName, paste(srcname, collapse=" "))
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
    return (dyn.unload(dllName(basename)))
}


dynReload <- function(basename, ...)
{
    tryCatch(dynUnload(basename), error=identity)
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
    return (min(100000, nrow(chain) / 4))
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


colMode <- function(x, na.rm=F)
{
    cols <- ncol(x)
    modes <- numeric(length=cols)
    names(modes) <- colnames(x)
    for (col in safefor(1:cols)) {
        dens <- density(x[, col], na.rm=na.rm)
        ind <- which.max(dens$y)
        modes[col] <- dens$x[ind]
    }

    return (modes)
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


thin_chain <- function(chain, nthin=10000)
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


acceptRate <- function(chain)
{
    rows   <- nrow(chain)
    accept <- 0
    for (i in safefor(2:rows)) {
        if (any(chain[ (i), ] != chain[ (i - 1), ])) {
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


# this depends on formLibPath() and/or .robPath existing in .Rprofile
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
