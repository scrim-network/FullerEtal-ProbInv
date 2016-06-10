# .Rprofile
# written by Robert W. Fuller in 2009

addLibPath <- function(...)
{
    args <- list(...)
    for (path in args) {

        # archdir is of the form x86_64-pc-linux-gnu-library
        archdir <- paste(R.version$platform, "-library", sep="")

        fqp     <- paste(path, "/", archdir, sep="")
        if (FALSE == file.exists(fqp)) {
            dir.create(fqp, recursive=TRUE)
        }
        .libPaths(c(fqp, .libPaths()))
    }
}


.First <- function()
{
    .robSeed <<- as.integer(Sys.time())
    set.seed(.robSeed)

    options(help_type="html")
}


#source("~/.Rprofile")
.Second <- .First

.First <- function()
{
    .Second()

    # this doesn't work too well -- breaks help
    #addLibPath(".")

    # but this is fine
    addLibPath(getwd())

    .robPath <<- .libPaths()[1]
}
