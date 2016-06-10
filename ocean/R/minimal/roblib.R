# roblib.R
# written by Robert W. Fuller in 2009


notDir <- function(filenames)
{
    fileinfo  <- file.info(filenames)
    filenames <- rownames(fileinfo[fileinfo[, "isdir"] == F, ])

    return (filenames)
}


last <- function(x)
{
    return (x[ length(x) ])
}
