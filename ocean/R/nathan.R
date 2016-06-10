# written by Robert W. Fuller on 20101003

source('data.R')
source('roblib.R')


# columns are models, rows are times

writeNathan <- function(path="../data/picntrl/run1")
{
    filenames <- notDir(list.files(path, full.names=T))
    nmodels <- length(filenames)

    mat <- matrix(ncol=nmodels, nrow=99)

    for (i in 1:nmodels) {
        filename <- basename(filenames[i])
        ts <- loadIpccTempDriftCorrectNathan(filename)
        mat[, i] <- ts[, "temp"]
    }

    modelnames <- sub("\\..*", "", basename(filenames))
    colnames(mat) <- modelnames
    rownames(mat) <- ts[, "time"]

    write.table(mat, file="ocean_temp.txt")
}


readNathan <- function()
{
    dataframe <- read.table(file="ocean_temp.txt")
    return (as.matrix(dataframe))
}
