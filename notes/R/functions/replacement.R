# replacement functions
#
# need a replacement function for appending to a list/vector?
#


# not sure if first() functions are useful
first <- function(x)
{
    return (x[1])
}


`first<-` <- function(x, value)
{
    x[1] <- value
    return (x)
}


`last<-` <- function(x, value)
{
    x[ length(x) ] <- value
    return (x)
}


last <- function(x)
{
    return (x[ length(x) ])
}


lrow <- function(x)
{
    return (x[ nrow(x), ])
}


`lrow<-` <- function(x, value)
{
    x[ nrow(x), ] <- value
    return (x)
}


lcol <- function(x)
{
    return (x[ , ncol(x) ])
}


`lcol<-` <- function(x, value)
{
    x[ , ncol(x) ] <- value
    return (x)
}
