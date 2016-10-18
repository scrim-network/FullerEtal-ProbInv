##################################################
#
# plot_sf_example.r   January 15, 2015
#
# Modified March 2015 by Kelsey Ruckert to add pdf and cdf function.
# Plus option to add mulitple lines to a single plot.
#
# Example use of the plot.sf() function in
# plot_sf.r script.
#
##################################################

# Source the function
source("plot_sf.r")

# Set the seed
set.seed(1234)

# Generate some data to use
my.norm <- rnorm(10000, 10, 2)
my.unif <- runif(10000)
my.weib <- rweibull(10000, 20, 5)
my.lnorm <- rlnorm(10000, 1, 0.5)


# Make the plots ----------------------
dev.new()
par(mfrow=c(2,2), mar=c(5,4,1,1)+0.1)

# Default plot settings
plot.sf(my.norm)

# Function wraps the standard "plot" function, so you can pass
# the standard "plot" parameters to the function
plot.sf(my.unif, type="l", lwd=2, col="blue", bty="l",
        ylab="Survival", xlab="Uniform Distribution")

# If the parameter "left.tail" is true, the plot turns into 
# a cumulative frequency plot (kind of like a CDF) that's plotted
# on a log scale.  This is good for when your data exhibits a left or
# negative skew.
plot.sf(my.weib, type="l", left.tail=TRUE, xlab="Left-tailed Weibull Dist.")

# The function invisibly returns the survival function value.
lnorm.sf <- plot.sf(my.lnorm, type="l")
points(my.lnorm, lnorm.sf$order.sf, col="red")
legend("topright", bty="n", 
       legend=c("Function Call", "Using returned values"), 
       lty=c(1,NA), pch=c(NA,1), col=c("black", "red") )

# The 'make.plot' parameter toggles plotting.
# Useful if you just want the survival function values. 
# Especially, if you want to plot multiple SF lines on 
# the same plot.
norm.sf <- plot.sf(my.norm, make.plot=FALSE)

# Example plot of multiple SF lines with plotting 
# the pdf and cdf.
library(RColorBrewer)
test.colors = brewer.pal(5, "BrBG")


dev.new()
par(mfrow=c(3,1), mgp=c(1.5,.5,0), mar=c(3.5,4,1,1))
plot(lnorm.sf$pdf, main="", col=test.colors[2], lwd=3, ylab="Probabilty density", 
     xlab="my.lnorm and my.norm", yaxt="n") 
lines(norm.sf$pdf, lwd=3, col=test.colors[4])

plot(lnorm.sf$cdf, main="", col=test.colors[2], lwd=3, ylab="Cumulative density", 
     xlab="my.lnorm and my.norm") 
lines(norm.sf$cdf, lwd=3, col=test.colors[4])

plot(lnorm.sf$sf.num, lnorm.sf$sf, log="y", type="l", ylab="Survival function 
     [1-cumulative frequency]", lwd=3, col=test.colors[2],
     xlab="my.lnorm and my.norm",main="",ylim=c(1e-04,1), yaxt="n") 
lines(norm.sf$sf.num, norm.sf$sf, lwd=3, col=test.colors[4]) 
axis(2, at=10^(-4:0), label=parse(text=paste("10^", -4:0, sep="")))



