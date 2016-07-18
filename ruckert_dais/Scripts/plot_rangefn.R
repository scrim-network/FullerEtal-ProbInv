#######################################################################
#
# plot_sf.r    15 Jan 2015
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# Function Name: plot.sf
# Parameters:
#   x - Data vector to be plotted
#   xlab - X-axis label (default = name of the data object passed
#          to the function)
#   left.tail - Should the plot highlight the left-tail instead of the 
#               right-tail? (default = F)
#   ylab - Y-axis label (default = SF  [1 - Cum. Freq.]) 
#   make.plot - Boolean value determining whether or not to make a plot
#               (default = T)
#   ... - Other parameters to be passed to the plot() function
#
#
#######################################################################
# Function to plot a range on an existing plot
plotrange <- function(low.number, median, high.number, year=T, height=F, color){
  
  range <- c(low.number, median, high.number)
  
  if(year){
    plot.position <- c(year, year, year)
    
    points(plot.position[1:2], c(range[1], range[3]), col=color, cex=1, pch="-")
    points(plot.position[1], range[2], col=color, cex=1, pch=8)
    segments(plot.position[1], range[1], plot.position[3], range[3], lwd=2, col=color)
  }
  
  if(height){
    plot.position <- c(height, height, height)
    
    points(c(range[1], range[3]), plot.position[1:2], col=color, cex=1, pch="|")
    points(range[2], plot.position[1], col=color, cex=1, pch=8)
    segments(range[1], plot.position[1], range[3], plot.position[3], lwd=2, col=color)
  }
}

#Function to plot a blox plot on an existing plot with adding lines showing select probabilities
add.hor.box <- function(data.numbers, probabilities, width.size, where.at, tick.length, line.width, color){
  quants <- quantile(data.numbers, probabilities)
  boxplot(data.numbers, add=TRUE, horizontal=TRUE, axes=FALSE, outline=FALSE, col=color, boxwex=width.size, at=where.at)
  segments(x0 = quants, y0 = rep(where.at - tick.length, 2),
           x1 = quants, y1 = rep(where.at + tick.length, 2), col = c("blue", "blue"), lwd = line.width)
}