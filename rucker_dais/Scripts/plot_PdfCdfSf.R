#######################################################################
#
# plot_PdfCdfSf.r    Aug 2015
#
# Author Gregory Garner (ggg121@psu.edu)
# Edits by: Kelsey Ruckert (klr324@psu.edu)
#
# Function that returns the survival function, probability density of SLE values, and
# the probability density of parameters for the DAIS model
#
# This is a modification to Gregory Garners "plot_sf.r" function:
#
# Version History:
#   1.0 - 15 Jan 2015 - Initial coding (G.G.)
#   1.1 - 03 Feb 2015 - Added support for left.tail (G.G.)
#   1.2 - Aug 2015 - added the addition of estimating probability for the DAIS model (K.R.)
#
# THIS CODE IS PROVIDED AS-IS WITH NO WARRANTY (NEITHER EXPLICIT
# NOR IMPLICIT).  I SHARE THIS CODE IN HOPES THAT IT IS USEFUL, 
# BUT I AM NOT LIABLE FOR THE BEHAVIOR OF THIS CODE IN YOUR OWN
# APPLICATION.  YOU ARE FREE TO SHARE THIS CODE SO LONG AS THE
# AUTHOR(S) AND VERSION HISTORY REMAIN INTACT.
#
# The function invisibly returns the survival function value for the
# passed data vector, even if make.plot is false.  This is useful when
# you want to customize your plot or use the survival function data
# in further analyses.
#
# See the accompanying "plot_sf_example.r" file for usage and examples.
#
#######################################################################
# Function to calculate the pdfs of each parameter
parameter.pdfs <- function(year.parameter) {
  
  gam.pdf = density(year.parameter[,1]); alp.pdf = density(year.parameter[,2])
  mu.pdf = density(year.parameter[,3]); eta.pdf = density(year.parameter[,4])
  po.pdf = density(year.parameter[,5]); kap.pdf = density(year.parameter[,6])
  fo.pdf = density(year.parameter[,7]); ho.pdf = density(year.parameter[,8])
  co.pdf = density(year.parameter[,9]); bo.pdf = density(year.parameter[,10])
  s.pdf = density(year.parameter[,11]); sigma.pdf = density(year.parameter[,12])
  
  return(list(gam = gam.pdf, alp = alp.pdf, mu = mu.pdf, eta = eta.pdf, po = po.pdf,
              kap = kap.pdf, fo = fo.pdf, ho = ho.pdf, co = co.pdf, bo = bo.pdf, s = s.pdf, sigma = sigma.pdf))
}

# Function to find SLE values in certain years
fn.prob.proj <- function(dais.pre.cali, year.pcs, survive.num, un.constr=F) {
  
  if(un.constr){
    prob_proj = mat.or.vec(survive.num, length(year.pcs)) #create an empty matrix
    for (i in 1:length(year.pcs)){
      prob_proj[,i] = dais.pre.cali[,year.pcs[i]]
    } #End for loop
    
  }  else {
    prob_proj = mat.or.vec(length(survive.num), length(year.pcs)) #create an empty matrix
    for (i in 1:length(year.pcs)){
      prob_proj[,i] = dais.pre.cali[survive.num[1:length(survive.num)],year.pcs[i]]
    } #End for loop
  } #End else bracket
  
  colnames(prob_proj, do.NULL = FALSE) # Name the columns in the matrix
  colnames(prob_proj) = c("-120000", "-20000", "-6000", "2002", "2050", "2100","2300")
  return(prob_proj)
}

# Function to calculate the pdf cdf and SF for certain years
plot.sf <- function(x, xlab=deparse(substitute(x)), left.tail=F,
  ylab=ifelse(left.tail, "SF [Cum. Freq.]", "SF  [1 - Cum. Freq.]"),
  make.plot=T, ...)
{
  num.x <- length(x)
  num.ytics <- floor(log10(num.x))
  sf <- seq(1,1/num.x,by=-1/num.x)
  pdf.x <- density(x)
  cdf.x <- ecdf(x)
  
  if(left.tail){
    order.x <- order(x, decreasing=T)
    order.sf <- sf[order(order.x)]
    
  }  else {
    order.x <- order(x)
    order.sf <- sf[order(order.x)]
  }
  
  if(make.plot) {
    plot(x[order.x], sf, log="y", xlab=xlab, ylab=ylab, yaxt="n", ...)
    axis(2, at=10^(-num.ytics:0), label=parse(text=paste("10^", -num.ytics:0, sep="")), las=1)
  }
  invisible(return(list(sf.num = x[order.x], pdf = pdf.x, cdf = cdf.x, sf = sf, num.ytics = num.ytics)))
}
