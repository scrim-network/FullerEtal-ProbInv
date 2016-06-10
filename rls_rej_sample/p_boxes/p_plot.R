############################################
##  file: p_plot.R
############################################
rm(list = ls(all = TRUE))
library(timeSeries)

t=seq(2000,2100,.25)
t.length=length(t)

c.red=14 + 0.3*(t-2010)
c.green=11 + 0.4*(t-2010)
c.blue=5 + 0.4*(t-2010)


##############################################
### uniform distribution
##############################################
array.u=read.table("../uniform/output/array_uniform.txt")

print(min(array.u[,5]))
print(max(array.u[,5]))
print(min(array.u[,6]))
print(max(array.u[,6]))

##############################################
### beta distribution
##############################################
array.b=read.table("../beta/output/array_beta.txt")

print(min(array.b[,5]))
print(max(array.b[,5]))
print(min(array.b[,6]))
print(max(array.b[,6]))

### plot the pairs  ###

pdf(file="./plots/all_constraints.pdf") #height=10, width=6)
par(mfcol=c(2,2))

plot(array.b[,5],array.b[,6],pch=19,col="black",cex=0.1,
     main="Beta",
    xlim=c(2010,2100),
     ylim=c(-10,35),
     xlab="t* [year]",
     ylab="c* [mm/year]")
 lines(t,c.red,lwd=5,lty=1,col="red")
 lines(t,c.green,lwd=5,lty=1,col="green")
 lines(t,c.blue,lwd=5,lty=1,col="blue")
mtext("p ~ .13",  side=3, adj=0, col="green", cex=1.1)
mtext("p ~ .14",side=3,  col="red", cex=1.1)
mtext("p ~ .35",side=3, adj=1,  col="blue", cex=1.1)

plot(array.u[,5],array.u[,6],pch=19,col="black",cex=0.1,
     main="Uniform",
    xlim=c(2010,2100),
     ylim=c(-10,35),
     xlab="t* [year]",
     ylab="c* [mm/year]")
 lines(t,c.red,lwd=5,lty=1,col="red")
 lines(t,c.green,lwd=5,lty=1,col="green")
 lines(t,c.blue,lwd=5,lty=1,col="blue")
mtext("p ~ .16",side=3, adj=0,  col="green", cex=1.1)
mtext("p ~.16",  side=3,col="red", cex=1.1)
mtext("p ~ .34" , side=3, adj=1,  col="blue", cex=1.1)


dev.off()
