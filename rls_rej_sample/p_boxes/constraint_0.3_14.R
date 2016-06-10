############################################
##  file: p_plot.R
############################################
rm(list = ls(all = TRUE))
library(timeSeries)

# c* > 14 mm/yr + 0.3 mm/yr (t* - 2010)

t=seq(2000,2100,.25)
t.length=length(t)

c.red=14 + 0.3*(t-2010)


##############################################
### uniform distribution
##############################################
array.u=read.table("../uniform/output/array_uniform.txt")
n.u=length(array.u[,1])
 print(n.u)
filter.red_unif=array(NA, dim=c(n.u,6))

for(n in 1:n.u)                                         
  {
    for(i in 1:t.length)
    {
    if (array.u[n,5] < t[i]  && array.u[n,6] > c.red[i])
      {
       filter.red_unif[n,1]=array.u[n,2]
       filter.red_unif[n,2]=array.u[n,3]
       filter.red_unif[n,3]=array.u[n,4]
       filter.red_unif[n,4]=array.u[n,5]
       filter.red_unif[n,5]=array.u[n,6]
       filter.red_unif[n,6]=array.u[n,7]
      }
    }
  }

filter.red_unif=removeNA(filter.red_unif)
#print(filter.red_unif)
  print(length(filter.red_unif[,1]))

##############################################
### beta distribution
##############################################
array.b=read.table("../beta/output/array_beta.txt")
n.b=length(array.b[,1])
 print(n.b)
filter.red_beta=array(NA, dim=c(n.b,6))

for(n in 1:n.b)
  {
    for(i in 1:t.length)
    {
    if (array.b[n,5] < t[i]  && array.b[n,6] > c.red[i])
      {
       filter.red_beta[n,1]=array.b[n,2]
       filter.red_beta[n,2]=array.b[n,3]
       filter.red_beta[n,3]=array.b[n,4]  
       filter.red_beta[n,4]=array.b[n,5]
       filter.red_beta[n,5]=array.b[n,6]
       filter.red_beta[n,6]=array.b[n,7]
      }
    }
  }
     
filter.red_beta=removeNA(filter.red_beta)
#print(filter.red_beta)
  print(length(filter.red_beta[,1]))


### calculate p-values
p_box.red_unif=length(filter.red_unif[,1])/n.u
 print(p_box.red_unif)

p_box.red_beta=length(filter.red_beta[,1])/n.b
 print(p_box.red_beta)


### plot the pairs  ###

pdf(file="./plots/p_box_0.3_14.pdf") #height=10, width=6)
par(mfcol=c(2,2))

plot(array.b[,5],array.b[,6],pch=19,col="black",cex=0.15,
     main="Beta",
    xlim=c(2010,2100),
     ylim=c(-10,35),
     xlab="t* [year]",
     ylab="c* [mm/year]")
 lines(t,c.red,lwd=5,lty=1,col="blue")
mtext(paste("p ~ ",  format(p_box.red_beta, digits=2)),side=3,  col="blue", cex=1.25)

plot(array.u[,5],array.u[,6],pch=19,col="black",cex=0.15,
     main="Uniform",
    xlim=c(2010,2100),
     ylim=c(-10,35),
     xlab="t* [year]",
     ylab="c* [mm/year]")
 lines(t,c.red,lwd=5,lty=1,col="blue")
mtext(paste("p ~ ",  format(p_box.red_unif, digits=2)),side=3,  col="blue", cex=1.25)


dev.off()
