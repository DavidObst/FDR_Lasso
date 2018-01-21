source("rsandwich.R")
source("lsandwich.R")
source("epsilonDT.R")
source("powermax.R")
source("fdrlasso.R")

library(latex2exp)

plot_qstar <- function(delta,epsilon,mode="empty",col="blue",ylim=c(0,1),main=NULL,add=FALSE)
{
## Plots the (TPP,FDP=q*(TPP))  boundary region
  
## Arguments:
#    delta - a float strictly greater than 0
#    epsilon - a float between 0 and 1
#    mode - a string, either "empty" or "full". The former plots only the frontier, the latter fills
# the impossible region
#    col - a string giving the color of the boundary
#    ylim - a size 2 vector giving the limits of the y axis
#    main - a string giving the title of the plot. May present latex language.
#    add - a boolean, if the generated plot should be added over an already existing one.
  
  tpp_seq = seq(1e-3,0.999,length.out=1000)
  fdp_seq = sapply(tpp_seq,fdrlasso,delta=delta,epsi=epsilon)

  if(mode=="empty")
  {
    
    if(add==TRUE) {par(new=TRUE)}
    
    plot(tpp_seq,fdp_seq,type="l",lwd=2,xlim=c(0,1),ylim=ylim,
         xlab="TPP",ylab="FDP",col=col,main=TeX(eval(main)))
  }
  
  else if(mode=="full")
  {
    plot(tpp_seq,fdp_seq,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),
         xlab="TPP",ylab="FDP")
    polygon(c(0,tpp_seq,1),c(0,fdp_seq,0),col='skyblue')
  }
  
  else{
    stop("Unknown mode !")
  }
}



