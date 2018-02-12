
plot_qstar <- function(delta,epsilon,mode="empty",col="blue",ylim=c(0,1),main=NULL,add=FALSE,add2=FALSE)
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
#    add2 - a boolen, if the generated plot should be added as a line plot to an existing one
# Genrally leave add = FALSE and if you want to add set add2 = TRUE
  
  tpp_seq = seq(1e-3,0.999,length.out=1000)
  fdp_seq = sapply(tpp_seq,fdrlasso,delta=delta,epsi=epsilon)

  if(mode=="empty")
  {
    if(add==FALSE & add2==FALSE)
    {
      plot(tpp_seq,fdp_seq,type="l",lwd=2,xlim=c(0,1),ylim=ylim,
           xlab="TPP",ylab="FDP",col=col,main=main)
    }
    
    if(add==TRUE) {
      par(new=TRUE)
      plot(tpp_seq,fdp_seq,type="l",lwd=2,xlim=c(0,1),ylim=ylim,
           xlab="TPP",ylab="FDP",col=col,main=main)
      
    }
    else if(add2==TRUE) {
      lines(tpp_seq,fdp_seq,lwd=2,col=col)
    }
    
    
  }
  
  else if(mode=="full")
  {
    plot(tpp_seq,fdp_seq,type="l",lwd=2,xlim=c(0,1),ylim=c(0,1),
         xlab="TPP",ylab="FDP",main=main)
    polygon(c(0,tpp_seq,1),c(0,fdp_seq,0),col='skyblue')
  }
  
  else{
    stop("Unknown mode !")
  }
}



