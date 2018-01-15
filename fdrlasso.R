fdrlasso <- function(tpp,delta,epsi)
{
#--------------------------------------------------------------------------
#This function calculates the Lasso trade-off curve given tpp (true
# positive proportion), delta = n/p (shape of the design matrix, or
# subsampling rate), and epsi = k/p (sparsity ratio).
# All tpp, delta, and epsi are between 0 and 1; if the
# pair (delta, epsi) is above the Donoho-Tanner phase transition, tpp
# should be no larger than u^\star = powermax(delta, epsi)
#  
# --------------------------------------------------------------------------
#  Copyright @ Weijie Su, Malgorzata Bogdan, and Emmanuel Candes, 2015
# Translated from Matlab to R by M. Guglielmino and D. Obst (ENSTA ParisTech)
# --------------------------------------------------------------------------
   
  if(tpp > powermax(delta,epsi))
  {
    stop("Invalid input !")
  }
  
  if(tpp == 0)
  {
    q = 0
  }
  
  ## Make stepsize smaller for higher acc
  stepsize = 0.1
  tmax = max(10,sqrt(delta/epsi/tpp)+1)
  tmin = tmax - stepsize
  
  while(tmin>0)
  {
    if(lsandwich(tmin,tpp,delta,epsi)<rsandwich(tmin,tpp))
    {
      break
    }
    
    tmax = tmin 
    tmin = tmax - stepsize
  }
  
  if(tmin <= 0)
  {
    stepsize = stepsize/100
    tmax = max(10,sqrt(delta/epsi/tpp)+1)
    tmin = tmax - stepsize
    
    while(tmin>0)
    {
      if(lsandwich(tmin,tpp,delta,epsi)<rsandwich(min,tpp))
      {
        break
      }
      
      tmax = tmin
      tmin = tmax - stepsize
    }
  }
  
  diff = tmax - tmin
  
  while(diff>1e-6)
  {
    tmid = 0.5*tmax + 0.5*tmin
    
    if(lsandwich(tmid,tpp,delta,epsi)>rsandwich(tmid,tpp))
    {
      tmax = tmid
    }
    
    else
    {
      tmin = tmid
    }
    
    diff = tmax - tmin
  }
  
  t = (tmax+tmin)/2
  
  q = 2*(1-epsi)*pnorm(-t)/(2*(1-epsi)*pnorm(-t)+epsi*tpp)
  
  return(q)
                                                               
}