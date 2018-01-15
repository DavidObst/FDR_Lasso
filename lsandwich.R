lsandwich <- function(t,tpp,delta,epsi)
{
  Lnume = (1-epsi)*(2*(1+t^2)*pnorm(-t)-2*t*dnorm(t))+epsi*(1+t^2) - delta
  Ldeno = epsi*((1+t^2)*(1-2*pnorm(-t))+ 2*t*dnorm(t))
  
  L = Lnume/Ldeno
  
  return(L)
}