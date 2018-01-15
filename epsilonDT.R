epsilonDT <- function(delta)
{
  minus_f <- function(x)
  {
    value = -(1+2/delta*x*dnorm(x) - 2/delta*(1+x^2)*pnorm(-x))/(1+x^2 -
          2*(1+x^2)*pnorm(-x)+2*x*dnorm(x))*delta
    
    return(value)
  }
  
  alpha_phase = optimize(minus_f,c(0,8))
  
  epsi = -alpha_phase$objective
  
  return(epsi)
}