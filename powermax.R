powermax <- function(delta,epsilon)
{
  if(delta>=1)
  {
    power = 1
    return(power)
  }
  
  epsilon_star = epsilonDT(delta)
  
  if(epsilon <= epsilon_star)
  {
    power = 1
    return(power)
  }
  
  power = (epsilon-epsilon_star)*(delta-epsilon_star)/epsilon/(1-epsilon_star) + epsilon_star/epsilon
  
  return(power)
}