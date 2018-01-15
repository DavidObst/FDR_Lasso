rsandwich <- function(t,tpp)
{
  R = (1-tpp)/(1-2*pnorm(-t))
  return(R)
}