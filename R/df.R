#' @title compute the df of cauchy distribution
#' @description using MC integration
#' @param x which the df take values in
#' @param theta the scale parameter
#' @param eta the location parameter
#' @return a number
#' @examples
#' \dontrun{
#' df(1,1,0)
#' } 
#' @importFrom stats integrate pcauchy
#' @export
df=function(x,theta,eta)
{
  set.seed(34567)
  ## the density function
  f0=function(x,theta,eta)
  {
    a=((x-eta)/theta)^2
    b=theta*pi*(1+a)
    c=1/b
    return (c)
  }
  ## compute one part of the df
  s1=integrate(f0,lower=-Inf,upper=0,rel.tol = .Machine$double.eps^0.25,theta=theta,eta=eta)$value
  ## compute another part of the df
  u=runif(1e6)
  s2=mean(x/(theta*pi*(1+((u*x-eta)/theta)^2)))
  
  ##the df
  if (x>0)  
    s=s1+s2
  else
    s=integrate(f0,lower=-Inf,upper=x,rel.tol = .Machine$double.eps^0.25,theta=theta,eta=eta)$value
  return(c(s,pcauchy(x,eta,theta)))
}