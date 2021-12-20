#' @title A esitimator for a parameter in a distribution
#' @description estimate a parameter and compute its sd and bias
#' @param sigma the parameter to estimate
#' @param n the number of samples
#' @return a random vector
#' @examples
#' \dontrun{
#' est(1,1e4)
#' } 
#' @importFrom stats optimize runif sd
#' @export
est=function(sigma,n)
{
  ##generate random sample from the distribution using inverse transform
  u <- runif(n)
  x <- sqrt((-2)*((sigma)^2)*log(1-u))
  
  ##bootstrap
  B=10
  sigma_=numeric(B)
  for (b in 1:B)
  {
    i=sample(x,size=n,replace=TRUE)
    s1=sum(log(i))
    s2=sum(i^2)
    
    ##the function that needs to optimize
    f0=function(sigma)
      s1-2*n*log(sigma)-s2/(2*(sigma^2))
    ##the MLE
    sigma_[b]=optimize(f0,lower = 0,upper=sigma+4,maximum = TRUE)$maximum
  } 
  estimator=mean(sigma_)
  ##the bias and se of bootstrap
  bias=mean(sigma_)-sigma
  se=sd(sigma_)
  val=numeric(3)
  val=c(estimator,bias,se)
  return(val)
}