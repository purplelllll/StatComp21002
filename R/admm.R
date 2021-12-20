#' @title admm to estimate some parameters
#' @description estimate some parameters in a regression model by the admm method
#' @param beta0 one parameter in a regression model
#' @param beta1 one parameter in a regression model
#' @param beta2 one parameter in a regression model
#' @param m the number of the sample
#' @param n the observation number for a sample
#' @return a table and a number
#' @examples
#' \dontrun{
#' admm(1.5,1,0,20,5)
#' } 
#' @import knitr
#' @importFrom stats rnorm
#' @export
admm=function(beta0,beta1,beta2,m,n)
{
  set.seed(78910)
  ##generate the sample
  X1=matrix(0,m,n)
  for (i in 1:m)
    for (j in 1:n)
      X1[i,j]=rnorm(1)
  
  X2=matrix(0,m,n)
  for (i in 1:m)
    for ( j in 1:n)
      X2[i,j]=rnorm(1)
  
  sigma=matrix(0,m,n)
  for (i in 1:m)
    for (j in 1:n)
      sigma[i,j]=rnorm(1,0,sd=0.8)
  
  #compute the yij
  y=matrix(0,m,n)
  for (i in 1:m)
    for (j in 1:n)
      y[i,j]=beta0+X1[i,j]*beta1+X2[i,j]*beta2+sigma[i,j]
  
  #set the initial value of the parameters
  beta0_=beta0
  beta1_=beta1
  beta2_=beta2
  
  theta_=matrix(0,m,n)
  for (i in 1:m)
    for ( j in 1:n)
      theta_[i,j]=rnorm(1)
  
  xi=1
  tau=0.5
  
  ##Aij Zij
  A=matrix(0,m,n)
  Z=matrix(0,m,n)
  
  N=1000
  tol <-1e-3
  r=0
  for (l in 1:N){
    for (i in 1:m)
    {
      for (j in 1:n)
      {
        A[i,j]=y[i,j]-(beta0_+X1[i,j]*beta1_+X2[i,j]*beta2_)
        if ((A[i,j]-theta_[i,j]/xi)>tau/xi)
          Z[i,j]=A[i,j]-tau/xi-theta_[i,j]/xi
        if ((A[i,j]-theta_[i,j]/xi)<=(tau-1)/xi)
          Z[i,j]=A[i,j]+(1-tau)/xi-theta_[i,j]/xi 
      }
    }
    
    ##iteration
    s1=sum(Z-y+X1*beta1_+X2*beta2_)
    s=sum(theta_)
    beta0_1=-(s1*xi+s)/(xi*m*n)
    e1=abs(beta0_1-beta0_)
    beta0_=beta0_1
    s2=sum((Z-y+beta0_+X2*beta2_)*X1)
    s20=sum(theta_*X1)
    beta1_1=-(s2*xi+s20)/(xi*sum(X1^2))
    e2=abs(beta1_1-beta1_)
    beta1_=beta1_1
    s3=sum((Z-y+beta0_+X1*beta1_)*X2)
    s30=sum(theta_*X2)
    beta2_1=-(s3*xi+s30)/(xi*sum(X2^2))
    e3=abs(beta2_1-beta2_)
    beta2_=beta2_1
    theta_1=theta_+xi*(Z-(y-(beta0_+X1*beta1_+X2*beta2_)))
    theta_=theta_1
    
    r=r+1
    if ((e1<tol)&&(e2<tol)&&(e3<tol))  break  
  }
  parameters=c(beta0,beta1,beta2)
  estimators=c(beta0_,beta1_,beta2_)
  r1=data.frame(parameters,estimators)
  
  return( list( (knitr::kable(t(r1),align="c")),r ) )
}

