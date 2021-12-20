#include <Rcpp.h>
using namespace Rcpp;
//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param n the number of samples
//' @param a related to the distribution
//' @param b related to the distribution
//' @return a random matrix
//' @examples
//' \dontrun{
//' gibbsC1(1,2,10) 
//' } 
//' @importFrom stats rbinom rbeta
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC1(int a, int b, int n) {
  
  int N=2000; 
  NumericMatrix X(N, 2);

  double x1=0;
  double x2=0.6;
  double m1=0;
  double m2=0;
  for (int i=1 ; i<N ; i++) {
    
    x1 = rbinom(1, n, x2)[0];
    m1=x1+a;
    m2=n-x1+b;
    x2 = rbeta(1, m1, m2)[0];
    X(i,0)=x1;
    X(i,1)=x2;
  }
 
    return (X);
}