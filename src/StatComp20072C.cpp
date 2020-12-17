#include <Rcpp.h>
using namespace Rcpp;

//' @title Generate the standard Laplace distribution with a normal increment using Rcpp
//' @description  Generate the standard Laplace distribution with a normal increment using Rcpp
//' @param N the simulation times
//' @param sigma the standard error
//' @param x0 the intial value
//' @return a random samples \code{n}
//' @examples
//' \dontrun{
//' sigma<-1
//' x0<-25
//' N<-3000
//' MC<-MetropolisC(sigma,x0,N)
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix MetropolisC(double sigma,double x0,int N) {
  NumericMatrix mat(N, 2);
  mat(0,0)=x0;
  mat(0,1)=0;
  double y=0,u=0;
  for(int i = 2; i < N+1; i++) {
    y=rnorm(1,mat(i-2,0),sigma)[0];
    u=runif(1,0,1)[0];
    if(u<=exp(-abs(y))/exp(-abs(mat(i-2,0)))){
      mat(i-1,0)=y;
      mat(i-1,1)=mat(i-2,1);
    }
    else{
      mat(i-1,0)=mat(i-2,0);
      mat(i-1,1)=mat(i-2,1)+1;
    }
  }
  return(mat);
}
