#' @title A  dataset
#' @name data
#' @description Data sets for presentation only
#' @examples
#' \dontrun{
#' data(data)
#' attach(data)
#' print(summary(s1)
#' }
NULL



#'@title Compute skewness coeff
#'@name skewness
#'@description compute the skewness coeff of a set of data
#'@param x A set of data
#'@return skewness coeff 
#'@examples
#'\dontrun{p.reject <- numeric(length(n)) 
#'  m <- 10000
#'  for (i in 1:length(n)) { sktests <- numeric(m)
#'  for (j in 1:m) { x <- rnorm(n[i])
#'  sktests[j] <- as.integer(abs(skewness(x)) >= cv[i] ) 
#'  }
#'  p.reject[i] <- mean(sktests)}
#'}
#'@import microbenchmark
#'@importFrom Rcpp evalCpp
#'@importFrom stats rnorm runif var var.test 
#'@export
#'@useDynLib StatComp20072
skewness<-function(x){
  xbar<-mean(x)
  m3<-mean((x-xbar)^3)
  m2<-mean((x-xbar)^2)
  return(m3/m2^1.5)
}


#'@title Hypothesis Test with CountFive
#'@name count5test
#'@description Use CountFive method to decide whether or not to accept the null hypothesis
#'@param x A set of data
#'@param y A set of data that is different from \code{x}
#'@return 0 or 1
#'@examples
#'\dontrun{
#'sigma1<-1 
#'sigma2<-1.5
#'  m<-1000
#'  n<-20
#'  power<-numeric(m)
#'  for(j in 1:length(n)){
#'      for(i in 1:m){
#'    p.count5test<-mean(replicate(1000,expr={
#'         x<-rnorm(n[j],0,sigma1)
#'          y<-rnorm(n[j],0,sigma2)
#'          count5test(x,y)
#'      }))
#'          }
#'      power[j]<-p.count5test
#'    }
#'    power
#'}
#'@export
count5test <- function(x, y) {
      X <- x - mean(x) 
      Y <- y - mean(y) 
      outx <- sum(X > max(Y)) + sum(X < min(Y)) 
      outy <- sum(Y > max(X)) + sum(Y < min(X))
      return(as.integer(max(c(outx, outy)) > 5))
}



#'@title F Test to Compare Two Variances
#'@name Ftest
#'@description F Test to Compare Two Variances and compute the P value
#'@param x A set of data
#'@param y A set of data that is different from \code{x}
#'@return p-value of F Test to Compare Two Variances
#'@examples
#'\dontrun{
#'sigma1 <- 1 
#'sigma2 <- 1.5
#'m<-1000
#'n<-c(10,50,100)
#' power<-numeric(m)
#'for(j in 1:length(n)){
#'      for(i in 1:m){
#'       x<-rnorm(n[j],0,sigma1)
#'        y<-rnorm(n[j],0,sigma2)
#'        p.Ftest<-Ftest(x,y)
#'      }
#'     power[j]<-p.Ftest
#'    }
#'    power

#'}
#'@export
Ftest<-function(x,y){
  P<-var.test(x,y,ratio = 1,alternative = c("two.sided","less","greater"),
              conf.level = 0.945)
  return(P$p.value)
}  


#'@title Number of extreme points
#'@name maxout
#'@description Counts the number of extreme points of each sample relative to the range of the other sample.
#'@param x A set of data
#'@param y A set of data compared to \code{x}
#'@return the number of extreme points
#'@examples
#'\dontrun{
#'set.seed(124)
#'x<-rnorm(20,4,2)
#'y<-rnorm(10,4,2)
#'test<-maxout(x,y)
#'}
#'@export
maxout<-function(x,y){
      X <- x - mean(x) 
      Y <- y - mean(y) 
      outx <- sum(X > max(Y)) + sum(X < min(Y))
      outy <- sum(Y > max(X)) + sum(Y < min(X)) 
      return(max(c(outx, outy)))
}




#'@title Generate the standard Laplace distribution with a normal increment
#'@name MetropolisR
#'@description Generate the standard Laplace distribution with a normal increment
#'@param sigma the standard error
#'@param x0 the intial value
#'@param N the simulation times
#'@return a random sample of size \code{n}
#'@examples
#'\dontrun{
#'N<-3000
#'sigma<-1
#'x0<-25
#'MR1<-MetropolisR(sigma,x0,N)
#'}
#'@export
MetropolisR<-function(sigma,x0,N){
  x<-numeric(N)
  x[1]<-x0
  u<-runif(N)
  k<-0
  for(i in 2:N){
    y<-rnorm(1,x[i-1],sigma)
    if(u[i]<=(exp(-abs(y))/exp(-abs(x[i-1]))))
      x[i]<-y else{
        x[i]<-x[i-1]
        k<-k+1
      }
  }
  return(list(x=x,k=k))
}




#'@title G-R statistic 
#'@name Rubin
#'@description Compute G-R statistic
#'@param psi A matrix with n rows and 2 columns
#'@return G-R statistic
#'@examples
#'\dontrun{
#'sigma1<-0.05
#'N<-20000
#'b<-1000
#'K<-4
#'x0<-c(-8,-4,4,8)
#'X<-matrix(0,K,N)
#'for(i in 1:K)
#'  X[i,]<-Metropolis(sigma1,x0[i],N)
#'psi1<-t(apply(X,1,cumsum))
#'for(i in 1:nrow(psi1))
#'  psi1[i,]<-psi1[i,]/(1:ncol(psi1))
#'print(Rubin(psi1))
#'}
#'@export
Rubin<-function(psi){
  psi<-as.matrix(psi)
  n<-ncol(psi)
  k<-nrow(psi)
  psi.means<-rowMeans(psi)
  B<-n*var(psi.means)
  psi.w<-apply(psi,1,"var")
  W<-mean(psi.w)
  v.hat<-W*(n-1)/n+(B/n)
  r.hat<-v.hat/W
  return(r.hat)
}



#'@title Observed data likelihood
#'@name observed
#'@description Observed data likelihood
#'@param x initial value of p and q
#'@return data likelihood
#'@examples
#'\dontrun{
#'theta<- optim(c(0.35,0.2),fn=observed)$par
#'}
#'@export
observed<- function(x){
  p <- x[1]
  q <- x[2]
  r <- 1-p-q
  f <- 444*log(p^2+2*p*r) + 132*log(q^2+2*q*r) + 2*361*log(r) + 63*(log(2)+log(p)+log(q))
  return(-f)
}

#'@title Simulations in VAR systems
#'@name rSim
#'@description Simulations in VAR systems(Explicit loop)
#'@param coeff the coefficient matrix
#'@param errors error terms
#'@return Data matrix
#'@examples
#'\dontrun{
#'a<-matrix(c(0.5,0.1,0.1,0.5),nrow=2)
#'u<-matrix(rnorm(10000),ncol=2)
#'rData<-rSim(a,u)
#'}
#'@export
rSim<-function(coeff,errors){
  simdata<-matrix(0,nrow(errors),ncol(errors))
 for (i in 2:nrow(errors)){
  simdata[i,]=coeff %*% simdata[(i-1),]+errors[i,]
  }
  return(simdata)
}



#'@title Mixed normal distribution
#'@name Mnorm
#'@description Generate a mixed normal distribution by weight
#'@param p The weight vector
#'@param N The number of random numbers
#'@param mu0 The mean vector
#'@param sigma0 The standard error vector
#'@return a random sample of size \code{n}
#'@examples 
#'\dontrun{
#'p<-c(0.2,0.2,0.1,0.5)
#'sigma0<-c(1,2,4,5)
#'mu0<-c(1,3,4,2)
#'N<-20
#'Mnorm(p,N,mu0,sigma0)
#'}
#'@export
Mnorm<-function(p,N,mu0,sigma0){
  sigma1<-sample(sigma0,size = N,replace = T,prob = p)
  mu1<-sample(mu0,size = N,replace = T,prob = p)
  x<-rnorm(N,mu1,sigma1) 
  return(x)
}

