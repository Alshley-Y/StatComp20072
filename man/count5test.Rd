% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StatComp20072R.R
\name{count5test}
\alias{count5test}
\title{Hypothesis Test with CountFive}
\usage{
count5test(x, y)
}
\arguments{
\item{x}{A set of data}

\item{y}{A set of data that is different from \code{x}}
}
\value{
0 or 1
}
\description{
Use CountFive method to decide whether or not to accept the null hypothesis
}
\examples{
\dontrun{
sigma1<-1 
sigma2<-1.5
 m<-1000
 n<-20
 power<-numeric(m)
 for(j in 1:length(n)){
     for(i in 1:m){
   p.count5test<-mean(replicate(1000,expr={
        x<-rnorm(n[j],0,sigma1)
         y<-rnorm(n[j],0,sigma2)
         count5test(x,y)
     }))
         }
     power[j]<-p.count5test
   }
   power
}
}
