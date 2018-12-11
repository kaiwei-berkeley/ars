# ARS
\name{ars}
\alias{ars}
\title{Adeptive Rejection Sampling}
\usage{
ars(h, start, end, N, k = 3, x1 = NULL, xk = NULL)
}
\arguments{
\item{h}{the original function we want to sample from}

\item{start}{lower bound of the domain of h(x)}

\item{end}{upper bound of the domain of h(x)}

\item{N}{sample size}

\item{k}{number of starting points, the default is 3}

\item{x1}{the right starting point, if NULL, the function will find one}

\item{xk}{the left starting point, if NULL, the function will find one}
}
\value{
a vector of N sampled value from the density h(x)
}
\description{
Adaptive Rejection Sampling from log-concave density functions h(x)
}
\examples{
library(ars)
h = function(x){dnorm(x)}
sample = ars(h = h,start = -Inf , end = Inf,N = 100)
hist(sample)
}
