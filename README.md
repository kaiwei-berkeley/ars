# Adeptive Rejection Sampling
## Description
Adaptive Rejection Sampling From Log-concave Density Functions h(x)
## Usage
ars(h, start, end, N, k = 3, x1 = NULL, xk = NULL)
## Arguments
*h*	
the original function we want to sample from

*start*	
lower bound of the domain of h(x)

*end*	
upper bound of the domain of h(x)

*N*	
sample size

*k*	
number of starting points, the default is 3

*x1*	
the right starting point, if NULL, the function will find one

*xk*	
the left starting point, if NULL, the function will find one
