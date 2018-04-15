## This file contains:
## Gradients of densities of common distribution functions on the form
## g[dist], where "dist" can be one of "logis", "norm", and
## "cauchy". These functions are used in Newton-Raphson algorithms
## when fitting CLMs and CLMMs in clm(), clm2(), clmm() and
## clmm2(). Similar gradients are implemented for the gumbel,
## log-gamma, and Aranda-Ordaz distributions.

glogis <- function(x)
### gradient of dlogis
    .C("glogis",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gnorm <- function(x)
### gradient of dnorm(x) wrt. x
    .C("gnorm",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gcauchy <- function(x)
### gradient of dcauchy(x) wrt. x
    .C("gcauchy",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

glogisR <- function(x) {
### glogis in R
  res <- rep(0, length(x))
  isFinite <- !is.infinite(x)

  x <- x[isFinite]
  isNegative <- x < 0
  q <- exp(-abs(x))
  q <- 2*q^2*(1 + q)^-3 - q*(1 + q)^-2
  q[isNegative] <- -q[isNegative]
  res[isFinite] <- q
  res
}

gnormR <- function(x)
### gnorm in R
    -x * dnorm(x)

gcauchyR <- function(x)
### gcauchy(x) in R
    -2*x/pi*(1+x^2)^-2
