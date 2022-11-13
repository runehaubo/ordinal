#############################################################################
##    Copyright (c) 2010-2022 Rune Haubo Bojesen Christensen
##
##    This file is part of the ordinal package for R (*ordinal*)
##
##    *ordinal* is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 2 of the License, or
##    (at your option) any later version.
##
##    *ordinal* is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    A copy of the GNU General Public License is available at
##    <https://www.r-project.org/Licenses/> and/or
##    <http://www.gnu.org/licenses/>.
#############################################################################
## This file contains:
## Gradients of densities of common distribution functions on the form
## g[dist], where "dist" can be one of "logis", "norm", and
## "cauchy". These functions are used in Newton-Raphson algorithms
## when fitting CLMs and CLMMs in clm(), clm2(), clmm() and
## clmm2(). Similar gradients are implemented for the gumbel,
## log-gamma, and Aranda-Ordaz distributions.

glogis <- function(x)
### gradient of dlogis
    .C("glogis_C",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gnorm <- function(x)
### gradient of dnorm(x) wrt. x
    .C("gnorm_C",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x

gcauchy <- function(x)
### gradient of dcauchy(x) wrt. x
    .C("gcauchy_C",
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
