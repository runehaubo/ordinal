#############################################################################
#    Copyright (c) 2010-2018 Rune Haubo Bojesen Christensen
#
#    This file is part of the ordinal package for R (*ordinal*)
#
#    *ordinal* is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    *ordinal* is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    <https://www.r-project.org/Licenses/> and/or
#    <http://www.gnu.org/licenses/>.
#############################################################################
## This file contains:
## [pdg]lgamma functions for the log-gamma distribution [lgamma].
## Here glgamma is the gradient of the density function, dlgamma.
## The log-gamma distribution is
## used as a flexible link function in clm2() and clmm2().

plgamma <- function(q, lambda, lower.tail = TRUE)
    .C("plgamma_C",
       q = as.double(q),
       length(q),
       as.double(lambda[1]),
       as.integer(lower.tail[1]),
       NAOK = TRUE)$q

plgammaR <- function(eta, lambda, lower.tail = TRUE) {
    q <- lambda
    v <- q^(-2) * exp(q * eta)
    if(q < 0)
        p <- 1 - pgamma(v, q^(-2))
    if(q > 0)
        p <- pgamma(v, q^(-2))
    if(isTRUE(all.equal(0, q, tolerance = 1e-6)))
        p <- pnorm(eta)
    if(!lower.tail) 1 - p else p
}

dlgamma <- function(x, lambda, log = FALSE) {
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  .C("dlgamma_C",
     x = as.double(x),
     length(x),
     as.double(lambda),
     as.integer(log),
     NAOK = TRUE)$x
}

dlgammaR <- function(x, lambda, log = FALSE) {
    q <- lambda
    q.2 <- q^(-2)
    qx <- q * x
    log.d <- log(abs(q)) + q.2 * log(q.2) -
        lgamma(q.2) + q.2 * (qx - exp(qx))
    if (!log) exp(log.d) else log.d
}

glgamma <- function(x, lambda) {
  stopifnot(length(lambda) == 1)
  .C("glgamma_C",
     x = as.double(x),
     length(x),
     as.double(lambda[1]),
     NAOK = TRUE)$x
}

glgammaR <- function(x, lambda) {
  stopifnot(length(lambda) == 1)
  (1 - exp(lambda * x))/lambda * dlgamma(x, lambda)
}

glgammaR2 <- function(x, lambda) {
  stopifnot(length(lambda == 1))
  if(lambda == 0)
    return(gnorm(x))
  y <- dlgamma(x, lambda)
  y[!is.na(y) && y > 0] <- y * (1 - exp(lambda * x))
  return(y)
}

