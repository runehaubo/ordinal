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
## [pdg]AO functions for the Aranda-Ordaz distribution. Here gAO is
## the gradient of the density function, dAO. The AO distribution is
## used as a flexible link function in clm2() and clmm2().

pAOR <- function(q, lambda, lower.tail = TRUE) {
    if(lambda < 1e-6)
        stop("'lambda' has to be positive. lambda = ", lambda, " was supplied")
    p <- 1 - (lambda * exp(q) + 1)^(-1/lambda)
    if(!lower.tail) 1 - p else p
}

pAO <- function(q, lambda, lower.tail = TRUE)
    .C("pAO_C",
       q = as.double(q),
       length(q),
       as.double(lambda[1]),
       as.integer(lower.tail),
       NAOK = TRUE)$q

dAOR <- function(eta, lambda, log = FALSE) {
### exp(eta) * (lambda * exp(eta) + 1)^(-1-1/lambda)
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  if(lambda < 1e-6)
    stop("'lambda' has to be positive. lambda = ", lambda,
         " was supplied")
  log.d <- eta - (1 + 1/lambda) * log(lambda * exp(eta) + 1)
  if(!log) exp(log.d) else log.d
}

dAO <- function(eta, lambda, log = FALSE) {
  stopifnot(length(lambda) == 1 &&
            length(log) == 1)
  .C("dAO_C",
     eta = as.double(eta),
     length(eta),
     as.double(lambda),
     as.integer(log),
     NAOK = TRUE)$eta
}

gAOR <- function(eta, lambda) {
  stopifnot(length(lambda) == 1)
  lex <- lambda * exp(eta)
  dAO(eta, lambda) * (1 - (1 + 1/lambda) * lex/(1 + lex))
}

gAO <- function(eta, lambda) {
  stopifnot(length(lambda) == 1)
  .C("gAO_C",
     eta = as.double(eta),
     length(eta),
     as.double(lambda[1]),
     NAOK = TRUE)$eta
}

