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
## [pdqrg]gumbel functions for the gumbel distribution.
## Here ggumbel is the gradient of the density function, dgumbel.

pgumbel <-
  function(q, location = 0, scale = 1, lower.tail = TRUE, max = TRUE)
### CDF for Gumbel max and min distributions
### Currently only unit length location and scale are supported.
{
  if(max)  ## right skew, loglog link
    .C("pgumbel_C",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail),
       NAOK = TRUE)$q
  else ## left skew, cloglog link
    .C("pgumbel2_C",
       q = as.double(q),
       length(q),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(lower.tail),
       NAOK = TRUE)$q
}

pgumbelR <- function(q, location = 0, scale = 1, lower.tail = TRUE)
### R equivalent of pgumbel()
{
    q <- (q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) 1 - p else p
}

pgumbel2R <- function(q, location = 0, scale = 1, lower.tail = TRUE)
{
    q <- (-q - location)/scale
    p <- exp(-exp(-q))
    if (!lower.tail) p else 1 - p
}

dgumbel <-
  function(x, location = 0, scale = 1, log = FALSE, max = TRUE)
### PDF for the Gumbel max and mon distributions
{
  if(max)  ## right skew, loglog link
    .C("dgumbel_C",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log),
       NAOK = TRUE)$x
  else ## left skew, cloglog link
    .C("dgumbel2_C",
       x = as.double(x),
       length(x),
       as.double(location)[1],
       as.double(scale)[1],
       as.integer(log),
       NAOK = TRUE)$x
}

dgumbelR <- function(x, location = 0, scale = 1, log = FALSE)
### dgumbel in R
{
    q <- (x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

dgumbel2R <- function(x, location = 0, scale = 1, log = FALSE)
{
    q <- (-x - location)/scale
    log.d <- -exp(-q) - q - log(scale)
    if (!log) exp(log.d) else log.d
}

ggumbel <- function(x, max = TRUE) {
### gradient of dgumbel(x) wrt. x
  if(max) ## right skew, loglog link
    .C("ggumbel_C",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x
  else ## left skew, cloglog link
    .C("ggumbel2_C",
       x = as.double(x),
       length(x),
       NAOK = TRUE)$x
}

ggumbelR <- function(x){
### ggumbel in R
  q <- exp(-x)
  ifelse(q == Inf, 0, {
    eq <- exp(-q)
    -eq*q + eq*q*q
  })
}

ggumbel2R <- function(x) -ggumbelR(-x)


rgumbel <- function(n, location = 0, scale = 1, max = TRUE) {
  if(max)
    location - scale * log(-log(runif(n)))
  else
    location + scale * log(-log(runif(n)))
}

qgumbel <- function(p, location = 0, scale = 1, lower.tail = TRUE, max = TRUE) {
  if(!lower.tail) p <- 1 - p
  if(max)  ## right skew, loglog link
    location - scale * log(-log(p))
  else ## left skew, cloglog link
    location + scale * log(-log(1 - p))
}
