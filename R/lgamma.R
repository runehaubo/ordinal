## This file contains:
## [pdg]lgamma functions for the log-gamma distribution [lgamma].
## Here glgamma is the gradient of the density function, dlgamma.
## The log-gamma distribution is
## used as a flexible link function in clm2() and clmm2().

plgamma <- function(q, lambda, lower.tail = TRUE)
    .C("plgamma",
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
  .C("dlgamma",
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
  .C("glgamma",
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

