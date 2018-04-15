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
    .C("pAO",
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
  .C("dAO",
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
  .C("gAO",
     eta = as.double(eta),
     length(eta),
     as.double(lambda[1]),
     NAOK = TRUE)$eta
}

