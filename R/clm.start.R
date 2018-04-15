## This file contains:
## Functions to compute starting values for CLMs in clm().

set.start <-
  function(rho, start=NULL, get.start=TRUE, threshold, link, frames)
{
  ## set starting values for the parameters:
  if(get.start) {
    start <- ## not 'starting' scale effects:
        clm.start(y.levels=frames$y.levels, threshold=threshold, X=frames$X,
                  NOM=frames$NOM, has.intercept=TRUE)
    if(NCOL(frames$S) > 1 || link == "cauchit") {
### NOTE: only special start if NCOL(frames$S) > 1 (no reason for
### special start if scale is only offset and no predictors).
### NOTE: start cauchit models at the probit estimates if start is not
### supplied:
      rho$par <- start
      if(link == "cauchit") setLinks(rho, link="probit")
      else setLinks(rho, link)
      tempk <- rho$k
      rho$k <- 0
      ## increased gradTol and relTol:
      fit <- try(clm.fit.NR(rho, control=list(gradTol=1e-3, relTol=1e-3)),
                 silent=TRUE)
      if(class(fit) == "try-error")
        stop("Failed to find suitable starting values: please supply some",
             call.=FALSE)
      start <- c(fit$par, rep(0, NCOL(frames$S) - 1))
      attr(start, "start.iter") <- fit$niter
      rho$k <- tempk
    }
  }
  ## test start:
  stopifnot(is.numeric(start))
  length.start <- ncol(rho$B1) + NCOL(frames$S) - 1 #- length(rho$alised)
  if(length(start) != length.start)
    stop(gettextf("length of start is %d should equal %d",
                  length(start), length.start), call.=FALSE)

  return(start)
}

start.threshold <-
  function(y.levels, threshold = c("flexible", "symmetric", "symmetric2", "equidistant"))
### args:
### y.levels - levels of the model response, at least of length two
### threshold - threshold structure, character.
{
  ## match and test arguments:
  threshold <- match.arg(threshold)
  ny.levels <- length(y.levels)
  ntheta <- ny.levels - 1L
  if(threshold %in% c("symmetric", "symmetric2", "equidistant") && ny.levels < 3)
    stop(gettextf("symmetric and equidistant thresholds are only
meaningful for responses with 3 or more levels"))

  ## default starting values:
  start <- qlogis((1:ntheta) / (ntheta + 1) ) # just a guess

  ## adjusting for threshold functions:
  if(threshold == "symmetric" && ntheta %% 2) { ## ntheta odd >= 3
    nalpha <- (ntheta + 1) / 2
    start <- c(start[nalpha], diff(start[nalpha:ntheta])) ## works for
    ## ntheta >= 1
  }
  if(threshold == "symmetric" && !ntheta %% 2) {## ntheta even >= 4
    nalpha <- (ntheta + 2) / 2
    start <- c(start[c(nalpha - 1, nalpha)],
               diff(start[nalpha:ntheta])) ## works for ntheta >= 2
  }
  if(threshold == "symmetric2" && ntheta %% 2) { ## ntheta odd >= 3
    nalpha <- (ntheta + 3) / 2
    start <- start[nalpha:ntheta] ## works for ntheta >= 3
  }
  if(threshold == "symmetric2" && !ntheta %% 2) {## ntheta even >= 4
    nalpha <- (ntheta + 2) / 2
    start <- start[nalpha:ntheta] ## works for ntheta >= 2
  }
  if(threshold == "equidistant")
    start <- c(start[1], mean(diff(start)))

  ## return starting values for the threshold parameters:
  return(as.vector(start))
}

start.beta <- function(X, has.intercept = TRUE)
  return(rep(0, NCOL(X) - has.intercept))

## clm.start <- function(y.levels, threshold, X, has.intercept = TRUE)
##   return(c(start.threshold(y.levels, threshold),
##            start.beta(X, has.intercept)))

clm.start <- function(y.levels, threshold, X, NOM=NULL, S=NULL,
                       has.intercept=TRUE)
{
  st <- start.threshold(y.levels, threshold)
  if(!is.null(NOM) && NCOL(NOM) > 1)
    st <- c(st, rep(rep(0, length(st)), ncol(NOM)-1))
  start <- c(st, start.beta(X, has.intercept))
  if(!is.null(S) && NCOL(S) > 1)
    start <- c(start, rep(0, ncol(S) - 1))
  start
}

