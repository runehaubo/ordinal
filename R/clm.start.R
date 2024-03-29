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
## Functions to compute starting values for CLMs in clm().

set.start <-
  function(rho, start=NULL, get.start=TRUE, threshold, link, frames)
{
  ## set starting values for the parameters:
  nScol <- if(is.null(frames[["S"]])) 0 else ncol(frames[["S"]]) # no cols in S 
  nSpar <- pmax(0, nScol - 1) # no Scale parameters
  if(get.start) {
    start <- ## not 'starting' scale effects:
        clm.start(y.levels=frames$y.levels, threshold=threshold, X=frames$X,
                  NOM=frames$NOM, has.intercept=TRUE)
    if(nSpar > 0 || # NCOL(frames[["S"]]) > 1 
       link == "cauchit" || length(rho$lambda)) {
### NOTE: only special start if nSpar > 0 (no reason for
### special start if scale is only offset and no predictors).
### NOTE: start cauchit models at the probit estimates if start is not
### supplied:
### NOTE: start models with lambda at model with probit link
      rho$par <- start
      if(link %in% c("Aranda-Ordaz", "log-gamma", "cauchit")) {
        setLinks(rho, link="probit")
      } else { 
        setLinks(rho, link)
      }
      tempk <- rho$k
      rho$k <- 0
      ## increased gradTol and relTol:
      fit <- try(clm_fit_NR(rho, control=list(gradTol=1e-3, relTol=1e-3)),
                 silent=TRUE)
      if(inherits(fit, "try-error"))
        stop("Failed to find suitable starting values: please supply some",
             call.=FALSE)
      start <- c(fit$par, rep(0, nSpar))
      if(length(rho$lambda) > 0) start <- c(start, rho$lambda)
      attr(start, "start.iter") <- fit$niter
      rho$k <- tempk
      setLinks(rho, link) # reset link in rho
    }
  }
  ## test start:
  stopifnot(is.numeric(start))
  length.start <- ncol(rho$B1) + nSpar + length(rho$lambda) 
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
  return(rep(0, ncol(X) - has.intercept))

## clm.start <- function(y.levels, threshold, X, has.intercept = TRUE)
##   return(c(start.threshold(y.levels, threshold),
##            start.beta(X, has.intercept)))

clm.start <- function(y.levels, threshold, X, NOM=NULL, S=NULL,
                       has.intercept=TRUE)
{
  st <- start.threshold(y.levels, threshold)
  if(!is.null(NOM) && ncol(NOM) > 1)
    st <- c(st, rep(rep(0, length(st)), ncol(NOM)-1))
  start <- c(st, start.beta(X, has.intercept))
  if(!is.null(S) && ncol(S) > 1)
    start <- c(start, rep(0, ncol(S) - 1))
  start
}

