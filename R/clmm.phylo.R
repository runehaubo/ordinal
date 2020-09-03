#############################################################################
##    Copyright (c) 2010-2020 Rune Haubo Bojesen Christensen
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
## Implementation of Cumulative Link Mixed Models in clmm().

if(getRversion() >= '2.15.1')
    utils::globalVariables(c("ths", "link", "threshold", "optRes",
                             "neval", "Niter", "tJac", "y.levels"))

clmmPhylo <-
  function(formula, data, weights, start, subset,
           na.action, contrasts, Hess = TRUE, model = TRUE,
           link = c("logit", "probit", "cloglog", "loglog",
             "cauchit"), ##, "Aranda-Ordaz", "log-gamma"), ## lambda,
           doFit = TRUE, control = list(), nAGQ = 1L, cormat=NULL,
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)
{
### Extract the matched call and initial testing:
  mc <- match.call(expand.dots = FALSE)
### OPTION: Possibly call clm() when there are no random effects?
  link <- match.arg(link)
  threshold <- match.arg(threshold)
  if(missing(formula))  stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL
  ## set control parameters:
  control <- getCtrlArgs(control, list(...))
  nAGQ <- as.integer(round(nAGQ))
  formulae <- clmm.formulae(formula=formula)
  ## mf, y, X, wts, off, terms:
  frames <- clmm.frames(modelcall=mc, formulae=formulae, contrasts)
### QUEST: What should 'method="model.frame"' return? Do we want Zt
### included here as well?
  if(control$method == "model.frame") return(frames)
  ## Test rank deficiency and possibly drop some parameters:
  ## X is guarantied to have an intercept at this point.
  frames$X <- drop.coef(frames$X, silent=FALSE)
  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  ths <- makeThresholds(levels(frames$y), threshold)
  ## Set rho environment:
  rho <- with(frames, {
    clm.newRho(parent.frame(), y=y, X=X, weights=wts,
               offset=off, tJac=ths$tJac) })
  ## compute grouping factor list, and Zt and ST matrices:
  retrms <- getREterms(frames=frames, formulae$formula)
  ## For each r.e. term, test if Z has more columns than rows to detect
  ## unidentifiability:
  test_no_ranef(Zt_list=retrms$retrms, frames=frames, checkRanef=control$checkRanef)
### OPTION: save (the evaluated) formula in frames, so we only need the
### frames argument to getREterms() ?
  use.ssr <- (retrms$ssr && !control$useMatrix)
  if(!is.null(cormat)) use.ssr <- FALSE

  ## Set inverse link function and its two first derivatives (pfun,
  ## dfun and gfun) in rho:
  setLinks(rho, link)

  ## Compute list of dimensions for the model fit:
  rho$dims <- getDims(frames=frames, ths=ths, retrms=retrms)

  ## Update model environment with r.e. information:
  if(use.ssr) {
      rho.clm2clmm.ssr(rho=rho, retrms = retrms, ctrl=control$ctrl)
      ## Set starting values for the parameters:
      if(missing(start)) start <- c(fe.start(frames, link, threshold), 0)
      rho$par <- start
      nbeta <- rho$nbeta <- ncol(frames$X) - 1 ## no. fixef parameters
      nalpha <- rho$nalpha <- ths$nalpha ## no. threshold parameters
      ntau <- rho$ntau <- length(retrms$gfList) ## no. variance parameters
      stopifnot(is.numeric(start) &&
                length(start) == (nalpha + nbeta + ntau))
  } else {
      rho.clm2clmm(rho=rho, retrms=retrms, ctrl=control$ctrl)
      if(missing(start)) {
          rho$fepar <- fe.start(frames, link, threshold)
          rho$ST <- STstart(rho$ST)
          start <- c(rho$fepar, ST2par(rho$ST))
      } else {
          stopifnot(is.list(start) && length(start) == 2)
          stopifnot(length(start[[1]]) == rho$dims$nfepar)
          stopifnot(length(start[[2]]) == rho$dims$nSTpar)
          rho$fepar <- as.vector(start[[1]])
          rho$ST <- par2ST(as.vector(start[[2]]), rho$ST)
      }
  }
### OPTION: set starting values in a more elegant way.
  
  #############
  ## Adjust rho$Zt for cholesky factor of known correlation matrix:
  # stopifnot(all(rownames(Zt) == rownames(CORR)))
  if (!is.null(cormat)) {
    if(length(retrms$retrms) > 1)
      stop("Only 1 random-effects term supported with non-null 'cormat'")
    ranef_names <- levels(retrms$gfList[[1]])
    Zt <- getZt(retrms$retrms)
    stopifnot(nrow(cormat) == ncol(cormat),
              nrow(Zt) == nrow(cormat),
              all(colnames(cormat) == rownames(cormat)),
              all(rownames(cormat) == ranef_names))
    chol_cormat <- try(chol(cormat))
    if(inherits(chol_cormat, "try-error")) 
      stop("Computing cholesky factor of cormat failed")
    R <- as(t(chol_cormat), class(Zt))
    rho$Zt <- as(t(crossprod(Zt, R)), class(Zt))
    rho$L <- Cholesky(tcrossprod(crossprod(getLambda(rho$ST, rho$dims$nlev.re), rho$Zt)),
                      LDL = TRUE, super = FALSE, Imult = 1)
  }
  
  # adjust_Zt_for_CORR <- function(retrms, CORR) {
  #   ## does Zt have meaningful rownames?
  #   if(!retrms$ssr) stop("CORR argument is only allowed with a single scalar RE term")
  #   stopifnot(all(rownames(Zt) == rownames(CORR)))
  #   cCORR <- try(chol(CORR))
  #   if(inherits(cCORR, "try-error")) stop("Computing cholesky factor of CORR failed")
  #   R <- as(t(cCORR), class(Zt))
  #   retrms$Zt <- t(crossprod(Zt, R))
  #   retrms
  # }
  # if (!is.null(CORR)) {
  #   retrmtadjust_Zt_for_CORR(retrms, CORR)
  # }
  #############
  ## Set AGQ parameters:
  set.AGQ(rho, nAGQ, control, use.ssr)

  ## Possibly return the environment, rho without fitting:
  if(!doFit)  return(rho)

  ## Fit the clmm:
  fit <-
      if(use.ssr) clmm.fit.ssr(rho, control = control$optCtrl,
                               method=control$method, Hess)
      else clmm.fit.env(rho, control = control$optCtrl,
                        method=control$method, Hess)

  ## Modify and return results:
  fit$nAGQ <- nAGQ
  fit$link <- link
  fit$start <- start
  fit$threshold <- threshold
  fit$call <- match.call()
  fit$formula <- formulae$formula
  fit$gfList <- retrms$gfList
  fit$control <- control
  res <- clmm.finalize(fit=fit, frames=frames, ths=ths, use.ssr)

  ## add model.frame to results list?
  if(model) res$model <- frames$mf

  return(res)
  }

