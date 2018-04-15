## This file contains:
## Implementation of of nominal_test.clm() and scale_test.clm() for
## automatic testing of nominal and scale effects in clm()s. These
## functions work in a fashion similar to add1().

nominal_test <- function(object, ...) {
  UseMethod("nominal_test")
}

scale_test <- function(object, ...) {
  UseMethod("scale_test")
}

nominal_test.clm <-
  function(object, scope, trace=FALSE, ...)
### Test nominal effects for all (or selected) terms in location
### and scale formulas.
{
    ## get scope: vector of terms names which to add to nominal:
    termsnm <- attr(object$terms, "term.labels")
    if(!is.null(object$S.terms))
        termsnm <- union(termsnm, attr(object$S.terms, "term.labels"))
    if(!missing(scope) && !is.null(scope)) {
        if(!is.character(scope))
            scope <- attr(terms(update.formula(object, scope)),
                          "term.labels")
        if(!all(match(scope, termsnm, 0L) > 0L))
            stop("scope is not a subset of term labels")
    } else {
        scope <- termsnm
    }
    if(!is.null(object$nom.terms)) {
        scope <- scope[!scope %in% attr(object$nom.terms,
                                        "term.labels")]
    }
    if(!length(scope))
        message("\nno additional terms to add to nominal\n")
    env <- environment(formula(object))
    ## get list of (updated) nominal formulas:
    nomforms <- if(!is.null(object$call$nominal))
        lapply(scope, function(tm) {
            update.formula(old=formula(object$nom.terms),
                           new=as.formula(paste("~. + ", tm)))
        })  else lapply(scope, function(tm) {
            as.formula(paste("~", tm), env=env) })
    ns <- length(scope)
    ## results matrix:
    ans <- matrix(nrow = ns + 1L, ncol = 3L,
                  dimnames = list(c("<none>", scope),
                  c("df", "logLik", "AIC")))
    ans[1L, ] <- c(object$edf, object$logLik, AIC(object))
    n0 <- nobs(object)
    ## for all terms in scope:
    i <- 1
    for(i in seq(ns)) {
        if(trace) {
            cat("trying +", scope[i], "\n", sep = " ")
            utils::flush.console()
        }
        ## update and fit model with nominal effect added:
        nfit <- try(update(object, nominal=nomforms[[i]],
                           convergence="silent"), silent=TRUE)
        ## model may not be identifiable or converge:
        if(!inherits(nfit, "try-error") &&
### NOTE: non-negative convergence codes indicate that the likelihood
### is correctly determined:
           nfit$convergence$code >= 0) {
            ans[i + 1L, ] <- c(nfit$edf, nfit$logLik, AIC(nfit))
            nnew <- nobs(nfit)
            if(all(is.finite(c(n0, nnew))) && nnew != n0)
                stop("number of rows in use has changed: remove missing values?")
        }
    }
    dfs <- ans[, 1L] - ans[1L, 1L]
    dfs[1L] <- NA
    aod <- data.frame(Df = dfs, logLik = ans[, 2L], AIC = ans[, 3L])
    rownames(aod) <- rownames(ans)
    ## compute likelihood ratio statistic and p-values:
    LR <- 2*(ans[, 2L] - ans[1L, 2L])
    LR[1L] <- NA
    nas <- !is.na(LR)
    P <- LR
    P[nas] <- pchisq(LR[nas], dfs[nas], lower.tail = FALSE)
    aod[, c("LRT", "Pr(>Chi)")] <- list(LR, P)
    head <- c("Tests of nominal effects",
              paste("\nformula:", Deparse(formula(object$terms))))
    if(!is.null(object$call$scale))
        head <- c(head, paste("scale:  ",
                              Deparse(formula(object$S.terms))))
    if(!is.null(object$call$nominal))
        head <- c(head, paste("nominal:",
                              Deparse(formula(object$nom.terms))))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

scale_test.clm <-
  function(object, scope, trace=FALSE, ...)
### Test scale effects for all (or selected) terms in formula
{
  ## get scope: vector of terms names which to add to scale:
  termsnm <- attr(object$terms, "term.labels")
  if(!missing(scope) && !is.null(scope)) {
    if(!is.character(scope))
      scope <- attr(terms(update.formula(object, scope)),
                    "term.labels")
    if(!all(match(scope, termsnm, 0L) > 0L))
      stop("scope is not a subset of term labels")
  } else {
    scope <- termsnm
  }
  ## if(!is.null(object$nom.terms)) {
  ##     scope <- scope[!scope %in% attr(object$nom.terms,
  ##                                     "term.labels")]
  ## }
  if(!is.null(object$S.terms)) {
      scope <- scope[!scope %in% attr(object$S.terms,
                                      "term.labels")]
  }
  if(!length(scope))
      message("\nno relevant terms to add to scale\n")
  env <- environment(formula(object))
  ## get list of (updated) scale formulas:
  scaleforms <-
    if(!is.null(object$call$scale))
      lapply(scope, function(tm) {
        update.formula(old=formula(object$S.terms),
                       new=as.formula(paste("~. + ", tm)))
      })
    else
      lapply(scope, function(tm) as.formula(paste("~", tm), env=env))
  ns <- length(scope)
  ## results matrix:
  ans <- matrix(nrow = ns + 1L, ncol = 3L,
                dimnames = list(c("<none>", scope),
                  c("df", "logLik", "AIC")))
  ans[1L, ] <- c(object$edf, object$logLik, AIC(object))
  n0 <- nobs(object)
  ## for all terms in scope:
  for(i in seq(ns)) {
      if(trace) {
          cat("trying +", scope[i], "\n", sep = " ")
          utils::flush.console()
      }
      ## update and fit model with scale effect added:
      nfit <- try(update(object, scale=scaleforms[[i]]), silent=TRUE)
      ## model may not be identifiable or converge:
      if(!inherits(nfit, "try-error") &&
         nfit$convergence$code >= 0) {
          ans[i + 1L, ] <- c(nfit$edf, nfit$logLik, AIC(nfit))
          nnew <- nobs(nfit)
          if (all(is.finite(c(n0, nnew))) && nnew != n0)
              stop("number of rows in use has changed: remove missing values?")
      }
  }
  dfs <- ans[, 1L] - ans[1L, 1L]
  dfs[1L] <- NA
  aod <- data.frame(Df = dfs, logLik = ans[, 2L], AIC = ans[, 3L])
  rownames(aod) <- rownames(ans)
  ## compute likelihood ratio statistic and p-values:
  LR <- 2*(ans[, 2L] - ans[1L, 2L])
  LR[1L] <- NA
  nas <- !is.na(LR)
  P <- LR
  P[nas] <- pchisq(LR[nas], dfs[nas], lower.tail = FALSE)
  aod[, c("LRT", "Pr(>Chi)")] <- list(LR, P)
  head <- c("Tests of scale effects",
            paste("\nformula:", Deparse(formula(object$terms))))
  if(!is.null(object$call$scale))
    head <- c(head, paste("scale:  ",
                          Deparse(formula(object$S.terms))))
  if(!is.null(object$call$nominal))
    head <- c(head, paste("nominal:",
                          Deparse(formula(object$nom.terms))))
  class(aod) <- c("anova", "data.frame")
  attr(aod, "heading") <- head
  aod
}

