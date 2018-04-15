## This file contains:
## Functions for fitting CLMMs with a single simple random-effects
## term (ssr).

rho.clm2clmm.ssr <- function(rho, retrms, ctrl)
### Version of rho.clm2clmm that is set up to use the C
### implementations of Laplace, AGQ and GHQ for a single random
### effect.
{
    gfList <- retrms$gfList
    rho$grFac <- gfList[[1]]
    rho$ctrl <- ctrl
    rho$sigma <- rep(1, nrow(rho$B1))
    rho$lambda <- 0
    rho$nlev <- as.vector(sapply(gfList, nlevels))
    rho$random.names <- sapply(gfList, levels)
    rho$tau.names <- names(gfList)
    rho$nrandom <- sum(rho$nlev) ## no. random effects
    rho$Niter <- 0L
    rho$neval <- 0L
    rho$u <- rho$uStart <- rep(0, rho$nrandom)
    rho$linkInt <- switch(rho$link,
                          logit = 1L,
                          probit = 2L,
                          cloglog = 3L,
                          loglog = 4L,
                          cauchit = 5L)
    rho$ST <- lapply(retrms$retrms, `[[`, "ST")
}

## set.AGQ <- function(rho, nAGQ) {
##   rho$nAGQ <- nAGQ
##   if(nAGQ %in% c(0L, 1L)) return(invisible())
##   ghq <- gauss.hermite(abs(nAGQ))
##   rho$ghqns <- ghq$nodes
##   rho$ghqws <-
##     if(nAGQ > 0) ghq$weights ## AGQ
##     else log(ghq$weights) + (ghq$nodes^2)/2 ## GHQ
## }

clmm.fit.ssr <-
  function(rho, control = list(), method=c("nlminb", "ucminf"),
           Hess = FALSE)
### Fit a clmm with a single simple random effects term using AGQ, GHQ
### or Laplace.
{
    optim.error <- function(fit, method)
        if(inherits(fit, "try-error"))
            stop("optimizer ", method, " terminated with an error", call.=FALSE)
### FIXME: Could have an argument c(warn, fail, ignore) to optionally
### return the fitted model despite the optimizer failing.

    method <- match.arg(method)
    ## Set appropriate objective function:
    obj.fun <-
        if(rho$nAGQ < 0) getNGHQ.ssr
        else if(rho$nAGQ > 1) getNAGQ.ssr
        else getNLA.ssr ## nAGQ %in% c(0, 1)

    init.val <- obj.fun(rho, rho$par)
    if(!is.finite(init.val))
        stop(gettextf("non-finite likelihood at starting value (%g)",
                      init.val), call.=FALSE)

    ## Fit the model:
    if(method == "ucminf") {
        fit <- try(ucminf(rho$par, function(par) obj.fun(rho, par),
                          control = control), silent=TRUE)
        ## Check if optimizer converged without error:
        optim.error(fit, method)
        ## Save return value:
        value <- fit$value
    } else if(method == "nlminb") {
        ## hack to remove ucminf control settings:
        keep <- !names(control) %in% c("grad", "grtol")
        control <- if(length(keep)) control[keep] else list()
        fit <- try(nlminb(rho$par, function(par) obj.fun(rho, par),
                          control = control), silent=TRUE)
        ## Check if optimizer converged without error:
        optim.error(fit, method)
        ## Save return value:
        value <- fit$objective
    }
    else stop("unkown optimization method: ", method)
    ## Extract parameters from optimizer results:
    rho$par <- fit$par

    ## Ensure random mode estimation at optimum:
    nllBase.uC(rho)
    update.uC(rho)

    rho$ST <- par2ST(rho$tau, rho$ST)
    names(rho$ST) <- names(rho$dims$nlev.re)

    ## Format ranef modes and condVar:
    ranef <- rho$u * rho$tau
    condVar <- 1/rho$D * rho$tau^2
    ## names(ranef) <- names(condVar) <- rho$random.names
    ## ranef <- list(ranef)
    ## condVar <- list(condVar)
    ## names(ranef) <- names(condVar) <- rho$tau.names

    ## Prepare list of results:
    res <- list(coefficients = fit$par[1:rho$dims$nfepar],
                ST = rho$ST,
                optRes = fit,
                logLik = -value,
                fitted.values = rho$fitted,
                ranef = ranef,
                condVar = condVar,
                dims = rho$dims,
                u = rho$u)
    ## Add gradient vector and optionally Hessian matrix:
    ## bound <- as.logical(paratBoundary2(rho))
    ## optpar <- fit$par[!bound]
    if(Hess) {
        ## gH <- deriv12(function(par) obj.fun(rho, par, which=!bound),
        gH <- deriv12(function(par) obj.fun(rho, par),
                      x=fit$par)
        res$gradient <- gH$gradient
        res$Hessian <- gH$Hessian
    } else {
        ## res$gradient <- grad.ctr(function(par) getNLA(rho, par, which=!bound),
        res$gradient <- grad.ctr(function(par) obj.fun(rho, par),
                                 x=fit$par)
    }
    ## Setting Niter and neval after gradient and Hessian evaluations:
    res$Niter <- rho$Niter
    res$neval <- rho$neval
    return(res)
}


getNGHQ.ssr <- function(rho, par) {
### negative log-likelihood by standard Gauss-Hermite quadrature
### implemented in C:
    if(!missing(par)) {
        rho$par <- par
        if(any(!is.finite(par)))
            stop(gettextf(paste(c("Non-finite parameters occured:",
                                  formatC(par, format="g")), collapse=" ")))
    }
    rho$neval <- rho$neval + 1L
  nllBase.uC(rho) ## Update tau, eta1Fix and eta2Fix
  with(rho, {
    .C("getNGHQ",
       nll = double(1),
       as.integer(grFac),
       as.double(tau),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(o1),
       as.double(o2),
       as.double(sigma),
       as.double(wts),
       length(sigma),
       length(uStart),
       as.double(ghqns),
       as.double(ghqws),
       as.integer(abs(nAGQ)),
       as.integer(linkInt),
       as.double(ghqns * tau),
       as.double(lambda))$nll
  })
}

getNAGQ.ssr <- function(rho, par) {
### negative log-likelihood by adaptive Gauss-Hermite quadrature
### implemented in C:
    if(!missing(par)) {
        rho$par <- par
        if(any(!is.finite(par)))
            stop(gettextf(paste(c("Non-finite parameters occured:",
                                  formatC(par, format="g")), collapse=" ")))
    }
    rho$neval <- rho$neval + 1L
    if(!update.uC(rho)) return(Inf)
  if(any(rho$D < 0)) return(Inf)
  with(rho, {
    .C("getNAGQ",
       nll = double(1),
       as.integer(grFac),
       as.double(tau),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(o1),
       as.double(o2),
       as.double(sigma),
       as.double(wts),
       length(sigma),
       length(uStart),
       as.double(ghqns),
       as.double(log(ghqws)),
       as.double(ghqns^2),
       as.double(u),
       as.double(D),
       as.integer(abs(nAGQ)),
       as.integer(linkInt),
       as.double(lambda))$nll
  })
}

getNLA.ssr <- function(rho, par) {
### negative log-likelihood by the Laplace approximation
### (with update.u2 in C or R):
  if(!missing(par)) {
      rho$par <- par
      if(any(!is.finite(par)))
          stop(gettextf(paste(c("Non-finite parameters occured:",
                                formatC(par, format="g")), collapse=" ")))
  }
    rho$neval <- rho$neval + 1L
  if(!update.uC(rho)) return(Inf)
  if(any(rho$D <= 0)) return(Inf)
  logDetD <- sum(log(rho$D))
  rho$negLogLik - rho$nrandom*log(2*pi)/2 + logDetD/2
}

nllBase.uC <- function(rho) {
### updates tau, eta1Fix and eta2Fix given new parameter values
  with(rho, {
    tau <- exp(par[nalpha + nbeta + 1:ntau])
    eta1Fix <- drop(B1 %*% par[1:(nalpha + nbeta)])
    eta2Fix <- drop(B2 %*% par[1:(nalpha + nbeta)])
  })
  return(invisible())
}

update.uC <- function(rho) {
### C-implementation of NR-algorithm.
  nllBase.uC(rho) ## update: tau, eta1Fix, eta2Fix
  fit <- with(rho, {
    .C("NRalgv3",
       as.integer(ctrl$trace),
       as.integer(ctrl$maxIter),
       as.double(ctrl$gradTol),
       as.integer(ctrl$maxLineIter),
       as.integer(grFac), ## OBS
       as.double(tau), # stDev
       as.double(o1),
       as.double(o2),
       as.double(eta1Fix),
       as.double(eta2Fix),
       as.double(sigma), ## rep(1, n)
       as.integer(linkInt), ##
       as.double(wts), ## pre. weights
       u = as.double(uStart),
       fitted = as.double(fitted), ## pre. pr
       funValue = double(1),
       gradValues = as.double(uStart),
       hessValues = as.double(rep(1, length(uStart))),
       length(fitted),
       length(uStart),
       maxGrad = double(1),
       conv = 0L,
       as.double(lambda), ##
       Niter = as.integer(Niter) ## OBS
       )[c("u", "fitted", "funValue", "gradValues",
           "hessValues", "maxGrad", "conv", "Niter")] })
    ## Get message:
  message <- switch(as.character(fit$conv),
                    "1" = "max|gradient| < tol, so current iterate is probably solution",
                    "0" = "Non finite negative log-likelihood",
                    "-1" = "iteration limit reached when updating the random effects",
                    "-2" = "step factor reduced below minimum when updating the random effects")
  ## Check for convergence and report warning/error:
  if(rho$ctrl$trace > 0 && fit$conv == 1)
    cat("\nOptimizer converged! ", "max|grad|:",
        fit$maxGrad, message, fill = TRUE)
  if(fit$conv != 1 && rho$ctrl$innerCtrl == "warnOnly")
    warning(message, "\n  at iteration ", rho$Niter)
  else if(fit$conv != 1 && rho$ctrl$innerCtrl == "giveError")
    stop(message, "\n  at iteration ", rho$Niter)
  ## Store values and return:
  rho$Niter <- fit$Niter
  rho$fitted <- fit$fitted
  rho$u <- fit$u
  rho$D <- fit$hessValue
  rho$gradient <- fit$gradValue
  if(!is.finite(rho$negLogLik <- fit$funValue))
    return(FALSE)
  return(TRUE)
}

