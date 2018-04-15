## This file contains:
## Functions to fit/estimate CLMs (clm.fit.NR, clm.fit.optim) and
## functions implementing the negative log-likelihood, its gradient
## and hessian (.nll, .grad, .hess). These functions are rarely to be
## called directly from outside the package.

clm.fit.NR <-
  function(rho, control = list())
### The main work horse: Where the actual fitting of the clm goes on.
### Fitting the clm via modified Newton-Raphson with step halving.

### -------- Assumes the existence of the following functions:
### clm.nll - negative log-likelihood
### clm.grad - gradient of nll wrt. par
### clm.hess - hessian of nll wrt. par
### Trace - for trace information
{
    control <- do.call(clm.control, control)
    stepFactor <- 1
    innerIter <- modif.iter <- abs.iter <- 0L
    conv <- 2L  ## Convergence flag (iteration limit reached)
    nll <- rho$clm.nll(rho)
    if(!is.finite(nll))
        stop("Non-finite log-likelihood at starting value")
    ## do.newton <-
    ## rel.conv <- FALSE
    ## stephalf <- TRUE

    ## Newton-Raphson algorithm:
    for(i in 1:(control$maxIter + 1L)) {
        gradient <- rho$clm.grad(rho)
        maxGrad <- max(abs(gradient))
        if(control$trace > 0) {
            Trace(iter=i+innerIter-1, stepFactor, nll, maxGrad,
                  rho$par, first=(i==1))
            if(control$trace > 1 && i > 1) {
                cat("\tgrad: ")
                cat(paste(formatC(gradient, digits=3, format="e")))
                cat("\n\tstep: ")
                cat(paste(formatC(-step, digits=3, format="e")))
                cat("\n\teigen: ")
                cat(paste(formatC(eigen(hessian, symmetric=TRUE,
                                        only.values=TRUE)$values, digits=3,
                                  format="e")))
                cat("\n")
            }
        }
        abs.conv <- (maxGrad < control$gradTol)
        if(abs.conv) abs.iter <- abs.iter + 1L
        hessian <- rho$clm.hess(rho)
        ## Compute cholesky factor of Hessian: ch = Ut U
        ch <- try(chol(hessian), silent=TRUE)
### NOTE: solve(hessian, gradient) is not good enough because it will
### compute step for negative-definite Hessians and we don't want
### that.

### FIXME: What if Hessian is closely singular but slightly positive?
### Could we do something better in that case?
        if(inherits(ch, "try-error")) {
            if(abs.conv) { ## step.ok not true.
                conv <- 1L
                break ## cannot meet relative criterion.
            }
            ## If Hessian is non-positive definite:
            min.ev <- min(eigen(hessian, symmetric=TRUE,
                                only.values=TRUE)$values)
            inflation.factor <- 1
            ## Inflate diagonal of Hessian to make it positive definite:
            inflate <- abs(min.ev) + inflation.factor
            hessian <- hessian + diag(inflate, nrow(hessian))
            if(control$trace > 0)
                cat(paste("Hessian is singular at iteration", i-1, "inflating diagonal with",
                          formatC(inflate, digits=5, format="f"), "\n"))
            ch <- try(chol(hessian), silent=TRUE)
            if(inherits(ch, "try-error"))
                stop(gettextf("Cannot compute Newton step at iteration %d",
                              i-1), call.=FALSE)
            modif.iter <- modif.iter + 1L
            ## do.newton <- FALSE
        } else
            modif.iter <- 0L
        if(modif.iter >= control$maxModIter) {
            conv <- 4L
            break
        }

        ## solve U'y = g for y, then
        ## solve U step = y for step:
        step <- c(backsolve(ch, backsolve(ch, gradient, transpose=TRUE)))
        rel.conv <- (max(abs(step)) < control$relTol)
        ## Test if step is in a descent direction -
        ## otherwise use step <- grad / max|grad|:
        ## if(crossprod(gradient, step) < 0) {
        ##     if(control$trace > 0)
        ##         cat("Newton step is not in descent direction; using gradient instead\n")
        ##     step <- c(gradient / max(abs(gradient)))
        ## } else
        if(abs.conv && rel.conv) {
            conv <- 0L
            ## no need to step back as stephalf was false so the new
            ## par are just better.
            break
        }

        ## update parameters:
        rho$par <- rho$par - stepFactor * step
        nllTry <- rho$clm.nll(rho)
        lineIter <- 0
        stephalf <- (nllTry > nll)
### NOTE: sometimes nllTry > nll just due to noise, so we also check
### reduction in gradient for small diffs:
        if(stephalf && abs(nll - nllTry) < 1e-10)
            stephalf <- maxGrad < max(abs(rho$clm.grad(rho)))

        ## Assess convergence:
        ## (only attempt to sattisfy rel.conv if abs.conv is true and
        ## it is possible to take the full newton step)
### OPTION: And if 'step' is not close to 1 or 1/2, but
### small. Otherwise this just indicates that the parameter is
### infinite.
        ## if(abs.conv && !step.ok) {
        if(abs.conv && stephalf) {
            conv <- 1L
            ## we need to step back to the par for which abs.conv
            ## was true:
            rho$par <- rho$par + stepFactor * step
            rho$clm.nll(rho)
            break
        }
        ## if(abs.conv && rel.conv) {
        ##     conv <- 0L
        ##     rho$par <- rho$par + stepFactor * step
        ##     rho$clm.nll(rho)
        ##     ## no need to step back as stephalf was false so the new
        ##     ## par are just better.
        ##     break
        ## }
        if(abs.conv && abs.iter >= 5L) {
            ## Cannot satisy rel.conv in 5 iterations after satisfying
            ## abs.conv. Probably some parameters are unbounded.
            conv <- 1L
            break
        }
        ## Step halving if nll increases:
        while(stephalf) {
            stepFactor <- stepFactor/2
            rho$par <- rho$par + stepFactor * step
            nllTry <- rho$clm.nll(rho)
            lineIter <- lineIter + 1
            if(control$trace > 0) {
                cat("step halving:\n")
                cat("nll reduction: ", formatC(nll - nllTry, digits=5, format="e"), "\n")
                Trace(i+innerIter-1, stepFactor, nll, maxGrad,
                      rho$par, first = FALSE)
            }
            if(lineIter > control$maxLineIter){
                conv <- 3L
                break
            }
            innerIter <- innerIter + 1
            stephalf <- (nllTry > nll)
            if(stephalf && abs(nll - nllTry) < 1e-10)
                stephalf <- (maxGrad < max(abs(rho$clm.grad(rho))))
        } ## end step halving
        if(conv == 3L) break

        if(control$trace > 0)
            cat("nll reduction: ", formatC(nll - nllTry, digits=5, format="e"), "\n")
        nll <- nllTry
        ## Double stepFactor if needed:
        stepFactor <- min(1, 2 * stepFactor)
    } ## end Newton iterations
    message <- switch(as.character(conv),
                      "0" = "Absolute and relative convergence criteria were met",
                      "1" = "Absolute convergence criterion was met, but relative criterion was not met",
                      "2" = "iteration limit reached",
                      "3" = "step factor reduced below minimum",
                      "4" = "maximum number of consecutive Newton modifications reached")
    if(conv <= 1L && control$trace > 0) {
        cat("\nOptimizer converged! ", message, fill = TRUE)
    }
    if(conv > 1 && control$trace > 0) {
        cat("\nOptimization failed ", message, fill = TRUE)
    }
    ## return results:
    gradient <- c(rho$clm.grad(rho))
    res <- list(par = rho$par,
                gradient = gradient, ##as.vector(gradient),
                ## Hessian = hessian,
                Hessian = rho$clm.hess(rho), ## ensure hessian is evaluated
                ## at optimum
                logLik = -nll,
                convergence = conv,
                ## 0: abs and rel criteria meet
                ## 1: abs criteria meet, rel criteria not meet
                ## 2: iteration limit reached
                ## 3: step factor reduced below minium
                message = message,
                maxGradient = max(abs(gradient)),
                niter = c(outer = i-1, inner = innerIter),
                fitted = rho$fitted)
    return(res)
}


clm.fit.optim <-
  function(rho, method = c("ucminf", "nlminb", "optim"), control=list())
{
  method <- match.arg(method)
  ## optimize the likelihood:
  optRes <-
    switch(method,
           "nlminb" = nlminb(rho$par,
             function(par) clm.nll(rho, par),
             function(par) clm.grad_direct(rho, par),
             control=control),
           "ucminf" = ucminf(rho$par,
             function(par) clm.nll(rho, par),
             function(par) clm.grad_direct(rho, par),
             control=control),
           "optim" = optim(rho$par,
             function(par) clm.nll(rho, par),
             function(par) clm.grad_direct(rho, par),
             method="BFGS",
             control=control),
           )
  ## save results:
  rho$par <- optRes[[1]]
  res <- list(par = rho$par,
              logLik = -clm.nll(rho),
              gradient = clm.grad(rho),
              Hessian = clm.hess(rho),
              fitted = rho$fitted)
  res$maxGradient = max(abs(res$gradient))
  res$optRes <- optRes
  res$niter <- switch(method, "nlminb" = optRes$evaluations,
                      "ucminf" = c(optRes$info["neval"], 0),
                      "optim" = optRes$counts)
  res$convergence <-
    switch(method, "nlminb" = optRes$convergence,
           "ucminf" = optRes$convergence,
           "optim" = optRes$convergence)

  return(res)
}


clm.nll <- function(rho, par) {
  if(!missing(par)) rho$par <- par
  with(rho, {
      if(k > 0)
      sigma <- Soff * exp(drop(S %*% par[n.psi + 1:k]))
### NOTE: we have to divide by sigma even if k=0 since there may be an
### offset but no predictors in the scale model:
    eta1 <- (drop(B1 %*% par[1:n.psi]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:n.psi]) + o2)/sigma
  })
### NOTE: getFitted is not found from within rho, so we have to
### evalueate it outside of rho
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(all(is.finite(rho$fitted)) && all(rho$fitted > 0))
### NOTE: Need test here because some fitted <= 0 if thresholds are
### not ordered increasingly.
    -sum(rho$wts * log(rho$fitted))
  else Inf
}

## clm.nll <- function(rho) { ## negative log-likelihood
## ### For linear models
##   with(rho, {
##     eta1 <- drop(B1 %*% par) + o1
##     eta2 <- drop(B2 %*% par) + o2
##   })
## ### NOTE: getFitted is not found from within rho, so we have to
## ### evalueate it outside of rho
##   rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
##   if(all(rho$fitted > 0))
## ### NOTE: Need test here because some fitted <= 0 if thresholds are
## ### not ordered increasingly.
## ### It is assumed that 'all(is.finite(pr)) == TRUE'
##     -sum(rho$wts * log(rho$fitted))
##   else Inf
## }

## clm.grad <- function(rho) { ## gradient of the negative log-likelihood
## ### return: vector of gradients
## ### For linear models
##   with(rho, {
##     p1 <- dfun(eta1)
##     p2 <- dfun(eta2)
##     wtpr <- wts/fitted
##     dpi.psi <- B1 * p1 - B2 * p2
##     -crossprod(dpi.psi, wtpr)
## ### NOTE: It is assumed that all(fitted > 0) == TRUE and that
## ### all(is.finite(c(p1, p2))) == TRUE
##   })
## }

clm.grad <- function(rho) {
### requires that clm.nll has been called prior to
### clm.grad.
  with(rho, {
    p1 <- dfun(eta1)
    p2 <- dfun(eta2)
    wtpr <- wts/fitted
    C2 <- B1*p1/sigma - B2*p2/sigma
    if(k <= 0) return(-crossprod(C2, wtpr))
    C3 <- -(eta1 * p1 - eta2 * p2) * S
    return(-crossprod(cbind(C2, C3), wtpr))
### NOTE: C2 and C3 are used by clm.hess
  })
}

clm.grad_direct <- function(rho, par) {
### does not require that clm.nll has been called prior to
### clm.grad.
    clm.nll(rho, par)
    clm.grad(rho)
}

## clm.hess <- function(rho) { ## hessian of the negative log-likelihood
## ### return Hessian matrix
## ### For linear models
##   with(rho, {
##     dg.psi <- crossprod(B1 * gfun(eta1) * wtpr, B1) -
##       crossprod(B2 * gfun(eta2) * wtpr, B2)
##     -dg.psi + crossprod(dpi.psi, (dpi.psi * wtpr / fitted))
## ### NOTE: It is assumed that all(fitted > 0) == TRUE and that
## ### all(is.finite(c(g1, g2))) == TRUE
##   })
## }

clm.hess <- function(rho) {
### requires that clm.grad has been called prior to this.
    with(rho, {
        g1 <- gfun(eta1)
        g2 <- gfun(eta2)
        wtprpr <- wtpr/fitted ## Phi3
        dg.psi <- crossprod(B1 * gfun(eta1) * wtpr / sigma^2, B1) -
            crossprod(B2 * gfun(eta2) * wtpr / sigma^2, B2)
        ## upper left:
        D <- dg.psi - crossprod(C2, (C2 * wtprpr))
        if(k <= 0) return(-D) ## no scale predictors
        ## upper right (lower left transpose):
        wtprsig <- wtpr/sigma
        epg1 <- p1 + g1*eta1
        epg2 <- p2 + g2*eta2
        Et <- crossprod(B1, -wtprsig * epg1 * S) -
            crossprod(B2, -wtprsig * epg2 * S) -
                crossprod(C2, wtprpr * C3)
        ## lower right:
        F <- -crossprod(S, wtpr * ((eta1*p1 - eta2*p2)^2 / fitted -
                                   (eta1*epg1 - eta2*epg2)) * S)
        ## combine and return hessian:
        H <- rbind(cbind(D    , Et),
                   cbind(t(Et), F))
        return(-H)
    })
}

