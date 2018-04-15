## This file contains:
## Utility functions for fitting CLMMs with clmm2().

### FIXME: Could make use of getFittedC throughout this file...

.negLogLikMfastR <- function(rho) { ## negative log-likelihood
### .negLogLikMfast in R
    with(rho, {
        o21 <- u[grFac] * stDev
        o11 <- o1 - o21
        o21 <- o2 - o21
        eta1 <- (eta1Fix + o11)/sigma
        eta2 <- (eta2Fix + o21)/sigma
        pr <-
            if(nlambda) pfun(eta1, lambda) - pfun(eta2, lambda)
            else pfun(eta1) - pfun(eta2)
        if(any(is.na(pr)) || any(pr <= 0))
            nll <- Inf
        else
            nll <- -sum(weights * log(pr)) -
                sum(dnorm(x=u, mean=0, sd=1, log=TRUE))
        nll
    })
}

.negLogLikM <- function(rho) { ## negative log-likelihood
    with(rho, {
        if(estimStDev)
            stDev <- exp(par[p+nxi+k+estimLambda+ 1:s])
        o21 <- u[grFac] * stDev
        o11 <- o1 - o21
        o21 <- o2 - o21
        if(estimLambda > 0)
            lambda <- par[nxi + p + k + 1:estimLambda]
        sigma <-
            if(k > 0) expSoffset * exp(drop(Z %*% par[nxi+p + 1:k]))
            else expSoffset
        eta1Fix <- drop(B1 %*% par[1:(nxi + p)])
        eta2Fix <- drop(B2 %*% par[1:(nxi + p)])
        eta1 <- (eta1Fix + o11)/sigma
        eta2 <- (eta2Fix + o21)/sigma
        pr <-
            if(nlambda) pfun(eta1, lambda) - pfun(eta2, lambda)
            else pfun(eta1) - pfun(eta2)
        if(any(is.na(pr)) || any(pr <= 0))
            nll <- Inf
        else
            nll <- -sum(weights * log(pr)) -
                sum(dnorm(x=u, mean=0, sd=1, log=TRUE))
        nll
    })
}


.gradM <- function(rho) { ## gradient of the negative log-likelihood
    with(rho, {
        if(nlambda) {
            p1 <- dfun(eta1, lambda)
            p2 <- dfun(eta2, lambda)
        }
        else {
            p1 <- dfun(eta1)
            p2 <- dfun(eta2)
        }
        wtprSig <- weights/pr/sigma
        .C("gradC",
           as.double(stDev),
           as.double(p1),
           as.double(p2),
           as.double(wtprSig),
           as.integer(grFac),
           length(wtprSig),
           u = as.double(u),
           length(u))$u
        ## tapply(stDev * wtprSig * (p1 - p2), grFac, sum) + u
    })
}

.gradC <- function(rho) {
    tmp <- with(rho, {
        .C("grad",
           as.double(stDev),
           p1 = double(length(pr)),
           p2 = double(length(pr)),
           as.double(pr),
           as.double(weights),
           as.double(sigma),
           wtprSig = double(length(pr)),
           as.double(eta1),
           as.double(eta2),
           gradValues = double(length(u)),
           as.double(u),
           as.integer(grFac),
           length(pr),
           length(u),
           as.double(lambda),
           as.integer(linkInt))[c("p1", "p2", "wtprSig", "gradValues")]
    })
    rho$wtprSig <- tmp$wtprSig
    rho$p1 <- tmp$p1
    rho$p2 <- tmp$p2
    tmp$gradValues
}

.hessC <- function(rho) {
    with(rho, {
        .C("hess",
           as.double(stDev),
           as.double(p1),
           as.double(p2),
           as.double(pr),
           as.double(wtprSig),
           as.double(eta1),
           as.double(eta2),
           as.integer(linkInt),
           as.integer(grFac),
           length(pr),
           hessValues = double(length(u)),
           as.double(lambda),
           length(u))$hessValues
    })
}

.hessianM <- function(rho)  ## hessian of the negative log-likelihood
    with(rho,{
        if(nlambda) {
            g1 <- gfun(eta1, lambda)
            g2 <- gfun(eta2, lambda)
        }
        else {
            g1 <- gfun(eta1)
            g2 <- gfun(eta2)
        }
        .C("hessC",
           as.double(stDev),
           as.double(p1),
           as.double(p2),
           as.double(pr),
           as.double(g1),
           as.double(g2),
           as.double(wtprSig),
           as.integer(grFac),
           length(pr),
           z = double(length(u)),
           length(u))$z
        ## tapply(((p1 - p2)^2 / pr - g1 + g2) * wtprSig, grFac, sum) *
        ##     stDev^2 + 1
    })


update.u2.v2 <- function(rho) {
### second version: C-implementation of NR-algorithm.
    .negLogLikBase(rho) ## update: par, stDev, eta1Fix, eta2Fix eta2Fix, sigma
    fit <- with(rho,
                .C("NRalg",
                   as.integer(ctrl$trace),
                   as.integer(ctrl$maxIter),
                   as.double(ctrl$gradTol),
                   as.integer(ctrl$maxLineIter),
                   as.integer(grFac),
                   as.double(stDev),
                   as.double(o1),
                   as.double(o2),
                   as.double(eta1Fix),
                   as.double(eta2Fix),
                   as.double(eta1),
                   as.double(eta2),
                   as.double(sigma),
                   as.integer(linkInt),
                   as.double(weights),
                   u = as.double(uStart),
                   pr = as.double(pr),
                   funValue = as.double(nll),
                   gradValues = as.double(uStart),
                   hessValues = as.double(uStart),
                   length(pr),
                   length(uStart),
                   maxGrad = double(1),
                   conv = 0L,
                   double(length(pr)), # p1
                   double(length(pr)), # p2
                   double(length(pr)), # wtprSig
                   as.double(lambda),
                   Niter = as.integer(Niter)
                   )[c("u", "pr", "funValue", "gradValues",
                       "hessValues", "maxGrad", "conv", "Niter")] )
    ## Get message:
    message <- switch(as.character(fit$conv),
                      "1" = "max|gradient| < tol, so current iterate is probably solution",
                      "0" = "Non finite negative log-likelihood",
                      "-1" = "iteration limit reached when updating the random effects",
                      "-2" = "step factor reduced below minimum when updating the random effects")
    ## check for convergence and report warning/error:
    if(rho$ctrl$trace > 0 && fit$conv == 1)
        cat("\nOptimizer converged! ", "max|grad|:",
            fit$maxGrad, message, fill = TRUE)
    if(fit$conv != 1 && rho$ctrl$innerCtrl == "warnOnly")
        warning(message, "\n  at iteration ", rho$Niter)
    else if(fit$conv != 1 && rho$ctrl$innerCtrl == "giveError")
        stop(message, "\n  at iteration ", rho$Niter)
    ## Store values and return:
    rho$Niter <- fit$Niter
    rho$u <- fit$u
    rho$D <- fit$hessValue
    rho$gradient <- fit$gradValue
    if(!is.finite(rho$negLogLik <- fit$funValue))
        return(FALSE)

    return(TRUE)
}

update.u2 <- function(rho)
{
    stepFactor <- 1
    innerIter <- 0
    rho$u <- rho$uStart
    rho$negLogLik <- .negLogLikM(rho)
    if(!is.finite(rho$negLogLik)) return(FALSE)
    rho$gradient <- .gradC(rho)
    maxGrad <- max(abs(rho$gradient))
    conv <- -1  ## Convergence flag
        message <- "iteration limit reached when updating the random effects"
    if(rho$ctrl$trace > 0)
        Trace(iter=0, stepFactor, rho$negLogLik, maxGrad, rho$u, first=TRUE)

    ## Newton-Raphson algorithm:
    for(i in 1:rho$ctrl$maxIter) {
        if(maxGrad < rho$ctrl$gradTol) {
            message <- "max|gradient| < tol, so current iterate is probably solution"
            if(rho$ctrl$trace > 0)
                cat("\nOptimizer converged! ", "max|grad|:",
                    maxGrad, message, fill = TRUE)
            conv <- 0
            break
        }
        rho$D <- .hessC(rho)
        ## rho$D <- .hessianM(rho)
        step <- rho$gradient / rho$D
        rho$u <- rho$u - stepFactor * step
        negLogLikTry <- .negLogLikMfast(rho)
        lineIter <- 0
        ## simple line search, i.e. step halfing:
        while(negLogLikTry > rho$negLogLik) {
            stepFactor <- stepFactor/2
            rho$u <- rho$u + stepFactor * step
            negLogLikTry <- .negLogLikMfast(rho)
            lineIter <- lineIter + 1
            if(rho$ctrl$trace > 0)
                Trace(i+innerIter, stepFactor, rho$negLogLik, maxGrad,
                      rho$u, first=FALSE)
            if(lineIter > rho$ctrl$maxLineIter){
                message <- "step factor reduced below minimum when updating the random effects"
                conv <- 1
                break
            }
            innerIter <- innerIter + 1
        }
        rho$negLogLik <- negLogLikTry
        rho$gradient <- .gradC(rho)
        maxGrad <- max(abs(rho$gradient))
        if(rho$ctrl$trace > 0)
            Trace(i+innerIter, stepFactor, rho$negLogLik, maxGrad, rho$u, first=FALSE)
        stepFactor <- min(1, 2 * stepFactor)
    }
    if(conv != 0 && rho$ctrl$innerCtrl == "warnOnly") {
        warning(message, "\n  at iteration ", rho$Niter)
        utils::flush.console()
    }
    else if(conv != 0 && rho$ctrl$innerCtrl == "giveError")
        stop(message, "\n  at iteration ", rho$Niter)
    rho$Niter <- rho$Niter + i
    rho$D <- .hessC(rho)
    if(!is.finite(rho$negLogLik))
        return(FALSE)
    else
        return(TRUE)
}

.hessMinC <- function(rho) {
    with(rho,{
        if(nlambda) {
            g1 <- gfun(eta1, lambda)
            g2 <- gfun(eta2, lambda)
        }
        else {
            g1 <- gfun(eta1)
            g2 <- gfun(eta2)
        }
        .C("hessC",
           as.double(stDev),
           as.double(p1),
           as.double(p2),
           as.double(pr),
           as.double(g1),
           as.double(g2),
           as.double(wtprSig),
           as.integer(grFac),
           length(pr),
           z = double(length(u)),
           length(u))$z
    })
}

.gradMinC <- function(stDev, p1, p2, wtprSig, grFac, u)
    .C("gradC",
       as.double(stDev),
       as.double(p1),
       as.double(p2),
       as.double(wtprSig),
       as.integer(unclass(grFac)),
       as.integer(length(wtprSig)),
       u = as.double(u),
       as.integer(length(u)))$u

.gradMinC <- function(rho) {
    with(rho, {
        if(nlambda) {
            p1 <- dfun(eta1, lambda)
            p2 <- dfun(eta2, lambda)
        }
        else {
            p1 <- dfun(eta1)
            p2 <- dfun(eta2)
        }
        wtprSig <- weights/pr/sigma
        .C("gradC",
           as.double(stDev),
           as.double(p1),
           as.double(p2),
           as.double(wtprSig),
           as.integer(grFac),
           length(wtprSig),
           u = as.double(u),
           length(u))$u
    })
}

grFacSumC <- function(x, grFac, u)
    .C("grFacSum",
       as.double(x),
       as.integer(grFac),
       as.integer(length(x)),
       u = as.double(u),
       as.integer(length(u)))$u


grFacSum <- function(x, grFac, n.x, u, n.u) {
    ## i, j, z
    z <- 0
    for (i in 1:n.u) {
        for (j in 1:n.x)
            if(grFac[j] == i)
                z <- z + x[j]
        u[i] <- z + u[i]
        z <- 0
    }
    u
}

getNAGQ2 <- function(rho, par) {
### Not in use
    if(!missing(par))
        rho$par <- par
    if(!update.u2(rho)) return(Inf)
    if(any(rho$D < 0)) return(Inf)
    with(rho, {
        K <- sqrt(2/D)
        agqws <- K %*% t(ghqws)
        agqns <- apply(K %*% t(ghqns), 2, function(x) x + u)
        ranNew <- apply(agqns, 2, function(x) x[grFac] * stDev)

        eta1Tmp <- (eta1Fix + o1 - ranNew) / sigma
        eta2Tmp <- (eta2Fix + o2 - ranNew) / sigma
        if(nlambda)
            ## PRnn <- (pfun(eta1Tmp, lambda) - pfun(eta2Tmp, lambda))^weights
            ## This is likely a computationally more safe solution:
            PRnn <- exp(weights * log(pfun(eta1Tmp, lambda) - pfun(eta2Tmp, lambda)))
        else
            ## PRnn <- (pfun(eta1Tmp) - pfun(eta2Tmp))^weights
            PRnn <- exp(weights * log(pfun(eta1Tmp) - pfun(eta2Tmp)))
        for(i in 1:r)
            ## PRrn[i,] <- apply(PRnn[grFac == i, ], 2, prod)
### FIXME: Should this be: ???
            PRrn[i,] <- apply(PRnn[grFac == i, ,drop = FALSE], 2, prod)
        PRrn <- PRrn * agqws * dnorm(x=agqns, mean=0, sd=1)
### FIXME: Could this be optimized by essentially computing dnorm 'by hand'?
    })
    -sum(log(rowSums(rho$PRrn)))
}

getNGHQ <- function(rho, par) {
### Not in use
    if(!missing(par))
        rho$par <- par
    .negLogLikM(rho) ## Update lambda, stDev, sigma and eta*Fix
    with(rho, {
        eta1Tmp <- (eta1Fix + o1 - ranNew * stDev) / sigma
        eta2Tmp <- (eta2Fix + o2 - ranNew * stDev) / sigma
        if(nlambda)
            PRnn <- exp(weights * log(pfun(eta1Tmp, lambda) - pfun(eta2Tmp, lambda)))
        else
            PRnn <- exp(weights * log(pfun(eta1Tmp) - pfun(eta2Tmp)))
        for(i in 1:r)
            PRrn[i,] <- apply(PRnn[grFac == i, ,drop = FALSE], 2, prod)
        PRrn <- PRrn * agqws * dnorm(x=agqns, mean=0, sd=1)
    })
    -sum(log(rowSums(rho$PRrn)))
}
