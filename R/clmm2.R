## This file contains:
## The main clmm2 function and some related auxiliary functions.

clmm2.control <-
    function(method = c("ucminf", "nlminb", "model.frame"),
             ..., trace = 0, maxIter = 50, gradTol = 1e-4,
             maxLineIter = 50,
             innerCtrl = c("warnOnly", "noWarn", "giveError"))
{
    method <- match.arg(method)
    innerCtrl <- match.arg(innerCtrl)
    ctrl <- list(trace=ifelse(trace < 0, 1, 0),
                 maxIter=maxIter,
                 gradTol=gradTol,
                 maxLineIter=maxLineIter,
                 innerCtrl=innerCtrl)
    optCtrl <- list(trace = abs(trace), ...)

    if(!is.numeric(unlist(ctrl[-5])))
        stop("maxIter, gradTol, maxLineIter and trace should all be numeric")
    if(any(ctrl[-c(1, 5)] <= 0))
       stop("maxIter, gradTol and maxLineIter have to be > 0")
    if(method == "ucminf" && !"grtol" %in% names(optCtrl))
        optCtrl$grtol <- 1e-5
    if(method == "ucminf" && !"grad" %in% names(optCtrl))
        optCtrl$grad <- "central"

    list(method = method, ctrl = ctrl, optCtrl = optCtrl)
}


.negLogLikBase <- function(rho) {
### Update stDev, sigma, eta1Fix, and eta2Fix given new par:
    with(rho, {
        if(estimLambda > 0)
            lambda <- par[nxi + p + k + 1:estimLambda]
        if(estimStDev)
            stDev <- exp(par[p+nxi+k+estimLambda+ 1:s])
        sigma <-
            if(k > 0) expSoffset * exp(drop(Z %*% par[nxi+p + 1:k]))
            else expSoffset
        eta1Fix <- drop(B1 %*% par[1:(nxi + p)])
        eta2Fix <- drop(B2 %*% par[1:(nxi + p)])
    })
    return(invisible())
}

.negLogLikMfast <- function(rho) { ## negative log-likelihood
    fit <- with(rho, {
        .C("nll",
           as.double(u),
           length(u),
           as.integer(grFac),
           as.double(stDev),
           as.double(o1),
           as.double(o2),
           length(o1),
           eta1 = as.double(eta1),
           eta2 = as.double(eta2),
           as.double(eta1Fix),
           as.double(eta2Fix),
           as.double(sigma),
           pr = as.double(pr),
           as.double(weights),
           as.double(lambda),
           as.integer(linkInt),
           nll = double(1)
           )[c("eta1", "eta2", "pr", "nll")]
    })
    rho$eta1 <- fit$eta1
    rho$eta2 <- fit$eta2
    rho$pr <- fit$pr
    fit$nll
}

update.u2.v3 <- function(rho) {
### third version: C-implementation of NR-algorithm.
    .negLogLikBase(rho) ## update: par, stDev, eta1Fix, eta2Fix eta2Fix, sigma
    fit <- with(rho,
                .C("NRalgv3",
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
                   as.double(sigma),
                   as.integer(linkInt),
                   as.double(weights),
                   u = as.double(uStart),
                   pr = as.double(pr),
                   funValue = double(1),
                   gradValues = as.double(uStart),
                   hessValues = as.double(rep(1, length(uStart))),
                   length(pr),
                   length(uStart),
                   maxGrad = double(1),
                   conv = 0L,
                   as.double(lambda),
                   Niter = as.integer(Niter)
                   )[c("u", "funValue", "gradValues",
                       "hessValues", "maxGrad", "conv", "Niter")] )
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
    rho$u <- fit$u
    rho$D <- fit$hessValue
    rho$gradient <- fit$gradValue
    if(!is.finite(rho$negLogLik <- fit$funValue))
        return(FALSE)
    return(TRUE)
}

clmm2 <-
  function(location, scale, nominal, random, data, weights, start, subset,
           na.action, contrasts, Hess = FALSE, model = TRUE, sdFixed,
           link = c("logistic", "probit", "cloglog", "loglog",
           "cauchit", "Aranda-Ordaz", "log-gamma"), lambda,
           doFit = TRUE, control, nAGQ = 1,
           threshold = c("flexible", "symmetric", "equidistant"), ...)
    ## Handle if model = FALSE
### Marginal fitted values? (pr | u = 0) or (pr | u = u.tilde) ?
### How can we (should we?) get u.tilde and var(u.tilde) with GHQ?
### Make safeStart function if !is.finite(negLogLik)
### Write test suite for doFit argument
{
    R <- match.call(expand.dots = FALSE)
    Call <- match.call()
    if(missing(random)) {
        Call[[1]] <- as.name("clm2")
        return(eval.parent(Call))
    }
    if(missing(lambda)) lambda <- NULL
    if(missing(contrasts)) contrasts <- NULL
    if(missing(control)) control <- clmm2.control(...)
    if(!setequal(names(control), c("method", "ctrl", "optCtrl")))
       stop("specify 'control' via clmm2.control()")
    if (missing(data)) data <- environment(location)
    if (is.matrix(eval.parent(R$data)))
        R$data <- as.data.frame(data)
### Collect all variables in a single formula and evaluate to handle
### missing values correctly:
    m <- match(c("location", "scale", "nominal"), names(R), 0)
    F <- lapply(as.list(R[m]), eval.parent) ## evaluate in parent
    varNames <- unique(unlist(lapply(F, all.vars)))
    longFormula <-
        eval(parse(text = paste("~", paste(varNames, collapse = "+")))[1])
    m <- match(c("location", "data", "subset", "weights", "random",
                 "na.action"), names(R), 0)
    R <- R[c(1, m)]
    R$location <- longFormula
    R$drop.unused.levels <- TRUE
    R[[1]] <- as.name("model.frame")
    names(R)[names(R) == "location"] <- "formula"
    R <- eval.parent(R)
    nonNA <- rownames(R)
### Append nonNA index to Call$subset to get the right design matrices
### from clm2:
    Call$subset <-
        if(is.null(Call$subset)) nonNA
        else c(paste(deparse(Call$subset), "&"), nonNA)
    Call$start <-
        if(is.null(Call$start) || !is.null(Call$sdFixed)) Call$start
        else start[-length(start)]
    Call$random <- Call$control <- Call$nAGQ <- Call$sdFixed <-
        Call$innerCtrl <- NULL
    Call$method <- control$method
    Call$doFit <- Call$Hess <- FALSE
    Call[[1]] <- as.name("clm2")
    rhoM <- eval.parent(Call)
    if(control$method == "model.frame")
        return(rhoM)
    rhoM$call <- match.call()
    rhoM$randomName <- deparse(rhoM$call$random)
### Set grouping factor and stDev parameter:
    rhoM$grFac <- R[,"(random)"]
    if(!missing(sdFixed) && !is.null(sdFixed)) {
        stopifnot(length(sdFixed) == 1 && sdFixed > 0)
        rhoM$estimStDev <- FALSE
        rhoM$stDev <- sdFixed
    }
    else
        rhoM$estimStDev <- TRUE
    with(rhoM, {
        r <- nlevels(grFac) ## no. random effects
        grFac <- as.integer(unclass(grFac))
        if(r <= 2) stop("Grouping factor must have 3 or more levels")
        s <- ifelse(estimStDev, 1L, 0L) ## no. variance parameters
        Niter <- 0L
    })
### set starting values:
    if(missing(start)) {
        suppressWarnings(fitCLM(rhoM))
        if(rhoM$estimStDev) rhoM$start <- rhoM$par <- c(rhoM$par, log(1))
        else rhoM$start <- rhoM$par
    } else
        rhoM$start <- rhoM$par <- start
    rhoM$uStart <- rhoM$u <- rep(0, rhoM$r)
### Test starting values:
    if(length(rhoM$start) !=
       with(rhoM, nxi + p + k + estimLambda + estimStDev))
        stop("'start' is ", length(rhoM$start),
             " long, but should be ", with(rhoM, nxi + p + k + estimLambda + estimStDev),
             " long")
    if(rhoM$ncolXX == 0) {
        if(!all(diff(c(rhoM$tJac %*% rhoM$start[1:rhoM$nalpha])) > 0))
            stop("Threshold starting values are not of increasing size")
    }
### Change the lower limit if lambda is estimated with the
### Aranda-Ordaz link and sdFixed is not supplied:
    if(rhoM$estimLambda > 0 && rhoM$link == "Aranda-Ordaz" &&
       is.null(rhoM$call$sdFixed))
        rhoM$limitLow <- c(rep(-Inf, length(rhoM$par)-2), 1e-5, -Inf)
### This should hardly ever be the case:
    .negLogLikBase(rhoM) ## set lambda, stDev, sigma, eta1Fix and eta2Fix
    if(!is.finite(.negLogLikMfast(rhoM)))
        stop("Non-finite integrand at starting values")
    rhoM$ctrl <- control$ctrl
    rhoM$optCtrl <- control$optCtrl
    if(rhoM$method == "nlminb") {
        m <- match(names(rhoM$optCtrl), c("grad","grtol"), 0)
        rhoM$optCtrl <- rhoM$optCtrl[!m]
    }
### Match doFit:
    if(is.logical(doFit) || is.numeric(doFit)) {
        if(doFit) doFit <- "C"
        else doFit <- "no"
    }
    else if(!is.character(doFit) || !doFit %in% c("no", "R", "C"))
        stop("argument 'doFit' not recognized. 'doFit' should be\n
numeric, logical or one of c('no', 'R', 'C')")

### Set ObjFun parameters:
    ObjFun <- getNLA2 ## same for "R" and "C"
    rhoM$updateU <-
        if(doFit == "R") update.u2
        else update.u2.v3
    rhoM$nAGQ <- as.integer(nAGQ)
    if(rhoM$nAGQ >= 2) {
        ghq <- gauss.hermite(rhoM$nAGQ)
        rhoM$ghqns <- ghq$nodes
        rhoM$ghqws <- ghq$weights
        if(doFit == "R") {
            ObjFun <- getNAGQinR
            rhoM$PRnn <- array(0, dim=c(rhoM$n, rhoM$nAGQ))
            rhoM$PRrn <- array(0, dim=c(rhoM$r, rhoM$nAGQ))
            rhoM$ghqws <- ghq$weights * exp(rhoM$ghqns^2)
        }
        else
            ObjFun <- getNAGQinC
    }
    if(rhoM$nAGQ <= -1) {
        ghq <- gauss.hermite(abs(rhoM$nAGQ))
        rhoM$ghqns <- ghq$nodes
        rhoM$ghqws <- ghq$weights * exp((ghq$nodes^2)/2)
        if(doFit == "R"){
            ObjFun <- getNGHQinR
        }
        else {
            ObjFun <- getNGHQinC
            rhoM$ghqws <- log(ghq$weights) + (ghq$nodes^2)/2
        }
    }
    stopifnot(rhoM$nAGQ != 0) ## test needed?

### Fit the model:
    if(!doFit %in% c("C", "R"))
        return(rhoM)
    if(rhoM$nAGQ > -1)
        rhoM$updateU(rhoM) # Try updating the random effects
    rhoM$optRes <- switch(rhoM$method,
                       "ucminf" = ucminf(rhoM$start, function(x)
                       ObjFun(rhoM, x), control=rhoM$optCtrl),
                       "nlminb" = nlminb(rhoM$start, function(x)
                       ObjFun(rhoM, x), control=rhoM$optCtrl,
                       lower = rhoM$limitLow, upper = rhoM$limitUp))
    rhoM$par <- rhoM$optRes[[1]]
    if(Hess) {
        if(rhoM$link == "Aranda-Ordaz" &&
           rhoM$estimLambda > 0 && rhoM$lambda < 1e-3)
            message("Cannot get Hessian because lambda = ",rhoM$lambda
                    ," is too close to boundary.\n",
                    " Fit model with link == 'logistic' to get Hessian")
        else {
            rhoM$Hessian <- myhess(function(x) ObjFun(rhoM, x),
                                    rhoM$par)
            rhoM$par <- rhoM$optRes[[1]]
        }
    }
    .negLogLikMfast(rhoM) ## update pr
    ## if(rhoM$nAGQ > -1)
    rhoM$updateU(rhoM) # Makes sure ranef's are evaluated at the optimum
### Post processing:
    res <- finalizeRhoM(rhoM)
    res$call <- match.call()
    res$na.action <- attr(R, "na.action")
    res$contrasts <- contrasts
    class(res) <- c("clmm2", "clm2")
    res
}

getNLA2 <- function(rho, par) {
### negative log-likelihood by the Laplace approximation
### (with update.u2 in C or R):
    if(!missing(par)) rho$par <- par
    if(!rho$updateU(rho)) return(Inf)
    if(any(rho$D < 0)) return(Inf)
    ## logDetD <- sum(log(rho$D/(2*pi)))
    logDetD <- sum(log(rho$D)) - rho$r * log(2*pi)
    rho$negLogLik + logDetD / 2
}

getNAGQinR <- function(rho, par) {
### negative log-likelihood by adaptive Gauss-Hermite quadrature
### implemented in R:
    if(!missing(par))
        rho$par <- par
    if(!rho$updateU(rho)) return(Inf)
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
          PRnn <- exp(weights * log(pfun(eta1Tmp, lambda) -
                                    pfun(eta2Tmp, lambda)))
        else
            ## PRnn <- (pfun(eta1Tmp) - pfun(eta2Tmp))^weights
            PRnn <- exp(weights * log(pfun(eta1Tmp) - pfun(eta2Tmp)))
### FIXME: The fitted values could be evaluated with getFittedC for
### better precision.
        for(i in 1:r)
            ## PRrn[i,] <- apply(PRnn[grFac == i, ], 2, prod)
### FIXME: Should this be: ???
            PRrn[i,] <- apply(PRnn[grFac == i, ,drop = FALSE], 2, prod)
        PRrn <- PRrn * agqws * dnorm(x=agqns, mean=0, sd=1)
### FIXME: Could this be optimized by essentially computing dnorm 'by hand'?
    })
    -sum(log(rowSums(rho$PRrn)))
}

## tmpAGQ(rho)

tmpAGQ <- function(rho, par) {
    if(!missing(par))
        rho$par <- par
    with(rho, {
        ls()
        stDev <- exp(ST[[1]][1, 1])
        nlambda <- 0
        K <- sqrt(2/D)
        agqws <- K %*% t(ghqws)
        agqns <- apply(K %*% t(ghqns), 2, function(x) x + u)
        grFac <- unclass(grFac)
        ranNew <- apply(agqns, 2, function(x) x[grFac] * stDev)
        eta1Tmp <- (eta1Fix + o1 - ranNew) / sigma
        eta2Tmp <- (eta2Fix + o2 - ranNew) / sigma
        if(nlambda)
            PRnn <- exp(weights * log(pfun(eta1Tmp, lambda) -
                                      pfun(eta2Tmp, lambda)))
        else
            PRnn <- exp(wts * log(pfun(eta1Tmp) - pfun(eta2Tmp)))

        dim(eta1Tmp)


        exp(wts[IND] * log(pfun(eta1Tmp[IND, ]) - pfun(eta2Tmp[IND, ])))

        PRrn <- do.call(rbind, lapply(1:dims$q, function(i) {
            apply(PRnn[grFac == i, ,drop = FALSE], 2, prod)
        }))
        head(PRrn)

        PRrn <- do.call(rbind, lapply(1:dims$q, function(i) {
            apply(PRnn[grFac == i, ,drop = FALSE], 2, function(x) sum(log(x)))
        }))
        head(PRrn)
        ## Could we do something like
        PRnn <- wts * log(pfun(eta1Tmp) - pfun(eta2Tmp))
        PRrn <- do.call(rbind, lapply(1:dims$q, function(i) {
            apply(PRnn[grFac == i, ,drop = FALSE], 2, function(x) sum(x))
        }))
        head(PRrn, 20)
        ## to avoid first exp()ing and then log()ing?
        head(exp(PRrn), 20)
        range(PRrn)
        exp(range(PRrn))

        out <- PRrn + log(agqws) + log(dnorm(x=agqns, mean=0, sd=1))


        log(2 * 3)
        log(2) + log(3)

        PRnn[grFac == 12, , drop=FALSE]
        IND <- which(grFac == 12)
        cbind(IND, wts[IND], PRnn[IND, ])

        dim(PRrn)
        ## There seems to be underfloow allready in the computations
        ## in PRnn, which propagates to PRrn
        PRrn <- PRrn * agqws * dnorm(x=agqns, mean=0, sd=1)
    })
    -sum(log(rowSums(rho$PRrn)))
}

getNAGQinC <- function(rho, par) {
### negative log-likelihood by adaptive Gauss-Hermite quadrature
### implemented in C:
    if(!missing(par))
        rho$par <- par
    if(!rho$updateU(rho)) return(Inf)
    if(any(rho$D < 0)) return(Inf)
    with(rho, {
        .C("getNAGQ",
           nll = double(1), ## nll
           as.integer(grFac), ## grFac
           as.double(stDev), ## stDev
           as.double(eta1Fix),
           as.double(eta2Fix),
           as.double(o1),
           as.double(o2),
           as.double(sigma), ## Sigma
           as.double(weights),
           length(sigma), ## nx - no. obs
           length(uStart), ## nu - no. re
           as.double(ghqns),
           as.double(log(ghqws)), ## lghqws
           as.double(ghqns^2), ## ghqns2
           as.double(u),
           as.double(D),
           as.integer(abs(nAGQ)),
           as.integer(linkInt),
           as.double(lambda))$nll
    })
}

getNGHQinR <- function(rho, par) {
### negative log-likelihood by standard Gauss-Hermite quadrature
### implemented in R:
  if(!missing(par))
    rho$par <- par
  .negLogLikBase(rho) ## Update lambda, stDev, sigma and eta*Fix
  with(rho, {
    ns <- ghqns * stDev
    SS <- numeric(r) ## summed likelihood
    for(i in 1:r) {
      ind <- grFac == i
      eta1Fi <- eta1Fix[ind]
      eta2Fi <- eta2Fix[ind]
      o1i <- o1[ind]
      o2i <- o2[ind]
      si <- sigma[ind]
      wt <- weights[ind]
      for(h in 1:abs(nAGQ)) {
        eta1s <- (eta1Fi + o1i - ns[h]) / si
        eta2s <- (eta2Fi + o2i - ns[h]) / si
        ## SS[i] <- exp(sum(wt * log(pfun(eta1s) - pfun(eta2s)))) *
        ##     ghqws[h] * exp(ghqns[h]^2) * dnorm(x=ghqns[h]) + SS[i]
        SS[i] <- exp(sum(wt * log(pfun(eta1s) - pfun(eta2s)))) *
          ghqws[h] + SS[i]
### FIXME: The fitted values could be evaluated with getFittedC for
### better precision.
      }
    }
    -sum(log(SS)) + r * log(2*pi)/2
  })
}

getNGHQinC <- function(rho, par) {
### negative log-likelihood by standard Gauss-Hermite quadrature
### implemented in C:
    if(!missing(par))
        rho$par <- par
    .negLogLikBase(rho) ## Update lambda, stDev, sigma and eta*Fix
    with(rho, {
        .C("getNGHQ",
           nll = double(1),
           as.integer(grFac),
           as.double(stDev),
           as.double(eta1Fix),
           as.double(eta2Fix),
           as.double(o1),
           as.double(o2),
           as.double(sigma),
           as.double(weights),
           length(sigma),
           length(uStart),
           as.double(ghqns),
           as.double(ghqws),
           as.integer(abs(nAGQ)),
           as.integer(linkInt),
           as.double(ghqns * stDev),
           as.double(lambda))$nll
    })
}

finalizeRhoM <- function(rhoM) {
    if(rhoM$method == "ucminf") {
        if(rhoM$optRes$info[1] > rhoM$optCtrl[["grtol"]])
            warning("clmm2 may not have converged:\n  optimizer 'ucminf' terminated with max|gradient|: ",
                    rhoM$optRes$info[1], call.=FALSE)
        rhoM$convergence <-
            ifelse(rhoM$optRes$info[1] > rhoM$optCtrl[["grtol"]], FALSE, TRUE)
    }
    if(rhoM$method == "nlminb") {
        rhoM$convergence <- ifelse(rhoM$optRes$convergence == 0, TRUE, FALSE)
        if(!rhoM$convergence)
            warning("clmm2 may not have converged:\n  optimizer 'nlminb' terminated with message: ",
                    rhoM$optRes$message, call.=FALSE)
    }
    if(rhoM$ctrl$gradTol < max(abs(rhoM$gradient)))
        warning("Inner loop did not converge at termination:\n  max|gradient| = ",
                max(abs(rhoM$gradient)))
    with(rhoM, {
        if(nxi > 0) {
            xi <- par[1:nxi]
            names(xi) <- xiNames
            thetaNames <- paste(lev[-length(lev)], lev[-1], sep="|")
            Alpha <- Theta <- matrix(par[1:nxi], nrow=ncolXX, byrow=TRUE)
            Theta <- t(apply(Theta, 1, function(x) c(tJac %*% x)))
            if(ncolXX > 1){
                dimnames(Theta) <- list(dnXX[[2]], thetaNames)
                dimnames(Alpha) <- list(dnXX[[2]], alphaNames)
            }
            else {
                Theta <- c(Theta)
                Alpha <- c(Alpha)
                names(Theta) <- thetaNames
                names(Alpha) <- alphaNames
            }
            coefficients <- xi
        }
        else coefficients <- numeric(0)
        if(p > 0) {
            beta <- par[nxi + 1:p]
            names(beta) <- dnX[[2]]
            coefficients <- c(coefficients, beta)
        }
        if(k > 0) {
            zeta <- par[nxi+p + 1:k]
            names(zeta) <- dnZ[[2]]
            coefficients <- c(coefficients, zeta)
        }
        if(estimLambda > 0) {
            names(lambda) <- "lambda"
            coefficients <- c(coefficients, lambda)
        }
        if(s > 0) {
            stDev <- exp(par[nxi+p+k + estimLambda + 1:s])
            coefficients <- c(coefficients, stDev)
        }
        names(stDev) <- randomName
        if(exists("Hessian", inherits=FALSE))
            dimnames(Hessian) <- list(names(coefficients),
                                      names(coefficients))
        edf <- p + nxi + k + estimLambda + s
        nobs <- sum(weights)
        fitted.values <- pr
        df.residual = nobs - edf
        ranef <- u * stDev
        condVar <- 1/D * stDev^2
        logLik <- -optRes[[2]]
    })
    res <- as.list(rhoM)
    keepNames <-
        c("ranef", "df.residual", "fitted.values", "edf", "start",
          "stDev", "beta", "coefficients", "zeta", "Alpha", "Theta",
          "xi", "lambda", "convergence", "Hessian",
          "gradient", "optRes", "logLik", "Niter", "grFac", "call",
          "scale", "location", "nominal", "method", "y", "lev",
          "nobs", "threshold", "estimLambda", "link", "nAGQ",
          "condVar", "contrasts", "na.action")
    m <- match(keepNames, names(res), 0)
    res <- res[m]
    res
}

anova.clmm2 <- function (object, ..., test = c("Chisq", "none"))
{
    anova.clm2(object, ..., test = c("Chisq", "none"))
}

print.clmm2 <- function(x, ...)
{
  if(x$nAGQ >= 2)
    cat(paste("Cumulative Link Mixed Model fitted with the adaptive",
              "Gauss-Hermite \nquadrature approximation with",
              x$nAGQ ,"quadrature points"), "\n\n")
  else if(x$nAGQ <= -1)
    cat(paste("Cumulative Link Mixed Model fitted with the",
              "Gauss-Hermite \nquadrature approximation with",
              abs(x$nAGQ) ,"quadrature points"), "\n\n")
  else
    cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
        fill=TRUE)
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$stDev)) {
        cat("\nRandom effects:\n")
        varMat <- matrix(c(x$stDev^2, x$stDev), nrow =
                         length(x$stDev), ncol=2)
        rownames(varMat) <- names(x$stDev)
        colnames(varMat) <- c("Var", "Std.Dev")
        print(varMat, ...)
    } else {
        cat("\nNo random effects\n")
    }
    if(length(x$beta)) {
        cat("\nLocation coefficients:\n")
        print(x$beta, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(length(x$zeta)) {
        cat("\nScale coefficients:\n")
        print(x$zeta, ...)
    } else {
        cat("\nNo Scale coefficients\n")
    }
    if(x$estimLambda > 0) {
        cat("\nLink coefficient:\n")
        print(x$lambda)
    }
    if(length(x$xi) > 0) {
        cat("\nThreshold coefficients:\n")
        print(x$Alpha, ...)
        if(x$threshold != "flexible") {
            cat("\nThresholds:\n")
            print(x$Theta, ...)
        }
    }
    cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
    cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    invisible(x)
}

vcov.clmm2 <- function(object, ...)
{
    if(is.null(object$Hessian)) {
        stop("Model needs to be fitted with Hess = TRUE")
    }
    dn <- names(object$coefficients)
    structure(solve(object$Hessian), dimnames = list(dn, dn))
}

summary.clmm2 <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
    estimStDev <- !("sdFixed" %in% names(as.list(object$call)))
    edf <- object$edf
    coef <- with(object,
                 matrix(0, edf-estimStDev, 4,
                        dimnames =
                         list(names(coefficients[seq_len(edf-estimStDev)]),
                        c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))))
    coef[, 1] <- object$coefficients[seq_len(edf-estimStDev)]
    if(is.null(object$Hessian)) {
      stop("Model needs to be fitted with Hess = TRUE")
    }
    vc <- try(vcov(object), silent = TRUE)
    if(class(vc) == "try-error") {
        warning("Variance-covariance matrix of the parameters is not defined")
        coef[, 2:4] <- NaN
        if(correlation) warning("Correlation matrix is unavailable")
        object$condHess <- NaN
    } else {
        sd <- sqrt(diag(vc))
        coef[, 2] <- sd[seq_len(edf - estimStDev)]
        object$condHess <-
            with(eigen(object$Hessian, only.values = TRUE),
                 abs(max(values) / min(values)))
        coef[, 3] <- coef[, 1]/coef[, 2]
        coef[, 4] <- 2*pnorm(abs(coef[, 3]), lower.tail=FALSE)
        if(correlation)
            object$correlation <- (vc/sd)/rep(sd, rep(object$edf, object$edf))
    }
    object$coefficients <- coef
    object$digits <- digits
    varMat <- matrix(c(object$stDev^2, object$stDev),
                     nrow = length(object$stDev), ncol=2)
    rownames(varMat) <- names(object$stDev)
    colnames(varMat) <- c("Var", "Std.Dev")
    object$varMat <- varMat
    class(object) <- "summary.clmm2"
    object
}

print.summary.clmm2 <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
    if(x$nAGQ >=2)
        cat(paste("Cumulative Link Mixed Model fitted with the adaptive",
                  "Gauss-Hermite \nquadrature approximation with",
                  x$nAGQ ,"quadrature points\n\n"))
    else if(x$nAGQ <= -1)
      cat(paste("Cumulative Link Mixed Model fitted with the",
                "Gauss-Hermite \nquadrature approximation with",
                abs(x$nAGQ) ,"quadrature points"), "\n\n")
    else
        cat("Cumulative Link Mixed Model fitted with the Laplace approximation\n",
            fill=TRUE)
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
    if(length(x$stDev)) {
        cat("\nRandom effects:\n")
        print(x$varMat, ...)
    } else {
        cat("\nNo random effects\n")
    }
### FIXME: Should the number of obs. and the number of groups be
### displayed as in lmer?
    coef <- format(round(x$coefficients, digits=digits))
    coef[,4] <- format.pval(x$coefficients[, 4])
    p <- length(x$beta); nxi <- length(x$xi)
    k <- length(x$zeta); u <- x$estimLambda
    if(p > 0) {
        cat("\nLocation coefficients:\n")
        print(coef[nxi + 1:p, , drop=FALSE],
              quote = FALSE, ...)
    } else {
        cat("\nNo location coefficients\n")
    }
    if(k > 0) {
      cat("\nScale coefficients:\n")
      print(coef[(nxi+p+1):(nxi+p+k), , drop=FALSE],
            quote = FALSE, ...)
    } else {
      cat("\nNo scale coefficients\n")
    }
    if(x$estimLambda > 0) {
        cat("\nLink coefficients:\n")
        print(coef[(nxi+p+k+1):(nxi+p+k+u), , drop=FALSE],
              quote = FALSE, ...)
    }
    if(nxi > 0) {
        cat("\nThreshold coefficients:\n")
        print(coef[1:nxi, -4, drop=FALSE], quote = FALSE, ...)
    }

    cat("\nlog-likelihood:", format(x$logLik, nsmall=2), "\n")
    cat("AIC:", format(-2*x$logLik + 2*x$edf, nsmall=2), "\n")
    cat("Condition number of Hessian:", format(x$condHess, nsmall=2), "\n")
    if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
    if(!is.null(correl <- x$correlation)) {
        cat("\nCorrelation of Coefficients:\n")
        ll <- lower.tri(correl)
        correl[ll] <- format(round(correl[ll], digits))
        correl[!ll] <- ""
        print(correl[-1, -ncol(correl)], quote = FALSE, ...)
    }
    invisible(x)
}

## ranef.clmm2 <- function(x) {
##     x$ranef
## }

## Trace <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
##     t1 <- sprintf(" %3d:  %-5e:    %.3f:   %1.3e:  ",
##                   iter, stepFactor, val, maxGrad)
##     t2 <- formatC(par)
##     if(first)
##         cat("iter:  step factor:     Value:     max|grad|:   Parameters:\n")
##     cat(t1, t2, "\n")
## }

gauss.hermite <- function (n)
{
    n <- as.integer(n)
    if (n < 0)
        stop("need non-negative number of nodes")
    if (n == 0)
        return(list(nodes = numeric(0), weights = numeric(0)))
    i <- 1:n
    i1 <- i[-n]
    muzero <- sqrt(pi)
    a <- rep(0, n)
    b <- sqrt(i1/2)

    A <- rep(0, n * n)
    A[(n + 1) * (i1 - 1) + 2] <- b
    A[(n + 1) * i1] <- b
    dim(A) <- c(n, n)
    vd <- eigen(A, symmetric = TRUE)
    w <- rev(as.vector(vd$vectors[1, ]))
    w <- muzero * w^2
    x <- rev(vd$values)
    list(nodes = x, weights = w)
}

profile.clmm2 <-
    function(fitted, alpha = 0.01, range, nSteps = 20, trace = 1, ...)
{
    estimStDev <- !("sdFixed" %in% names(as.list(fitted$call)))
    if(!estimStDev) ##  || is.null(fitted$Hessian))
        fitted <- update(fitted, Hess = TRUE, sdFixed = NULL)
    MLogLik <- fitted$logLik
    MLstDev <- fitted$stDev
    if(missing(range) && is.null(fitted$Hessian))
        stop("'range' should be specified or model fitted with 'Hess = TRUE'")
    if(missing(range) && !is.null(fitted$Hessian)) {
        range <- log(fitted$stDev) + qnorm(1 - alpha/2) *
            c(-1, 1) * sqrt(vcov(fitted)[fitted$edf, fitted$edf])
        range <- exp(range)
        pct <- paste(round(100*c(alpha/2, 1-alpha/2), 1), "%")
        ci <- array(NA, dim = c(1, 2),
                    dimnames = list("stDev", pct))
        ci[] <- range
    }
    stopifnot(all(range > 0))
    logLik <- numeric(nSteps)
    stDevSeq <- seq(min(range), max(range), length.out = nSteps)
    if(trace) message("Now profiling stDev with ", nSteps,
                      " steps: i =")
    if(trace) cat(1, "")
    rho <- update(fitted, Hess = FALSE, sdFixed = min(range))
    logLik[1] <- rho$logLik
    start <- as.vector(rho$coefficients)

    for(i in 2:nSteps){
        if(trace) cat(i, "")
        rho <- update(rho, sdFixed = stDevSeq[i], start = start)
        logLik[i] <- rho$logLik
        start <- as.vector(rho$coefficients)
    }
    if(trace) cat("\n")

    if(any(logLik > fitted$logLik))
        warning("Profiling found a better optimum,",
                "  so original fit had not converged")
    sgn <- 2*(stDevSeq > MLstDev) -1
    Lroot <- sgn * sqrt(2) * sqrt(-logLik + MLogLik)
    res <- data.frame("Lroot" = c(0, Lroot),
                      "stDev" = c(MLstDev, stDevSeq))
    res <- res[order(res[,1]),]
    if(!all(diff(res[,2]) > 0))
        warning("likelihood is not monotonically decreasing from maximum,\n",
                "  so profile may be unreliable for stDev")
    val <- structure(list(stDev = res), original.fit = fitted)
    if(exists("ci", inherits=FALSE)) attr(val, "WaldCI") <- ci
    class(val) <- c("profile.clmm2", "profile")
    val
}

confint.profile.clmm2 <-
    function(object, parm = seq_along(Pnames), level = 0.95, ...)
{
    Pnames <- names(object)
    confint.profile.clm2(object, parm = parm, level = level, ...)
}

plot.profile.clmm2 <-
    function(x, parm = seq_along(Pnames), level = c(0.95, 0.99),
             Log = FALSE, relative = TRUE, fig = TRUE, n = 1e3, ...,
             ylim = NULL)
{
    Pnames <- names(x)
    plot.profile.clm2(x, parm = parm, level = level, Log = Log,
                     relative = relative, fig = fig,
                     n = n, ...,  ylim = ylim)
}

update.clmm2 <-
    function(object, formula., location, scale, nominal, ...,
             evaluate = TRUE)
{
    call <- object$call
    if (is.null(call))
        stop("need an object with call component")
    extras <- match.call(expand.dots = FALSE)$...
    if (!missing(location))
        call$location <-
            update.formula(formula(attr(object$location, "terms")),
                                   location)
    if (!missing(scale))
        call$scale <-
            if(!is.null(object$scale))
                update.formula(formula(attr(object$scale, "terms")), scale)
            else
                scale

    if (!missing(nominal))
        call$nominal <-
            if(!is.null(object$nominal))
                update.formula(formula(attr(object$nominal, "terms")), nominal)
            else
                nominal

    if (length(extras)) {
        existing <- !is.na(match(names(extras), names(call)))
        for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
        if (any(!existing)) {
            call <- c(as.list(call), extras[!existing])
            call <- as.call(call)
        }
    }
    if (evaluate)
        eval(call, parent.frame())
    else call
}

