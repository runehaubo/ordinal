## This file contains:
## An alternate (and older) implementation of CLMs in clm2(). The new
## and recommended implementation is available in clm(), cf. ./R/clm.R

clm2.control <-
    function(method = c("ucminf", "Newton", "nlminb", "optim",
             "model.frame"), ..., convTol = 1e-4,
             trace = 0, maxIter = 100, gradTol = 1e-5,
             maxLineIter = 10)
{
    method <- match.arg(method)
    ctrl <-
        if(method == "Newton")
            list(convTol=convTol, trace=trace, maxIter=maxIter,
                 gradTol=gradTol, maxLineIter=maxLineIter)
        else
            list(trace = abs(trace), ...)

    if(!all(is.numeric(c(maxIter, gradTol, maxLineIter, convTol))))
        stop("maxIter, gradTol, maxLineIter, convTol should all be numeric")
    if(convTol <= 0)
        stop("convTol should be > 0")
    if(method == "ucminf" && !"grtol" %in% names(ctrl))
        ctrl$grtol <- gradTol
##     if(method == "ucminf" && convTol > ctrl$grtol)
##         stop("convTol should be <= grtol/gradTol")
##     if(method == "Newton" && convTol > gradTol)
##         stop("convTol should be <= gradTol")

    list(method = method, convTol = convTol, ctrl = ctrl)
}

newRho <- function(parent, XX, X, Z, y, weights, Loffset, Soffset, ## OK
                   link, lambda, theta, threshold, Hess, control)
### FIXME: Could remove theta argument?
{
    rho <- new.env(parent = parent)
    rho$X <- X
    rho$dnX <- dimnames(X)
    dimnames(rho$X) <- NULL
    rho$Z <- Z
    rho$dnZ <- dimnames(Z)
    dimnames(rho$Z) <- NULL
    rho$weights <- weights
    rho$Loffset <- Loffset
    rho$expSoffset <- rho$sigma <- exp(Soffset)
    rho$Hess <- ifelse(Hess, 1L, 0L)
    rho$method <- control$method
    rho$convTol <- control$convTol
    rho$ctrl <- control$ctrl

    rho$pfun <- switch(link,
                       logistic = plogis,
                       probit = pnorm,
                       cloglog = function(x) pgumbel(x, max=FALSE),
                       cauchit = pcauchy,
                       loglog = pgumbel,
                       "Aranda-Ordaz" = function(x, lambda) pAO(x, lambda),
                       "log-gamma" = function(x, lambda) plgamma(x, lambda))
    rho$dfun <- switch(link,
                       logistic = dlogis,
                       probit = dnorm,
                       cloglog = function(x) dgumbel(x, max=FALSE),
                       cauchit = dcauchy,
                       loglog = dgumbel,
                       "Aranda-Ordaz" = function(x, lambda) dAO(x, lambda),
                       "log-gamma" = function(x, lambda) dlgamma(x, lambda))
    rho$gfun <- switch(link,
                       logistic = glogis,
                       probit = function(x) -x * dnorm(x),
                       cloglog = function(x) ggumbel(x, max=FALSE),
                       cloglog = ggumbel,
                       cauchit = gcauchy,
                       "Aranda-Ordaz" = function(x, lambda) gAO(x, lambda), ## shouldn't happen
                       "log-gamma" = function(x, lambda) glgamma(x, lambda)
                       )
    rho$link <- link
    rho$linkInt <- switch(link,
                          logistic = 1L,
                          probit = 2L,
                          cloglog = 3L,
                          loglog = 4L,
                          cauchit = 5L,
                          "Aranda-Ordaz" = 6L,
                          "log-gamma" = 7L)
    rho$estimLambda <- ifelse(link %in% c("Aranda-Ordaz", "log-gamma") &&
                              is.null(lambda), 1L, 0L)
    rho$nlambda <- 0L
    rho$lambda <-
        if(!is.null(lambda)) lambda
        else 1
    if(link %in% c("Aranda-Ordaz", "log-gamma"))
        rho$nlambda <- 1L

    if(rho$estimLambda > 0 & rho$link == "Aranda-Ordaz" &
       rho$method != "nlminb"){
        message("Changing to nlminb optimizer to accommodate optimization with bounds")
        m <- match( names(rho$ctrl), "grtol", 0)
        rho$ctrl <- rho$ctrl[!m]
        rho$method <- "nlminb"
    }
    if(rho$method == "nlminb") {
        rho$limitUp <- Inf
        rho$limitLow <- -Inf
    }

    rho$n <- n <- length(y)
    rho$p <- ifelse(missing(X), 0, ncol(X))
    rho$k <- ifelse(missing(Z), 0, ncol(Z))
    rho$y <- y
    rho$threshold <- threshold
    rho$ncolXX <- ncol(XX)
    rho$dnXX <- dimnames(XX)
    rho$lev <- levels(y)
    rho$ntheta <- nlevels(y) - 1
    rho$B2 <- 1 * (col(matrix(0, n, rho$ntheta + 1)) == c(unclass(y)))
### Setting elements of o[12] to [+-]Inf cause problems in
### getGnll and clmm-related functions because 1) 0*Inf = NaN, while
### 0*large.value = 0, so several computations have to be handled
### specially and 2) Inf-values are not by default allowed in .C calls
### and all specials would have to be handled separately.
    ## o1 <- B2[, rho$ntheta + 1, drop = TRUE]
    ## o1[o1 == 1] <- Inf
    ## rho$o1 <- o1 - rho$Loffset
    ## o2 <- B2[,1, drop = TRUE]
    ## o2[o2 == 1] <- -Inf
    ## rho$o2 <- o2 - rho$Loffset
    inf.value <- 1e5
    rho$o1 <- c(inf.value * rho$B2[, rho$ntheta + 1]) - rho$Loffset
    rho$o2 <- c(-inf.value * rho$B2[,1]) - rho$Loffset
    rho$B1 <- rho$B2[,-(rho$ntheta + 1), drop = FALSE]
    rho$B2 <- rho$B2[,-1, drop = FALSE]
    makeThresholds2(rho, threshold)
    rho$B1 <- rho$B1 %*% rho$tJac
    rho$B2 <- rho$B2 %*% rho$tJac
    rho$xiNames <- rho$alphaNames
    rho$nxi <- rho$nalpha * rho$ncolXX
    if(rho$ncolXX > 1) { ## test actually not needed
        rho$xiNames <- paste(rep(rho$alphaNames, rho$ncolXX), ".",
                             rep(colnames(XX), each=rho$nalpha), sep="")
        LL1 <- lapply(1:rho$ncolXX, function(x) rho$B1 * XX[,x])
        rho$B1 <- do.call(cbind, LL1)
        LL2 <- lapply(1:rho$ncolXX, function(x) rho$B2 * XX[,x])
        rho$B2 <- do.call(cbind, LL2)
    }
    if(rho$p > 0) {
        rho$B1 <- cbind(rho$B1, -X)
        rho$B2 <- cbind(rho$B2, -X)
    }
    dimnames(rho$B1) <- NULL
    dimnames(rho$B2) <- NULL
    return(rho)
} # populates the rho environment

setStart <- function(rho) ## Ok
{ ## set starting values in the rho environment
    ## try logistic/probit regression on 'middle' cut
    q1 <- max(1, rho$ntheta %/% 2)
    y1 <- (c(unclass(rho$y)) > q1)
    x <- cbind(Intercept = rep(1, rho$n), rho$X)
    fit <-
        switch(rho$link,
               "logistic"= glm.fit(x, y1, rho$weights, family = binomial(), offset = rho$Loffset),
               "probit" = glm.fit(x, y1, rho$weights, family = binomial("probit"), offset = rho$Loffset),
               ## this is deliberate, a better starting point
               "cloglog" = glm.fit(x, y1, rho$weights, family = binomial("probit"), offset = rho$Loffset),
               "loglog" = glm.fit(x, y1, rho$weights, family = binomial("probit"), offset = rho$Loffset),
               "cauchit" = glm.fit(x, y1, rho$weights, family = binomial("cauchit"), offset = rho$Loffset),
               "Aranda-Ordaz" = glm.fit(x, y1, rho$weights, family = binomial("probit"), offset = rho$Loffset),
               "log-gamma" = glm.fit(x, y1, rho$weights, family = binomial("probit"), offset = rho$Loffset))
    if(!fit$converged)
        stop("attempt to find suitable starting values failed")
    coefs <- fit$coefficients
    if(any(is.na(coefs))) {
        warning("design appears to be rank-deficient, so dropping some coefs")
        keep <- !is.na(coefs)
        coefs <- coefs[keep]
        rho$X <- rho$X[, keep[-1], drop = FALSE]
        rho$dnX[[2]] <- rho$dnX[[2]][keep[-1]]
        rho$B1 <- rho$B1[, c(rep(TRUE, rho$nxi), keep[-1]), drop = FALSE]
        rho$B2 <- rho$B2[, c(rep(TRUE, rho$nxi), keep[-1]), drop = FALSE]
        rho$p <- ncol(rho$X)
    }
    ## Intercepts:
    spacing <- qlogis((1:rho$ntheta)/(rho$ntheta+1)) # just a guess
    if(rho$link != "logit") spacing <- spacing/1.7
    ## if(rho$threshold == "flexible") # default
    alphas <- -coefs[1] + spacing - spacing[q1]
    if(rho$threshold == "symmetric" && rho$ntheta %% 2) ## ntheta odd
        alphas <- c(alphas[q1+1],cumsum(rep(spacing[q1+2], rho$nalpha-1)))
    if(rho$threshold == "symmetric" && !rho$ntheta %% 2) ## ntheta even
        alphas <- c(alphas[q1:(q1+1)], cumsum(rep(spacing[q1+1], rho$nalpha-2)))
    if(rho$threshold == "symmetric2" && rho$ntheta %% 2) ## ntheta odd
        alphas <- cumsum(rep(spacing[q1+2], rho$nalpha-1))
    if(rho$threshold == "symmetric2" && !rho$ntheta %% 2) ## ntheta even
        alphas <- cumsum(rep(spacing[q1+1], rho$nalpha-2))
    if(rho$threshold == "equidistant")
        alphas <- c(alphas[1], mean(diff(spacing)))
    ## initialize nominal effects to zero:
    if(rho$ncolXX > 1) {
        xi <- c(alphas, rep(rep(0, rho$nalpha), rho$ncolXX-1))
        stopifnot(length(xi) == rho$nalpha * rho$ncolXX)}
    else xi <- alphas
    if(rho$estimLambda > 0){
        rho$lambda <- 1
        names(rho$lambda) <- "lambda"
    }
    start <- c(xi, coefs[-1], rep(0, rho$k), rep(1, rho$estimLambda))
    names(start) <- NULL
    rho$start <- rho$par <- start
}

getPar <- function(rho) rho$par ## OK

getNll <- function(rho, par)  { ## ok
  if(!missing(par))
    rho$par <- par
  with(rho, {
    if(estimLambda > 0)
      lambda <- par[nxi + p + k + 1:estimLambda]
    sigma <-
      if(k > 0) expSoffset * exp(drop(Z %*% par[nxi+p + 1:k]))
      else expSoffset
    eta1 <- (drop(B1 %*% par[1:(nxi + p)]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:(nxi + p)]) + o2)/sigma
    pr <-
      if(nlambda) pfun(eta1, lambda) - pfun(eta2, lambda)
      else pfun(eta1) - pfun(eta2)
    if(all(is.finite(pr)) && all(pr > 0)) -sum(weights * log(pr))
    else Inf
  })
}

getGnll <- function(rho, par)  { ## ok
  if(!missing(par))
    rho$par <- par
  with(rho, {
    if(estimLambda > 0)
      lambda <- par[nxi + p + k + 1:estimLambda]
    sigma <-
      if(k > 0) expSoffset * exp(drop(Z %*% par[nxi+p + 1:k]))
      else expSoffset
    eta1 <- (drop(B1 %*% par[1:(nxi + p)]) + o1)/sigma
    eta2 <- (drop(B2 %*% par[1:(nxi + p)]) + o2)/sigma
    if(nlambda) {
      pr <- pfun(eta1, lambda) - pfun(eta2, lambda)
      p1 <- dfun(eta1, lambda)
      p2 <- dfun(eta2, lambda)
    }
    else {
      pr <- pfun(eta1) - pfun(eta2)
      p1 <- dfun(eta1)
      p2 <- dfun(eta2)
    }
    prSig <- pr * sigma
    ## eta1 * p1 is complicated because in theory eta1 contains
    ## Inf(-Inf) where p1 contains 0 and 0 * Inf = NaN...
    ## eta.p1 <- ifelse(p1 == 0, 0, eta1 * p1)
    ## eta.p2 <- ifelse(p2 == 0, 0, eta2 * p2)
    gradSigma <-
      ## if(k > 0) crossprod(Z, weights * (eta.p1 - eta.p2)/pr)
      if(k > 0) crossprod(Z, weights * (eta1 * p1 - eta2 * p2)/pr)
      else numeric(0)
    gradThetaBeta <-
      if(nxi > 0) -crossprod((B1*p1 - B2*p2), weights/prSig)
      else -crossprod((X * (p2 - p1)), weights/prSig)
    grad <-
    ##   if (all(is.finite(pr)) && all(pr > 0))
    ##     c(gradThetaBeta, gradSigma)
    ##   else rep(Inf, nxi + p + k)
      c(gradThetaBeta, gradSigma)
  })
  if(rho$estimLambda > 0)
    c(rho$grad, grad.lambda(rho, rho$lambda, rho$link))
  else
    rho$grad
}

getHnll <- function(rho, par)  { ## ok
    if(!missing(par))
        rho$par <- par
    with(rho, {
        eta1 <- drop(B1 %*% par[1:(nxi + p)]) + o1
        eta2 <- drop(B2 %*% par[1:(nxi + p)]) + o2
        pr <- pfun(eta1) - pfun(eta2)
        p1 <- dfun(eta1)
        p2 <- dfun(eta2)
        g1 <- gfun(eta1)
        g2 <- gfun(eta2)
        wtpr <- weights/pr

        dS.psi <- -crossprod(B1 * g1*wtpr, B1) +
            crossprod(B2 * g2*wtpr, B2)
        dpi.psi <- B1 * p1 - B2 * p2
###    dS.pi <- dpi.psi * wtpr/pr
        if (all(pr > 0))
            dS.psi + crossprod(dpi.psi, (dpi.psi * wtpr/pr))
        else array(NA, dim = c(nxi + p, nxi + p))
    })
}

.negLogLik <- function(rho) { ## negative log-likelihood ## OK
    with(rho, {
        eta1 <- drop(B1 %*% par[1:(nxi + p)]) + o1
        eta2 <- drop(B2 %*% par[1:(nxi + p)]) + o2
        pr <- pfun(eta1) - pfun(eta2)
        if (all(pr > 0))
            -sum(weights * log(pr))
        else Inf
    })
}

.grad <- function(rho) { ## gradient of the negative log-likelihood ## OK
    with(rho, {
        p1 <- dfun(eta1)
        p2 <- dfun(eta2)
        wtpr <- weights/pr
        if (all(pr > 0))
            -crossprod((B1 * p1 - B2 * p2), wtpr)
        else rep(NA, nalpha + p)
    })
}

.hessian <- function(rho) { ## hessian of the negative log-likelihood ## OK
    with(rho, {
        dS.psi <- crossprod(B1 * gfun(eta1)*wtpr, B1) -
            crossprod(B2 * gfun(eta2)*wtpr, B2)
        dpi.psi <- B1 * p1 - B2 * p2
        if (all(pr > 0))
            -dS.psi + crossprod(dpi.psi, (dpi.psi * wtpr/pr))
        else array(NA, dim = c(nxi+p, nxi+p))
    })
}

fitNR <- function(rho) ## OK
{
    ctrl <- rho$ctrl
    stepFactor <- 1
    innerIter <- 0
    conv <- 1  ## Convergence flag
    message <- "iteration limit reached"
    rho$negLogLik <- .negLogLik(rho)
    if(rho$negLogLik == Inf)
        stop("Non-finite log-likelihood at starting value")
    rho$gradient <- .grad(rho)
    maxGrad <- max(abs(rho$gradient))
    if(ctrl$trace > 0)
        Trace(iter=0, stepFactor, rho$negLogLik, maxGrad, rho$par, first=TRUE)

    ## Newton-Raphson algorithm:
    for(i in 1:ctrl$maxIter) {
        if(maxGrad < ctrl$gradTol) {
            message <- "max|gradient| < tol, so current iterate is probably solution"
            if(ctrl$trace > 0)
                cat("\nOptimizer converged! ", "max|grad|:",
                    maxGrad, message, fill = TRUE)
            conv <- 0
            break
        }
        rho$Hessian <- .hessian(rho)
        ## step <- .Call("La_dgesv", rho$Hessian, rho$gradient, .Machine$double.eps,
        ##               PACKAGE = "base") ## solve H*step = g for 'step'
        step <- as.vector(solve(rho$Hessian, rho$gradient))
        rho$par <- rho$par - stepFactor * step
        negLogLikTry <- .negLogLik(rho)
        lineIter <- 0
        ## simple line search, i.e. step halfing:
        while(negLogLikTry > rho$negLogLik) {
            stepFactor <- stepFactor/2
            rho$par <- rho$par + stepFactor * step
            negLogLikTry <- .negLogLik(rho)
            lineIter <- lineIter + 1
            if(ctrl$trace > 0)
                Trace(i+innerIter, stepFactor, rho$negLogLik, maxGrad,
                      rho$par, first=FALSE)
            if(lineIter > ctrl$maxLineIter){
                message <- "step factor reduced below minimum"
                conv <- 2
                break
            }
            innerIter <- innerIter + 1
        }
        rho$negLogLik <- negLogLikTry
        rho$gradient <- .grad(rho)
        maxGrad <- max(abs(rho$gradient))
        if(ctrl$trace > 0)
            Trace(iter=i+innerIter, stepFactor, rho$negLogLik,
                  maxGrad, rho$par, first=FALSE)
        stepFactor <- min(1, 2 * stepFactor)
    }
    if(conv > 0)
        if(ctrl$trace > 0) cat(message, fill = TRUE)
    ## Save info
    rho$optRes$niter <- c(outer = i, inner = innerIter)
    rho$logLik <- -rho$negLogLik
    rho$maxGradient <- maxGrad
    rho$gradient <- as.vector(rho$gradient)
    rho$Hessian <- .hessian(rho)
    rho$optRes$message <- message
    rho$optRes$convergence <- conv
}


fitCLM <- function(rho) { ## OK
    if(rho$method == "Newton") {
        if(rho$k != 0)
            stop("Newton scheme not implemented for models with scale")
        if(rho$ncolXX > 1)
            stop("Newton scheme not implemented for models with nominal effects")
        if(rho$link %in% c("Aranda-Ordaz", "log-gamma"))
            stop("Newton scheme not implemented for models with",
                 rho$link, "link function")
        fitNR(rho)
        return(invisible())
    }
    optRes <- switch(rho$method,
                  "nlminb" = nlminb(getPar(rho), function(par)
                  getNll(rho, par), function(par) getGnll(rho, par),
                  control=rho$ctrl, lower = rho$limitLow,
                  upper = rho$limitUp),
                  "ucminf" = ucminf(getPar(rho), function(par)
                  getNll(rho, par), function(par) getGnll(rho, par),
                  control=rho$ctrl),
                  "optim" = optim(getPar(rho), function(par)
                  getNll(rho, par), function(par) getGnll(rho, par),
                  method="BFGS", control=rho$ctrl),
                  )
    rho$par <- optRes[[1]]
    rho$logLik <- - getNll(rho, optRes[[1]])
    rho$optRes <- optRes

    rho$gradient <- c(getGnll(rho))
    rho$maxGradient <- max(abs(rho$gradient))
    if(rho$maxGradient > rho$convTol)
        warning("clm2 may not have converged:\n  optimizer ", rho$method,
                " terminated with max|gradient|: ", rho$maxGradient,
                call.=FALSE)
    return(invisible())
}

finalizeRho <- function(rho) { ## OK
    if(rho$method != "Newton") {
        rho$gradient <- c(getGnll(rho))
        rho$maxGradient <- max(abs(rho$gradient))
        rho$par <- rho$optRes[[1]]
        if(rho$Hess) {
            if(rho$k > 0 || rho$threshold != "flexible" ||
               rho$ncolXX > 1 || rho$nlambda > 0) {
                if(rho$link == "Aranda-Ordaz" &&
                   rho$estimLambda > 0 && rho$lambda < 1e-3)
                    message("Cannot get Hessian because lambda = ",rho$lambda
                            ," is too close to boundary.\n",
                            " Fit model with link == 'logistic' to get Hessian")
                else {
                    rho$Hessian <- myhess(function(par) getNll(rho, par),
                                           rho$par)
                    getNll(rho, rho$optRes[[1]]) # to reset the variables:
                                        # (par, pr)
                }
            }
            else
                rho$Hessian <- getHnll(rho, rho$optRes[[1]])
        }
    }
    rho$convergence <-
        ifelse(rho$maxGradient > rho$convTol, FALSE, TRUE)

    with(rho, {
        if(nxi > 0) {
            xi <- par[seq_len(nxi)]
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
        names(gradient) <- names(coefficients)
        edf <- p + nxi + k + estimLambda
        nobs <- sum(weights)

        fitted.values <- pr
        df.residual <- nobs - edf
        if(exists("Hessian", inherits=FALSE)) {
            dimnames(Hessian) <- list(names(coefficients),
                                      names(coefficients))
        }
    })
    res <- as.list(rho)
    keepNames <-
        c("df.residual", "fitted.values", "edf", "start",
          "beta", "coefficients", "zeta", "Alpha", "Theta",
          "xi", "lambda", "convergence", "Hessian", "convTol",
          "gradient", "optRes", "logLik", "call",
          "scale", "location", "nominal", "method", "y", "lev",
          "nobs", "threshold", "estimLambda", "link",
          "contrasts", "na.action")
    m <- match(keepNames, names(res), 0)
    res <- res[m]
    res
}

clm2 <- ## OK
  function(location, scale, nominal, data, weights, start, subset,
           na.action, contrasts, Hess = TRUE, model = TRUE,
           link = c("logistic", "probit", "cloglog", "loglog",
           "cauchit", "Aranda-Ordaz", "log-gamma"), lambda,
           doFit = TRUE, control,
           threshold = c("flexible", "symmetric", "equidistant"), ...)
{
    L <- match.call(expand.dots = FALSE)
    if(missing(location))
        stop("Model needs a specification of the location")
    if(missing(lambda)) lambda <- NULL
    if(missing(contrasts)) contrasts <- NULL
    link <- match.arg(link)
    if(!(link %in% c("Aranda-Ordaz", "log-gamma")) & !is.null(lambda)){
        warning("lambda ignored with link ", link)
        lambda <- NULL
    }
    if(!is.null(lambda) & length(lambda) > 1) {
        lambda <- lambda[1]
        warning("lambda is ", length(lambda),
                " long. Only the first element ", lambda[1], " is used")
    }
    if(!is.null(lambda) & link == "Aranda-Ordaz")
        if(lambda < 1e-6)
            stop("lambda has to be positive and lambda < 1e-6 not allowed for numerical reasons. lambda = ",
                 lambda, " was supplied.")
    if (missing(control)) control <- clm2.control(...)
    if(!setequal(names(control), c("method", "convTol", "ctrl")))
        stop("specify 'control' via clm2.control()")

    if (missing(data)) L$data <- environment(location)
    if (is.matrix(eval.parent(L$data)))
        L$data <- as.data.frame(L$data)

### Collect variables in location, scale and nominal formulae in a
### single formula, evaluate the model.frame and get index of row
### names for the rows to keep in the individual model.frames:
    m <- match(c("location", "scale", "nominal"), names(L), 0)
    F <- lapply(as.list(L[m]), eval.parent) ## evaluate in parent
    ## frame to allow 'f <- formula(sureness ~ prod); clm2(f, ...)'
    varNames <- unique(unlist(lapply(F, all.vars)))
    longFormula <-
        eval(parse(text = paste("~", paste(varNames, collapse = "+")))[1])
    m <- match(c("location", "data", "subset", "weights",
                  "na.action"), names(L), 0)
    L0 <- L[c(1, m)]
    if(!missing(scale) || !missing(nominal))
        L0$location <- longFormula
    L0$drop.unused.levels <- TRUE
    L0[[1]] <- as.name("model.frame")
    names(L0)[names(L0) == "location"] <- "formula"
    L0 <- eval.parent(L0)
    m <- match(c("location", "scale", "nominal", "data", "subset",
                 "weights", "na.action"), names(L), 0)
    L <- L[c(1, m)]
    L$drop.unused.levels <- TRUE
    L[[1]] <- as.name("model.frame")
    S <- L ## L: Location, S: Scale
    L$scale <- L$nominal <- NULL
    names(L)[names(L) == "location"] <- "formula"
    L <- eval.parent(L)
    keep <- match(rownames(L0), rownames(L))
    L <- L[keep, , drop = FALSE]
    TermsL <- attr(L, "terms")

### format response:
    y <- model.response(L)
    if(!is.factor(y))
        stop("response needs to be a factor")

### format thresholds:
    threshold <- match.arg(threshold)

### format location:
    X <- model.matrix(TermsL, L, contrasts)
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if (Xint > 0) X <- X[, -Xint, drop = FALSE]
    else warning("an intercept is needed and assumed in the location")
    n <- nrow(X)
    if(is.null(wt <- model.weights(L))) wt <- rep(1, n)
    if(is.null(Loffset <- model.offset(L))) Loffset <- rep(0, n)

### Format nominal:
    if(!missing(nominal)) {
        Nom <- S
        Nom$location <- Nom$scale <- NULL
        names(Nom)[names(Nom) == "nominal"] <- "formula"
        Nom <- eval.parent(Nom)
        Nom <- Nom[match(rownames(L0), rownames(Nom)), ,drop=FALSE]
        TermsNom <- attr(Nom, "terms")
        XX <- model.matrix(TermsNom, Nom)## , contrasts)
### Not allowing other than treatment contrasts in nominal
        if(is.null(Noffset <- model.offset(Nom))) Noffset <- rep(0, n)
        Nint <- match("(Intercept)", colnames(XX), nomatch = 0)
        if(Nint != 1)
            stop("An intercept is needed in the nominal formula")
### Are there any requirements about the presence of an
### intercept in the nominal formula?
    }
    else
        XX <- array(1, dim=c(n, 1))

### format scale:
    if(!missing(scale)) {
        S$location <- S$nominal <- NULL
        names(S)[names(S) == "scale"] <- "formula"
        S <- eval.parent(S)
        S <- S[match(rownames(L0), rownames(S)), ,drop=FALSE]
        TermsS <- attr(S, "terms")
### Should contrasts be allowed for the scale?
        Z <- model.matrix(TermsS, S, contrasts)
        Zint <- match("(Intercept)", colnames(Z), nomatch = 0)
        if(Zint > 0) Z <- Z[, -Zint, drop = FALSE]
        else warning("an intercept is needed and assumed in the scale")
        if(is.null(Soffset <- model.offset(S))) Soffset <- rep(0, n)
        if(ncol(Z) > 0 && n != nrow(Z)) # This shouldn't happen
            stop("Model needs same dataset in location and scale")
    } else if(missing(scale) && !is.factor(y)){
        Z <- array(1, dim = c(n, 1))
        Soffset <- rep(0, n)
    } else {
        Z <- array(dim = c(n, 0))
        Soffset <- rep(0, n)
    }

### return model.frame?
    if(control$method == "model.frame") {
        mf <- list(location = L)
        if(!missing(scale)) mf$scale <- S
        if(!missing(nominal)) mf$nominal <- Nom
        return(mf)
    }

### initialize and populate rho environment:
    rho <- newRho(parent.frame(), XX = XX, X=X, Z=Z, y=y, weights=wt,
                  Loffset=Loffset, Soffset=Soffset, link=link,
                  lambda = lambda, threshold=threshold,
                  Hess = Hess, control = control)

### get starting values:
    if(missing(start))
        setStart(rho)
    else
        rho$start <- rho$par <- start
    if(rho$estimLambda > 0 & rho$link == "Aranda-Ordaz")
        rho$limitLow <- c(rep(-Inf, length(rho$par)-1), 1e-5)
    if(length(rho$start) != with(rho, nxi + p + k + estimLambda))
        stop("'start' is not of the correct length")
### FIXME: Better check of increasing thresholds when ncol(XX) > 0
    if(ncol(XX) == 0) {
        if(!all(diff(c(rho$tJac %*% rho$start[1:rho$nalpha])) > 0))
            stop("Threshold starting values are not of increasing size")
    }
    if(!getNll(rho) < Inf)
        stop("Non-finite log-likelihood at starting values")
    if(model) {
        rho$location <- L
        if(!missing(scale)) rho$scale <- S
        if(!missing(nominal)) rho$nominal <- Nom
    }

### fit the model:
    if(!doFit)
        return(rho)
    fitCLM(rho)
    res <- finalizeRho(rho)

### add to output:
    res$call <- match.call()
    res$na.action <- attr(L0, "na.action")
    res$contrasts <- contrasts
    class(res) <- "clm2"
    res
}

print.clm2 <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
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

vcov.clm2 <- function(object, ...)
{
  if(is.null(object$Hessian)) {
    message("\nRe-fitting to get Hessian\n")
    utils::flush.console()
    object <- update(object, Hess=TRUE, start=object$coefficients)
  }
  dn <- names(object$coefficients)
  H <- object$Hessian
  ## To handle NaNs in the Hessian resulting from parameter
  ## unidentifiability:
  if(any(His.na <- !is.finite(H))) {
    H[His.na] <- 0
    VCOV <- ginv(H)
    VCOV[His.na] <- NaN
  }
  else
    VCOV <- ginv(H)
  structure(VCOV, dimnames = list(dn, dn))
}

summary.clm2 <- function(object, digits = max(3, .Options$digits - 3),
                        correlation = FALSE, ...)
{
  if(is.null(object$Hessian))
    stop("Model needs to be fitted with Hess = TRUE")
  coef <- matrix(0, object$edf, 4,
                 dimnames = list(names(object$coefficients),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coef[, 1] <- object$coefficients
  vc <- try(vcov(object), silent = TRUE)
  if(class(vc) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coef[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    coef[, 2] <- sd <- sqrt(diag(vc))
    ## Cond is Inf if Hessian contains NaNs:
    object$condHess <-
      if(any(is.na(object$Hessian))) Inf
      else with(eigen(object$Hessian, only.values = TRUE),
                abs(max(values) / min(values)))
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2*pnorm(abs(coef[, 3]), lower.tail=FALSE)
    if(correlation)
        object$correlation <- (vc/sd)/rep(sd, rep(object$edf, object$edf))
  }
  object$coefficients <- coef
  object$digits <- digits
    class(object) <- "summary.clm2"
  object
}

print.summary.clm2 <- function(x, digits = x$digits, signif.stars =
                              getOption("show.signif.stars"), ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl, control=NULL)
    }
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
        print(coef[seq_len(nxi), -4, drop=FALSE], quote = FALSE, ...)
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

anova.clm2 <- function (object, ..., test = c("Chisq", "none"))
{
  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0)
    stop('anova is not implemented for a single "clm2" object')
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df.residual)
  s <- order(dflis, decreasing = TRUE)
  mlist <- mlist[s]
  if (any(!sapply(mlist, inherits, "clm2")))
    stop('not all objects are of class "clm2"')
  ns <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(ns != ns[1]))
    stop("models were not all fitted to the same size of dataset")
  rsp <- unique(sapply(mlist, function(x) {
                       tmp <- attr(x$location, "terms")
                       class(tmp) <- "formula"
                       paste(tmp[2]) } ))
  mds <- sapply(mlist, function(x) {
      tmp1 <- attr(x$location, "terms")
      class(tmp1) <- "formula"
      if(!is.null(x$scale)) {
          tmp2 <- attr(x$scale, "terms")
          class(tmp2) <- "formula"
          tmp2 <- tmp2[2]
      }
      else tmp2 <- ""
      if(!is.null(x$nominal)) {
          tmp3 <- attr(x$nominal, "terms")
          class(tmp3) <- "formula"
          tmp3 <- tmp3[2]
      }
      else tmp3 <- ""
      paste(tmp1[3], "|", tmp2, "|", tmp3) } )
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) -2*x$logLik)
  tss <- c("", paste(1:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1], df[-1]))
  out <- data.frame(Model = mds, Resid.df = dfs, '-2logLik' = lls,
                    Test = tss, Df = df, LRtest = x2, Prob = pr)
  names(out) <- c("Model", "Resid. df", "-2logLik", "Test",
                  "   Df", "LR stat.", "Pr(Chi)")
  if (test == "none") out <- out[, 1:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of cumulative link models\n",
      paste("Response:", rsp))
  out
}

predict.clm2 <- function(object, newdata, ...)
{
    if(!inherits(object, "clm2")) stop("not a \"clm2\" object")
    if(missing(newdata)) pr <- object$fitted
    else {
        newdata <- as.data.frame(newdata)
        Terms <- attr(object$location, "terms")
        m <- model.frame(Terms, newdata, na.action = function(x) x)#,
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        X <- model.matrix(Terms, m, contrasts = object$contrasts)
        xint <- match("(Intercept)", colnames(X), nomatch=0)
        if(xint > 0) X <- X[, -xint, drop=FALSE]
        n <- nrow(X)
        y <- m[,names(cl)[attr(Terms, "response")]]
        if(length(object$zeta) > 0) {
            Terms <- attr(object$scale, "terms")
            m <- model.frame(Terms, newdata, na.action = function(x) x)#,
            if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, m)
            Z <- model.matrix(Terms, m, contrasts = object$contrasts)
            zint <- match("(Intercept)", colnames(Z), nomatch=0)
            if(zint > 0) Z <- Z[, -zint, drop=FALSE]
        }
        if(!is.null(object$nominal)) {
            Terms <- attr(object$nominal, "terms")
            m <- model.frame(Terms, newdata, na.action = function(x) x)#,
            if (!is.null(cl <- attr(Terms, "dataClasses")))
                .checkMFClasses(cl, m)
            XX <- model.matrix(Terms, m, contrasts = object$contrasts)
            namC <- colnames(XX)
        }
        B2 <- 1 * (col(matrix(0, n, nlevels(y))) == unclass(y))
        o1 <- c(100 * B2[, nlevels(y)])
        o2 <- c(-100 * B2[,1])
        B1 <- B2[,-nlevels(y), drop=FALSE]
        B2 <- B2[,-1, drop=FALSE]
        locationPar <- c(t(object$Theta))
        if(!is.null(object$nominal)) {
          ncolXX <- ncol(XX)
          LL1 <- lapply(1:ncolXX, function(x) B1 * XX[,x])
          B1 <- do.call(cbind, LL1)
          LL2 <- lapply(1:ncolXX, function(x) B2 * XX[,x])
          B2 <- do.call(cbind, LL2)
        }
        if(ncol(X) > 0) {
            B1 <- cbind(B1, -X)
            B2 <- cbind(B2, -X)
            locationPar <- c(locationPar, object$beta)
        }
        pfun <- switch(object$link,
                       logistic = plogis,
                       probit = pnorm,
                       cloglog = function(x) pgumbel(x, max=FALSE),
                       ## cloglog = pgumbel,
                       cauchit = pcauchy,
                       loglog = pgumbel,
                       "Aranda-Ordaz" = function(x, lambda) pAO(x, lambda),
                       "log-gamma" = function(x, lambda) plgamma(x, lambda))
        sigma <- 1
        if(length(object$zeta) > 0)
            sigma <- sigma * exp(drop(Z %*% object$zeta))
        eta1 <- (drop(B1 %*% locationPar) + o1) / sigma
        eta2 <- (drop(B2 %*% locationPar) + o2) / sigma
        if(object$link %in% c("Aranda-Ordaz", "log-gamma"))
          pr <- pfun(eta1, object$lambda) - pfun(eta2, object$lambda)
        else
          pr <- pfun(eta1) - pfun(eta2)
    }
    if(missing(newdata) && !is.null(object$na.action))
        pr <- napredict(object$na.action, pr)
    as.vector(pr)
}

profile.clm2 <- function(fitted, whichL = seq_len(p),
                        whichS = seq_len(k), lambda = TRUE, alpha = 0.01,
                        maxSteps = 50, delta = LrootMax/10, trace = 0,
                        stepWarn = 8, ...)
{
    rho <- update(fitted, doFit=FALSE)
    if(rho$estimLambda > 0 & rho$link == "Aranda-Ordaz")
        rho$limitLow <- c(rep(-Inf, length(rho$par)-2), 1e-5)
    nxi <- rho$nxi; k <- rho$k; p <- rho$p; X <- rho$X; Z <- rho$Z
    B1 <- rho$B1; B2 <- rho$B2
    sO <- rho$expSoffset; O1 <- rho$o1; O2 <- rho$o2
    beta0 <- with(fitted, coefficients[nxi + seq_len(p+k)])
    Lnames <- names(beta0[seq_len(p)])
    Snames <- names(beta0[p + seq_len(k)])
    Pnames <- c(Lnames, Snames)
    if(is.character(whichL)) whichL <- match(whichL, Lnames)
    if(is.character(whichS)) whichS <- match(whichS, Snames)
    nL <- length(whichL); nS <- length(whichS)
    summ <- summary(fitted)
    std.err <- summ$coefficients[nxi + seq_len(p+k), "Std. Error"]
    if(trace < 0) rho$ctrl$trace <- trace <- 1
    origLogLik <- fitted$logLik
    LrootMax <- qnorm(1 - alpha/2)
    prof <- vector("list", length = nL + nS)
    names(prof) <-
        c(paste("loc", Lnames[whichL], sep=".")[seq_len(nL)],
          paste("scale", Snames[whichS], sep=".")[seq_len(nS)])
    for(where in c("loc", "scale")[c(nL>0, nS>0)]) {
        if(where == "loc") {
            rho$p <- max(0, p - 1)
            which <- whichL }
        if(where == "scale") {
            which <- whichS
            rho$o1 <- O1
            rho$o2 <- O2
            rho$p <- p
            rho$k <- max(0, k - 1)
            rho$X <- X
            if(rho$nxi > 0) {
                rho$B1 <- B1
                rho$B2 <- B2 } }
        for(i in which) {
            if(where == "loc") {
                rho$X <- X[, -i, drop=FALSE]
                if(nxi > 0) {
                    rho$B1 <- B1[, -(nxi+i), drop=FALSE]
                    rho$B2 <- B2[, -(nxi+i), drop=FALSE] } }
            else {
                rho$Z <- Z[, -i, drop=FALSE]
                i <- i + p }
            res.i <- c(0, beta0[i])
            for(sgn in c(-1, 1)) {
                if(trace) {
                    message("\nParameter: ", where, ".", c(Lnames, Snames)[i],
                            c(" down", " up")[(sgn + 1)/2 + 1])
                    utils::flush.console() }
                rho$par <- fitted$coefficients[-(nxi+i)]
                step <- 0; Lroot <- 0
                while((step <- step + 1) < maxSteps && abs(Lroot) < LrootMax) {
                    beta.i <- beta0[i] + sgn * step * delta * std.err[i]
                    if(where=="loc") {
                        rho$o1 <- O1 - X[, i] * beta.i
                        rho$o2 <- O2 - X[, i] * beta.i }
                    else rho$expSoffset <- exp(sO + Z[, (i - p)] * beta.i)
                    fitCLM(rho)
                    Lroot <- sgn * sqrt(2*(-rho$logLik + origLogLik))
                    res.i <- rbind(res.i, c(Lroot, beta.i)) }
                if(step - 1 < stepWarn)
                    warning("profile may be unreliable for ",
                            where, ".", c(Lnames, Snames)[i],
                            " because only ", step - 1, "\n  steps were taken ",
                            c("downwards", "upwards")[(sgn + 1)/2 + 1])
            }
            rownames(res.i) <- NULL
            prof[[paste(where, c(Lnames, Snames)[i], sep=".")]] <- # -p+nL
                structure(data.frame(res.i[order(res.i[,1]),]),
                          names = c("Lroot",  c(Lnames, Snames)[i]))
            if(!all(diff(prof[[length(prof)]][,2]) > 0))
                warning("likelihood is not monotonically decreasing from maximum,\n",
                        "  so profile may be unreliable for ",
                        names(prof)[length(prof)])
        }
    }
    if(lambda & rho$nlambda)
        prof$lambda <- profileLambda(fitted, trace = trace, ...)
    val <- structure(prof, original.fit = fitted, summary = summ)
    class(val) <- c("profile.clm2", "profile")
    val
}

profileLambda <-
    function(fitted, link = fitted$link, range,
             nSteps = 20, trace = 0, ...)
{
    if(link == "log-gamma" & missing(range))
        range <- c(-4, 4)
    if(link == "Aranda-Ordaz" & missing(range))
        range <- c(1e-4, 4)
    if(!link %in% c("log-gamma", "Aranda-Ordaz"))
        stop("link needs to be 'log-gamma' or 'Aranda-Ordaz';", link,
             "not recognized")
    if(link == "Aranda-Ordaz" & min(range) <= 0)
        stop("range should be > 0 for the 'Aranda-Ordaz' link")
    if(fitted$estimLambda == 0)
        fitted <- update(fitted, Hess = FALSE, link = link,
                          lambda = NULL)
    MLogLik <- fitted$logLik
    MLlambda <- fitted$lambda
    logLik <- numeric(nSteps)
    rho <- update(fitted, Hess = FALSE, link = link,
                  lambda = min(range))
    logLik[1] <- rho$logLik
    rho <- update(rho, doFit = FALSE)
    lambdaSeq <- seq(min(range), max(range), length.out = nSteps)
    if(trace)  message("\nNow profiling lambda with ", nSteps - 1,
                       " steps: i =")
    for(i in 2:nSteps){
        if(trace) cat(i-1, " ")
        rho$lambda <- lambdaSeq[i]
        fitCLM(rho)
        logLik[i] <- rho$logLik
    }
    if(trace) cat("\n\n")
    if(any(logLik > fitted$logLik))
        warning("Profiling found a better optimum,",
                "  so original fit had not converged")

    sgn <- 2*(lambdaSeq > MLlambda) -1
    Lroot <- sgn * sqrt(2) * sqrt(-logLik + MLogLik)
    res <- data.frame("Lroot" = c(0, Lroot),
                      "lambda" = c(MLlambda, lambdaSeq))
    res <- res[order(res[,1]),]
    if(!all(diff(res[,2]) > 0))
        warning("likelihood is not monotonically decreasing from maximum,\n",
                "  so profile may be unreliable for lambda")
    res
}

confint.clm2 <-
    function(object, parm, level = 0.95,  whichL = seq_len(p),
             whichS = seq_len(k), lambda = TRUE, trace = 0, ...)
{
    p <- length(object$beta); k <- length(object$zeta)
    if(trace) {
        message("Waiting for profiling to be done...")
        utils::flush.console() }
    object <- profile(object, whichL = whichL, whichS = whichS,
                      alpha = (1. - level)/4., lambda = lambda,
                      trace = trace)
    confint(object, level=level, ...)
}

confint.profile.clm2 <-
  function(object, parm = seq_along(Pnames), level = 0.95, ...)
{
    of <- attr(object, "original.fit")
    Pnames <- names(object)
    if(is.character(parm))  parm <- match(parm, Pnames, nomatch = 0)
    a <- (1-level)/2
    a <- c(a, 1-a)
    pct <- paste(round(100*a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2),
                dimnames = list(Pnames[parm], pct))
    cutoff <- qnorm(a)
    for(pm in parm) {
        pro <- object[[ Pnames[pm] ]]
        sp <- spline(x = pro[, 2], y = pro[, 1])
        ci[Pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
    }
    ci
}

plot.profile.clm2 <-
    function(x, parm = seq_along(Pnames), level = c(0.95, 0.99),
             Log = FALSE, relative = TRUE, fig = TRUE, n = 1e3, ...,
             ylim = NULL)
### Should this function have a 'root' argument to display the
### likelihood root statistic (approximate straight line)?
{
    Pnames <- names(x)
    ML <- attr(x, "original.fit")$logLik
    for(pm in parm) {
        lim <- sapply(level, function(x)
                      exp(-qchisq(x, df=1)/2) )
        pro <- x[[ Pnames[pm] ]]
        sp <- spline(x = pro[, 2], y = pro[, 1], n=n)
        sp$y <- -sp$y^2/2
        if(relative & !Log) {
            sp$y <- exp(sp$y)
            ylab <- "Relative likelihood"
            dots <- list(...)
            if(missing(ylim))
                ylim <- c(0, 1)
        }
        if(relative & Log) {
            ylab <- "Relative log-likelihood"
            lim <- log(lim)
        }
        if(!relative & Log) {
            sp$y <- sp$y + ML
            ylab <- "Log-likelihood"
            lim <- ML + log(lim)
        }
        if(!relative & !Log) {
            stop("Not supported: at least one of 'Log' and 'relative' ",
                 "have to be TRUE")
            sp$y <- exp(sp$y + ML)
            ylab <- "Likelihood"
            lim <- exp(ML + log(lim))
        }
        x[[ Pnames[pm] ]] <- sp
        if(fig) {
            plot(sp$x, sp$y, type = "l", ylim = ylim,
                 xlab = Pnames[pm], ylab = ylab, ...)
            abline(h = lim)
        }
    }
    attr(x, "limits") <- lim
    invisible(x)
}

logLik.clm2 <- function(object, ...)
    structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clm2 <- function(fit, scale = 0, k = 2, ...)
{
    edf <- fit$edf
    c(edf, -2*fit$logLik + k * edf)
}

update.clm2 <-
    function(object, formula., location, scale, nominal, ..., evaluate = TRUE)
### This method makes it possible to use the update.formula features
### for location and scale formulas in clm2 objects. This includes the
### possibility of using e.g.
### update(obj, loc = ~ . - var1, sca = ~ . + var2)
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

dropterm.clm2 <-
  function(object, scope, scale = 0, test = c("none", "Chisq"),
           k = 2, sorted = FALSE, trace = FALSE,
           which = c("location", "scale"), ...)
### Most of this is lifted from MASS::dropterm.default, but adapted to
### the two formulas (location and scale) in the model.
{
    which <- match.arg(which)
    Terms <-
        if(which == "location") attr(object$location, "terms")
        else attr(object$scale, "terms")
    tl <- attr(Terms, "term.labels")
    if(missing(scope)) scope <- drop.scope(Terms)
    else {
        if(!is.character(scope))
            scope <- attr(terms(update.formula(Terms, scope)),
                          "term.labels")
        if(!all(match(scope, tl, FALSE)))
            stop("scope is not a subset of term labels")
    }
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames =  list(c("<none>", scope), c("df", "AIC")))
    ans[1,  ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$fitted)
    for(i in seq(ns)) {
        tt <- scope[i]
        if(trace) {
	    message("trying -", tt)
	    utils::flush.console()
	}
        Call <- as.list(object$call)
        Call[[which]] <-
            update.formula(Terms, as.formula(paste("~ . -", tt)))
        nfit <- eval.parent(as.call(Call))
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$fitted) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[1,1] - ans[,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    o <- if(sorted) order(aod$AIC) else seq_along(aod$AIC)
    test <- match.arg(test)
    if(test == "Chisq") {
        dev <- ans[, 2] - k*ans[, 1]
        dev <- dev - dev[1] ; dev[1] <- NA
        nas <- !is.na(dev)
        P <- dev
        P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail = FALSE)
        aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    Call <- as.list(object$call)
    Call <- Call[names(Call) %in% c("location", "scale")]
    head <- c("Single term deletions", "\nModel:",
              paste(names(Call), ":",  Call))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

addterm.clm2 <-
    function(object, scope, scale = 0, test = c("none", "Chisq"),
             k = 2, sorted = FALSE, trace = FALSE,
             which = c("location", "scale"), ...)
### Most of this is lifted from MASS::addterm.default, but adapted to
### the two formulas (location and scale) in the model.
{
    which <- match.arg(which)
    if (which == "location")
        Terms <- attr(object$location, "terms")
    else if(!is.null(object$scale))
        Terms <- attr(object$scale, "terms")
    else
        Terms <- as.formula(" ~ 1")
    if(missing(scope) || is.null(scope)) stop("no terms in scope")
    if(!is.character(scope))
        scope <- add.scope(Terms, update.formula(Terms, scope))
    if(!length(scope))
        stop("no terms in scope for adding to object")
    ns <- length(scope)
    ans <- matrix(nrow = ns + 1, ncol = 2,
                  dimnames = list(c("<none>", scope), c("df", "AIC")))
    ans[1,  ] <- extractAIC(object, scale, k = k, ...)
    n0 <- length(object$fitted)
    for(i in seq(ns)) {
        tt <- scope[i]
        if(trace) {
	    message("trying +", tt)
	    utils::flush.console()
        }
        Call <- as.list(object$call)
        Call[[which]] <-
            update.formula(Terms, as.formula(paste("~ . +", tt)))
        nfit <- eval.parent(as.call(Call))
	ans[i+1, ] <- extractAIC(nfit, scale, k = k, ...)
        if(length(nfit$fitted) != n0)
            stop("number of rows in use has changed: remove missing values?")
    }
    dfs <- ans[,1] - ans[1,1]
    dfs[1] <- NA
    aod <- data.frame(Df = dfs, AIC = ans[,2])
    o <- if(sorted) order(aod$AIC) else seq_along(aod$AIC)
    test <- match.arg(test)
    if(test == "Chisq") {
	dev <- ans[,2] - k*ans[, 1]
	dev <- dev[1] - dev; dev[1] <- NA
	nas <- !is.na(dev)
	P <- dev
	P[nas] <- pchisq(dev[nas], dfs[nas], lower.tail=FALSE)
	aod[, c("LRT", "Pr(Chi)")] <- list(dev, P)
    }
    aod <- aod[o, ]
    Call <- as.list(object$call)
    Call <- Call[names(Call) %in% c("location", "scale")]
    head <- c("Single term additions", "\nModel:",
              paste(names(Call), ":",  Call))
    if(scale > 0)
        head <- c(head, paste("\nscale: ", format(scale), "\n"))
    class(aod) <- c("anova", "data.frame")
    attr(aod, "heading") <- head
    aod
}

## addterm <- function(object, ...) UseMethod("addterm")
## dropterm <- function(object, ...) UseMethod("dropterm")

##################################################################
## Additional utility functions:

grad.lambda <- function(rho, lambda, link, delta = 1e-6) {
    ll <- lambda + c(-delta, delta)
    if(link == "Aranda-Ordaz") ll[ll < 0] <- 0
    par <- rho$par
    f <- sapply(ll, function(x) getNll(rho, c(par[-length(par)], x)))
    rho$lambda <- lambda
    rho$par <- par
    diff(f) /  diff(ll)
}

TraceC <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    .C("trace",
       as.integer(iter[1]),
       as.double(stepFactor[1]),
       as.double(val[1]),
       as.double(maxGrad[1]),
       as.double(par),
       length(par),
       as.integer(first[1]))

    return(invisible())
}

TraceR <- function(iter, stepFactor, val, maxGrad, par, first=FALSE) {
    t1 <- sprintf(" %3d:     %.2e:   %.3f:   %1.3e:  ",
                  iter, stepFactor, val, maxGrad)
    t2 <- formatC(par)
    if(first)
        cat("iter:  step factor:      Value:   max|grad|:   Parameters:\n")
    cat(t1, t2, "\n")
}

print.Anova <- function (x, ...)
## Lifted from package MASS:
{
    heading <- attr(x, "heading")
    if (!is.null(heading))
        cat(heading, sep = "\n")
    attr(x, "heading") <- NULL
    res <- format.data.frame(x, ...)
    nas <- is.na(x)
    res[] <- sapply(seq_len(ncol(res)), function(i) {
        x <- as.character(res[[i]])
        x[nas[, i]] <- ""
        x
    })
    print.data.frame(res)
    invisible(x)
}

fixed <- function(theta, eps = 1e-3) {
    res <- vector("list")
    res$name <- "fixed"
    if(!missing(theta) && length(theta) > 1) {
        if(length(theta) < 3)
            stop("'length(theta) = ", length(theta),
                 ", but has to be 1 or >= 3")
        res$eps <- NULL
        res$theta <- theta
        res$getTheta <- function(y, theta, eps) theta
    }
    else if(!missing(theta) && length(theta) == 1) {
        if(as.integer(theta) < 3)
            stop("'as.integer(theta)' was ", as.integer(theta),
                 ", but has to be > 2")
        res$eps <- NULL
        res$theta <- theta
        res$getTheta <- function(y, theta, eps) {
            eps <- diff(range(y)) / (theta - 1)
            seq(min(y) - eps/2, max(y) + eps/2, len = theta + 1)
        }
    }
    else if(missing(theta) && length(eps) == 1) {
        res$eps <- eps
        res$theta <- NULL
        res$getTheta <- function(y, theta, eps) {
            J <- diff(range(y))/eps + 1
            seq(min(y) - eps/2, max(y) + eps/2, len = J)
        }
    }
    else
        stop("inappropriate arguments")
    class(res) <- "threshold"
    res
}

makeThresholds2 <- function(rho, threshold, ...) {
    if(threshold == "flexible") {
        rho$tJac <- diag(rho$ntheta)
        rho$nalpha <- rho$ntheta
        rho$alphaNames <-
            paste(rho$lev[-length(rho$lev)], rho$lev[-1], sep="|")
    }
    if(threshold == "symmetric") {
        if(!rho$ntheta >=2)
            stop("symmetric thresholds are only meaningful for responses with 3 or more levels")
        if(rho$ntheta %% 2) { ## ntheta is odd
            rho$nalpha <- (rho$ntheta + 1)/2 ## No. threshold parameters
            rho$tJac <- t(cbind(diag(-1, rho$nalpha)[rho$nalpha:1, 1:(rho$nalpha-1)],
                                diag(rho$nalpha)))
            rho$tJac[,1] <- 1
            rho$alphaNames <-
                c("central", paste("spacing.", 1:(rho$nalpha-1), sep=""))
        }
        else { ## ntheta is even
            rho$nalpha <- (rho$ntheta + 2)/2
            rho$tJac <- cbind(rep(1:0, each=rho$ntheta/2),
                              rbind(diag(-1, rho$ntheta/2)[(rho$ntheta/2):1,],
                                    diag(rho$ntheta/2)))
            rho$tJac[,2] <- rep(0:1, each=rho$ntheta/2)
            rho$alphaNames <- c("central.1", "central.2",
                                paste("spacing.", 1:(rho$nalpha-2), sep=""))
        }
    }
    if(threshold == "equidistant") {
        if(!rho$ntheta >=2)
            stop("symmetric thresholds are only meaningful for responses with 3 or more levels")
        rho$tJac <- cbind(1, 0:(rho$ntheta-1))
        rho$nalpha <- 2
        rho$alphaNames <- c("threshold.1", "spacing")
    }
}
