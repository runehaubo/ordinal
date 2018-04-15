## This file contains:
## The function clm.fit() - an lm.fit or glm.fit equivalent for CLMs.

clm.fit <- function(y, ...) {
    UseMethod("clm.fit")
}

clm.fit.factor <-
  function(y, X, S, N, weights = rep(1, nrow(X)),
           offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
           control = list(), start, doFit=TRUE,
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"),
           ...)
### This function basically does the same as clm, but without setting
### up the model matrices from formulae, and with minimal post
### processing after parameter estimation.
{
    ## Initial argument matching and testing:
    threshold <- match.arg(threshold)
    link <- match.arg(link)
    control <- do.call(clm.control, control)
    if(missing(y)) stop("please specify y")
    if(missing(X)) X <- cbind("(Intercept)" = rep(1, length(y)))
    stopifnot(is.factor(y), is.matrix(X))
    if(missing(weights) || is.null(weights))
        weights <- rep(1, length(y))
    if(missing(offset) || is.null(offset))
        offset <- rep(0, length(y))
    if(missing(S.offset) || is.null(S.offset))
        S.offset <- rep(0, length(y))
    stopifnot(length(y) == nrow(X) &&
              length(y) == length(weights) &&
              length(y) == length(offset) &&
              length(y) == length(S.offset))
    frames <- list(y=y, X=X)
    y[weights <= 0] <- NA
    y.levels <- levels(droplevels(y))
    struct <- namedList(y, X, weights, offset, S.offset, y.levels,
                        threshold, link, control, doFit)
    ## S and N are optional:
    if(!missing(S) && !is.null(S)) {
        struct$S <- S
        stopifnot(is.matrix(S),
                  length(y) == nrow(S))
    }
    if(!missing(N) && !is.null(N)) {
        struct$NOM <- N
        stopifnot(is.matrix(N),
                  length(y) == nrow(N))
    }
    clm.fit.default(struct)
}

clm.fit.default <- function(y, ...)
### y: design object with the following components: ...
### (tJac=NULL), (y.levels=NULL), threshold, (aliased=NULL),
### (start=NULL), link, control, weights, (coef.names=NULL), y, X,
### (S=NULL), (NOM=NULL), doFit=TRUE, S.offset=NULL
{
    ## check args:
    stopifnot(is.list(y))
    y <- c(y, list(...))
    stopifnot(all(
        c("y", "X", "offset", "weights", "link", "threshold",
          "control", "doFit") %in% names(y) ))
    ## preprocess design objects if needed:
    if(is.null(y$y.levels)) y$y.levels <- levels(y$y)
    if(is.null(y$tJac)) {
        y <- c(y, makeThresholds(y$y.levels, y$threshold))
    }
    if(is.null(y$aliased))
        y <- drop.cols(y, silent=TRUE, drop.scale=FALSE)
    ## Make model environment:
    rho <- do.call(clm.newRho, y)
    start <- set.start(rho, start=y$start, get.start=is.null(y$start),
                       threshold=y$threshold, link=y$link,
                       frames=y)
    rho$par <- as.vector(start) ## remove attributes
    setLinks(rho, y$link)
    if(y$doFit == FALSE) return(rho)
    ## Fit the model:
    fit <- if(y$control$method == "Newton") {
        clm.fit.NR(rho, y$control)
    } else {
        clm.fit.optim(rho, y$control$method, y$control$ctrl) }
    ## Adjust iteration count:
    if(y$control$method == "Newton" &&
       !is.null(start.iter <- attr(start, "start.iter")))
        fit$niter <- fit$niter + start.iter
    ## Update coefficients, gradient, Hessian, edf, nobs, n,
    ## fitted.values, df.residual:
    fit <- clm.finalize(fit, y$weights, y$coef.names, y$aliased)
    fit$tJac <- format_tJac(y$tJac, y$y.levels, y$alpha.names)
    th.res <- formatTheta(fit$alpha, fit$tJac, y)
    ## Check convergence:
    conv <- conv.check(fit, control=y$control, Theta.ok=th.res$Theta.ok,
                       tol=y$control$tol)
    print.conv.check(conv, action=y$control$convergence) ## print convergence message
    th.res$Theta.ok <- NULL
    fit <- c(fit, conv[c("vcov", "cond.H")], th.res)
    fit$convergence <- conv[!names(conv) %in% c("vcov", "cond.H")]
    fit <- fit[sort(names(fit))]
    class(fit) <- "clm.fit"
    fit
}

clm.finalize <- function(fit, weights, coef.names, aliased)
### extracFromFit
###
### distinguishing between par and coef where the former does not
### contain aliased coefficients.
{
    nalpha <- length(aliased$alpha)
    nbeta <- length(aliased$beta)
    nzeta <- length(aliased$zeta)
    ncoef <- nalpha + nbeta + nzeta ## including aliased coef
    npar <- sum(!unlist(aliased)) ## excluding aliased coef
    stopifnot(length(fit$par) == npar)
    fit <- within(fit, {
        coefficients <- rep(NA, nalpha + nbeta + nzeta)
        ## ensure correct order of alpha, beta and zeta:
        keep <- match(c("alpha", "beta", "zeta"), names(aliased),
                      nomatch=0)
        aliased <- lapply(aliased[keep], as.logical)
        for(i in names(aliased))
            names(aliased[[i]]) <- coef.names[keep][[i]]
        names(coefficients) <- unlist(coef.names[keep])
        par.names <- names(coefficients)[!unlist(aliased)]
        coefficients[!unlist(aliased)] <- par
        alpha <- coefficients[1:nalpha]
        if(nbeta) beta <- coefficients[nalpha + 1:nbeta]
        if(nzeta) zeta <- coefficients[nalpha + nbeta + 1:nzeta]
        names(gradient) <- par.names
        dimnames(Hessian) <- list(par.names, par.names)
        edf <- npar ## estimated degrees of freedom
        nobs <- sum(weights)
        n <- length(weights)
        fitted.values <- fitted
        df.residual = nobs - edf
        ## keep <- i <- fitted <- par.names <- par <- coef.names <- NULL
    })
    notkeep <- c("keep", "i", "fitted", "par.names", "par",
                 "coef.names")
    fit[!names(fit) %in% notkeep]
}
