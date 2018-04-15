## This file contains:
## Functions to assess and check convergence of CLMs. Some
## functions/methods are exported and some are used internally in
## clm().

convergence <- function(object, ...) {
  UseMethod("convergence")
}

convergence.clm <-
  function(object, digits = max(3, getOption("digits") - 3),
           tol = sqrt(.Machine$double.eps), ...)
### Results: data.frame with columns:
### Estimate
### Std. Error
### Gradient - gradient of the coefficients at optimizer termination
### Error - the signed error in the coefficients at termination
### Rel. Error - the relative error in the coefficeints at termination
###
### The (signed) Error is determined as the Newton step, so this is
### only valid close to the optimum where the likelihood function is
### quadratic.
###
### The relative error equals step/Estimate.
{
    ## get info table and coef-table:
    info <- object$info[c("nobs", "logLik", "niter", "max.grad", "cond.H")]
    ## Initialize coef-table with NAs:
    coefs <- coef(object, na.rm=TRUE)
    g <- object$gradient
    H <- object$Hessian
    tab <- matrix(NA_real_, nrow=length(coefs), ncol=6L,
                  dimnames=list(names(coef(object, na.rm=TRUE)),
                  c("Estimate", "Std.Err", "Gradient",
                    "Error", "Cor.Dec", "Sig.Dig")))
    tab[, c(1L, 3L)] <- cbind(coefs, g)
    res <- list(info=info, coefficients=tab, original.fit=object)
    class(res) <- "convergence.clm"
    if(!all(is.finite(H))) {
        warning("non-finite values in Hessian: illegitimate model fit")
        return(res)
    }
    ## Get eigen values of Hessian:
    res$eigen.values <- e.val <-
        eigen(H, symmetric=TRUE, only.values=TRUE)$values
    ## Compute Cholesky factor of Hessian:
    ch <- try(chol(H), silent=TRUE)
    if(any(abs(e.val) <= tol) || inherits(ch, "try-error")) {
        return(res)
    }
    ## Hessian is positive definite:
    ## Compute approximate error in the coefficients:
    step <- c(backsolve(ch, backsolve(ch, g, transpose=TRUE)))
    if(max(abs(step)) > 1e-2)
        warning("convergence assessment may be unreliable ",
                "due to large numerical error")
    ## Compute approximate error in the log-likelihood function:
    env <- get_clmRho(object)
    ## Note: safer to get env this way.
    ## env <- update(object, doFit=FALSE)
    env$par <- coef(object, na.rm=TRUE) - step
    new.logLik <- -env$clm.nll(env)
    new.max.grad <- max(abs(env$clm.grad(env)))
    if(new.max.grad > max(abs(g)) && max(abs(step)) > tol)
        warning("Convergence assessment may be unreliable: ",
                "please assess the likelihood with slice()")
### NOTE: we only warn if step is larger than a tolerance, since if
### step \sim 1e-16, the max(abs(grad)) may increase though stay
### essentially zero.
    logLik.err <- object$logLik - new.logLik
    err <- format.pval(logLik.err, digits=2, eps=1e-10)
    if(!length(grep("<", err)))
        err <- formatC(as.numeric(err), digits=2, format="e")
    res$info$logLik.Error <- err
    ## Fill in the coef-table:
    se <- sqrt(diag(chol2inv(ch)))
    res$coefficients[, c(2, 4:6)] <-
        cbind(se, step, cor.dec(step),
              signif.digits(coefs, step))
    res
}

print.convergence.clm <-
    function(x, digits = max(3, getOption("digits") - 3), ...)
{
    ## Prepare for printing:
    print(x$info, row.names=FALSE, right=FALSE)
    cat("\n")
    tab.print <- coef(x)
    for(i in 1:2)
        tab.print[,i] <- format(c(coef(x)[,i]), digits=digits)
    for(i in 3:4) tab.print[,i] <-
        format(c(coef(x)[,i]), digits=max(1, digits - 1))
    print(tab.print, quote=FALSE, right=TRUE, ...)
    ## Print eigen values:
    cat("\nEigen values of Hessian:\n")
    cat(format(x$eigen.values, digits=digits), "\n")
    conv <- x$original.fit$convergence
    cat("\nConvergence message from clm:\n")
    for(i in seq_along(conv$code)) {
        Text <- paste("(", conv$code[i], ") ", conv$messages[i], sep="")
        cat(Text, "\n")
    }
    if(!is.null(alg.text <- conv$alg.message))
        cat(paste("In addition:", alg.text), "\n")
    cat("\n")
    ## for(i in seq_along(conv$code)) {
    ##     cat("Code:  Message:\n", fill=TRUE)
    ##     cat(conv$code[i], "      ", conv$message[i], "\n", fill=TRUE)
    ## }
    ## if(!is.null(alg.text <- conv$alg.message)) {
    ##     cat("\nIn addition: ", alg.text, "\n\n", fill=TRUE)
    ## }
    return(invisible(x))
}


cor.dec <- function(error) {
### computes the no. correct decimals in a number if 'error' is the
### error in the number.
### The function is vectorized.
  xx <- -log10(abs(error))
  lead <- floor(xx)
  res <- ifelse(xx < lead - log10(.5), lead-1, lead)
  res[abs(error) >= .05] <- 0
  as.integer(round(res))
}

signif.digits <- function(value, error) {
### Determines the number of significant digits in 'value' if the
### absolute error in 'value' is 'error'.
### The function is vectorized.
  res <- cor.dec(error) + ceiling(log10(abs(value)))
  res[res < 0] <- 0
  as.integer(round(res))
}

conv.check <-
    function(fit, control=NULL, Theta.ok=NULL, tol=sqrt(.Machine$double.eps), ...)
    ## function(gr, Hess, conv, method, gradTol, relTol,
    ##          tol=sqrt(.Machine$double.eps), ...)
### Compute variance-covariance matrix and check convergence along the
### way.
### fit: clm-object or the result of clm.fit.NR() | gradient, Hessian,
### (control), convergence
### control: (tol), (method), gradTol, relTol
###
### Return: list with elements
### vcov, conv, cond.H, messages and
{
    if(missing(control))
        control <- fit$control
    if(is.null(control))
        stop("'control' not supplied - cannot check convergence")

    if(!is.null(control$tol))
        tol <- control$tol
    if(tol < 0)
        stop(gettextf("numerical tolerance is %g, expecting non-negative value",
                      tol), call.=FALSE)
### FIXME: test this.
    H <- fit$Hessian
    g <- fit$gradient
    max.grad <- max(abs(g))
    cov <- array(NA_real_, dim=dim(H), dimnames=dimnames(H))
    cond.H <- NA_real_
    res <- list(vcov=cov, code=integer(0L), cond.H=cond.H,
                messages=character(0L))
    class(res) <- "conv.check"
    if(is.list(code <- fit$convergence))
        code <- code[[1L]]
    mess <-
        switch(as.character(code),
               "0" = "Absolute and relative convergence criteria were met",
               "1" = "Absolute convergence criterion was met, but relative criterion was not met",
               "2" = "iteration limit reached",
               "3" = "step factor reduced below minimum",
               "4" = "maximum number of consecutive Newton modifications reached")
    if(control$method != "Newton") mess <- NULL
### FIXME: get proper convergence message from optim, nlminb, ucminf etc.
    res <- c(res, alg.message=mess)
    ## }
    evd <- eigen(H, symmetric=TRUE, only.values=TRUE)$values
    negative <- sum(evd < -tol)
    if(negative) {
        res$code <- -2L
        res$messages <-
            gettextf(paste("Model failed to converge:",
                           "degenerate Hessian with %d negative eigenvalues"),
                     negative)
        return(res)
    }
    ## Add condition number to res:
    res$cond.H <- max(evd) / min(evd)
    ## Compute Newton step:
    ch <- try(chol(H), silent=TRUE)
    if(max.grad > control$gradTol) {
        res$code <- -1L
        res$messages <-
            gettextf("Model failed to converge with max|grad| = %g (tol = %g)",
                     max.grad, control$gradTol)
        return(res)
    }
    if(!is.null(Theta.ok) && !Theta.ok) {
        res$code <- -3L
        res$messages <-
            "not all thresholds are increasing: fit is invalid"
        return(res)
    }
    zero <- sum(abs(evd) < tol)
    if(zero || inherits(ch, "try-error")) {
        res$code <- 1L
        res$messages <-
            "Hessian is numerically singular: parameters are not uniquely determined"
        return(res)
    }
### NOTE: Only do the following if 'ch <- try(chol(H), silent=TRUE)'
### actually succedded:
    step <- c(backsolve(ch, backsolve(ch, g, transpose=TRUE)))
    ## Compute var-cov:
    res$vcov[] <- chol2inv(ch)
### NOTE: we want res$vcov to be present in all of the situations
### below.
    if(max(abs(step)) > control$relTol) {
        res$code <- c(res$code, 1L)
        corDec <- as.integer(min(cor.dec(step)))
        res$messages <-
            c(res$messages,
              gettextf("some parameters may have only %d correct decimals",
                       corDec))
    }
    if(max(evd) * tol > 1) {
        res$code <- c(res$code, 2L)
        res$messages <-
            c(res$messages,
              paste("Model is nearly unidentifiable: ",
                    "very large eigenvalue",
                    "\n - Rescale variables?", sep=""))
    }
    if((min(evd) / max(evd)) < tol) {
        res$code <- c(res$code, 3L)
        if(!5L %in% res$code) {
            res$messages <-
                c(res$messages,
                  paste("Model is nearly unidentifiable: ",
                        "large eigenvalue ratio",
                        "\n - Rescale variables?", sep=""))
        }
    }
    if(!length(res$code)) {
        res$code <- 0L
        res$messages <- "successful convergence"
    }
    res
}

cov.conv <- conv.check

### FIXME: let convergence() print convergence info from clm using
### print.conv.check

print.conv.check <-
    function(x, action=c("warn", "silent", "stop", "message"), ...)
{
    action <- match.arg(action)
    if(x$code == 0L || action == "silent") return(invisible())

    Text <- paste("(", x$code[1L], ") ", x$messages[1L], sep="")
    if(!is.null(alg.text <- x$alg.message))
        Text <- paste(Text, "\nIn addition:", alg.text)
    switch(action,
           "stop" = stop(Text, call.=FALSE),
           "warn" = warning(Text, call.=FALSE),
           "message" = message(Text))
}
