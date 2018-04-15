## This file contains:
## Implementation of various methods for clm objects.

print.clm <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
{
  cat("formula:", Deparse(formula(x$terms)), fill=TRUE)
### NOTE: deparse(x$call$formula) will not always work since this may
### not always be appropriately evaluated.
  if(!is.null(x$call$scale))
    cat("scale:  ", Deparse(formula(x$S.terms)), fill=TRUE)
  if(!is.null(x$call$nominal))
    cat("nominal:", Deparse(formula(x$nom.terms)), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", Deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", Deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  if(length(x$beta)) {
    if(sum(x$aliased$beta) > 0) {
      cat("\nCoefficients: (", sum(x$aliased$beta),
          " not defined because of singularities)\n", sep = "")
    }
    else cat("\nCoefficients:\n")
    print.default(format(x$beta, digits = digits),
                  quote = FALSE)
  }
    if(length(x$zeta)) {
    if(sum(x$aliased$zeta) > 0)
      cat("\nlog-scale coefficients: (", sum(x$aliased$zeta),
          " not defined because of singularities)\n", sep = "")
    else cat("\nlog-scale coefficients:\n")
    print.default(format(x$zeta, digits = digits),
                  quote = FALSE)
  }
  if(length(x$alpha) > 0) {
    if(sum(x$aliased$alpha) > 0)
      cat("\nThreshold coefficients: (", sum(x$aliased$alpha),
          " not defined because of singularities)\n", sep = "")
    else cat("\nThreshold coefficients:\n")
    if(!is.null(x$call$nominal))
      print.default(format(x$alpha.mat, digits = digits),
                    quote = FALSE)
    else
      print.default(format(x$alpha, digits = digits),
                    quote = FALSE)
  }

  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  return(invisible(x))
}

vcov.clm <-
    function(object, tol = sqrt(.Machine$double.eps),
             method = c("clm", "Cholesky", "svd", "eigen", "qr"), ...)
{
    method <- match.arg(method)
    if(method == "clm")
        return(object$vcov)
    if(is.null(object$Hessian))
        stop("Model needs to be fitted with Hess = TRUE")
    dn <- dimnames(object$Hessian)
    H <- object$Hessian
    if(!all(is.finite(H)))
        stop("cannot compute vcov: non-finite values in Hessian")
    if(method == "svd") {
        Hsvd <- svd(H)
        ## positive <- Hsvd$d > max(tol * Hsvd$d[1L], tol)
        positive <- Hsvd$d > tol
        if(!all(positive))
            stop(gettextf("Cannot compute vcov: \nHessian is numerically singular with min singular value = %g",
                          min(Hsvd$d)))
        cov <- Hsvd$v %*% (1/Hsvd$d * t(Hsvd$u))
    }
    else if(method == "eigen") {
        evd <- eigen(H, symmetric=TRUE)
        ## tol <- max(tol * evd$values[1L], tol) ## if evd$values[1L] < 0
        if(any(evd$values < tol))
            stop(gettextf("Cannot compute vcov: \nHessian is not positive definite with min eigenvalue = %g",
                          min(evd$values)))
        cov <- with(evd, vectors %*% diag(1/values) %*% t(vectors))
    }
    else if(method == "Cholesky") {
        cholH <- try(chol(H), silent=TRUE)
        if(inherits(cholH, "try-error"))
            stop("Cannot compute vcov: \nHessian is not positive definite")
        cov <- chol2inv(cholH)
    }
    else if(method == "qr") {
        qrH <- qr(H, tol=sqrt(.Machine$double.eps))
        if(qrH$rank < nrow(H))
            stop("Cannot compute vcov: \nHessian is numerically singular")
        cov <- solve.qr(qrH)
    }
    else
        stop("method not recognized")
    ## Need to test for negative variances, since some methods (svd,
    ## qr) may produce a vcov-matrix if the Hessian is *negative*
    ## definite:
    if(any(diag(cov) < 0)) {
        stop("Cannot compute vcov: \nHessian is not positive definite")
        }
    structure(cov, dimnames=dn)
}

summary.clm <- function(object, correlation = FALSE, ...)
{
    vcov <- object$vcov
    coefs <- matrix(NA, length(object$coefficients), 4,
                    dimnames = list(names(object$coefficients),
                    c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
    coefs[, 1] <- object$coefficients
    if(!all(is.finite(vcov))) {
        ## warning("Variance-covariance matrix of the parameters is not defined")
        coefs[, 2:4] <- NA
        if(correlation) warning("Correlation matrix is unavailable")
    }
    else {
        alias <- unlist(object$aliased)
        coefs[!alias, 2] <- sd <- sqrt(diag(vcov))
        ## Cond is Inf if Hessian contains NaNs:
        object$cond.H <-
            if(any(is.na(object$Hessian))) Inf
            else with(eigen(object$Hessian, symmetric=TRUE, only.values = TRUE),
                      abs(max(values) / min(values)))
        coefs[!alias, 3] <- coefs[!alias, 1]/coefs[!alias, 2]
        coefs[!alias, 4] <- 2 * pnorm(abs(coefs[!alias, 3]),
                                      lower.tail=FALSE)
        if(correlation)
            object$correlation <- cov2cor(vcov)
    }
    object$coefficients <- coefs
    class(object) <- "summary.clm"
    return(object)
}

print.summary.clm <-
    function(x, digits = max(3, getOption("digits") - 3),
             signif.stars = getOption("show.signif.stars"), ...)
{
  cat("formula:", Deparse(formula(x$terms)), fill=TRUE)
### NOTE: deparse(x$call$formula) will not always work since this may
### not always be appropriately evaluated.
  if(!is.null(x$call$scale))
    cat("scale:  ", Deparse(formula(x$S.terms)), fill=TRUE)
  if(!is.null(x$call$nominal))
    cat("nominal:", Deparse(formula(x$nom.terms)), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", Deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", Deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  nalpha <- length(x$alpha)
  nbeta <- length(x$beta)
  nzeta <- length(x$zeta)
  if(nbeta > 0) {
    if(sum(x$aliased$beta) > 0)
      cat("\nCoefficients: (", sum(x$aliased$beta),
          " not defined because of singularities)\n", sep = "")
    else cat("\nCoefficients:\n")
    printCoefmat(x$coefficients[nalpha + 1:nbeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...)
  } ## else  cat("\nNo Coefficients\n")
  if(nzeta > 0) {
    if(sum(x$aliased$zeta) > 0)
      cat("\nlog-scale coefficients: (", sum(x$aliased$zeta),
          " not defined because of singularities)\n", sep = "")
    else cat("\nlog-scale coefficients:\n")
    printCoefmat(x$coefficients[nalpha + nbeta + 1:nzeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...)
  }
  if(nalpha > 0) { ## always true
    if(sum(x$aliased$alpha) > 0)
      cat("\nThreshold coefficients: (", sum(x$aliased$alpha),
          " not defined because of singularities)\n", sep = "")
    else cat("\nThreshold coefficients:\n")
    printCoefmat(x$coefficients[seq_len(nalpha), -4, drop=FALSE],
                 digits=digits, has.Pvalue=FALSE, signif.stars=FALSE,
                 ...)
  }

  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1, -ncol(correl)], quote = FALSE, ...)
  }
  return(invisible(x))
}

logLik.clm <- function(object, ...)
  structure(object$logLik, df = object$edf, nobs=object$nobs,
            class = "logLik")

extractAIC.clm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}

### NOTE: AIC.clm implicitly defined via logLik.clm

anova.clm <- function(object, ...)
### requires that clm objects have components:
###  edf: no. parameters used
###  call$formula
###  link (character)
###  threshold (character)
###  logLik
###
{
  mc <- match.call()
  dots <- list(...)
  ## remove 'test' and 'type' arguments from dots-list:
  not.keep <- which(names(dots) %in% c("test", "type"))
  if(length(not.keep)) {
    message("'test' and 'type' arguments ignored in anova.clm\n")
    dots <- dots[-not.keep]
  }
  if(length(dots) == 0)
    stop('anova is not implemented for a single "clm" object')
  mlist <- c(list(object), dots)
  if(!all(sapply(mlist, function(model)
                 inherits(model, c("clm", "clmm")))))
    stop("only 'clm' and 'clmm' objects are allowed")
  nfitted <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(nfitted != nfitted[1L]))
    stop("models were not all fitted to the same dataset")
### FIXME: consider comparing y returned by the models for a better
### check?
  no.par <- sapply(mlist, function(x) x$edf)
  ## order list with increasing no. par:
  ord <- order(no.par, decreasing=FALSE)
  mlist <- mlist[ord]
  no.par <- no.par[ord]
  no.tests <- length(mlist)
  ## extract formulas, links, thresholds, scale formulas, nominal
  ## formulas:
  forms <- sapply(mlist, function(x) Deparse(x$call$formula))
  links <- sapply(mlist, function(x) x$link)
  thres <- sapply(mlist, function(x) x$threshold)
  nominal <- sapply(mlist, function(x) Deparse(x$call$nominal))
  scale <- sapply(mlist, function(x) Deparse(x$call$scale))
  models <- data.frame(forms)
  models.names <- 'formula:'
  if(any(!nominal %in% c("~1", "NULL"))) {
    nominal[nominal == "NULL"] <- "~1"
    models$nominal <- nominal
    models.names <- c(models.names, "nominal:")
  }
  if(any(!scale %in% c("~1", "NULL"))) {
    scale[scale == "NULL"] <- "~1"
    models$scale <- scale
    models.names <- c(models.names, "scale:")
  }
  models.names <- c(models.names, "link:", "threshold:")
  models <- cbind(models, data.frame(links, thres))
  ## extract AIC, logLik, statistics, df, p-values:
  AIC <- sapply(mlist, function(x) -2*x$logLik + 2*x$edf)
  logLiks <- sapply(mlist, function(x) x$logLik)
  statistic <- c(NA, 2*diff(sapply(mlist, function(x) x$logLik)))
  df <- c(NA, diff(no.par))
  pval <- c(NA, pchisq(statistic[-1], df[-1], lower.tail=FALSE))
  pval[!is.na(df) & df==0] <- NA
  ## collect results in data.frames:
  tab <- data.frame(no.par, AIC, logLiks, statistic, df, pval)
  tab.names <- c("no.par", "AIC", "logLik", "LR.stat", "df",
                 "Pr(>Chisq)")
  mnames <- sapply(as.list(mc), Deparse)[-1]
  colnames(tab) <- tab.names
  rownames(tab) <- rownames(models) <- mnames[ord]
  colnames(models) <- models.names
  attr(tab, "models") <- models
  attr(tab, "heading") <-
    "Likelihood ratio tests of cumulative link models:\n"
  class(tab) <- c("anova.clm", "data.frame")
  tab
}

print.anova.clm <-
  function(x, digits=max(getOption("digits") - 2, 3),
           signif.stars=getOption("show.signif.stars"), ...)
{
  if (!is.null(heading <- attr(x, "heading")))
    cat(heading, "\n")
  models <- attr(x, "models")
  print(models, right=FALSE)
  cat("\n")
  printCoefmat(x, digits=digits, signif.stars=signif.stars,
               tst.ind=4, cs.ind=NULL, # zap.ind=2, #c(1,5),
               P.values=TRUE, has.Pvalue=TRUE, na.print="", ...)
  return(invisible(x))
}

model.matrix.clm <- function(object, type = c("design", "B"), ...) {
    type <- match.arg(type)
    mf <- try(model.frame(object), silent=TRUE)
    if(inherits(mf, "try-error"))
        stop("Cannot extract model.matrix: refit model with 'model=TRUE'?")
### NOTE: we want to stop even if type="B" since the fullmf is needed
### in get_clmRho also and this way the error message is better.
    if(type == "design") {
        contr <- c(object$contrasts, object$S.contrasts,
                   object$nom.contrasts)
        design <- get_clmDesign(fullmf=object$model,
                                terms.list=terms(object, "all"),
                                contrasts=contr)
        keep <- c("X", "NOM", "S")
        select <- match(keep, names(design), nomatch=0)
        ans <- design[select]
    } else { ## if type == "B":
        env <- get_clmRho.clm(object)
        ans <- list(B1 = env$B1, B2 = env$B2)
        ans$S <- env$S ## may not exist
    }
    return(ans)
}

model.frame.clm <- function(formula, ...) {
### returns a model frame with *all* variables used for fitting.
    if(is.null(mod <- formula$model))
        stop("Cannot extract model.frame: refit model with 'model=TRUE'")
    else
        mod
}

coef.clm <- function(object, na.rm = FALSE, ...) {
  if(na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs)]
  }
  else
    object$coefficients
}

coef.summary.clm <- function(object, na.rm = FALSE, ...) {
  if(na.rm) {
    coefs <- object$coefficients
    coefs[!is.na(coefs[,1]), , drop=FALSE]
  }
  else
    object$coefficients
}

nobs.clm <- function(object, ...) object$nobs

terms.clm <-
    function(x, type=c("formula", "scale", "nominal", "all"), ...)
{
    type <- match.arg(type)
    term.nm <- c("terms", "S.terms", "nom.terms")
    Terms <- x[names(x) %in% term.nm]
    ind <- match(term.nm, names(Terms), 0L)
    Terms <- Terms[ind]
    names(Terms) <- c("formula", "scale", "nominal")[ind != 0]
    if(type == "all") return(Terms)
    if(!type %in% names(Terms))
        stop(gettextf("no terms object for '%s'", type))
    Terms[[type]]
}
