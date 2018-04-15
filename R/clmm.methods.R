## This file contains:
## Implementation of various methods for clmm objects.

formatVC <- function(varc, digits = max(3, getOption("digits") - 2))
### "format()" the 'VarCorr' matrix of the random effects -- for
### show()ing
### Borrowed from lme4/R/lmer.R with minor modifications.
{
    recorr <- lapply(varc, attr, "correlation")
    reStdDev <- lapply(varc, attr, "stddev")
    reLens <- unlist(lapply(reStdDev, length))
    nr <- sum(reLens)
    reMat <- array('', c(nr, 4),
		   list(rep.int('', nr),
			c("Groups", "Name", "Variance", "Std.Dev.")))
    reMat[1+cumsum(reLens)-reLens, 1] <- names(reLens)
    reMat[,2] <- unlist(lapply(varc, colnames))
    reMat[,3] <- format(unlist(reStdDev)^2, digits = digits)
    reMat[,4] <- format(unlist(reStdDev), digits = digits)
    if(any(reLens > 1)) {
	maxlen <- max(reLens)
	corr <- do.call("rBind", lapply(recorr, function(x, maxlen)
                                    {
                                        if(is.null(x)) return("")
                                        x <- as(x, "matrix")
                                        cc <- format(round(x, 3), nsmall = 3)
                                        cc[!lower.tri(cc)] <- ""
                                        nr <- dim(cc)[1]
                                        if (nr >= maxlen) return(cc)
                                        cbind(cc, matrix("", nr, maxlen-nr))
                                    },
                                        maxlen))
	colnames(corr) <- c("Corr", rep.int("", maxlen - 1))
        cbind(reMat, corr)
    } else reMat
}


varcov <-
    function(object, format=FALSE,
             digits=max(3, getOption("digits") - 2), ...)
### VarCorr method for model environments - should be the same for
### fitted model objects.
{
    ## Compute variance-covariance matrices of the random effects.
    res <- lapply(object$ST, function(st) {
        ## Variance-covariance matrix for the random effects:
        VC <- tcrossprod(st)
        ## Standard deviations:
        stddev <- sqrt(diag(VC))
        corr <- t(VC / stddev)/stddev
        attr(VC, "stddev") <- stddev
        ## correlation:
        if(NCOL(st) > 1) {
            diag(corr) <- 1
            attr(VC, "correlation") <- corr
        }
        VC
    })
    names(res) <- names(object$dims$nlev.re)
    if(format) noquote(formatVC(res, digits=digits)) else res
}

VarCorr <- function(x, ...) UseMethod("VarCorr")
VarCorr.clmm <- function(x, ...) varcov(x, ...)

print.clmm <-
  function(x, digits = max(3, getOption("digits") - 3), ...)
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
  cat("formula:", deparse(x$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  cat("\nRandom effects:\n")
  print(formatVC(varcov(x), digits=digits), quote=FALSE, ...)
  nlev.char <- paste(names(x$dims$nlev.gf), " ", x$dims$nlev.gf, sep="", collapse=",  ")
  cat("Number of groups: ", nlev.char, "\n")

  if(length(x$beta)) {
    cat("\nCoefficients:\n")
    print(x$beta, digits=digits, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(length(x$alpha) > 0) {
    cat("\nThresholds:\n")
    print(x$alpha, digits=digits, ...)
  }

  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  return(invisible(x))
}

vcov.clmm <- function(object, ...) vcov.clm(object, method="Cholesky")

summary.clmm <- function(object, correlation = FALSE, ...)
{
  if(is.null(object$Hessian))
    stop("Model needs to be fitted with Hess = TRUE")

  nfepar <- object$dims$nfepar
  coef <- matrix(0, nfepar, 4,
                 dimnames = list(names(object$coefficients[1:nfepar]),
                   c("Estimate", "Std. Error", "z value", "Pr(>|z|)")))
  coef[, 1] <- object$coefficients[1:nfepar]
  vc <- try(vcov(object), silent = TRUE)
  if(class(vc) == "try-error") {
    warning("Variance-covariance matrix of the parameters is not defined")
    coef[, 2:4] <- NaN
    if(correlation) warning("Correlation matrix is unavailable")
    object$condHess <- NaN
  }
  else {
    coef[, 2] <- sd <- sqrt(diag(vc)[1:nfepar])
    ## Cond is Inf if Hessian contains NaNs:
    object$condHess <-
      if(any(is.na(object$Hessian))) Inf
      else with(eigen(object$Hessian, only.values = TRUE),
                abs(max(values) / min(values)))
    coef[, 3] <- coef[, 1]/coef[, 2]
    coef[, 4] <- 2 * pnorm(abs(coef[, 3]), lower.tail=FALSE)
    if(correlation) ## {
      ## sd <- sqrt(diag(vc))
      object$correlation <- cov2cor(vc)
    ## (vc / sd) / rep(sd, rep(object$edf, object$edf))
  }
  object$info$cond.H <- formatC(object$condHess, digits=1, format="e")
  object$coefficients <- coef
  class(object) <- "summary.clmm"
  return(object)
}

print.summary.clmm <-
  function(x, digits = max(3, getOption("digits") - 3),
           signif.stars = getOption("show.signif.stars"), ...)
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
  cat("formula:", deparse(x$formula), fill=TRUE)
  if(!is.null(data.name <- x$call$data))
    cat("data:   ", deparse(data.name), fill=TRUE)
  if(!is.null(x$call$subset))
    cat("subset: ", deparse(x$call$subset), fill=TRUE)
  cat("\n")

  print(x$info, row.names=FALSE, right=FALSE)

  cat("\nRandom effects:\n")
  print(formatVC(varcov(x), digits=digits), quote=FALSE, ...)
  nlev.char <- paste(names(x$dims$nlev.gf), " ", x$dims$nlev.gf, sep="", collapse=",  ")
  cat("Number of groups: ", nlev.char, "\n")

  nbeta <- length(x$beta)
  nalpha <- length(x$alpha)
  if(nbeta > 0) {
    cat("\nCoefficients:\n")
    printCoefmat(x$coefficients[nalpha + 1:nbeta, , drop=FALSE],
                 digits=digits, signif.stars=signif.stars,
                 has.Pvalue=TRUE, ...)
  } else {
    cat("\nNo Coefficients\n")
  }
  if(nalpha > 0) { ## always true
    cat("\nThreshold coefficients:\n")
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

## anova.clmm <- function(object, ...)
##   anova.clm(object, ...)

anova.clmm <- function(object, ...) {
### This essentially calls anova.clm(object, ...), but the names of
### the models were not displayed correctly in the printed output
### unless the following dodge is enforced.
  mc <- match.call()
  arg.list <- as.list(mc)
  arg.list[[1]] <- NULL
  return(do.call(anova.clm, arg.list))
}

logLik.clmm <- function(object, ...)
  structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clmm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}

nobs.clmm <- function(object, ...) object$dims$nobs

### FIXME: define edf method
edf.clmm <- function(object, ...) object$dims$edf

## anova.clmm <- function(object, ...)
##   anova.clm(object, ...)

anova.clmm <- function(object, ...) {
### This essentially calls anova.clm(object, ...), but the names of
### the models were not displayed correctly in the printed output
### unless the following dodge is enforced.
  mc <- match.call()
  arg.list <- as.list(mc)
  arg.list[[1]] <- NULL
  return(do.call(anova.clm, arg.list))
}

logLik.clmm <- function(object, ...)
  structure(object$logLik, df = object$edf, class = "logLik")

extractAIC.clmm <- function(fit, scale = 0, k = 2, ...) {
  edf <- fit$edf
  c(edf, -2*fit$logLik + k * edf)
}

