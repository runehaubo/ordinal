## This file contains:
## profile and confint methods for clm objects.

profile.clm <-
  function(fitted, which.beta = seq_len(nbeta),
           which.zeta = seq_len(nzeta), alpha = 0.001,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### match and tests arguments and dispatch to .zeta and .beta
### functions for the actual profiling.

### which.[beta, zeta] - numeric or character vectors.

### Works for models with nominal and scale effects and for any number
### of aliased coefs.
{
  ## match and test arguments:
    if(any(is.na(diag(vcov(fitted)))))
        stop("Cannot get profile when vcov(fitted) contains NAs", call.=FALSE)
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
### BETA:
  beta.names <- names(fitted$beta) ## possible beta
  nbeta <- length(fitted$beta)
  if(is.character(which.beta))
    which.beta <- match(which.beta, beta.names, nomatch = 0)
  ## which.beta is a numeric vector
  if(!all(which.beta %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
### ZETA:
  zeta.names <- names(fitted$zeta) ## possible zeta
  nzeta <- length(fitted$zeta)
  if(is.character(which.zeta))
    which.zeta <- match(which.zeta, zeta.names, nomatch = 0)
  ## which.zeta is a numeric vector
  if(!all(which.zeta %in% seq_len(nzeta)))
    stop("invalid 'parm' argument")
  ## the actual profiling for beta and zeta par:
  prof.beta <- if(nbeta)
    profile.clm.beta(fitted, which.beta, alpha, max.steps, nsteps,
                     trace, step.warn, control, ...)
  else NULL
  prof.zeta <- if(nzeta)
    profile.clm.zeta(fitted, which.zeta, alpha, max.steps, nsteps,
                     trace, step.warn, control, ...)
  else NULL
  ## collect and return results:
  val <- structure(c(prof.beta, prof.zeta), original.fit = fitted)
  class(val) <- c("profile.clm")
  return(val)
}

profile.clm.beta <-
  function(fitted, which.beta, alpha = 0.001,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### which.beta is assumed to be a numeric vector
{
    lroot.max <- qnorm(1 - alpha/2)
    delta = lroot.max/nsteps
    nbeta <- length(fitted$beta)
    beta.names <- names(fitted$beta)
    nalpha <- length(fitted$alpha)
    orig.par <- c(fitted$alpha, fitted$beta)
    if(!is.null(zeta <- fitted$zeta)) {
        names(zeta) <- paste("sca", names(fitted$zeta), sep=".")
        orig.par <- c(orig.par, zeta)
    }
### NOTE: we need to update zeta.names to make names(orig.par)
### unique. This is needed to correctly construct the resulting
### par.vals matrix and to extract from it again.
    std.err <- coef(summary(fitted))[nalpha + 1:nbeta, "Std. Error"]
    if(any(is.na(std.err)))
        stop("Cannot profile model where standard errors are NA",
             call.=FALSE)
    ## results list:
    prof.list <- vector("list", length = length(which.beta))
    names(prof.list) <- beta.names[which.beta]
    ## get model matrices and model environment:
### NOTE: Fixing the fragile update approach:
    ## mf <- update(fitted, method = "model.frame")
    ## Need to subset by wts to make nrow(X) == nrow(B1)
    ## X <- with(mf, X[wts > 0, , drop=FALSE]) ## containing alias cols
    wts <- getWeights(model.frame(fitted))
    X <- model.matrix(fitted)$X[wts > 0, , drop=FALSE]
    rho <- get_clmRho(fitted)
    ## rho <- update(fitted, doFit = FALSE)
    orig <- as.list(rho)[c("B1", "B2", "o1", "o2")]
    rho$n.psi <- rho$n.psi - 1 ## needed for models with scale
    nalpha.clean <- sum(!fitted$aliased$alpha)
    par.clean <- orig.par[!is.na(orig.par)]
    ## which of which.beta are NA:
    alias.wb <- fitted$aliased$beta[which.beta]
    ## For each which.beta move up or down, fit the model and store the
    ## signed likelihood root statistic and parameter values:
    for(wb in which.beta) {
        if(alias.wb[wb == which.beta]) next ## ignore aliased coef
        rem <- nalpha.clean +
            (which.beta - cumsum(alias.wb))[wb == which.beta]
        par.wb <- matrix(coef(fitted), nrow = 1) ## MLE
        wb.name <- beta.names[wb]
        lroot.wb <- 0 ## lroot at MLE
        ## set variables in fitting environment:
        rho$B1 <- orig$B1[, -rem, drop=FALSE]
        rho$B2 <- orig$B2[, -rem, drop=FALSE]
        for(direction in c(-1, 1)) { ## move down or up
            if(trace) {
                message("\nParameter: ", wb.name,
                        c(" down", " up")[(direction + 1)/2 + 1])
                utils::flush.console()
            }
            ## reset starting values:
            rho$par <- par.clean[-rem]
            for(step in seq_len(max.steps)) {
                ## increment beta.i, offset and refit model without wb parameter:
                beta.i <- fitted$beta[wb] +
                    direction * step * delta * std.err[wb]
                new.off <- X[, 1+wb, drop=TRUE] * beta.i
                rho$o1 <- orig$o1 - new.off
                rho$o2 <- orig$o2 - new.off
                fit <- clm.fit.NR(rho, control)
                ## save likelihood root statistic:
                lroot <- -direction * sqrt(2*(fitted$logLik - fit$logLik))
                ## save lroot and pararameter values:
                lroot.wb <- c(lroot.wb, lroot)
                temp.par <- orig.par
                temp.par[names(fit$par)] <- fit$par
                temp.par[wb.name] <- beta.i
                par.wb <- rbind(par.wb, temp.par)
                ## break for loop if profile is far enough:
                if(abs(lroot) > lroot.max) break
            } ## end 'step in seq_len(max.steps)'
            ## test that lroot.max is reached and enough steps are taken:
            if(abs(lroot) < lroot.max)
                warning("profile may be unreliable for ", wb.name,
                        " because lroot.max was not reached for ",
                        wb, c(" down", " up")[(direction + 1)/2 + 1])
            if(step <= step.warn)
                warning("profile may be unreliable for ", wb.name,
                        " because only ", step, "\n  steps were taken ",
                        c("down", "up")[(direction + 1)/2 + 1])
        } ## end 'direction in c(-1, 1)'
        ## order lroot and par values and collect in a data.frame:
        lroot.order <- order(lroot.wb, decreasing = TRUE)
        prof.list[[wb.name]] <-
            structure(data.frame(lroot.wb[lroot.order]), names = "lroot")
        prof.list[[wb.name]]$par.vals <- par.wb[lroot.order, ]

        if(!all(diff(par.wb[lroot.order, wb.name]) > 0))
            warning("likelihood is not monotonically decreasing from maximum,\n",
                    "  so profile may be unreliable for ", wb.name)
    } ## end 'wb in which.beta'
    prof.list
}

profile.clm.zeta <-
  function(fitted, which.zeta, alpha = 0.001,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### which.zeta is assumed to be a numeric vector
{
  lroot.max <- qnorm(1 - alpha/2)
  delta = lroot.max/nsteps
  nzeta <- length(fitted$zeta)
  nbeta <- length(fitted$beta)
  zeta <- fitted$zeta
  names(zeta) <- zeta.names <- paste("sca", names(fitted$zeta), sep=".")
### NOTE: we need to update zeta.names to make names(orig.par)
### unique. This is needed to correctly construct the resulting
### par.vals matrix and to extract from it again.
  orig.par <- c(fitted$alpha, fitted$beta, zeta)
  nalpha <- length(fitted$alpha)
  std.err <- coef(summary(fitted))[nalpha+nbeta+1:nzeta, "Std. Error"]
  if(any(is.na(std.err)))
      stop("Cannot profile model where standard errors are NA",
           call.=FALSE)
  ## results list:
  prof.list <- vector("list", length = length(which.zeta))
  names(prof.list) <- names(zeta)[which.zeta]
  ## get model environment:
  rho <- get_clmRho(fitted)
  ## rho <- update(fitted, doFit = FALSE)
  S <- rho$S ## S without intercept
  Soff <- rho$Soff
  rho$k <- max(0, rho$k - 1)
  ab <- c(fitted$alpha, fitted$beta)
  ab.clean <- ab[!is.na(ab)]
  zeta.clean <- zeta[!fitted$aliased$zeta]
  ## which of which.zeta are NA:
  alias.wz <- fitted$aliased$zeta[which.zeta]
  ## For each which.zeta move up or down, fit the model and store the
  ## signed likelihood root statistic and parameter values:
  for(wz in which.zeta) {
    if(alias.wz[wz]) next ## ignore aliased coef
    ## rem: which columns of S to remove
    rem <- (which.zeta - cumsum(alias.wz))[wz]
    par.wz <- matrix(coef(fitted), nrow = 1) ## MLE
    wz.name <- zeta.names[wz]
    lroot.wz <- 0 ## lroot at MLE
    ## set variables in fitting environment:
    rho$S <- S[, -rem, drop=FALSE]
    for(direction in c(-1, 1)) { ## move down or up
      if(trace) {
        message("\nParameter: ", wz.name,
                c(" down", " up")[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      ## reset starting values:
      rho$par <- c(ab.clean, zeta.clean[-rem])
      ## rho$par <- coef(fitted, na.rm = TRUE)[-rem]
      for(step in seq_len(max.steps)) {
        ## increment zeta.i, offset and refit model without wz parameter:
        zeta.i <- zeta[wz] +
          direction * step * delta * std.err[wz]
        rho$Soff <- rho$sigma <- Soff * exp(S[, wz, drop=TRUE] * zeta.i)
### NOTE: Need to update sigma in addition to Soff since otherwise
### sigma isn't updated when k=0 (single scale par)
        fit <- clm.fit.NR(rho, control)
        ## save likelihood root statistic:
        lroot <- -direction * sqrt(2*(fitted$logLik - fit$logLik))
        ## save lroot and pararameter values:
        lroot.wz <- c(lroot.wz, lroot)
        temp.par <- orig.par
        temp.par[names(fit$par)] <- fit$par
        temp.par[wz.name] <- zeta.i
        par.wz <- rbind(par.wz, temp.par)
        ## break for loop if profile is far enough:
        if(abs(lroot) > lroot.max) break
      } ## end 'step in seq_len(max.steps)'
      ## test that lroot.max is reached and enough steps are taken:
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for ", wz.name,
                " because qnorm(1 - alpha/2) was not reached when profiling ",
                c(" down", " up")[(direction + 1)/2 + 1])
      if(step <= step.warn)
        warning("profile may be unreliable for ", wz.name,
                " because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'direction in c(-1, 1)'
    ## order lroot and par values and collect in a data.frame:
    ## lroot.order <- order(lroot.wz, decreasing = TRUE)
    lroot.order <- order(par.wz[, wz.name], decreasing = FALSE)
### FIXME: Need to change how values are ordered here. We should order
### with par.wz[, wz.name] instead of lroot.wz since if lroot.wz is
### flat, the order may be incorrect.
    prof.list[[wz.name]] <-
      structure(data.frame(lroot.wz[lroot.order]), names = "lroot")
    prof.list[[wz.name]]$par.vals <- par.wz[lroot.order, ]

    if(!all(diff(lroot.wz[lroot.order]) <= sqrt(.Machine$double.eps)))
        warning("likelihood is not monotonically decreasing from maximum,\n",
                "  so profile may be unreliable for ", wz.name)
} ## end 'wz in which.zeta'
  prof.list
}

## profile.sclm <- ## using clm.fit.env()
##   function(fitted, which.beta = seq_len(nbeta), alpha = 0.001,
##            max.steps = 50, nsteps = 8, trace = FALSE,
##            step.warn = 5, control = list(), ...)
## ### NOTE: seq_len(nbeta) works for nbeta = 0: numeric(0), while
## ### 1:nbeta gives c(1, 0).
##
## ### This is almost a copy of profile.clm2, which use clm.fit rather
## ### than clm.fit.env. The current implementation is the fastest, but
## ### possibly less readable.
## {
##   ## match and test arguments:
##   stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
##             alpha > 0 && alpha < 1)
##   stopifnot(round(max.steps) > round(nsteps))
##   stopifnot(round(nsteps) > round(step.warn))
##   stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
##   max.steps <- round(max.steps)
##   nsteps <- round(nsteps)
##   step.warn <- round(step.warn)
##   trace <- as.logical(trace)[1]
##   ## possible parameters on which to profile (including aliased coef):
##   beta.names <- names(fitted$beta)
##   nbeta <- length(fitted$beta)
##   if(is.character(which.beta))
##     which.beta <- match(which.beta, beta.names, nomatch = 0)
##   ## which.beta is a numeric vector
##   if(!all(which.beta %in% seq_len(nbeta)))
##     stop("invalid 'parm' argument")
##   stopifnot(length(which.beta) > 0)
##   std.err <- coef(summary(fitted))[-(1:length(fitted$alpha)),
##                                    "Std. Error"]
##   ## profile limit:
##   lroot.max <- qnorm(1 - alpha/2)
##   ## profile step length:
##   delta <- lroot.max / nsteps
##   ## results list:
##   prof.list <- vector("list", length = length(which.beta))
##   names(prof.list) <- beta.names[which.beta]
##   ## get model.frame:
##   X <- update(fitted, method = "model.frame")$X ## containing alias cols
##   rho <- update(fitted, doFit = FALSE)
##   orig <- as.list(rho)[c("B1", "B2", "o1", "o2")]
##   rho$n.psi <- rho$n.psi - 1
##   nalpha.clean <- sum(!fitted$aliased$alpha)
##   ## which of which.beta are NA:
##   alias.wb <- fitted$aliased$beta[which.beta]
##   ## For each which.beta move up or down, fit the model and store the
##   ## signed likelihood root statistic and parameter values:
##   for(wb in which.beta) {
##     if(alias.wb[wb]) next ## ignore aliased coef
##     rem <- nalpha.clean + (which.beta - cumsum(alias.wb))[wb]
##     par.wb <- matrix(coef(fitted), nrow = 1) ## MLE
##     wb.name <- beta.names[wb]
##     lroot.wb <- 0 ## lroot at MLE
##     ## set variables in fitting environment:
##     rho$B1 <- orig$B1[, -rem, drop=FALSE]
##     rho$B2 <- orig$B2[, -rem, drop=FALSE]
##     for(direction in c(-1, 1)) { ## move down or up
##       if(trace) {
##         message("\nParameter: ", wb.name,
##                 c(" down", " up")[(direction + 1)/2 + 1])
##         utils::flush.console()
##       }
##       ## reset starting values:
##       rho$par <- coef(fitted, na.rm = TRUE)[-rem]
##       ## rho$par <- orig.par[-wb.name]
##       for(step in seq_len(max.steps)) {
##         ## increment beta.i, offset and refit model without wb parameter:
##         beta.i <- fitted$beta[wb] +
##           direction * step * delta * std.err[wb]
##         new.off <- X[, 1+wb, drop=TRUE] * beta.i
##         rho$o1 <- orig$o1 - new.off
##         rho$o2 <- orig$o2 - new.off
##         fit <- clm.fit.env(rho, control)
##         ## save likelihood root statistic:
##         lroot <- -direction * sqrt(2*(fitted$logLik - fit$logLik))
##         ## save lroot and pararameter values:
##         lroot.wb <- c(lroot.wb, lroot)
##         temp.par <- coef(fitted)
##         temp.par[names(fit$par)] <- fit$par
##         temp.par[wb.name] <- beta.i
##         par.wb <- rbind(par.wb, temp.par)
##         ## break for loop if profile is far enough:
##         if(abs(lroot) > lroot.max) break
##       } ## end 'step in seq_len(max.steps)'
##       ## test that lroot.max is reached and enough steps are taken:
##       if(abs(lroot) < lroot.max)
##         warning("profile may be unreliable for ", wb.name,
##                 " because lroot.max was not reached for ",
##                 wb, c(" down", " up")[(direction + 1)/2 + 1])
##       if(step <= step.warn)
##         warning("profile may be unreliable for ", wb.name,
##                 " because only ", step, "\n  steps were taken ",
##                 c("down", "up")[(direction + 1)/2 + 1])
##     } ## end 'direction in c(-1, 1)'
##     ## order lroot and par. values and collect in a data.frame:
##     lroot.order <- order(lroot.wb, decreasing = TRUE)
##     prof.list[[wb.name]] <-
##       structure(data.frame(lroot.wb[lroot.order]), names = "lroot")
##     prof.list[[wb.name]]$par.vals <- par.wb[lroot.order, ]
##
##     if(!all(diff(par.wb[lroot.order, wb.name]) > 0))
##       warning("likelihood is not monotonically decreasing from maximum,\n",
##               "  so profile may be unreliable for ", wb.name)
##   } ## end 'wb in which.beta'
##   val <- structure(prof.list, original.fit = fitted)
##   class(val) <- c("profile.clm")
##   return(val)
## }

format.perc <- function(probs, digits)
### function lifted from stats:::format.perc to avoid using ':::'
    paste(format(100 * probs, trim = TRUE, scientific = FALSE,
                 digits = digits), "%")


confint.clm <-
  function(object, parm, level = 0.95,
           type = c("profile", "Wald"), trace = FALSE, ...)
### parm argument is ignored - use confint.profile for finer control.
{
  ## match and test arguments
  type <- match.arg(type)
  stopifnot(is.numeric(level) && length(level) == 1 &&
            level > 0 && level < 1)
  trace <- as.logical(trace)[1]
  if(!(missing(parm) || is.null(parm)))
    message("argument 'parm' ignored")
  ## Wald CI:
  if(type == "Wald") {
    a <- (1 - level)/2
    a <- c(a, 1 - a)
    pct <- format.perc(a, 3)
    fac <- qnorm(a)
    coefs <- coef(object)
    ses <- coef(summary(object))[, 2]
    ci <- array(NA, dim = c(length(coefs), 2L),
                dimnames = list(names(coefs), pct))
    ci[] <- coefs + ses %o% fac
    return(ci)
  }
  ## profile likelhood CI:
  if(trace) {
    message("Wait for profiling to be done...")
    utils::flush.console()
  }
  ## get profile:
  object <- profile(object, alpha = (1 - level)/4, trace = trace, ...)
  ## get and return CIs:
  confint(object, level = level, ...)
}

## confint.clm <-
##   function(object, parm = seq_len(npar), level = 0.95,
##            type = c("profile", "Wald"), trace = FALSE, ...)
## ### parm: a 2-list with beta and zeta?
## ### or args which.beta, which.zeta while parm is redundant?
##
## ### make waldci.clm(object, which.alpha, which.beta, which.zeta, level
## ### = 0.95) ??
## {
##   ## match and test arguments
##   type <- match.arg(type)
##   stopifnot(is.numeric(level) && length(level) == 1 &&
##             level > 0 && level < 1)
##   trace <- as.logical(trace)[1]
##   mle <- object$beta
##   if(!is.null(zeta <- object$zeta)) {
##     names(zeta) <- paste("sca", names(zeta), sep=".")
##     mle <- c(mle, zeta)
##   }
##   npar <- length(mle)
##   beta.names <- names(mle)
##   if(is.character(parm)) stop("parm should be numeric")
##   ## parm <- match(parm, names(c(object$beta, object$zeta))), nomatch = 0)
##   if(!all(parm %in% seq_len(npar))) stop("invalid 'parm' argument")
##   stopifnot(length(parm) > 0)
##   ## Wald CI:
##   if(type == "Wald")
##     return(waldci.clm(object, parm, level))
##   ## return(confint.default(object = object, parm = beta.names[parm],
##   ##                        level = level))
##   ## profile likelhood CI:
##   if(trace) {
##     message("Waiting for profiling to be done...")
##     utils::flush.console()
##   }
##   ## get profile:
## ### Edit these calls:
##   object <- profile(object, which.beta = beta.names[parm],
##                     alpha = (1 - level)/4, trace = trace, ...)
##   ## get and return CIs:
##   confint(object, parm = beta.names[parm], level = level, ...)
## }

confint.profile.clm <-
  function(object, parm = seq_len(nprofiles), level = 0.95, ...)
### parm index elements of object (the list of profiles)
### each par.vals matrix of each profile will have
### sum(!unlist(of$aliased)) columns.
{
  ## match and test arguments:
  stopifnot(is.numeric(level) && length(level) == 1 &&
            level > 0 && level < 1)
  of <- attr(object, "original.fit")
  prof.names <- names(object)
  nprofiles <- length(prof.names)
  if(is.character(parm))
### Allow character here?
    parm <- match(parm, prof.names, nomatch = 0)
  if(!all(parm %in% seq_len(nprofiles)))
    stop("invalid 'parm' argument")
  stopifnot(length(parm) > 0)
  ## prepare CI:
  a <- (1-level)/2
  a <- c(a, 1-a)
  pct <- paste(round(100*a, 1), "%")
  ci <- array(NA, dim = c(length(parm), 2),
              dimnames = list(prof.names[parm], pct))
  cutoff <- qnorm(a)
  ## compute CI from spline interpolation of the likelihood profile:
  for(pr.name in prof.names[parm]) {
    if(is.null(pro <- object[[ pr.name ]])) next
    sp <- spline(x = pro[, "par.vals"][, pr.name], y = pro[, 1]) ## OBS
    ci[pr.name, ] <- approx(sp$y, sp$x, xout = rev(cutoff))$y
  }
  ## do not drop(ci) because rownames are lost for single coef cases:
  return(ci)
}

plot.profile.clm <-
    function(x, which.par = seq_len(nprofiles), level = c(0.95, 0.99),
             Log = FALSE, relative = TRUE, root = FALSE, fig = TRUE,
             approx = root, n = 1e3,
             ask = prod(par("mfcol")) < length(which.par) &&
             dev.interactive(), ..., ylim = NULL)
{
  ## match and test arguments:
  stopifnot(is.numeric(level) && all(level > 0) &&
            all(level < 1))
  stopifnot(n == round(n) && n > 0)
  Log <- as.logical(Log)[1]
  relative <- as.logical(relative)[1]
  root <- as.logical(root)[1]
  fig <- as.logical(fig)[1]
  approx <- as.logical(approx)[1]
  of <- attr(x, "original.fit")
  mle <- of$beta
  if(!is.null(zeta <- of$zeta)) {
    names(zeta) <- paste("sca", names(zeta), sep=".")
    mle <- c(mle, zeta)
  }
  prof.names <- names(x)
  nprofiles <- length(prof.names)
  if(is.character(which.par))
    which.par <- match(which.par, prof.names, nomatch = 0)
  if(!all(which.par %in% seq_len(nprofiles)))
    stop("invalid 'which.par' argument")
  stopifnot(length(which.par) > 0)
  ML <- of$logLik
  ## prepare return value:
  which.names <- prof.names[which.par]
  spline.list <- vector("list", length(which.par))
  names(spline.list) <- which.names
  if(approx) {
    std.err <- coef(summary(of))[-(1:length(of$alpha)), 2]
    names(std.err) <- names(mle)
  }
  ## aks before "over writing" the plot?
  if(ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  ## for each pm make the appropriate plot:
  for(pr.name in prof.names[which.par]) {
    ## confidence limits:
    lim <- sapply(level, function(x)
                  exp(-qchisq(x, df=1)/2) )
    if(is.null(pro <- x[[ pr.name ]])) next
    sp <- spline(x=pro[, "par.vals"][, pr.name], y=pro[, 1], n=n)
    if(approx) y.approx <- (mle[pr.name] - sp$x) / std.err[pr.name]
    if(root) {
      ylab <- "profile trace"
      lim <- c(-1, 1) %o% sqrt(-2 * log(lim))
      sp$y <- -sp$y
      if(approx) y.approx <- -y.approx
    } else { ## !root:
      sp$y <- -sp$y^2/2
      if(approx) y.approx <- -y.approx^2/2
      if(relative && !Log) {
        sp$y <- exp(sp$y)
        if(approx) y.approx <- exp(y.approx)
        ylab <- "Relative profile likelihood"
        if(missing(ylim)) ylim <- c(0, 1)
      }
      if(relative && Log) {
        ylab <- "Relative profile log-likelihood"
        lim <- log(lim)
      }
      if(!relative && Log) {
        sp$y <- sp$y + ML
        if(approx) y.approx <- y.approx + ML
        ylab <- "Profile log-likelihood"
        lim <- ML + log(lim)
      }
      if(!relative && !Log) {
        sp$y <- exp(sp$y + ML)
        if(approx) y.approx <- exp(y.approx + ML)
        ylab <- "Profile likelihood"
        lim <- exp(ML + log(lim))
      }
    }
    spline.list[[ pr.name ]] <- sp

    if(fig) { ## do the plotting:
      plot(sp$x, sp$y, type = "l", ylim = ylim,
           xlab = pr.name, ylab = ylab, ...)
      abline(h = lim)
      if(approx) lines(sp$x, y.approx, lty = 2)
      if(root)  points(mle[pr.name], 0, pch = 3)
    }
  }
  attr(spline.list, "limits") <- lim
  invisible(spline.list)
}

profileAlt.clm <- ## using clm.fit()
  function(fitted, which.beta = seq_len(nbeta), alpha = 0.01,
           max.steps = 50, nsteps = 8, trace = FALSE,
           step.warn = 5, control = list(), ...)
### NOTE: seq_len(nbeta) works for nbeta = 0: numeric(0), while
### 1:nbeta gives c(1, 0).

### args:
### alpha - The likelihood is profiled in the 100*(1-alpha)%
###   confidence region as determined by the profile likelihood
### max.steps - the maximum number of profile steps in each direction
### nsteps - the approximate no. steps determined by the quadratic
### approximation to the log-likelihood function
### trace - if trace > 0 information of progress is printed
### step.warn - a warning is issued if the profile in each direction
###   contains less than step.warn steps (due to lack of precision).
{
  ## match and test arguments:
  stopifnot(is.numeric(alpha) && length(alpha) == 1 &&
            alpha > 0 && alpha < 1)
  stopifnot(round(max.steps) > round(nsteps))
  stopifnot(round(nsteps) > round(step.warn))
  stopifnot(round(nsteps) > 0 && round(step.warn) >= 0)
  max.steps <- round(max.steps)
  nsteps <- round(nsteps)
  step.warn <- round(step.warn)
  trace <- as.logical(trace)[1]
  beta.names <- names(fitted$beta)
  nbeta <- length(fitted$beta)
  if(is.character(which.beta))
    which.beta <- match(which.beta, beta.names, nomatch = 0)
  if(!all(which.beta %in% seq_len(nbeta)))
    stop("invalid 'parm' argument")
   stopifnot(length(which.beta) > 0)
  ## Extract various things from the original fit:
  orig.par <- coef(fitted) ## c(alpha, beta)
  beta0 <- fitted$beta ## regression coef.
  nalpha <- length(fitted$alpha) ## no. threshold coef.
  nbeta <- length(beta0)
  beta.names <- names(beta0)
  orig.logLik <- fitted$logLik
  std.err <- coef(summary(fitted))[-(1:nalpha), "Std. Error"]
  link <- fitted$link
  threshold <- fitted$threshold
  ## profile limit:
  lroot.max <- qnorm(1 - alpha/2)
  ## profile step length:
  delta <- lroot.max / nsteps
  ## results list:
  prof.list <- vector("list", length = length(which.beta))
  names(prof.list) <- beta.names[which.beta]
  ## get model.frame:
### NOTE: Attempting the following fix for a safer extraction of
### model-design-objects:
  ## mf <- update(fitted, method = "model.frame")
  contr <- c(fitted$contrasts, fitted$S.contrasts, fitted$nom.contrasts)
  mf <- get_clmDesign(fitted$model, fitted$terms.list, contr)
  y <- mf$y
  X <- mf$X
  wts <- mf$wts
  orig.off <- mf$off
  ## For each which.beta move up or down, fit the model and store the
  ## signed likelihood root statistic and parameter values:
  for(wb in which.beta) {
    par.wb <- matrix(orig.par, nrow = 1) ## MLE
    wb.name <- beta.names[wb]
    lroot.wb <- 0 ## lroot at MLE
    X.wb <- X[, -(1+wb), drop=FALSE]
    for(direction in c(-1, 1)) { ## move down or up
      if(trace) {
        message("\nParameter: ", wb.name,
                c(" down", " up")[(direction + 1)/2 + 1])
        utils::flush.console()
      }
      ## (re)set starting values:
      start <- orig.par[-(nalpha + wb)]
      for(step in seq_len(max.steps)) {
        ## increment offset and refit model without wb parameter:
        beta.i <- beta0[wb] + direction * step * delta * std.err[wb]
        new.off <- orig.off + X[, 1+wb, drop=TRUE] * beta.i
        fit <- clm.fit(y=y, X=X.wb,
                       weights=wts, offset=new.off,
                       control=control, start=start, link=link,
                       threshold=threshold)
        ## save likelihood root statistic:
        lroot <- -direction * sqrt(2*(fitted$logLik - fit$logLik))
        ## save lroot and pararameter values:
        lroot.wb <- c(lroot.wb, lroot)
        temp.par <- orig.par
        temp.par[names(fit$par)] <- fit$par
        temp.par[wb.name] <- beta.i
        par.wb <- rbind(par.wb, temp.par)
        ## update starting values:
        start <- fit$par
        ## break for loop if profile is far enough:
        if(abs(lroot) > lroot.max) break
      } ## end 'step in seq_len(max.steps)'
      ## test that lroot.max is reached and enough steps are taken:
      if(abs(lroot) < lroot.max)
        warning("profile may be unreliable for ", wb.name,
                " because lroot.max was not reached for ",
                wb, c(" down", " up")[(direction + 1)/2 + 1])
      if(step <= step.warn)
        warning("profile may be unreliable for ", wb.name,
                " because only ", step, "\n  steps were taken ",
                c("down", "up")[(direction + 1)/2 + 1])
    } ## end 'direction in c(-1, 1)'
    ## order lroot and par. values and collect in a data.frame:
    lroot.order <- order(lroot.wb, decreasing = TRUE)
    prof.list[[wb.name]] <-
      structure(data.frame(lroot.wb[lroot.order]), names = "lroot")
    prof.list[[wb.name]]$par.vals <- par.wb[lroot.order, ]

    if(!all(diff(par.wb[lroot.order, wb.name]) > 0))
      warning("likelihood is not monotonically decreasing from maximum,\n",
              "  so profile may be unreliable for ", wb.name)
  } ## end 'wb in which.beta'
  val <- structure(prof.list, original.fit = fitted)
  class(val) <- c("profile.clm")
  return(val)
}
