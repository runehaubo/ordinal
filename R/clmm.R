## This file contains:
## Implementation of Cumulative Link Mixed Models in clmm().

if(getRversion() >= '2.15.1')
    utils::globalVariables(c("ths", "link", "threshold", "optRes",
                             "neval", "Niter", "tJac", "y.levels"))

clmm <-
  function(formula, data, weights, start, subset,
           na.action, contrasts, Hess = TRUE, model = TRUE,
           link = c("logit", "probit", "cloglog", "loglog",
             "cauchit"), ##, "Aranda-Ordaz", "log-gamma"), ## lambda,
           doFit = TRUE, control = list(), nAGQ = 1L,
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)
{
### Extract the matched call and initial testing:
  mc <- match.call(expand.dots = FALSE)
### FIXME: Possibly call clm() when there are no random effects?
  link <- match.arg(link)
  threshold <- match.arg(threshold)
  if(missing(formula))  stop("Model needs a formula")
  if(missing(contrasts)) contrasts <- NULL
  ## set control parameters:
  control <- getCtrlArgs(control, list(...))
  nAGQ <- as.integer(round(nAGQ))
  formulae <- clmm.formulae(formula=formula)
  ## mf, y, X, wts, off, terms:
  frames <- clmm.frames(modelcall=mc, formulae=formulae, contrasts)
### FIXME: What should 'method="model.frame"' return? Do we want Zt
### included here as well?
  if(control$method == "model.frame") return(frames)
  ## Test rank deficiency and possibly drop some parameters:
  ## X is guarantied to have an intercept at this point.
  frames$X <- drop.coef(frames$X, silent=FALSE)
  ## Compute the transpose of the Jacobian for the threshold function,
  ## tJac and the names of the threshold parameters, alpha.names:
  ths <- makeThresholds(levels(frames$y), threshold)
  ## Set rho environment:
  rho <- with(frames, {
    clm.newRho(parent.frame(), y=y, X=X, weights=wts,
               offset=off, tJac=ths$tJac) })
  ## compute grouping factor list, and Zt and ST matrices:
  retrms <- getREterms(frames = frames, formulae$formula)
### FIXME: save (the evaluated) formula in frames, so we only need the
### frames argument to getREterms() ?
  use.ssr <- (retrms$ssr && !control$useMatrix)

  ## Set inverse link function and its two first derivatives (pfun,
  ## dfun and gfun) in rho:
  setLinks(rho, link)

  ## Compute list of dimensions for the model fit:
  rho$dims <- getDims(frames=frames, ths=ths, retrms=retrms)

  ## Update model environment with r.e. information:
  if(use.ssr) {
      rho.clm2clmm.ssr(rho=rho, retrms = retrms, ctrl=control$ctrl)
      ## Set starting values for the parameters:
      if(missing(start)) start <- c(fe.start(frames, link, threshold), 0)
      rho$par <- start
      nbeta <- rho$nbeta <- ncol(frames$X) - 1 ## no. fixef parameters
      nalpha <- rho$nalpha <- ths$nalpha ## no. threshold parameters
      ntau <- rho$ntau <- length(retrms$gfList) ## no. variance parameters
      stopifnot(is.numeric(start) &&
                length(start) == (nalpha + nbeta + ntau))
  } else {
      rho.clm2clmm(rho=rho, retrms=retrms, ctrl=control$ctrl)
      if(missing(start)) {
          rho$fepar <- fe.start(frames, link, threshold)
          rho$ST <- STstart(rho$ST)
          start <- c(rho$fepar, ST2par(rho$ST))
      } else {
          stopifnot(is.list(start) && length(start) == 2)
          stopifnot(length(start[[1]]) == rho$dims$nfepar)
          stopifnot(length(start[[2]]) == rho$dims$nSTpar)
          rho$fepar <- as.vector(start[[1]])
          rho$ST <- par2ST(as.vector(start[[2]]), rho$ST)
      }
  }
### FIXME: set starting values in a more elegant way.

  ## Set AGQ parameters:
  set.AGQ(rho, nAGQ, control, use.ssr)

  ## Possibly return the environment, rho without fitting:
  if(!doFit)  return(rho)

  ## Fit the clmm:
  fit <-
      if(use.ssr) clmm.fit.ssr(rho, control = control$optCtrl,
                               method=control$method, Hess)
      else clmm.fit.env(rho, control = control$optCtrl,
                        method=control$method, Hess)

  ## Modify and return results:
  fit$nAGQ <- nAGQ
  fit$link <- link
  fit$start <- start
  fit$threshold <- threshold
  fit$call <- match.call()
  fit$formula <- formulae$formula
  fit$gfList <- retrms$gfList
  fit$control <- control
  res <- clmm.finalize(fit=fit, frames=frames, ths=ths, use.ssr)

  ## add model.frame to results list?
  if(model) res$model <- frames$mf

  return(res)
}

clmm.formulae <- function(formula) {
    ## Evaluate the formula in the enviroment in which clmm was called
    ## (parent.frame(2)) to get it evaluated properly:
    form <- eval.parent(formula, 2)
    ## get the environment of the formula. If this does not have an
    ## environment (it could be a character), then use the calling environment.
    form.envir <-
        if(!is.null(env <- environment(form))) env
        else parent.frame(2)
    ## ensure 'formula' is a formula-object:
    form <- try(formula(deparse(form), env = form.envir), silent=TRUE)
    ## report error if the formula cannot be interpreted
    if(class(form) == "try-error")
        stop("unable to interpret 'formula'")
    environment(form) <- form.envir
    ## Construct a formula with all (fixed and random) variables
    ## (fullForm) and a formula with only fixed-effects variables
    ## (fixedForm):
    fixedForm <- nobars(form) ## ignore terms with '|'
    fullForm <- subbars(form)      # substitute `+' for `|'
    ## Set the appropriate environments:
    environment(fullForm) <- environment(fixedForm) <-
        environment(form) <- form.envir
    list(formula = form, fullForm = fullForm, fixedForm = fixedForm)
}

clmm.frames <- function(modelcall, formulae, contrasts) {
    ## Extract full model.frame (fullmf):
    m <- match(c("data", "subset", "weights", "na.action", "offset"),
               names(modelcall), 0)
    mf <- modelcall[c(1, m)]
    mf$formula <- formulae$fullForm
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    fixedmf <- mf ## save call for later modification and evaluation
    fullmf <- eval(mf, envir = parent.frame(2)) ## '2' to get out of
    ## clmm.frames and clmm
### FIXME: What if data is a matrix?
    fixedmf$formula <- formulae$fixedForm
    fixedmf <- eval(fixedmf, envir = parent.frame(2))
    attr(fullmf, "terms") <- attr(fixedmf, "terms")
    ## return:
    list(mf = fullmf,
         y = getY(fullmf),
         X = getX(fullmf, fixedmf, contrasts),
         wts = getWeights(fullmf),
         off = getOffsetStd(fullmf),
         terms = attr(fixedmf, "terms")
         )
}

getY <- function(mf) {
### Extract model response:
    y <- model.response(mf)
    if(!is.factor(y)) stop("response needs to be a factor")
    y
}

getX <- function(fullmf, fixedmf, contrasts) {
    fixedTerms <- attr(fixedmf, "terms")
    X <- model.matrix(fixedTerms, fullmf, contrasts)
    n <- nrow(X)
    ## remove intercept from X:
    Xint <- match("(Intercept)", colnames(X), nomatch = 0)
    if(Xint <= 0) {
        X <- cbind("(Intercept)" = rep(1, n), X)
        warning("an intercept is needed and assumed")
    } ## intercept in X is garanteed.
    X
}

getZt <- function(retrms) {
    ZtList <- lapply(retrms, '[[', "Zt")
    Zt <- do.call(rBind, ZtList)
    Zt@Dimnames <- vector("list", 2)
    Zt
}

getREterms <- function(frames, formula) {
### NOTE: Need to parse mf - not just fullmf because we need the model
### fits for an identifiability check below.
    ## fullmf <- frames$mf
### FIXME: Should we have
    fullmf <- droplevels(with(frames, mf[wts > 0, ]))
    barlist <- expandSlash(findbars(formula[[3]]))
### FIXME: make sure 'formula' is appropriately evaluated and returned
### by clmm.formulae
    if(!length(barlist)) stop("No random effects terms specified in formula")
    term.names <- unlist(lapply(barlist, function(x) deparse(x)))
    names(barlist) <- unlist(lapply(barlist, function(x) deparse(x[[3]])))
### NOTE: Deliberately naming the barlist elements by grouping factors
### and not by r.e. terms.
    ## list of grouping factors for the random terms:
    rel <- lapply(barlist, function(x) {
        ff <- eval(substitute(as.factor(fac)[,drop = TRUE],
                              list(fac = x[[3]])), fullmf)
        ## per random term transpose indicator matrix:
        Zti <- as(ff, "sparseMatrix")
        ## per random term model matrix:
        mm <- model.matrix(eval(substitute(~ expr,
                                           list(expr = x[[2]]))), fullmf)
        Zt = do.call(rBind, lapply(seq_len(ncol(mm)), function(j) {
            Zti@x <- mm[,j]
            Zti } ))
### FIXME: can we drop rows from Zt when g has missing values in terms
### of the form (1 + g | f)?
        ST <- matrix(0, ncol(mm), ncol(mm),
                     dimnames = list(colnames(mm), colnames(mm)))
        list(f = ff, Zt = Zt, ST = ST)
### FIXME: return the i'th element of Lambda here.
    })
    ## For each r.e. term, test if Z has more columns than rows to detect
    ## unidentifiability:
    for(i in seq_along(barlist)) {
        Zti <- rel[[i]][["Zt"]]
        if(nrow(Zti) > ncol(Zti) ||
           (all(frames$wts == 1) && nrow(Zti) == ncol(Zti)))
            stop(gettextf("no. random effects (=%d) >= no. observations (=%d) for term: (%s)",
                          nrow(Zti), ncol(Zti), term.names[i]), call.=FALSE)
    }
    ## Test if total no. random effects >= total nobs:
    q <- sum(sapply(rel, function(x) nrow(x$Zt)))
    if(all(frames$wts == 1) && q >= nrow(fullmf))
        stop(gettextf("no. random effects (=%d) >= no. observations (=%d)",
                      q, nrow(fullmf)), call.=FALSE)
### NOTE: q > nrow(fullmf) is (sometimes) allowed if some frames$wts > 1
###
### NOTE: if all(frames$wts == 1) we cannot have observation-level
### random effects so we error if nrow(Zti) >= ncol(Zti)
###
### FIXME: Could probably also throw an error if q >= sum(frames$wts),
### but I am not sure about that.
###
### FIXME: It would be better to test the rank of the Zt matrix, but
### also computationally more intensive.
###
### FIXME: If the model is nested (all gr.factors are nested), then
### order the columns of Zt, such that they come in blocks
### corresponding to the levels of the coarsest grouping factor. Each
### block of Zt-columns contain first the j'th level of the 1st gr.fac.
### followed by columns for the 2nd gr.fac.
###
    ## single simple random effect on the intercept?
    ssr <- (length(barlist) == 1 && as.character(barlist[[1]][[2]]) == "1")
    ## order terms by decreasing number of levels in the factor but don't
    ## change the order if this is already true:
    nlev <- sapply(rel, function(re) nlevels(re$f))
    if (any(diff(nlev)) > 0) rel <- rel[rev(order(nlev))]
    nlev <- nlev[rev(order(nlev))]
    ## separate r.e. terms from the factor list:
    retrms <- lapply(rel, "[", -1)
    names(retrms) <- NULL
    ## list of grouping factors:
    gfl <- lapply(rel, "[[", "f")
    ## which r.e. terms are associated with which grouping factors:
    attr(gfl, "assign") <- seq_along(gfl)
    ## only save unique g.f. and update assign attribute:
    fnms <- names(gfl)
    ## check for repeated factors:
    if (length(fnms) > length(ufn <- unique(fnms))) {
        ## check that the lengths of the number of levels coincide
        gfl <- gfl[match(ufn, fnms)]
        attr(gfl, "assign") <- match(fnms, ufn)
        names(gfl) <- ufn
    }
    ## test that all variables for the random effects are factors and
    ## have at least 3 levels:
    stopifnot(all(sapply(gfl, is.factor)))
    stopifnot(all(sapply(gfl, nlevels) > 2))
    ## no. r.e. per level for each of the r.e. terms
    qi <- unlist(lapply(rel, function(re) ncol(re$ST)))
    stopifnot(q == sum(nlev * qi))
    dims <- list(n = nrow(fullmf), ## no. observations
                 nlev.re = nlev, ## no. levels for each r.e. term
                 nlev.gf = sapply(gfl, nlevels), ## no. levels for each grouping factor
                 qi = qi,
                 nretrms = length(rel), ## no. r.e. terms
                 ngf = length(gfl), ## no. unique grouping factors
                 ## total no. random effects:
                 q = sum(nlev * qi), ## = sum(sapply(rel, function(re) nrow(re$Zt)))
                 ## no. r.e. var-cov parameters:
                 nSTpar = sum(sapply(qi, function(q) q * (q + 1) / 2))
                 )
    ## c(retrms=retrms, list(gfList = gfl, dims = dims, ssr = ssr))
    list(retrms=retrms, gfList = gfl, dims = dims, ssr = ssr)
}

fe.start <- function(frames, link, threshold) {
    ## get starting values from clm:
    fit <- with(frames,
                clm.fit(y=y, X=X, weights=wts, offset=off, link=link,
                        threshold=threshold))
    unname(coef(fit))
}

getDims <- function(frames, ths, retrms)
### Collect and compute all relevant dimensions in a list
{
    dims <- retrms$dims ## n is also on retrms$dims
    dims$n <- sum(frames$wts > 0)
    dims$nbeta <- ncol(frames$X) - 1
    dims$nalpha <- ths$nalpha
    dims$nfepar <- dims$nalpha + dims$nbeta
    dims
}

rho.clm2clmm <- function(rho, retrms, ctrl)
### update environment, rho returned by clm.newRho().
{
### FIXME: write default list of control arguments?
    ## control arguments are used when calling update.u(rho)
    rho$ctrl = ctrl
    ## compute Zt design matrix:
    rho$Zt <- getZt(retrms$retrms)
    rho$ST <- lapply(retrms$retrms, `[[`, "ST")
    rho$allST1 <- all(sapply(rho$ST, ncol) == 1)
    ## Lambda <- getLambda(rho$ST, rho$dims$nlev.re)
    ## Vt <- crossprod(Lambda, rho$Zt)
    ## rho$L <- Cholesky(tcrossprod(Vt),
    ##                   LDL = TRUE, super = FALSE, Imult = 1)
    rho$L <- Cholesky(tcrossprod(crossprod(getLambda(rho$ST, rho$dims$nlev.re), rho$Zt)),
                      LDL = TRUE, super = FALSE, Imult = 1)
    rho$Niter <- 0L ## no. conditional mode updates
    rho$neval <- 0L ## no. evaluations of the log-likelihood function
    rho$u <- rho$uStart <- rep(0, rho$dims$q)
    rho$.f <- if(package_version(packageDescription("Matrix")$Version) >
                 "0.999375-30") 2 else 1
}

getLambda <- function(ST, nlev) {
### ST: a list of ST matrices
### nlev: a vector of no. random effects levels
    .local <- function(ST, nlev) {
        if(ncol(ST) == 1) .symDiagonal(n=nlev,
               x = rep(as.vector(ST[1, 1]), nlev)) else
        kronecker(as(ST, "sparseMatrix"), .symDiagonal(n=nlev))
        ## This would make sense if the columns in Z (rows in Zt) were ordered differently:
        ## kronecker(Diagonal(n=nlev), ST)
### NOTE: .symDiagonal() appears to be faster than Diagonal() here.
    }
    stopifnot(length(ST) == length(nlev))
    res <- if(length(ST) == 1) .local(ST[[1]], nlev) else
    .bdiag(lapply(seq_along(ST), function(i) .local(ST[[i]], nlev[i])))
    ## coerce to diagonal matrix if relevant:
    if(all(sapply(ST, ncol) == 1)) as(res, "diagonalMatrix") else
    as(res, "CsparseMatrix")
### QUESTION: Are there any speed gains by coerce'ing Lambda to
### 'diagonalMatrix' or 'CsparseMatrix'?
### QUESTION: What is the best way to form the kronecker product in .local()?
}

getNLA <- function(rho, par, which=rep(TRUE, length(par))) {
### negative log-likelihood by the Laplace approximation
    if(!missing(par)) {
        setPar.clmm(rho, par, which)
        if(any(!is.finite(par)))
            stop(gettextf(paste(c("Non-finite parameters not allowed:",
                                  formatC(par, format="g")), collapse=" ")))
    }
    rho$neval <- rho$neval + 1L
    if(!update.u(rho)) return(Inf)
    if(any(rho$D < 0)) return(Inf)
    logDetD <- c(suppressWarnings(determinant(rho$L)$modulus)) -
        rho$dims$q * log(2*pi) / 2
    rho$nll + logDetD
}

nll.u <- function(rho) { ## negative log-likelihood
    if(rho$allST1) { ## are all ST matrices scalars?
        rho$varVec <- rep.int(unlist(rho$ST), rho$dims$nlev.re)
        b.expanded <- as.vector(crossprod(rho$Zt, rho$varVec * rho$u))
### NOTE: Working with Lambda when it is diagonal will slow things
### down significantly.
    } else {
        rho$ZLt <- crossprod(getLambda(rho$ST, rho$dims$nlev.re), rho$Zt)
        b.expanded <- as.vector(crossprod(rho$ZLt, rho$u))
    }
    rho$eta1Fix <- drop(rho$B1 %*% rho$fepar)
    rho$eta2Fix <- drop(rho$B2 %*% rho$fepar)
    rho$eta1 <- as.vector(rho$eta1Fix - b.expanded + rho$o1)
    rho$eta2 <- as.vector(rho$eta2Fix - b.expanded + rho$o2)
    rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
    if(any(!is.finite(rho$fitted)) || any(rho$fitted <= 0))
        nll <- Inf
    else
        nll <- -sum(rho$wts * log(rho$fitted)) -
            sum(dnorm(x=rho$u, mean=0, sd=1, log=TRUE))
    nll
}

nllFast.u <- function(rho) { ## negative log-likelihood
  ## Does not update X %*% beta - fixed effect part.
    if(rho$allST1) {
        rho$varVec <- rep.int(unlist(rho$ST), rho$dims$nlev.re)
        b.expanded <- as.vector(crossprod(rho$Zt, rho$varVec * rho$u))
    } else {
        rho$ZLt <- crossprod(getLambda(rho$ST, rho$dims$nlev.re), rho$Zt)
        b.expanded <- as.vector(crossprod(rho$ZLt, rho$u))
    }
  rho$eta1 <- as.vector(rho$eta1Fix - b.expanded + rho$o1)
  rho$eta2 <- as.vector(rho$eta2Fix - b.expanded + rho$o2)
  rho$fitted <- getFittedC(rho$eta1, rho$eta2, rho$link)
  if(any(!is.finite(rho$fitted)) || any(rho$fitted <= 0))
    nll <- Inf
  else
    nll <- -sum(rho$wts * log(rho$fitted)) -
      sum(dnorm(x=rho$u, mean=0, sd=1, log=TRUE))
  nll
}

grad.u <- function(rho){ ## gradient of nll wrt. u (random effects)
### should only be called with up to date values of eta1, eta2, par
    ## compute phi1:
    rho$p1 <- rho$dfun(rho$eta1)
    rho$p2 <- rho$dfun(rho$eta2)
    rho$wtpr <- rho$wts/rho$fitted
    phi1 <- as.vector(rho$wtpr * (rho$p1 - rho$p2))
    if(rho$allST1)
        (rho$Zt %*% phi1) * rho$varVec + rho$u
    else
        rho$ZLt %*% phi1 + rho$u
}

hess.u <- function(rho) { ## Hessian of nll wrt. u (random effects)
### should only be called with up-to-date values of eta1, eta2, par,
### p1, p2
    g1 <- rho$gfun(rho$eta1) ## does not need to be saved in rho
    g2 <- rho$gfun(rho$eta2) ## does not need to be saved in rho
    phi2 <- rho$wts * ( ((rho$p1 - rho$p2) / rho$fitted)^2 -
                       ( (g1 - g2) / rho$fitted) )
    ## This may happen if the link function [pfun, dfun and gfun]
    ## evaluates its arguments inaccurately:
    if(any(phi2 < 0))  return(FALSE)
    if(rho$allST1)
        Vt <- crossprod(Diagonal(x = rho$varVec),
                        tcrossprod(rho$Zt, Diagonal(x = sqrt(phi2))))
    else
        Vt <- rho$ZLt %*% Diagonal(x = sqrt(phi2))
    rho$L <- update(rho$L, Vt, mult = 1)
    return(TRUE)
}

getPar.clmm <- function(rho)
### Extract vector of parameters from model-environment rho
  c(rho$fepar, ST2par(rho$ST))

setPar.clmm <- function(rho, par, which=rep(TRUE, length(par))) {
### Set parameters in model environment rho.
    which <- as.logical(as.vector(which))
    oldpar <- getPar.clmm(rho)
    stopifnot(length(which) == length(oldpar))
    stopifnot(sum(which) == length(par))
    ## over-wright selected elements of oldpar:
    oldpar[which] <- as.vector(par)
    ## assign oldpar to rho$fepar and rho$ST:
    rho$fepar <- oldpar[1:rho$dims$nfepar]
    rho$ST <- par2ST(oldpar[-(1:rho$dims$nfepar)], rho$ST)
}

ST2par <- function(STlist) {
### Compute parameter vector from list of ST matrices.
    unlist(lapply(STlist, function(ST) {
        ## if(ncol(ST) == 1) as.vector(ST) else
        as.vector(c(diag(ST), ST[lower.tri(ST)]))
    }))
}

par2ST <- function(STpar, STlist) {
### Fill in parameters in list of ST matrices. Reverse of ST2par().
    nc <- sapply(STlist, ncol)
    asgn <- rep(1:length(nc), sapply(nc, function(qi) qi * (qi + 1) / 2))
    STparList <- split(STpar, asgn)
    stopifnot(length(asgn) == length(ST2par(STlist)))

    for(i in 1:length(STlist)) {
        par <- STparList[[i]]
        if(nc[i] > 1) {
            diag(STlist[[i]]) <- par[1:nc[i]]
            STlist[[i]][lower.tri(STlist[[i]])] <- par[-(1:nc[i])]
        } else {
            STlist[[i]][] <- par
        }
    }
  STlist
}

STatBoundary <- function(STpar, STlist, tol=1e-3) {
### Compute dummy vector of which ST parameters are at the
### boundary of the parameters space (variance-parameters that are
### zero).
    STcon <- STconstraints(STlist)
    stopifnot(length(STpar) == length(STcon))
    as.integer(STcon == 1 & STpar <= tol)
}

paratBoundary <- function(rho, tol=1e-3)
### Compute dummy vector of which parameters are at the boundary of
### the parameter space.
    c(rep(0, rho$dims$nfepar),
      STatBoundary(ST2par(rho$ST), rho$ST, tol))

paratBoundary2 <- function(rho, tol=1e-3) {
    STcon <- STconstraints(rho$ST)
    c(rep(0L, rho$dims$nfepar),
      as.integer(STcon == 1 & ST2par(rho$ST) < tol))
}

STconstraints <- function(STlist) {
### Compute indicator vector of which variance parameters are constrained above zero. The
### variance parameters are non-negative, while the covariance parameters are not
### constrained.
###
### This function can also be used to generate starting values for the covar. parameters.
    nc <- sapply(STlist, ncol)
    unlist(lapply(nc, function(qi) {
        c(rep(1L, qi), rep(0L, qi * (qi - 1) / 2))
    } ))
}

parConstraints <- function(rho)
### Returns a dummy vector of the same length as getPar.clmm(rho)
### indicating which parameters are contrained to be non-negative.
    c(rep(0, rho$dims$nfepar), STconstraints(rho$ST))

STstart <- function(STlist) par2ST(STconstraints(STlist), STlist)

isNested <- function(f1, f2)
### Borrowed from lme4/R/lmer.R
### Checks if f1 is nested within f2.
{
    f1 <- as.factor(f1)
    f2 <- as.factor(f2)
    stopifnot(length(f1) == length(f2))
    sm <- as(new("ngTMatrix",
                 i = as.integer(f2) - 1L,
                 j = as.integer(f1) - 1L,
                 Dim = c(length(levels(f2)),
                 length(levels(f1)))),
             "CsparseMatrix")
    all(diff(sm@p) < 2)
}

set.AGQ <- function(rho, nAGQ, control, ssr) {
    ## Stop if arguments are incompatible:
    if(nAGQ != 1 && !ssr)
        stop("Quadrature methods are not available with more than one random effects term",
             call.=FALSE)
    if(nAGQ != 1 && control$useMatrix)
        stop("Quadrature methods are not available with 'useMatrix = TRUE'",
             call.=FALSE)
    rho$nAGQ <- nAGQ
    if(nAGQ %in% 0:1) return(invisible())
    ghq <- gauss.hermite(abs(nAGQ))
    rho$ghqns <- ghq$nodes
    rho$ghqws <-
        if(nAGQ > 0) ghq$weights ## AGQ
        else log(ghq$weights) + (ghq$nodes^2)/2 ## GHQ
}

clmm.fit.env <-
  function(rho, control = list(), method=c("nlminb", "ucminf"),
           Hess = FALSE)
### Fit the clmm by optimizing the Laplace likelihood.
### Returns a list with elements:
###
### coefficients
### ST
### logLik
### Niter
### dims
### u
### optRes
### fitted.values
### L
### Zt
### ranef
### condVar
### gradient
### (Hessian)
{
    method <- match.arg(method)
    if(method == "ucminf")
        warning("cannot use ucminf optimizer for this model, using nlminb instead")
    ## Compute lower bounds on the parameter vector
    lwr <- c(-Inf, 0)[parConstraints(rho) + 1]
    ## hack to remove ucminf control settings:
    keep <- !names(control) %in% c("grad", "grtol")
    control <- if(length(keep)) control[keep] else list()
    ## Fit the model with Laplace:
    fit <- try(nlminb(getPar.clmm(rho), function(par) getNLA(rho, par),
                      lower=lwr, control=control), silent=TRUE)
### FIXME: Make it possible to use the ucminf optimizer with
### log-transformed std-par instead.

    ## Check if optimizer converged without error:
    if(inherits(fit, "try-error"))
        stop("optimizer ", method, " failed to converge", call.=FALSE)
### FIXME: Could have an argument c(warn, fail, ignore) to optionally
### return the fitted model despite the optimizer failing.

    ## Ensure parameters in rho are set at the optimum:
    setPar.clmm(rho, fit$par)
    ## Ensure random mode estimation at optimum:
    nllFast.u(rho)
    update.u(rho)

    names(rho$ST) <- names(rho$dims$nlev.re)
    ## Prepare list of results:
    res <- list(coefficients = fit$par[1:rho$dims$nfepar],
                ST = rho$ST,
                logLik = -fit$objective,
                dims = rho$dims,
### FIXME: Should we evaluate hess.u(rho) to make sure rho$L contains
### the right values corresponding to the optimum?
                u = rho$u,
                optRes = fit,
                fitted.values = rho$fitted,
                L = rho$L,
                Zt = rho$Zt
                )
    ## save ranef and condVar in res:
    if(rho$allST1) {
        res$ranef <- rep.int(unlist(rho$ST), rho$dims$nlev.re) * rho$u
        res$condVar <- as.vector(diag(solve(rho$L)) *
                                 rep.int(unlist(rho$ST)^2, rho$dims$nlev.re))
    } else {
        Lambda <- getLambda(rho$ST, rho$dims$nlev.re)
        res$ranef <- Lambda %*% rho$u
        res$condVar <- tcrossprod(Lambda %*% solve(rho$L), Lambda)
    }
    ## Add gradient vector and optionally Hessian matrix:
    bound <- as.logical(paratBoundary2(rho))
    optpar <- fit$par[!bound]
    if(Hess) {
### NOTE: This is the Hessian evaluated for all parameters that are
### not at the boundary at the parameter space. The likelihood for
### models with boundary parameters is still defined as a function of
### all the parameters, so standard errors will differ whether or not
### boundary terms are included or not.
        gH <- deriv12(function(par) getNLA(rho, par, which=!bound),
                      x=optpar)
        res$gradient <- gH$gradient
        res$Hessian <- gH$Hessian
    } else {
        res$gradient <- grad.ctr(function(par) getNLA(rho, par, which=!bound),
                                 x=optpar)
    }
### FIXME: We should check that the (forward) gradient for variances at the
### boundary are not < -1e-5 (wrt. -logLik/nll/getNLA)
    ## Setting Niter and neval after gradient and Hessian evaluations:
    res$Niter <- rho$Niter
    res$neval <- rho$neval
    ## return value:
    res
}

update.u <- function(rho)
{
  stepFactor <- 1
  innerIter <- 0
  rho$u <- rho$uStart
  rho$nll <- nll.u(rho)
  if(!is.finite(rho$nll)) return(FALSE)
  rho$gradient <- grad.u(rho)
  maxGrad <- max(abs(rho$gradient))
  conv <- -1  ## Convergence flag
  message <- "iteration limit reached when updating the random effects"
  if(rho$ctrl$trace > 0)
    Trace(iter=0, stepFactor, rho$nll, maxGrad, rho$u, first=TRUE)
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
    if(!hess.u(rho)) return(FALSE)
    step <- as.vector(solve(rho$L, rho$gradient))
    rho$u <- rho$u - stepFactor * step
    nllTry <- nllFast.u(rho) ## no 'X %*% beta' update
    lineIter <- 0

    ## Step halfing:
    while(nllTry > rho$nll) {
      stepFactor <- stepFactor/2
      rho$u <- rho$u + stepFactor * step
      nllTry <- nllFast.u(rho) ## no 'X %*% beta' update
      lineIter <- lineIter + 1
      if(rho$ctrl$trace > 0)
        Trace(i+innerIter, stepFactor, rho$nll, maxGrad,
              rho$u, first=FALSE)
      if(lineIter > rho$ctrl$maxLineIter){
        message <- "step factor reduced below minimum when updating
the random effects"
        conv <- 1
        break
      }
      innerIter <- innerIter + 1
    }
    rho$nll <- nllTry
    rho$gradient <- grad.u(rho)
    maxGrad <- max(abs(rho$gradient))
    if(rho$ctrl$trace > 0)
      Trace(i+innerIter, stepFactor, rho$nll, maxGrad, rho$u, first=FALSE)
    stepFactor <- min(1, 2 * stepFactor)
  }
  if(conv != 0 && rho$ctrl$innerCtrl == "warnOnly") {
    warning(message, "\n  at iteration ", rho$Niter)
    utils::flush.console()
  }
  else if(conv != 0 && rho$ctrl$innerCtrl == "giveError")
        stop(message, "\n  at iteration ", rho$Niter)
  rho$Niter <- rho$Niter + i - 1
  if(!hess.u(rho)) return(FALSE)
  if(!is.finite(rho$nll))
    return(FALSE)
  else
    return(TRUE)
}

clmm.finalize <-
  function(fit, frames, ths, use.ssr)
{
    fit$tJac <- ths$tJac
    fit$contrasts <- attr(frames$X, "contrasts")
    fit$na.action <- attr(frames$mf, "na.action")
    fit$terms <- frames$terms
### FIXME: Should the terms object contain only the fixed effects
### terms?
    fit$xlevels <- .getXlevels(frames$terms, frames$mf)
    fit$y.levels <- levels(frames$y)
    fit <- within(fit, {
        ## extract coefficients from 'fit':
        names(coefficients) <- names(gradient) <-
            c(ths$alpha.names, colnames(frames$X)[-1])
        alpha <- coefficients[1:dims$nalpha]
        beta <- if(dims$nbeta > 0)
            coefficients[dims$nalpha + 1:dims$nbeta] else numeric(0)
### QUESTION: How does lmerTest get the Hessian of varcov-parameters
### when some are at the boundary?
        ## set various fit elements:
        edf <- dims$edf <- dims$nfepar + dims$nSTpar
        dims$nobs <- sum(frames$wts)
        dims$df.residual <- dims$nobs - dims$edf
        Theta <- alpha %*% t(tJac)
        nm <- paste(y.levels[-length(y.levels)], y.levels[-1], sep="|")
        dimnames(Theta) <- list("", nm)
        rm(nm)

        info <-
            data.frame("link" = link,
                       "threshold" = threshold,
                       "nobs" = dims$nobs,
                       "logLik" = formatC(logLik, digits=2, format="f"),
                       "AIC" = formatC(-2*logLik + 2*dims$edf, digits=2,
                       format="f"),
                       ## "niter" = paste(optRes$info["neval"], "(", Niter, ")",
                       ## sep=""),
                       "niter" = paste(neval, "(", Niter, ")",
                       sep=""),
                       "max.grad" = formatC(max(abs(gradient)), digits=2,
                       format="e")
                       ## BIC is not part of output since it is not clear what
                       ## the no. observations are.
                       )
    })
    bound <- if(use.ssr)  rep(FALSE, fit$dims$edf) else as.logical(paratBoundary2(fit))
    dn <- c(names(fit$coefficients),
            paste("ST", seq_len(fit$dims$nSTpar), sep=""))[!bound]
    names(fit$gradient) <- dn
    if(!is.null(fit$Hessian))
        dimnames(fit$Hessian) <- list(dn, dn)

    ## set class and return fit:
    class(fit) <- "clmm"
    return(fit)
}


