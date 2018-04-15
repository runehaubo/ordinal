## This file contains:
## Methods to compute and plot likelihood-slices for clm objects.

slice <- function(object, ...) {
  UseMethod("slice")
}

slice.clm <-
  function(object, parm = seq_along(par), lambda = 3, grid = 1e2,
           quad.approx = TRUE, ...)
{
  ## argument matching and testing:
  stopifnot(is.numeric(lambda) && lambda > 0)
  stopifnot(is.numeric(grid) && grid >= 1)
  grid <- as.integer(round(grid))
  par <- coef(object, na.rm=TRUE)
  par.names <- names(par)
  npar <- length(par)
  stopifnot(length(parm) == length(unique(parm)))
  if(is.character(parm))
    parm <- match(parm, par.names, nomatch = 0)
### disallow character argument due to ambiguity?
  if(!all(parm %in% seq_along(par)))
    stop("invalid 'parm' argument")
  stopifnot(length(parm) > 0)
  parm <- as.integer(round(parm))
  ## parm is an integer vector indexing non-aliased coef.
  ml <- object$logLik
  parm.names <- par.names[parm]

  ## get environment corresponding to object:
  rho <- get_clmRho(object)
  ## rho <- update(object, doFit = FALSE)
  names(par) <- NULL
  rho$par <- par ## set rho$par to mle
  stopifnot(isTRUE(all.equal(rho$clm.nll(rho), -object$logLik)))

  ## generate sequence of parameters at which to compute the
  ## log-likelihood:
  curv <- sqrt(1/diag(object$Hessian)) ## curvature in nll wrt. par
  par.range <- par + curv %o% c(-lambda, lambda)
  ## par.seq - list of length npar with a sequence of values for each
  ## parameter :
  par.seq <- lapply(parm, function(ind) {
    seq(par.range[ind, 1], par.range[ind, 2], length = grid) })
  ## compute relative logLik for all par.seq for each par:
  logLik <- lapply(seq_along(parm), function(i) { # for each par
    rho$par <- par ## reset par values to MLE
    sapply(par.seq[[ i ]], function(par.val) { # for each par.seq value
      rho$par[ parm[i] ] <- par.val
      -rho$clm.nll(rho) - ml ## relative logLik
    })
  })

  ## collect parameter sequences and relative logLik in a list of
  ## data.frames:
  res <- lapply(seq_along(parm), function(i) {
    structure(data.frame(par.seq[[ i ]], logLik[[ i ]]),
              names = c(parm.names[i], "logLik"))
  })

  ## set attributes:
  names(res) <- parm.names
  attr(res, "original.fit") <- object
  attr(res, "mle") <- par[parm]
  class(res) <- "slice.clm"

  if(!quad.approx) return(res)
  ## compute quadratic approx to *positive* logLik:
  Quad <- function(par, mle, curv)
    -((mle - par)^2 / curv^2 / 2)
  for(i in seq_along(parm))
    res[[ i ]]$quad <-
      Quad(par.seq[[ i ]], par[ parm[i] ], curv[ parm[i] ])

  return(res)
}

plot.slice.clm <-
  function(x, parm = seq_along(x), type = c("quadratic", "linear"),
           plot.mle = TRUE,
           ask = prod(par("mfcol")) < length(parm) && dev.interactive(),
           ...)
{
  ## Initiala argument matching and testing:
  type <- match.arg(type)
  stopifnot(is.numeric(parm))
  parm <- as.integer(round(parm))
  of <- attr(x, "original.fit")
  par <- coef(of)
  ml <- of$logLik

  ## take the signed sqrt of nll and quad:
  if(type == "linear") {
    sgn.sqrt <- function(par, mle, logLik)
      (2 * (par > mle) - 1) * sqrt(-logLik)
    mle <- coef(attr(x, "original.fit"))
    for(i in parm) {
      x[[i]]$logLik <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$logLik)
      if(!is.null(x[[i]]$quad))
        x[[i]]$quad <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$quad)
    }
    ylab <- "Signed log-likelihood root"
  }
  else
    ylab <- "Relative log-likelihood"

  if(ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }

  ## actual plotting:
  for(i in parm) {
    z <- x[[i]]
    plot(z[1:2], type = "l", ylab=ylab, ...)
    if(!is.null(z$quad))
      lines(z[[1]], z[[3]], lty = 2)
    if(plot.mle && type == "quadratic")
      ## abline(v = par[i])
      abline(v = attr(x, "mle")[i])
    ## abline(v = par[names(x)[i]])
  }

  return(invisible())
}

## slice.clm <-
##   function(object, parm = seq_along(par), lambda = 3, grid = 1e2,
##            quad.approx = TRUE, ...)
## {
##   ## argument matching and testing:
##   stopifnot(is.numeric(lambda) && lambda > 0)
##   stopifnot(is.numeric(grid) && grid >= 1)
##   grid <- as.integer(grid)
##   par <- coef(object)
##   par.names <- names(par)
##   npar <- length(par)
##   stopifnot(length(parm) == length(unique(parm)))
##   if(is.character(parm))
##     parm <- match(parm, par.names, nomatch = 0)
##   if(!all(parm %in% seq_along(par)))
##     stop("invalid 'parm' argument")
##   stopifnot(length(parm) > 0)
##   parm <- as.integer(parm)
##   ml <- object$logLik
##   parm.names <- par.names[parm]
##
##   ## get environment corresponding to object:
##   rho <- update(object, doFit = FALSE)
##   names(par) <- NULL
##   rho$par <- par ## set rho$par to mle
##   stopifnot(isTRUE(all.equal(rho$clm.nll(rho), -object$logLik)))
##
##   ## generate sequence of parameters at which to compute the
##   ## log-likelihood:
##   curv <- sqrt(1/diag(object$Hess)) ## curvature in nll wrt. par
##   par.range <- par + curv %o% c(-lambda, lambda)
##   ## par.seq - list of length npar:
##   par.seq <- sapply(parm, function(ind) {
##     seq(par.range[ind, 1], par.range[ind, 2], length = grid) },
##                     simplify = FALSE)
##   ## compute relative logLik for all par.seq for each par:
##   logLik <- lapply(seq_along(parm), function(i) { # for each par
##     rho$par <- par ## reset par values to MLE
##     sapply(par.seq[[ i ]], function(par.val) { # for each val
##       rho$par[ parm[i] ] <- par.val
##       -rho$clm.nll(rho) - ml ## relative logLik
##     })
##   })
##
##   ## collect results in a list of data.frames:
##   res <- lapply(seq_along(parm), function(i) {
##     structure(data.frame(par.seq[[ i ]], logLik[[ i ]]),
##               names = c(parm.names[i], "logLik"))
##   })
##
##   ## set attributes:
##   names(res) <- parm.names
##   attr(res, "original.fit") <- object
##   class(res) <- "slice.clm"
##
##   if(!quad.approx) return(res)
##   ## compute quadratic approx to *positive* logLik:
##   Quad <- function(par, mle, curv)
##     -((mle - par)^2 / curv^2 / 2)
##   for(i in seq_along(parm))
##     res[[ i ]]$quad <-
##       Quad(par.seq[[ i ]], par[ parm[i] ], curv[ parm[i] ])
##
##   return(res)
## }


