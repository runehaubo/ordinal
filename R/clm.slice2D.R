slice2D <- function(object, ...) {
  UseMethod("slice2D")
}

slice2D.clm <-
    function(object, parm=seq_along(par), lambda=3, grid=20,  ...)
{
    ## argument matching and testing:
    stopifnot(is.numeric(lambda) && lambda > 0)
    stopifnot(is.numeric(grid) && grid >= 1)
    grid <- as.integer(round(grid))
    par <- coef(object, na.rm=TRUE)
    par.names <- names(par)
    stopifnot(length(parm) == length(unique(parm)))
    if(is.character(parm))
        parm <- match(parm, par.names, nomatch = 0)
    if(!all(parm %in% seq_along(par)))
        stop("invalid 'parm' argument")
    stopifnot(length(parm) >= 2L)
    parm <- as.integer(parm)
    nparm <- length(parm)
    ## parm is an integer vector indexing non-aliased coef.
    ml <- object$logLik
    parm.names <- par.names[parm]
    mle <- par[parm]

    ## get environment corresponding to object:
    env <- get_clmRho(object)
    ## env <- update(object, doFit=FALSE)
    names(par) <- NULL
    env$par <- as.vector(par) ## set env$par to mle
    stopifnot(isTRUE(all.equal(env$clm.nll(env), -object$logLik)))

    ## generate sequence of parameters at which to compute the
    ## log-likelihood:
    curv <- sqrt(1/diag(object$Hessian)) ## curvature in nll wrt. par
    par.range <- par + curv %o% (c(-1, 1) * lambda)
    ## All pairwise combinations:
    pairs <- t(combn(seq_len(nparm), 2))
    ncombn <- nrow(pairs)
### Allow for sequential paired comparisons?
    par.seq <- lapply(parm, function(ind) {
        seq(par.range[ind, 1], par.range[ind, 2], length = grid) })
    names(par.seq) <- par.names
    zlist <- vector(mode="list", length=ncombn)
    names(zlist) <- paste(par.names[pairs[, 1]],
                          par.names[pairs[, 2]], sep=".")
    for(k in 1:ncombn) {
        i <- pairs[k, 1]
        j <- pairs[k, 2]
        xx <- expand.grid(x=par.seq[[i]], y=par.seq[[j]])
        ## Set parameter values to MLEs:
        env$par <- par
        ## Compute log-likelihood over entire grid:
        z <- apply(xx, 1, function(x) {
            env$par[c(i, j)] <- as.vector(x);
            env$clm.nll(env) })
        ## Store log-likelihood values in a matrix:
        zlist[[k]] <- matrix(z, ncol=grid)
    }
    res <-
        list(zlist=zlist, par.seq=par.seq, par.range=par.range,
             pairs=pairs,
             original.fit=object, mle=mle)
    class(res) <- c("slice2D.clm")
    res
}

safe.as.int <- function(x) as.integer(round(x))

plot.slice2D.clm <-
    function(x, parm = seq_along(orig.par), ## How to specify default values
             ## of parm?
             plot.mle = TRUE,
             ask = prod(par("mfcol")) < nrow(pairs) && dev.interactive(),
             ...)
### parm: a character vector or integer vector of length >= 2 with
### those par-combinations to make contour plots.
{
    ## stopifnot(all parm in names(par.seq))
    orig.par <- coef(x$original.fit, na.rm=TRUE)
### More parm stuff here...
    stopifnot(is.numeric(parm) && length(parm) >= 2L)
    parm <- as.integer(round(parm))
    par.names <- names(orig.par)
    ## of <- attr(x, "original.fit")
    ## par <- coef(of)
    ## ml <- of$logLik
    keep <- (x$pairs[, 1] %in% parm) & (x$pairs[, 2] %in% parm)
    pairs <- x$pairs[keep, , drop=FALSE]
    stopifnot(length(pairs) >= 2)

    if(ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }

    ## Plotting the contours:
    for(k in seq_len(nrow(pairs))) {
        i <- pairs[k, 1]
        j <- pairs[k, 2]
        contour(x$par.seq[[i]], x$par.seq[[j]], x$zlist[[k]],
                xlab = par.names[i], ylab = par.names[j])
        points(orig.par[i], orig.par[j], pch = 4, col = "red", lwd = 2)
    }
    return(invisible())
}


sliceg.clm <-
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
      rho$clm.nll(rho)
      rho$clm.grad(rho)[ parm[i] ]
    })
  })

  ## collect parameter sequences and relative logLik in a list of
  ## data.frames:
  res <- lapply(seq_along(parm), function(i) {
      structure(data.frame(par.seq[[ i ]], logLik[[ i ]]),
                ## names = c(parm.names[i], "logLik"))
                names = c(parm.names[i], "gradient"))
  })

  ## set attributes:
  names(res) <- parm.names
  attr(res, "original.fit") <- object
  attr(res, "mle") <- par[parm]
  ## class(res) <- "slice.clm"
  class(res) <- "sliceg.clm"

  ## if(!quad.approx) return(res)
  ## ## compute quadratic approx to *positive* logLik:
  ## Quad <- function(par, mle, curv)
  ##   -((mle - par)^2 / curv^2 / 2)
  ## for(i in seq_along(parm))
  ##   res[[ i ]]$quad <-
  ##     Quad(par.seq[[ i ]], par[ parm[i] ], curv[ parm[i] ])

  return(res)
}

plot.sliceg.clm <-
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
  ## if(type == "linear") {
  ##   sgn.sqrt <- function(par, mle, logLik)
  ##     (2 * (par > mle) - 1) * sqrt(-logLik)
  ##   mle <- coef(attr(x, "original.fit"))
  ##   for(i in parm) {
  ##     x[[i]]$logLik <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$logLik)
  ##     if(!is.null(x[[i]]$quad))
  ##       x[[i]]$quad <- sgn.sqrt(x[[i]][1], mle[i], x[[i]]$quad)
  ##   }
  ##   ylab <- "Signed log-likelihood root"
  ## }
  ## else
  ## ylab <- "Relative log-likelihood"
  ylab <- "Gradient"

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
