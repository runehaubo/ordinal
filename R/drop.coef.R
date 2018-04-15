## This file contains:
## Functions that can drop columns from rank-deficient design
## matrices. One is exported and others used internally.

drop.coef <- function(X, silent = FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0
{
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
  if(qr.X$rank == ncol(X))
    return(X) ## return X if X has full column rank
  if(!silent) ## message the no. dropped columns:
    message(gettextf("design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  ## did we succeed? stop-if-not:
  if(qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}

drop.coef2 <- function(X, tol = 1e-7, silent = FALSE, test.ans = FALSE)
### works if ncol(X) >= 0 and nrow(X) >= 0
{
  ## test and match arguments:
  stopifnot(is.matrix(X))
  silent <- as.logical(silent)[1]
  aliased <- rep.int(0, ncol(X))
  ## perform the qr-decomposition of X using LINPACK methods:
  qr.X <- qr(X, tol = tol, LAPACK = FALSE)
  if(qr.X$rank == ncol(X)) {
    ## return X if X has full column rank
    attr(X, "aliased") <- aliased
    attr(X, "orig.colnames") <- colnames(X)
    return(X)
  }
  if(!silent) ## message the no. dropped columns:
    message(gettextf("design is column rank deficient so dropping %d coef",
                     ncol(X) - qr.X$rank))
  ## return the columns correponding to the first qr.x$rank pivot
  ## elements of X:
  newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
  sel <- qr.X$pivot[-(1:qr.X$rank)]
  aliased[sel] <- 1
  attr(newX, "aliased") <- aliased
  attr(newX, "orig.colnames") <- colnames(X)
  ## Copy old attributes:
  attributes(newX)$contrasts <- attributes(X)$contrasts
  attr(newX, "assign") <- attr(X, "assign")[-sel]
  ## did we succeed? stop-if-not:
  if(test.ans && qr.X$rank != qr(newX)$rank)
    stop(gettextf("determination of full column rank design matrix failed"),
         call. = FALSE)
  return(newX)
}


drop.cols <- function(mf, silent = FALSE, drop.scale=TRUE)
### drop columns from X and possibly NOM and S to ensure full column
### rank.
### mf - list with X and possibly NOM and S design matrices. Includes
### alpha.names
###
### returns: updated version of mf.
{
    nalpha <- length(mf$alpha.names)
    ## X is assumed to contain an intercept at this point:
    Xint <- match("(Intercept)", colnames(mf$X), nomatch = 0)
    if(Xint <= 0) {
        mf$X <- cbind("(Intercept)" = rep(1, nrow(mf$X)), mf$X)
        warning("an intercept is needed and assumed")
    } ## intercept in X is guaranteed.
    if(!is.null(mf[["NOM"]])){
        ## store coef names:
        mf$coef.names <- list()
        mf$coef.names$alpha <-
            paste(rep(mf$alpha.names, ncol(mf$NOM)), ".",
                  rep(colnames(mf$NOM), each=nalpha), sep="")
        mf$coef.names$beta <- colnames(mf$X)[-1]
        ## drop columns from NOM:
        mf$NOM <- drop.coef2(mf$NOM, silent=silent)
        ## drop columns from X:
        NOMX <- drop.coef2(cbind(mf$NOM, mf$X[,-1, drop=FALSE]),
                           silent=silent)
        ## extract and store X:
        mf$X <- cbind("(Intercept)" = rep(1, nrow(mf$X)),
                      NOMX[,-seq_len(ncol(mf$NOM)), drop=FALSE])
        ## store alias information:
        mf$aliased <- list(alpha = rep(attr(mf$NOM, "aliased"),
                           each=nalpha))
        mf$aliased$beta <- attr(NOMX, "aliased")[-seq_len(ncol(mf$NOM))]
        if(drop.scale && !is.null(mf[["S"]])) {
            mf$coef.names$zeta <- colnames(mf$S)[-1]
            ## drop columns from S:
            NOMS <- drop.coef2(cbind(mf$NOM, mf$S[,-1, drop=FALSE]),
                               silent=silent)
            ## extract and store S:
            mf$S <- cbind("(Intercept)" = rep(1, nrow(mf$S)),
                          NOMS[,-seq_len(ncol(mf$NOM)), drop=FALSE])
            mf$aliased$zeta <- attr(NOMS,
                                    "aliased")[-seq_len(ncol(mf$NOM))]
        } else if(!is.null(mf[["S"]])) {
            Sint <- match("(Intercept)", colnames(mf$S), nomatch = 0)
            if(Sint <= 0) {
                mf$S <- cbind("(Intercept)" = rep(1, nrow(mf$S)), mf$S)
                warning("an intercept is needed and assumed in 'scale'",
                        call.=FALSE)
            } ## intercept in S is guaranteed.
            mf$coef.names$zeta <- colnames(mf$S)[-1]
            mf$S <- drop.coef2(mf$S, silent=silent)
            mf$aliased$zeta <- attr(mf$S, "aliased")[-1]
        }
        return(mf)
    } ## end !is.null(mf[["NOM"]])
    ## drop columns from X assuming an intercept:
    mf$coef.names <- list(alpha = mf$alpha.names,
                          beta = colnames(mf$X)[-1])
    mf$X <- drop.coef2(mf$X, silent=silent)
    mf$aliased <- list(alpha = rep(0, nalpha),
                       beta = attr(mf$X, "aliased")[-1])
    ## drop columns from S if relevant:
    if(!is.null(mf[["S"]])) {
        Sint <- match("(Intercept)", colnames(mf$S), nomatch = 0)
        if(Sint <= 0) {
            mf$S <- cbind("(Intercept)" = rep(1, nrow(mf$S)), mf$S)
            warning("an intercept is needed and assumed in 'scale'",
                    call.=FALSE)
        } ## intercept in S is guaranteed.
        mf$coef.names$zeta <- colnames(mf$S)[-1]
        mf$S <- drop.coef2(mf$S, silent=silent)
        mf$aliased$zeta <- attr(mf$S, "aliased")[-1]
    }
    return(mf)
}
