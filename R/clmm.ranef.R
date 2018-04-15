## This file contains:
## Implementation of ranef and condVar methods for clmm objects to
## extract the conditional model of the random-effects and their
## conditional variances.

## fixef.clmm <- function(object, ...)  coef(object, ...)
## object$coefficients
### FIXME: This function needs to change such that a *named* vector is
### returned.

ranef <- function(object, ...) UseMethod("ranef")
## fixef <- function(object, ...) UseMethod("fixef")

ranef.clmm <- function(object, condVar=FALSE, ...)
### This function...
### args...
### Returns....
{
    formatRanef <- function(relist, ST, gf.levels, assign, qi) {
        asgn <- split(seq_along(assign), assign)
        ## colnames of random effects:
        cn <- lapply(ST, colnames)
        cn <- lapply(asgn, function(ii) unlist(cn[ii]))
        ranefList <- lapply(seq_along(relist), function(i) {
            matrix(relist[[i]], ncol=qi[i])
        })
        ## Combine r.e. terms associated with the same grouping factors,
        ## set dimnames and coerce to data.frame:
        ranefList <- lapply(seq_along(asgn), function(i) {
            mat <- do.call(cbind, ranefList[ asgn[[i]] ])
            dimnames(mat) <- list(gf.levels[[i]], cn[[i]])
            as.data.frame(mat)
        })
        ## list of r.e. by grouping factors:
        names(ranefList) <- names(gflevs)
        ranefList
    }
    ## which r.e. terms are associated with which grouping factors:
    asgn <- attributes(object$gfList)$assign
    ## names of levels of grouping factors:
    gflevs <- lapply(object$gfList, levels)
    ## random effects indicator factor:
    reind <- with(object$dims, factor(rep.int(seq_len(nretrms),
                                              nlev.re * qi)))
    ## list of random effects by r.e. term:
    relist <- split(object$ranef, reind)
    ranefList <- formatRanef(relist, object$ST, gflevs, asgn,
                             object$dims$qi)
    if(condVar) {
### FIXME: Should we return matrices for vector-valued random effects
### as lmer does?
        ## Add conditional variances of the random effects:
        cond.var <- object$condVar
        if(NCOL(cond.var) > 1) cond.var <- diag(cond.var)
        cvlist <- split(cond.var, reind)
        cond.var <- formatRanef(cvlist, object$ST, gflevs, asgn,
                                object$dims$qi)
        for(i in seq_along(ranefList))
            attr(ranefList[[i]], "condVar") <- cond.var[[i]]
    }
    ranefList
}

condVar <- function(object, ...) UseMethod("condVar")
condVar.clmm <- function(object, ...)
    lapply(ranef.clmm(object, condVar=TRUE),
           function(y) attr(y, "condVar"))
