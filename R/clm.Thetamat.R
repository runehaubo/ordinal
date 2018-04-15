## This file contains:
## Functions (getThetamat) to compute a table of threshold
## coefficients from model fits (clm()s) with nominal effects.

getThetamat <-
  function(terms, alpha, assign, contrasts, tJac, xlevels)
### Compute matrix of thresholds for all combinations of levels of
### factors in the nominal formula.
###
### Input:
### terms: nominal terms object
### alpha: vector of threshold parameters
### assign: attr(NOM, "assign"), where NOM is the design matrix for
###   the nominal effects
### contrasts: list of contrasts for the nominal effects
### tJac: threshold Jacobian with appropriate dimnames.
### xlevels: names of levels of factors among the nominal effects.
###
### Output:
### Theta: data.frame of thresholds
### mf.basic: if nrow(Theta) > 1 a data.frame with factors in columns
###   and all combinations of the factor levels in rows.
{
    ## Make matrix of thresholds; Theta:
    Theta <- matrix(alpha, ncol=ncol(tJac), byrow=TRUE)
    ## Matrix with variables-by-terms:
    factor.table <- attr(terms, "factors")
    all.varnm <- rownames(factor.table)
### NOTE: need to index with all.varnm not to include (weights) and
### possibly other stuff.
    var.classes <- attr(terms, "dataClasses")[all.varnm]
    numeric.var <- which(var.classes != "factor")
### NOTE: Logical variables are treated as numeric variables.
    numeric.terms <- factor.terms <- numeric(0)
    if(length(factor.table)) {
        ## Terms associated with numerical variables:
        numeric.terms <-
            which(colSums(factor.table[numeric.var, , drop=FALSE]) > 0)
        ## Terms only involving factor variables:
        factor.terms <-
            which(colSums(factor.table[numeric.var, , drop=FALSE]) == 0)
    }
    ## Remove rows in Theta for numeric variables:
    if(length(numeric.terms)) {
### NOTE: ncol(NOM) == length(asgn) == nrow(Theta)
### length(attr(terms, "term.labels")) == ncol(factor.table)
### NOTE: length(var.classes) == nrow(factor.table)
        numeric.rows <- which(assign %in% numeric.terms)
        Theta <- Theta[-numeric.rows, , drop=FALSE]
        ## Drop terms so the design matrix, X for the factors does not
        ## include numeric variables:
        if(length(factor.terms))
            terms <- drop.terms(terms, dropx=numeric.terms,
                                keep.response=FALSE)
    }
    ## if some nominal effects are factors:
    if(length(factor.terms)) {
        ## get xlevels for factors, not ordered (factors)
        factor.var <- which(var.classes == "factor")
        factor.varnm <- names(var.classes)[factor.var]
        xlev <- xlevels[factor.varnm]
        ## minimal complete model frame:
        mf.basic <- do.call(expand.grid, xlev)
        ## minimal complete design matrix:
        X <- model.matrix(terms, data=mf.basic,
                          contrasts=contrasts[factor.varnm])
### NOTE: get_clmDesign adds an intercept if its not there, so we need
### to do that as well here. Otherwise 'X[, keep, drop=FALSE]' will
### fail:
        if(!"(Intercept)" %in% colnames(X))
            X <- cbind("(Intercept)" = rep(1, nrow(X)), X)
### NOTE: There are no contrasts for numerical variables, but there
### may be for ordered factors.
        ## From threshold parameters to thresholds:
### NOTE: some rows of Theta may contain NAs due to rank deficiency of
### the NOM design matrix.
        keep <- apply(Theta, 1, function(x) sum(is.na(x)) == 0)
        ## Theta <- apply(Theta, 2, function(th) X %*% th)
        tmp <- lapply(1:NCOL(Theta), function(i) {
            X[, keep, drop=FALSE] %*% Theta[keep, i]
        })
        Theta <- do.call(cbind, tmp)
    }
    ## Adjust each row in Theta for threshold functions:
    tmp <- lapply(seq_len(nrow(Theta)), function(i)
                  c(tJac %*% Theta[i, ]))
    Theta <- do.call(rbind, tmp)
### NOTE: apply returns a vector and not a matrix when ncol(Theta) ==
### 1, so we need to avoid it here.
    ## Theta <- t(apply(Theta, 1, function(th) tJac %*% th))
    colnames(Theta) <- rownames(tJac)
    res <- list(Theta = as.data.frame(Theta))
    ## add factor information if any:
    if(NROW(Theta) > 1)  res$mf.basic <- mf.basic
    ## return:
    res
}
