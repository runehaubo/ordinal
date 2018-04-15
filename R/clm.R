## This file contains:
## The main clm function and some auxiliary functions to generate the
## model frame and handle the model environment.

checkArgs.clm <- function(mc) {
    nm <- names(as.list(mc))
    if(!"formula" %in% nm) stop("Model needs a formula", call.=FALSE)
    if("offset" %in% nm)
        stop("offset argument not allowed: ",
             "specify 'offset' in formula or scale arguments instead")
    invisible()
}

clm <-
  function(formula, scale, nominal, data, weights, start, subset,
           doFit = TRUE, na.action, contrasts, model = TRUE,
           control = list(),
           link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
           threshold = c("flexible", "symmetric", "symmetric2", "equidistant"), ...)
{
    mc <- match.call(expand.dots = FALSE)
    link <- match.arg(link)
    threshold <- match.arg(threshold)
    if(missing(contrasts)) contrasts <- NULL
    if(missing(start)) start <- NULL
    checkArgs.clm(mc=match.call())
    ## set control parameters: ## getControl.clm
    control <- do.call(clm.control, c(control, list(...)))

    ## Extract and process formulas:
    call.envir <- parent.frame(n=1)
    formulas <- get_clmFormulas(mc, call.envir)
    ## Get full model.frame and terms.objects:
    fullmf <- get_clm.mf(mc, formulas$fullForm, attr(formulas, "envir"),
                         call.envir)
    if(control$method == "model.frame") return(fullmf)
    terms.list <-
        if(any(c("scale", "nominal") %in% names(formulas)))
            get_clmTerms(mc, formulas, call.envir) else list(formula=terms(fullmf))
    ## Get y, X, weights, off etc.:
    design <- get_clmDesign(fullmf, terms.list, contrasts)
    lst <- namedList(doFit, control, link, threshold, start, formulas)
    if(control$method == "design") return(c(design, lst))
    ## Get clm.struct:
    design <- c(design, makeThresholds(design$y.levels, threshold))
    design <- drop.cols(design, silent=TRUE, drop.scale=FALSE)
    clm.struct <- c(design, lst)
    if(control$method == "struct") return(clm.struct)
    ## Fit model, check convergence, or return a model environment:
    fit <- clm.fit.default(clm.struct)
    if(doFit == FALSE) return(fit)
    ## Format output, prepare result:
    keep <- c("terms", "contrasts", "xlevels", # formula
              "S.terms", "S.contrasts", "S.xlevels", # scale
              "nom.terms", "nom.contrasts", "nom.xlevels", # nominal
              "na.action", "y", "y.levels",
              "control", "link", "threshold", "start", "formulas")
    res <- c(fit, clm.struct[match(keep, names(clm.struct), 0L)],
             list(formula=lst$formulas$formula, call=match.call()))
    ## res$tJac <- format_tJac(res$tJac, res$y.levels, clm.struct$alpha.names)
    res$info=get_clmInfoTab(res)
    if(model) res$model <- fullmf
    res <- res[sort(names(res))]
    class(res) <- "clm"
    res
}

clm.newRho <-
    function(parent=parent.frame(), y, X, NOM=NULL, S=NULL, weights,
             offset, S.offset=NULL, tJac, ...)
### Setting variables in rho: B1, B2, o1, o2, weights.
{
    ## Make B1, B2, o1, o2 based on y, X and tJac:
    keep <- weights > 0
    y[!keep] <- NA
    y <- droplevels(y)
    ntheta <- nlevels(y) - 1
    y <- c(unclass(y))
    y[is.na(y)] <- 0
    n <- sum(keep)
    B2 <- 1 * (col(matrix(0, nrow(X), ntheta + 1)) == y)
    o1 <- c(1e5 * B2[keep, ntheta + 1]) - offset[keep]
    o2 <- c(-1e5 * B2[keep, 1]) - offset[keep]
    B1 <- B2[keep, -(ntheta + 1), drop = FALSE]
    B2 <- B2[keep, -1, drop = FALSE]
    ## adjust B1 and B2 for structured thresholds:
    B1 <- B1 %*% tJac
    B2 <- B2 %*% tJac
    ## update B1 and B2 with nominal effects:
    if(NCOL(NOM) > 1) { ## !is.null(NOM) && ncol(NOM) > 1) {
        ## if !is.null(NOM) and NOM is more than an intercept:
        LL1 <- lapply(1:ncol(NOM), function(x) B1 * NOM[keep, x])
        B1 <- do.call(cbind, LL1)
        LL2 <- lapply(1:ncol(NOM), function(x) B2 * NOM[keep, x])
        B2 <- do.call(cbind, LL2)
    }
    ## update B1 and B2 with location effects (X):
    nbeta <- NCOL(X) - 1
    if(nbeta > 0) {
        B1 <- cbind(B1, -X[keep, -1, drop = FALSE])
        B2 <- cbind(B2, -X[keep, -1, drop = FALSE])
    }
    dimnames(B1) <- NULL
    dimnames(B2) <- NULL
    n.psi <- ncol(B1) ## no. linear model parameters
    ## there may be scale offset without scale predictors:
    sigma <- Soff <-
        if(is.null(S.offset)) rep(1, n) else exp(S.offset[keep])
    ## save scale model matrix:
    k <- 0
    if(!is.null(S)) {
        S <- S[keep, -1, drop=FALSE]
        dimnames(S) <- NULL
        k <- ncol(S) ## no. scale parameters
    }
    has.scale <- ## TRUE if scale has to be considered.
        (!is.null(S) || any(S.offset != 0))
    ## initialize fitted values and weights:
    fitted <- numeric(length = n)
    wts <- weights[keep]
    lst <- namedList(B1, B2, o1, o2, n.psi, S, Soff, k, sigma, has.scale, fitted,
                     wts, clm.nll, clm.grad, clm.hess)
    list2env(x=lst, parent=parent)
}
