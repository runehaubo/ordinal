#############################################################################
#    Copyright (c) 2010-2018 Rune Haubo Bojesen Christensen
#
#    This file is part of the ordinal package for R (*ordinal*)
#
#    *ordinal* is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#
#    *ordinal* is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    A copy of the GNU General Public License is available at
#    <https://www.r-project.org/Licenses/> and/or
#    <http://www.gnu.org/licenses/>.
#############################################################################
# clm.anova.R

single_anova <- function(object,
                         type = c("III", "II", "I", "3", "2", "1", "marginal", "2b")) {
  type <- type[1L]
  if(!is.character(type)) type <- as.character(type)
  type <- match.arg(type)
  if(type %in% c("I", "II", "III"))
    type <- as.character(as.integer(as.roman(type)))
  if(any(is.na(vcov(object))))
    stop("anova table not available with non-finite values in vcov(object)")
  # Get list of contrast matrices (L) - one for each model term:
  L_list <- if(type == "1") {
    get_contrasts_type1(object)
  } else if(type == "2") {
    get_contrasts_type2_unfolded(object)
  } else if(type == "2b") {
    get_contrasts_type2(object)
  } else if(type == "3") {
    get_contrasts_type3(object)
  } else if(type == "marginal") {
    get_contrasts_marginal(object)
  } else {
    stop("'type' not recognized")
  }
  # Add cols to L for alpha, zeta and lambda params:
  L_list <- adjust_contrast_for_param(object, L_list)
  # Get F-test for each term and collect in table:
  table <- rbindall(lapply(L_list, function(L) contestMD(object, L)))
  # Format ANOVA table and return:
  if(length(nm <- setdiff(names(L_list), rownames(table)))) {
    tab <- array(NA_real_, dim=c(length(nm), ncol(table)),
                 dimnames = list(nm, colnames(table)))
    table <- rbind(table, tab)[names(L_list), ]
  }
  # Format 'type':
  type <- if(type == "marginal") {
    "Marginal"
  } else if(grepl("b|c", type)) {
    alph <- gsub("[0-9]", "", type)
    paste0("Type ", as.roman(as.integer(gsub("b|c", "", type))), alph)
  } else paste("Type", as.roman(as.integer(type)))
  attr(table, "heading") <-
    paste(type, "Analysis of Deviance Table with Wald chi-square tests\n")
  attr(table, "hypotheses") <- L_list
  class(table) <- c("anova", "data.frame")
  table
}

adjust_contrast_for_param <- function(model, L) {
  nalpha <- length(model$alpha)
  nzeta <- if(is.null(model$zeta)) 0L else length(model$zeta)
  nlambda <- if(is.null(model$lambda)) 0L else length(model$lambda)
  nextra <- nzeta + nlambda
  # pre and post add extra cols to L:
  add <- function(L) {
    pre <- array(0, dim=c(nrow(L), nalpha))
    post <- array(0, dim=c(nrow(L), nextra))
    cbind(pre, L[, -1L, drop=FALSE], post)
  }
  if(!is.list(L)) add(L) else lapply(L, add)
}

model_matrix <- function(object, ...) {
  if(!inherits(object, "clm")) return(model.matrix(object, ...))
  X <- model.matrix(object)$X
  if(!any(object$aliased$beta)) return(X)
  remove <- c(FALSE, object$aliased$beta)
  newX <- X[, !remove, drop=FALSE]
  attr(newX, "assign") <- attr(X, "assign")[!remove]
  contr <- attr(X, "contrasts")
  if(!is.null(contr)) attr(newX, "contrasts") <- contr
  newX
}

contestMD <- function(model, L, rhs=0, eps=sqrt(.Machine$double.eps), ...) {
  mk_Qtable <- function(Qvalue, df) {
    pvalue <- pchisq(q=Qvalue, df=df, lower.tail=FALSE)
    data.frame("Df"=df, "Chisq"=Qvalue, "Pr(>Chisq)"=pvalue, 
               check.names = FALSE)
  }
  if(!is.matrix(L)) L <- matrix(L, ncol=length(L))
  stopifnot(is.matrix(L), is.numeric(L),
            ncol(L) == length(coef(model, na.rm=TRUE)))
  if(length(rhs) == 1L) rhs <- rep(rhs, nrow(L))
  stopifnot(is.numeric(rhs), length(rhs) == nrow(L))
  if(nrow(L) == 0L) { # May happen if there are no fixed effects
    x <- numeric(0L)
    return(mk_Qtable(x, x))
  }
  if(any(is.na(L))) return(mk_Qtable(NA_real_, NA_real_))
  beta <- coef(model, na.rm=TRUE)
  vcov_beta <- vcov(model)
  # Adjust beta for rhs:
  if(!all(rhs == 0)) beta <- beta - drop(MASS::ginv(L) %*% rhs)
  # Compute Var(L beta) and eigen-decompose:
  VLbeta <- L %*% vcov_beta %*% t(L) # Var(contrast) = Var(Lbeta)
  eig_VLbeta <- eigen(VLbeta)
  P <- eig_VLbeta$vectors
  d <- eig_VLbeta$values
  tol <- max(eps * d[1], 0)
  pos <- d > tol
  q <- sum(pos) # rank(VLbeta)
  if(q < nrow(L) && !all(rhs == 0))
    warning("Contrast is rank deficient and test may be affected")
  if(q <= 0) { # shouldn't happen if L is a proper contrast
    x <- numeric(0L)
    return(mk_Qtable(x, x))
  }
  PtL <- crossprod(P, L)[1:q, ]
  # Compute t-squared values and Q-value:
  t2 <- drop(PtL %*% beta)^2 / d[1:q]
  Qvalue <- sum(t2) 
  mk_Qtable(Qvalue, df=q)
}

##############################################
######## get_contrasts_type3
##############################################
get_contrasts_type3 <- function(model, which=NULL) {
  term_names <- attr(terms(model), "term.labels")
  # Extract original design matrix:
  Xorig <- model_matrix(model)
  # Assumes Xorig is full (column) rank
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(Xorig) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))
  
  # Extract contrast coding in Xorig:
  codings <- unlist(attr(Xorig, "contrast"))
  # If only treatment contrasts are used we can just return the type 3
  # contrasts for contr.treatment coding:
  if(length(codings) > 0 &&
     all(is.character(codings)) && all(codings %in% c("contr.treatment")))
    return(extract_contrasts_type3(model, X=Xorig))
  # otherwise we need to map the type III contrasts to whatever contrast
  # coding was used:
  X <- get_model_matrix(model, type="remake", contrasts="contr.treatment")
  # Ensure that X is full (column) rank:
  X <- ensure_full_rank(X, silent=TRUE, test.ans=FALSE)
  # Extract contrasts assuming contr.treatment coding:
  type3ctr <- extract_contrasts_type3(model, X=X)
  map <- zapsmall(ginv(X) %*% Xorig) # Maps between contrast codings
  rownames(map) <- colnames(X)
  lapply(type3ctr[which], function(L) L %*% map)
}


##############################################
######## get_contrasts_type1
##############################################
get_contrasts_type1 <- function(model) {
  terms <- terms(model)
  X <- model_matrix(model)
  nalpha <- length(model$alpha)
  p <- ncol(X)
  if(p == 0L) return(list(matrix(numeric(0L), nrow=0L))) # no fixef
  if(p == 1L && attr(terms, "intercept")) # intercept-only model
    return(list(matrix(numeric(0L), ncol=nalpha)))
  # Compute 'normalized' doolittle factorization of XtX:
  L <- if(p == 1L) matrix(1L) else t(doolittle(crossprod(X))$L)
  dimnames(L) <- list(colnames(X), colnames(X))
  # Determine which rows of L belong to which term:
  ind.list <- term2colX(terms, X)[attr(terms, "term.labels")]
  lapply(ind.list, function(rows) L[rows, , drop=FALSE])
}

##############################################
######## get_contrasts_type2_unfolded
##############################################
get_contrasts_type2_unfolded <- function(model, which=NULL) {
  # Computes the 'genuine type II contrast' for all terms that are
  # contained in other terms. For all terms which are not contained in other
  # terms, the simple marginal contrast is computed.
  X <- model_matrix(model)
  Terms <- terms(model)
  term_names <- attr(Terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))
  
  is_contained <- containment(model)
  do_marginal <- names(is_contained)[sapply(is_contained, length) == 0L]
  do_type2 <- setdiff(term_names, do_marginal)
  
  if(!length(do_marginal)) list() else
    Llist <- get_contrasts_marginal(model, which=do_marginal)
  if(length(do_type2))
    Llist <- c(Llist, get_contrasts_type2(model, which=do_type2))
  Llist[term_names]
}

##############################################
######## get_contrasts_type2
##############################################
get_contrasts_type2 <- function(model, which=NULL) {
  # Computes the type 2 contrasts - either for all terms or for those
  # included in 'which' (a chr vector naming model terms).
  # returns a list
  X <- model_matrix(model)
  nalpha <- length(model$alpha)
  terms <- terms(model)
  data_classes <- attr(terms(model), "dataClasses")
  if(is.null(asgn <- attr(X, "assign")))
    stop("design matrix 'X' should have a non-null 'assign' attribute")
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))
  which <- setNames(as.list(which), which)
  
  # Compute containment:
  is_contained <- containment(model)
  
  # Compute term asignment list: map from terms to columns in X
  has_intercept <- attr(terms, "intercept") > 0
  col_terms <- if(has_intercept) c("(Intercept)", term_names)[asgn + 1] else
    term_names[asgn[asgn > 0]]
  if(!length(col_terms) == ncol(X)) # should never happen.
    stop("An error happended when computing Type II contrasts")
  term2colX <- split(seq_along(col_terms), col_terms)[unique(col_terms)]
  
  # Compute contrast for each term - return as named list:
  lapply(which, function(term) {
    # Reorder the cols in X to [, unrelated_to_term, term, contained_in_term]
    cols_term <- unlist(term2colX[c(term, is_contained[[term]])])
    Xnew <- cbind(X[, -cols_term, drop=FALSE], X[, cols_term, drop=FALSE])
    # Compute order of terms in Xnew:
    newXcol_terms <- c(col_terms[-cols_term], col_terms[cols_term])
    # Compute Type I contrasts for the reordered X:
    Lc <- t(doolittle(crossprod(Xnew))$L)
    dimnames(Lc) <- list(colnames(Xnew), colnames(Xnew))
    # Extract rows for term and get original order of columns:
    Lc[newXcol_terms == term, colnames(X), drop=FALSE]
  })
}

##############################################
######## get_contrasts_marginal
##############################################
#' @importFrom stats model.matrix terms
get_contrasts_marginal <- function(model, which=NULL) {
  # Computes marginal contrasts.
  #
  # No tests of conformity with coefficients are implemented
  #
  # returns a list
  X <- model_matrix(model)
  terms <- terms(model)
  term_names <- attr(terms, "term.labels")
  if(is.null(which)) {
    which <- term_names
    # If model has at most one term return Type I contrasts:
    if(ncol(X) <= 1L || length(term_names) <= 1L)
      return(get_contrasts_type1(model))
  } else stopifnot(is.character(which), all(which %in% term_names))

  # Compute map from terms to columns in X and contrasts matrix
  term2colX <- term2colX(terms, X)
  L <- structure(diag(ncol(X)), dimnames = list(colnames(X), colnames(X)))

  # Extract contrast for each term - return as named list:
  which <- setNames(as.list(which), which)
  lapply(which, function(term) {
    L[term2colX[[term]], , drop=FALSE]
  })
}

##############################################
######## rbindall
##############################################
rbindall <- function(...) do.call(rbind, ...)

cbindall <- function(...) do.call(cbind, ...)

