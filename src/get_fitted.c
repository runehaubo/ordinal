#include<R.h>
#include<Rmath.h>
#include <Rinternals.h>
#include "links.h"

SEXP get_fitted(SEXP, SEXP, SEXP, SEXP);

// -------------------------------------------------------


SEXP get_fitted(SEXP eta1p, SEXP eta2p, SEXP linkp, SEXP lambdap) {
  /* Compute fitted values (probabilities) from vectors of linear
     predictors (eta1 and eta2) given the link function (linkp) and an
     optional lambda parameter.

     eta1 and eta2 are required to be equal length  numeric vectors,
     linkp a character vector and lambdap a numeric scalar.

     return: vector of fittec values of same length as eta1 and eta2. 
    */
  SEXP ans = PROTECT(duplicate(coerceVector(eta1p, REALSXP)));
  eta2p = coerceVector(eta2p, REALSXP);
  linkp = coerceVector(linkp, STRSXP);
  const char *linkc = CHAR(asChar(linkp));
  double *eta1 = REAL(ans), *eta2 = REAL(eta2p), 
    lambda = asReal(lambdap);
  int i, nans = LENGTH(ans);

  if(LENGTH(eta2p) != nans) {
    UNPROTECT(1);
    error("'eta1' and 'eta2' should have the same length");
  }
  
  if(strcmp(linkc, "probit") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// pnorm(x, mu, sigma, lower_tail, give_log);
	eta1[i] = pnorm(eta1[i], 0.0, 1.0, 1, 0) -  
	  pnorm(eta2[i], 0.0, 1.0, 1, 0);
      else 
	eta1[i] = pnorm(eta2[i], 0.0, 1.0, 0, 0) - 
	  pnorm(eta1[i], 0.0, 1.0, 0, 0);
    }
  } 
  else if(strcmp(linkc, "logit") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// plogis(x, mu, sigma, lower_tail, give_log);
	eta1[i] = plogis(eta1[i], 0.0, 1.0, 1, 0) -  
	  plogis(eta2[i], 0.0, 1.0, 1, 0);
      else 
	eta1[i] = plogis(eta2[i], 0.0, 1.0, 0, 0) - 
	  plogis(eta1[i], 0.0, 1.0, 0, 0);
    }
  } 
  else if(strcmp(linkc, "loglog") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// d_pgumbel(double q, double loc, double scale, int lower_tail)
	eta1[i] = d_pgumbel(eta1[i], 0., 1., 1) -  
	  d_pgumbel(eta2[i], 0., 1., 1);
      else 
	eta1[i] = d_pgumbel(eta2[i], 0., 1., 0) - 
	  d_pgumbel(eta1[i], 0., 1., 0);
    }
  } 
  else if(strcmp(linkc, "cloglog") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// d_pgumbel2(double q, double loc, double scale, int lower_tail)
	eta1[i] = d_pgumbel2(eta1[i], 0., 1., 1) -  
	  d_pgumbel2(eta2[i], 0., 1., 1);
      else 
	eta1[i] = d_pgumbel2(eta2[i], 0., 1., 0) - 
	  d_pgumbel2(eta1[i], 0., 1., 0);
    }
  } 
  else if(strcmp(linkc, "cauchit") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// pcauchy(q, loc, scale, lower_tail, give_log)
	eta1[i] = pcauchy(eta1[i], 0., 1., 1, 0) -  
	  pcauchy(eta2[i], 0., 1., 1, 0);
      else 
	eta1[i] = pcauchy(eta2[i], 0., 1., 0, 0) - 
	  pcauchy(eta1[i], 0., 1., 0, 0);
    }
  } 
  else if(strcmp(linkc, "Aranda-Ordaz") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// d_pAO(q, lambda, lower_tail)
	eta1[i] = d_pAO(eta1[i], lambda, 1) -  
	  d_pAO(eta2[i], lambda, 1);
      else 
	eta1[i] = d_pAO(eta2[i], lambda, 0) - 
	  d_pAO(eta1[i], lambda, 0);
    }
  } 
  else if(strcmp(linkc, "log-gamma") == 0) {
    for(i = 0; i < nans; i++) {
      if(eta2[i] <= 0) 
	// d_plgamma(double eta, double lambda, int lower_tail)
	eta1[i] = d_plgamma(eta1[i], lambda, 1) - 
	  d_plgamma(eta2[i], lambda, 1);
      else 
	eta1[i] = d_plgamma(eta2[i], lambda, 0) - 
	  d_plgamma(eta1[i], lambda, 0);
    }
  } 
  else {
    UNPROTECT(1); // unprotecting before exiting with an error
    error("link not recognized");
  }
  UNPROTECT(1);
  return ans;
}
