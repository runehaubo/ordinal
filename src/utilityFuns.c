#include<R.h>
#include<Rmath.h>
#include "links.h"

double mu = 0, sigma = 1;
int give_log = 0, lower_tail = 1;

//---------------------------------

double d_pfun();
double d_pfun2(); // with lower_tail arg
double d_dfun();
double d_gfun();
double d_gAO();

//--- negative log-likelihood:

double d_nll();

//--- Utilities:

double mmax();
double maxAbs();
void Trace();

//---------------------------------

//------------------------------------------------------------------
// CDFs:

void pgumbel(double *q, int *nq, double *loc, double *scale,
	     int *lower_tail)
{
// pgumbel()
    int i;
// How can I handle if loc and scale are not of unit length?
    for(i = 0; i < *nq; i++)
	q[i] = d_pgumbel(q[i], *loc, *scale, *lower_tail);
}

void pgumbel2(double *q, int *nq, double *loc, double *scale,
	      int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = 1 - d_pgumbel(-q[i], *loc, *scale, *lower_tail);
}

void pAO(double *q, int *nq, double *lambda, int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = d_pAO(q[i], *lambda, *lower_tail);
}

void plgamma(double *q, int *nq, double *lambda, int *lower_tail)
{
    int i;
    for(i = 0; i < *nq; i++)
	q[i] = d_plgamma(q[i], *lambda, *lower_tail);
}

//------------------------------------------------------------------
// PDFs:

void dgumbel(double *x, int *nx, double *loc, double *scale,
	     int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dgumbel(x[i], *loc, *scale, *give_log);
}

void dgumbel2(double *x, int *nx, double *loc, double *scale,
	     int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dgumbel2(x[i], *loc, *scale, *give_log);
}

void dAO(double *x, int *nx, double *lambda, int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dAO(x[i], *lambda, *give_log);
}

void dlgamma(double *x, int *nx, double *lambda, int *give_log)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dlgamma(x[i], *lambda, *give_log);
}

//------------------------------------------------------------------
// gradients of PDFs:


void glogis(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_glogis(x[i]);
}

void gnorm(double *x, int *nx)
{
// Gradient of dnorm(x) wrt. x
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gnorm(x[i]);
}

void gcauchy(double *x, int *n)
{
// Gradient of dcauchy(x) wrt. x
    int i;
    for(i = 0; i < *n; i++)
	x[i] = d_gcauchy(x[i]);
}

void ggumbel(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_ggumbel(x[i]);
}

void ggumbel2(double *x, int *nx)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = -d_ggumbel(-x[i]);
    // or x[i] = d_ggumbel2(x[i]);
}

void gAO(double *x, int *nx, double *lambda)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gAO(x[i], *lambda);
}

void glgamma(double *x, int *nx, double *lambda)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_glgamma(x[i], *lambda);
}

//------------------------------------------------------------------
// link utility functions:

/* Link functions::
  1: logistic
  2: probit
  3: cloglog
  4: loglog
  5: cauchit
  6: Aranda-Ordaz
  7: log-gamma
 */

double d_pfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return plogis(x, mu, sigma, lower_tail, give_log);
    case 2: // probit
	return pnorm(x, mu, sigma, lower_tail, give_log);
    case 3: // cloglog
	return d_pgumbel(x, mu, sigma, lower_tail);
    case 4: // loglog
	return d_pgumbel2(x, mu, sigma, lower_tail);
    case 5: // cauchit
	return pcauchy(x, mu, sigma, lower_tail, give_log);
    case 6: // Aranda-Ordaz
	return d_pAO(x, lambda, lower_tail);
    case 7: // log-gamma
	return d_plgamma(x, lambda, lower_tail);
    default : // all other
	// if(link == 6)
	//     error("the Aranda-Ordaz link is not available");
	// if(link == 7)
	//     error("the log-gamma link is not available");
	// else
	error("link not recognized\n");
	return NA_REAL;
    }
}

double d_pfun2(double x, double lambda, int link, int lower_tail)
// 2nd version of d_pfun with a lower_tail arg
{
    switch(link) {
    case 1: // logistic
	return plogis(x, mu, sigma, lower_tail, give_log);
    case 2: // probit
	return pnorm(x, mu, sigma, lower_tail, give_log);
    case 3: // cloglog
	return d_pgumbel(x, mu, sigma, lower_tail);
    case 4: // loglog
	return d_pgumbel2(x, mu, sigma, lower_tail);
    case 5: // cauchit
	return pcauchy(x, mu, sigma, lower_tail, give_log);
    case 6: // Aranda-Ordaz
	return d_pAO(x, lambda, lower_tail);
    case 7: // log-gamma
	return d_plgamma(x, lambda, lower_tail);
    default : // all other
	// if(link == 6)
	//     error("the Aranda-Ordaz link is not available");
	// if(link == 7)
	//     error("the log-gamma link is not available");
	// else
	error("link not recognized\n");
	return NA_REAL;
    }
}

void pfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_pfun(x[i], *lambda, *link);
}

double d_dfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return dlogis(x, mu, sigma, give_log);
    case 2: // probit
	return dnorm(x, mu, sigma, give_log);
    case 3: // cloglog
	return d_dgumbel(x, mu, sigma, give_log);
    case 4: // loglog
	return d_dgumbel2(x, mu, sigma, give_log);
    case 5: // cauchit
	return dcauchy(x, mu, sigma, give_log);
    case 6:
	return d_dAO(x, lambda, give_log);
    case 7:
	return d_dlgamma(x, lambda, give_log);
    default : // all other
	error("link not recognized\n");
	return NA_REAL;
    }
}

void dfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_dfun(x[i], *lambda, *link);
}

double d_gfun(double x, double lambda, int link)
{
    switch(link) {
    case 1: // logistic
	return d_glogis(x);
    case 2: // probit
	return d_gnorm(x);
    case 3: // cloglog
	return d_ggumbel(x);
    case 4: // loglog
	return d_ggumbel2(x);
    case 5: // cauchit
	return d_gcauchy(x);
    case 6:
	return d_gAO(x, lambda);
    case 7:
	return d_glgamma(x, lambda);
    default : // all other
	error("link not recognized\n");
	return NA_REAL;
    }
}

void gfun(double *x, int *nx, double *lambda, int *link)
{
    int i;
    for(i = 0; i < *nx; i++)
	x[i] = d_gfun(x[i], *lambda, *link);
}

//------------------------------------------------------------------

void getFitted(double *eta1, double *eta2, int *neta) 
{
  // adjust for NA and NaN values?
  int i;
  for(i = 0; i < *neta; i++) {
    if(eta2[i] <= 0) 
      // pnorm(x, mu, sigma, lower_tail, give_log);
      eta1[i] = pnorm(eta1[i], 0.0, 1.0, 1, 0) -  
	pnorm(eta2[i], 0.0, 1.0, 1, 0);
    else 
      eta1[i] = pnorm(eta2[i], 0.0, 1.0, 0, 0) - 
	pnorm(eta1[i], 0.0, 1.0, 0, 0);
  }
}

void getFitted2(double *eta1, double *eta2, int *neta, double *lambda,
		int *link) 
// 2nd version now including a link arg
{
  // adjust for NA and NaN values?
  int i;
  for(i = 0; i < *neta; i++) {
    if(eta2[i] <= 0) 
      // d_pfun2(x, lambda, link, lower_tail)
      eta1[i] = d_pfun2(eta1[i], *lambda, *link, 1) - 
	d_pfun2(eta2[i], *lambda, *link, 1);
    else 
      eta1[i] = d_pfun2(eta2[i], *lambda, *link, 0) - 
	d_pfun2(eta1[i], *lambda, *link, 0);
  }
}


//------------------------------------------------------------------
// Gradients and Hessians for update.b in clmm2():

void grFacSum(double *x, int *grFac, int *nx, double *u, int *nu)
// compute tapply(x, grFac, sum) + u
{
    int i, j;
    double z = 0;

    for(i = 0; i < *nu; i++) {
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                z = z + x[j];
	}
	u[i] = u[i] + z;
        z = 0;
    }
}

// FIXME: grFacSum such that it can be used by gradC and hessC - this
// should simplify the code

double d_nll(double *u, int nu, int *grFac, double stDev,
	     double *o1, double *o2, int no, double *eta1,
	     double *eta2, double *eta1Fix, double *eta2Fix,
	     double *sigma, double *pr, double *weights,
	     double lambda, int *link)
/*
  Returns:          nll
  Updates:          eta1, eta2, pr  given the new value of u
  Leaves unchanged: u, grFac, stDev, o1, o2, eta1Fix, eta2Fix, sigma, weights
*/
{
  int i, j;
  double o, nll = 0.0;
  
  for(i = 0; i < no; i++) {
    o = u[grFac[i] - 1] * stDev;
    eta1[i] = (eta1Fix[i] + o1[i] - o) / sigma[i];
    eta2[i] = (eta2Fix[i] + o2[i] - o) / sigma[i];
    /* Accurate evaluation of pr (fitted values) even if eta1 and eta2
       are both large: */
    if(eta2[i] <= 0)
      pr[i] = d_pfun2(eta1[i], lambda, *link, 1) -
	d_pfun2(eta2[i], lambda, *link, 1);
    else
      pr[i] = d_pfun2(eta2[i], lambda, *link, 0) -
	d_pfun2(eta1[i], lambda, *link, 0);
    if(!R_FINITE(pr[i]) || pr[i] <= 0.) {
      return INFINITY;
    }
    nll -= weights[i] * log(pr[i]);
  }
  for(j = 0; j < nu; j++)
    nll -= dnorm(u[j], 0., 1., 1);
  return nll;
}

void nll(double *u, int *nu, int *grFac, double *stDev,
	 double *o1, double *o2, int *no, double *eta1,
	 double *eta2, double *eta1Fix, double *eta2Fix,
	 double *sigma, double *pr, double *weights,
	 double *lambda, int *link, double *nll)
{
    *nll = d_nll(u, *nu, grFac, *stDev, o1, o2, *no, eta1, eta2,
		 eta1Fix, eta2Fix, sigma, pr, weights, *lambda, link);
}

void grad(double *stDev, double *p1, double *p2, double *pr,
	  double *weights, double *sigma, double *wtprSig,
	  double *eta1, double *eta2, double *gradValues,
	  double *u, int *grFac, int *nx, int *ngv,
	  double *lambda, int *link)
/*
  Returns:          void
  Updates:          gradValues, p1, p2, wtprSig  given the new values of eta1, eta2
  Leaves unchanged: grFac, stDev, eta1, eta2, pr, sigma, weights, link, nx, ngv
  Assumes:
  nx: length of grFac, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2
  ngv: length of gradValues
 */
{
    int i, j;
    // double tmp[*nx], z = 0;

    // update p1, p2, wtprSig:
    for(i = 0; i < *nx; i++) {
	p1[i] = d_dfun(eta1[i], *lambda, *link);
	p2[i] = d_dfun(eta2[i], *lambda, *link);
	wtprSig[i] = weights[i] / pr[i] / sigma[i];
    }

    // sum for each level of the grouping factor:
    for(i = 0; i < *ngv; i++) {
	gradValues[i] = 0; // Could set these to
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
		gradValues[i] += *stDev * wtprSig[j] *
		    (p1[j] - p2[j]);
	}
	gradValues[i] += u[i];
    }
}

void gradC(double *stDev, double *p1, double *p2, double *wtprSig,
	   int *grFac, int *nx, double *u, int *nu)
{
    // gradient for update.b
    int i, j;
    double z = 0;

    for(i = 0; i < *nx; i++) {
	wtprSig[i] = *stDev * wtprSig[i] * (p1[i] - p2[i]);
    }

    for(i = 0; i < *nu; i++) {
        for (j = 0; j < *nx; j++) {
            if(grFac[j] == i + 1)
                z += wtprSig[j];
	}
	u[i] += z;
        z = 0;
    }
}

void hess(double *stDev, double *p1, double *p2, double *pr,
	  double *wtprSig, double *eta1, double *eta2, int *link,
	  int *grFac, int *nx, double *hessValues, double *lambda,
	  int *nhv)
/*
  Returns:          void
  Updates:          hessValues  given the new values of eta1, eta2
  Leaves unchanged: grFac, stDev, eta1, eta2, p1, p2, pr, sigma, weights, link, nx, ngv
  Assumes:
  nx: length of grFac, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2
  nhv: length of hessValues
 */
{
  int i, j;
  
  // sum for each level of the grouping factor:
  for(i = 0; i < *nhv; i++) {
    hessValues[i] = 0;
    for (j = 0; j < *nx; j++) {
      if(grFac[j] == i + 1)
	hessValues[i] +=
	  (R_pow_di(p1[j] - p2[j], 2) / pr[j] -
	   (d_gfun(eta1[j], *lambda, *link) -
	    d_gfun(eta2[j], *lambda, *link))) * wtprSig[j];
    }
    hessValues[i] = (hessValues[i] * *stDev * *stDev) + 1;
  }
}

void hessC(double *stDev, double *p1, double *p2, double *pr,
	   double *g1, double *g2, double *wtprSig,
	   int *grFac, int *nx, double *z, int *nz)
{
  // hessian for update.b
  int i, j;
  double sigma2;
  
  sigma2 = R_pow_di(*stDev, 2);
  
  for(i = 0; i < *nx; i++)
    pr[i] = (R_pow_di(p1[i] - p2[i], 2) / pr[i] -
	     (g1[i] - g2[i])) * wtprSig[i];
  
  for(i = 0; i < *nz; i++) {
    for (j = 0; j < *nx; j++) {
      if(grFac[j] == i + 1)
	z[i] = z[i] + pr[j];
    }
    z[i] = z[i] * sigma2 + 1;
  }
}

//------------------------------------------------------------------
// Trace function:

void Trace(int iter, double stepFactor, double val, double maxGrad,
	   double *par, int npar, int first)
{
    int i;

    if(first)
	Rprintf("iter:  step factor:     Value:     max|grad|:   Parameters:\n");
    Rprintf(" %3d:    %1.3e:   %.3f:     %1.3e:  ", iter, stepFactor, val, maxGrad);
    for(i = 0; i < npar; i++)
	Rprintf(" %.4f", par[i]);
    Rprintf("\n");
}

//------------------------------------------------------------------

void NRalg(int *trace, int *maxIter, double *gradTol,
	   int *maxLineIter, int *grFac,
	   double *stDev, double *o1, double *o2,
	   double *eta1Fix, double *eta2Fix, double *eta1,
	   double *eta2, double *sigma, int *link,
	   double *weights, double *u,
	   double *pr, double *funValue,
	   double *gradValues, double *hessValues,
	   int *nx, int *nu, double *maxGrad, int *conv,
	   double *p1, double *p2, double *wtprSig,
	   double *lambda, int *Niter)
{
/*
  nx: length(pr)
  r:  length(start) = length(u)

  updates: u, funValue, gradValues, hessValues, maxGrad,

  correct vector input:
  eta1, eta2, pr, funValue (grad is called before d_nll), u = 0,
  grFac, o1, o2, eta1Fix, eta2Fix, sigma, weights

  arbitrary input:
  p1, p2, wtprSig, gradValues, hessValues,

  needed output:
  u, funValue, gradValues, hessValues, conv, Niter,
*/
    int lineIter, innerIter = 0, i, j;
    double stepFactor = 1, funValueTry, step[*nu];

    *funValue = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1, eta2,
		      eta1Fix, eta2Fix, sigma, pr, weights,
		      *lambda, link);
    if(!R_FINITE(*funValue)) {
	*conv = 0;
	return ;
    }
    grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	 gradValues, u, grFac, nx, nu, lambda, link);
    *maxGrad = maxAbs(gradValues, *nu);
    *conv = -1; // Convergence flag
    if(*trace)
	Trace(0, stepFactor, *funValue, *maxGrad, u, *nu, 1);

    // Newton-Raphson algorithm:
    for(i = 0; i < *maxIter; i++) {
        if(*maxGrad < *gradTol) {
            *conv = 1;
            return ;
	}
	hess(stDev, p1, p2, pr, wtprSig, eta1, eta2, link,
	     grFac, nx, hessValues, lambda, nu);
	for(j = 0; j < *nu; j++) {
	    step[j] = gradValues[j] / hessValues[j];
	    u[j] -= stepFactor * step[j];
	}
	funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
			    eta2, eta1Fix, eta2Fix, sigma, pr,
			    weights, *lambda, link);
	lineIter = 0;
	//  simple line search, i.e. step halfing:
	while(funValueTry > *funValue) {
	    stepFactor *= 0.5;
	    for(j = 0; j < *nu; j++)
		u[j] += stepFactor * step[j];
	    funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
				eta2, eta1Fix, eta2Fix, sigma, pr,
				weights, *lambda, link);
	    lineIter++;
	    if(*trace)
		Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad,
		      u, *nu, 0);
	    if(lineIter > *maxLineIter){
		*conv = -2;
		return ;
	    }
	    innerIter++;
        }
        *funValue = funValueTry;
	grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	     gradValues, u, grFac, nx, nu, lambda, link);
	*maxGrad = maxAbs(gradValues, *nu);
	if(*trace)
	    Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad, u,
		  *nu, 0);
	stepFactor = fmin2(1., stepFactor * 2.);
	(*Niter)++;
    }
}

void NRalgv3(int *trace, int *maxIter, double *gradTol,
	     int *maxLineIter, int *grFac,
	     double *stDev, double *o1, double *o2,
	     double *eta1Fix, double *eta2Fix, double *sigma,
	     int *link, double *weights, double *u,
	     double *pr, double *funValue,
	     double *gradValues, double *hessValues,
	     int *nx, int *nu, double *maxGrad, int *conv,
	     double *lambda, int *Niter)
// Less input and slightly faster than NRalg().
{
/*
  control arguments from clmm - see ?clmm.control:
  trace, maxIter, gradTol, maxLineIter all of length 1

  length = nx: grFac, o1, o2, eta1Fix, eta2Fix, sigma, weights
  length = 1: stDev, funValue, nx, nu, maxGrad, conv, lambda, Niter
  length = nu: gradValues, hessValues, u

  updates: u, funValue, gradValues, hessValues, maxGrad, conv, Niter,
  pr,

  correct vector input:
  eta1, eta2, pr, u = 0, grFac, o1, o2, eta1Fix, eta2Fix, sigma,
  weights

  arbitrary input:
  gradValues, hessValues,

  needed output:
  u, funValue, gradValues, hessValues, conv, Niter,
*/
    int lineIter, innerIter = 0, i, j;
    double stepFactor = 1, funValueTry, step[*nu];
    double eta1[*nx], eta2[*nx], p1[*nx], p2[*nx], wtprSig[*nx];

    *funValue = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1, eta2,
		      eta1Fix, eta2Fix, sigma, pr, weights,
		      *lambda, link);
    if(!R_FINITE(*funValue)) {
	*conv = 0;
	return ;
    }
    grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	 gradValues, u, grFac, nx, nu, lambda, link);
    *maxGrad = maxAbs(gradValues, *nu);
    *conv = -1; // Convergence flag
    if(*trace)
	Trace(0, stepFactor, *funValue, *maxGrad, u, *nu, 1);

    // Newton-Raphson algorithm:
    for(i = 0; i < *maxIter; i++) {
        if(*maxGrad < *gradTol) {
            *conv = 1;
            return ;
	}
	hess(stDev, p1, p2, pr, wtprSig, eta1, eta2, link, grFac, nx,
	     hessValues, lambda, nu);
	for(j = 0; j < *nu; j++) {
	    /* Actually there is no need to store 'step' since
	       'gradValues' could hold the step values (maintained
	       here for code clarity) */
	    step[j] = gradValues[j] / hessValues[j];
	    u[j] -= stepFactor * step[j];
	}
	funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
			    eta2, eta1Fix, eta2Fix, sigma, pr,
			    weights, *lambda, link);
	lineIter = 0;
	//  simple line search, i.e. step halfing:
	while(funValueTry > *funValue) {
	    stepFactor *= 0.5;
	    for(j = 0; j < *nu; j++)
		u[j] += stepFactor * step[j];
	    funValueTry = d_nll(u, *nu, grFac, *stDev, o1, o2, *nx, eta1,
				eta2, eta1Fix, eta2Fix, sigma, pr,
				weights, *lambda, link);
	    lineIter++;
	    if(*trace)
		Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad,
		      u, *nu, 0);
	    if(lineIter > *maxLineIter){
		*conv = -2;
		return ;
	    }
	    innerIter++;
        }
        *funValue = funValueTry;
	grad(stDev, p1, p2, pr, weights, sigma, wtprSig, eta1, eta2,
	     gradValues, u, grFac, nx, nu, lambda, link);
	*maxGrad = maxAbs(gradValues, *nu);
	if(*trace)
	    Trace(i+1+innerIter, stepFactor, *funValue, *maxGrad, u,
		  *nu, 0);
	stepFactor = fmin2(1.0, stepFactor * 2.0);
	(*Niter)++;
    }
    (*Niter)--;
}

//------------------------------------------------------------------

void getNGHQ(double *nll, int *grFac, double *stDev,
	     double *eta1Fix, double *eta2Fix, double *o1, double *o2,
	     double *Sigma, double *weights, int *nx, int *nu,
	     double *ghqns, /* double *ghqws,*/ double *lghqws,
	     int *nGHQ, int *link, double *ns, double *lambda)
{
  int i, j, h;
  double SS = 0, SS1 = 0, SS2 = 0, eta1tmp, eta2tmp, pr_tmp;
  
  for(i = 0; i < *nu; i++) {
    for(h = 0; h < *nGHQ; h++) {
      for(j = 0; j < *nx; j++) {
	if(grFac[j] == i + 1) {
	  eta1tmp = (eta1Fix[j] + o1[j] - ns[h]) / Sigma[j];
	  eta2tmp = (eta2Fix[j] + o2[j] - ns[h]) / Sigma[j];
	  /* Accurate evaluation of differences of probabilities even
	     if eta1tmp and eta2tmp are large: */
	  if(eta2tmp <= 0)
	    pr_tmp = d_pfun2(eta1tmp, *lambda, *link, 1) -
	      d_pfun2(eta2tmp, *lambda, *link, 1);
	  else
	    pr_tmp = d_pfun2(eta2tmp, *lambda, *link, 0) -
	      d_pfun2(eta1tmp, *lambda, *link, 0);
	  // sum up contributions:
	  SS1 += weights[j] * log(pr_tmp);
	}
      }
      // SS2 += exp(SS1) * ghqws[h];
      // SS2 += exp(SS1 + log(ghqws[h]));
      SS2 += exp(SS1 + lghqws[h]);
      SS1 = 0;
    }
    SS += log(SS2);
    SS2 = 0;
  }
  *nll = -SS + *nu * log(M_PI * 2) * 0.5;
}

void getNAGQ(double *nll, int *grFac, double *stDev,
	     double *eta1Fix, double *eta2Fix, double *o1, double *o2,
	     double *Sigma, double *weights, int *nx, int *nu,
	     double *ghqns, double *lghqws, /* double *lghqws, */
	     double *ghqns2, double *u, double *D,
	     int *nAGQ, int *link, double *lambda)
/*
  nll: negative log-likelihood (return value)

  length = nx: grFac, o1, o2, eta1Fix, eta2Fix, Sigma, weights
  length = 1: stDev, nll, nx, nu, nAGQ, lambda, link
  length = nu: D, u
  length = nAGQ: ghqns, lghqws (log ghqws) / ghqws
 */
{
  int i, j, h;
  double SS1 = 0, SS2 = 0, eta1tmp, eta2tmp, K, ranNew, pr_tmp;
  *nll = 0;
  
  for(i = 0; i < *nu; i++) {
    K = sqrt(2. / D[i]);
    for(h = 0; h < *nAGQ; h++) {
      for(j = 0; j < *nx; j++) {
	if(grFac[j] == i + 1) {
	  ranNew = *stDev * (u[i] + K * ghqns[h]);
	  eta1tmp = (eta1Fix[j] + o1[j] - ranNew) / Sigma[j];
	  eta2tmp = (eta2Fix[j] + o2[j] - ranNew) / Sigma[j];
	  /* Accurate evaluation of differences of probabilities even
	     if eta1tmp and eta2tmp are large: */
	  if(eta2tmp <= 0)
	    pr_tmp = d_pfun2(eta1tmp, *lambda, *link, 1) -
	      d_pfun2(eta2tmp, *lambda, *link, 1);
	  else
	    pr_tmp = d_pfun2(eta2tmp, *lambda, *link, 0) -
	      d_pfun2(eta1tmp, *lambda, *link, 0);
	  // sum up contributions:
	  SS1 += weights[j] * log(pr_tmp);
	}
      }
      // SS2 += exp(SS1) * K * ghqws[h] *
      // 	dnorm(u[i] + K * ghqns[h], mu, sigma, give_log);
      //  	    SS2 += exp(SS1 + lghqws[h] + ghqns2[h] - //R_pow_di(ghqns[h], 2) +
      //  		       0.5 * R_pow_di(u[i] + K * ghqns[h], 2)) * K;
      SS2 += exp(SS1 + lghqws[h] + ghqns2[h] - //R_pow_di(ghqns[h], 2) +
		 0.5 * R_pow_di(u[i] + K * ghqns[h], 2));
      SS1 = 0;
    }
    // *nll -= log(SS2);
    *nll -= log(SS2) + log(K);
    SS2 = 0;
  }
  *nll += *nu * log(M_PI * 2) * 0.5;
}


//------------------------------------------------------------------

double mmax(double *x, int nx)
/*
   Return the maximum of the elements in x
   nx: length of x ( >= 1)
 */
{
    int i;
    double cmax; // current max

    cmax = x[0];
    if(nx == 1)
	return cmax;
    for(i = 1; i < nx; i++) {
	if(x[i] > cmax)
	    cmax = x[i];
    }
    return cmax;
}

double maxAbs(double *x, int nx)
/*
  Return max(abs(x))
  nx: length of x ( >= 1 )
 */
{
    int i;
    double cmax; // current max

    cmax = fabs(x[0]);
    if(nx == 1)
	return cmax;
    for(i = 1; i < nx; i++) {
	if(fabs(x[i]) > cmax)
	    cmax = fabs(x[i]);
    }
    return cmax;
}
