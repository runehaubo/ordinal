#include "links.h"

/* This file implements scalar distribution, density and gradient
   function */

/*-------------------------------------------------------*/
/* Scalar cumulative distribution functions (CDFs) */
/*-------------------------------------------------------*/

double d_pgumbel(double q, double loc, double scale, int lower_tail)
// Consider implementing 'int give_log' to follow the convention from
// pnorm etc. 
{
  if(ISNAN(q)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1.;
  else if(q == R_NegInf)
    q = 0.;
  else {
    q = (q - loc) / scale;
    q = exp( -exp( -q));
  }
  return !lower_tail ? 1 - q : q;
}

double d_pgumbel2(double q, double loc, double scale, int lower_tail)
// this is (partly) redundant since d_pgumbel2(q) = 1 - d_pgumbel(-q)
{
  if(ISNAN(q)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1;
  else if(q == R_NegInf)
    q = 0;
  else {
    q = (-q - loc) / scale;
    q = exp(-exp(-q));
  }
  return !lower_tail ? q : 1 - q;
}

double d_pAO(double q, double lambda, int lower_tail)
{
  if(ISNAN(q) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(q == R_PosInf) 
    q = 1;
  else if(q == R_NegInf)
    q = 0;
  else {
    if(lambda < 1.0e-6)
      error("'lambda' has to be positive. lambda = %e was supplied\n",
	    lambda);
    q = 1 - R_pow(lambda * exp(q) + 1, -1/lambda);
  }
  return !lower_tail ? 1 - q : q;
}

double d_plgamma(double eta, double lambda, int lower_tail)
{
  double v; 
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf) 
    v = 1;
  else if(eta == R_NegInf)
    v = 0;
  else {
    v = R_pow_di(lambda, -2) * exp(lambda * eta);
    if(lambda < 1.0e-6)
      v = 1 - pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		     lower_tail, 0 /*give_log*/);
    if(lambda > -1.0e-6)
      v = pgamma(v, R_pow_di(lambda, -2), /*scale = */ 1,
		 lower_tail, 0/*give_log*/);
    if(lambda >= -1.0e-6 && lambda <= 1.0e-6)
      // pnorm(x, mu, sigma, lower_tail, give_log);
      v = pnorm(eta, 0., 1., 1, 0);
  }
  return !lower_tail ? 1 - v : v;
}

/*-------------------------------------------------------*/
/* Scalar probability density functions (PDFs) */
/*-------------------------------------------------------*/

double d_dgumbel(double x, double loc, double scale, int give_log)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    // if(x == INFINITE || x == -INFINITE) // seems to work as well.
    return 0; // this special case needs to be handled separately 
  x = (x - loc) / scale;
  x = -exp(-x) - x - log(scale);
  return give_log ? x : exp(x);
}

double d_dgumbel2(double x, double loc, double scale, int give_log)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0;
  x = (-x - loc) / scale;
  x = -exp(-x) - x - log(scale);
  return give_log ? x : exp(x);
}

double d_dAO(double eta, double lambda, int give_log)
{
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf || eta == R_NegInf)
    return 0;
  if(lambda < 1.0e-6)
    error("'lambda' has to be positive. lambda = %e was supplied\n",
	  lambda);
  eta -= (1 + 1 / lambda) * log(lambda * exp(eta) + 1);
  return give_log ? eta : exp(eta);
}

double d_dlgamma(double x, double lambda, int give_log)
{
  if(ISNAN(x) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 
  if(lambda < 1.0e-5 && lambda > -1.0e-5) // lambda close to zero
    return dnorm(x, 0. , 1., give_log); 

  double q_2 = R_pow_di(lambda, -2);
  x *= lambda;
  x = log(fabs(lambda)) + q_2 * log(q_2) -
    lgammafn(q_2) + q_2 * (x - exp(x));
  return !give_log ? exp(x) : x;
}

/*-------------------------------------------------------*/
/* Scalar gradients of probability density functions */
/*-------------------------------------------------------*/

double d_glogis(double x)
{
  // Gradient of dlogis(x) wrt. x
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    // if(x == INFINITE || x == -INFINITE) // seems to work as well.
    return 0; // this special case needs to be handled separately 

  /* Store the sign of x, compute the gradient for the absolute value
     and restore the sign. This is needed to avoid exp(LARGE) to blow
     up and the function to return NaN.
  */
  int sign = x > 0; //could use fsign() instead...
  x = exp(-fabs(x));
  x = 2 * x * x * R_pow_di(1 + x, -3) - x *
    R_pow_di(1 + x, -2);
  return sign ? x : -x;
}

double d_gnorm(double x) 
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == INFINITY || x == -INFINITY)
    return 0;
  else
    return -x * dnorm(x, 0., 1., 0);
}

double d_gcauchy(double x)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 
  
  return x = -2 * x / M_PI * R_pow_di(1 + x * x, -2);
}

double d_ggumbel(double x)
{
  if(ISNAN(x)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0; 

  x = exp(-x);
  if(x == INFINITY)
    return 0;
  
  double eq = exp(-x);
  return -eq * x + eq * x * x;
}

double d_ggumbel2(double x)
// redundant function...
{
    return -d_ggumbel(-x);
}

double d_gAO(double eta, double lambda)
{
  if(ISNAN(eta) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(eta == R_PosInf || eta == R_NegInf)
    return 0; 

  double lex = lambda * exp(eta);
  if(lex == R_PosInf || lex == 0)
    return 0.;
  double y = d_dAO(eta, lambda, 0/*give_log*/);
  
  return y == 0. ? 0. : y * (1 - (1 + 1/lambda) * lex / (1 + lex)); 
}

double d_glgamma(double x, double lambda)
{
  if(ISNAN(x) || ISNAN(lambda)) // true for NA and NaN
    return NA_REAL;
  if(x == R_PosInf || x == R_NegInf)
    return 0.; 
  if(lambda < 1.0e-5 && lambda > -1.0e-5) // lambda close to zero
    return -x * dnorm(x, 0., 1., 0/*give_log*/);
  
  double z = exp(lambda * x);
  if(z == R_PosInf || z == 0.)
    return 0.;
  double y = d_dlgamma(x, lambda, 0/*give_log*/);
  return y <= 0. ? 0.0 : y * (1 - exp(lambda * x)) / lambda; 
  // Equivalent to:
  /* if(y <= 0) 
     return 0.0;
     else 
     return y * (1 - exp(lambda * x)) / lambda;
  */
}

