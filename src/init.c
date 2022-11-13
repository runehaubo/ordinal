/////////////////////////////////////////////////////////////////////////////
//    Copyright (c) 2010-2022 Rune Haubo Bojesen Christensen
//
//    This file is part of the ordinal package for R (*ordinal*)
//
//    *ordinal* is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 2 of the License, or
//    (at your option) any later version.
//
//    *ordinal* is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    A copy of the GNU General Public License is available at
//    <https://www.r-project.org/Licenses/> and/or
//    <http://www.gnu.org/licenses/>.
/////////////////////////////////////////////////////////////////////////////
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void dAO_C(void *, void *, void *, void *);
extern void dgumbel_C(void *, void *, void *, void *, void *);
extern void dgumbel2_C(void *, void *, void *, void *, void *);
extern void dlgamma_C(void *, void *, void *, void *);
extern void gAO_C(void *, void *, void *);
extern void gcauchy_C(void *, void *);
extern void getNAGQ(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void getNGHQ_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void ggumbel_C(void *, void *);
extern void ggumbel2_C(void *, void *);
extern void glgamma_C(void *, void *, void *);
extern void glogis_C(void *, void *);
extern void gnorm_C(void *, void *);
extern void grad_C(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gradC(void *, void *, void *, void *, void *, void *, void *, void *);
extern void grFacSum_C(void *, void *, void *, void *, void *);
extern void hess(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hessC(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void nll(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void NRalg(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void NRalgv3(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void pAO_C(void *, void *, void *, void *);
extern void pgumbel_C(void *, void *, void *, void *, void *);
extern void pgumbel2_C(void *, void *, void *, void *, void *);
extern void plgamma_C(void *, void *, void *, void *);

/* .Call calls */
extern SEXP get_fitted(SEXP, SEXP, SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
  {"dAO_C",      (DL_FUNC) &dAO_C,       4},
  {"dgumbel_C",  (DL_FUNC) &dgumbel_C,   5},
  {"dgumbel2_C", (DL_FUNC) &dgumbel2_C,  5},
  {"dlgamma_C",  (DL_FUNC) &dlgamma_C,   4},
  {"gAO_C",      (DL_FUNC) &gAO_C,       3},
  {"gcauchy_C",  (DL_FUNC) &gcauchy_C,   2},
  {"getNAGQ",    (DL_FUNC) &getNAGQ,    19},
  {"getNGHQ_C",  (DL_FUNC) &getNGHQ_C,  17},
  {"ggumbel_C",  (DL_FUNC) &ggumbel_C,   2},
  {"ggumbel2_C", (DL_FUNC) &ggumbel2_C,  2},
  {"glgamma_C",  (DL_FUNC) &glgamma_C,   3},
  {"glogis_C",   (DL_FUNC) &glogis_C,    2},
  {"gnorm_C",    (DL_FUNC) &gnorm_C,     2},
  {"grad_C",     (DL_FUNC) &grad_C,     16},
  {"gradC",      (DL_FUNC) &gradC,       8},
  {"grFacSum_C", (DL_FUNC) &grFacSum_C,  5},
  {"hess",       (DL_FUNC) &hess,       13},
  {"hessC",      (DL_FUNC) &hessC,      11},
  {"nll",        (DL_FUNC) &nll,        17},
  {"NRalg",      (DL_FUNC) &NRalg,      29},
  {"NRalgv3",    (DL_FUNC) &NRalgv3,    24},
  {"pAO_C",      (DL_FUNC) &pAO_C,       4},
  {"pgumbel_C",  (DL_FUNC) &pgumbel_C,   5},
  {"pgumbel2_C", (DL_FUNC) &pgumbel2_C,  5},
  {"plgamma_C",  (DL_FUNC) &plgamma_C,   4},
  {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
  {"get_fitted", (DL_FUNC) &get_fitted, 4},
  {NULL, NULL, 0}
};

void R_init_ordinal(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
