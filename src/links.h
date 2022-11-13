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
#ifndef _ORDINAL_LINKS_H_
#define _ORDINAL_LINKS_H_
/* That ifndef, etc. is an idiom to prevent the body of the header
 * being read more than once.
 */

#include <R.h>
#include <Rmath.h>

#ifdef	__cplusplus
extern "C" {
#endif
  /* That stanza allows the same header file to be used by C and C++
   * programs. There is a matching stanza at the end of this header
   * file.
   */
  
  /* Additional scalar cumulative probability functions */
  double d_pgumbel  (double,double,double,int);
  double d_pgumbel2 (double,double,double,int);
  double d_pAO      (double,double,int);
  double d_plgamma  (double,double,int);
  
  /* Additional scalar density functions */
  double d_dgumbel  (double,double,double,int);
  double d_dgumbel2 (double,double,double,int);
  double d_dAO      (double,double,int);
  double d_dlgamma  (double,double,int);
  
  /* Scalar density gradients */
  double d_glogis   (double);
  double d_gnorm    (double);
  double d_gcauchy  (double);
  double d_ggumbel  (double);
  double d_ggumbel2 (double);
  double d_gAO      (double,double);
  double d_glgamma  (double,double);
  
#ifdef	__cplusplus
}
#endif

#endif
