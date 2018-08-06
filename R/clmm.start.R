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
## This file contains:
## Functions to compute starting values for clmm()s.

clmm.start <- function(frames, link, threshold) {
  ## get starting values from clm:
  fit <- with(frames,
              clm.fit(y=y, X=X, weights=wts, offset=off, link=link,
                      threshold=threshold))

  ## initialize variance parameters to zero:
  start <- c(fit$par, rep(0, length(frames$grList)))
  return(start)
}

