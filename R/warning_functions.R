#############################################################################
##    Copyright (c) 2010-2022 Rune Haubo Bojesen Christensen
##
##    This file is part of the ordinal package for R (*ordinal*)
##
##    *ordinal* is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 2 of the License, or
##    (at your option) any later version.
##
##    *ordinal* is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    A copy of the GNU General Public License is available at
##    <https://www.r-project.org/Licenses/> and/or
##    <http://www.gnu.org/licenses/>.
#############################################################################
givesWarnings <- function(expr) countWarnings(expr) > 0L

countWarnings <- function(expr)
{
    .number_of_warnings <- 0L
    frame_number <- sys.nframe()
    ans <- withCallingHandlers(expr, warning = function(w)  {
        assign(".number_of_warnings", .number_of_warnings + 1L,
               envir = sys.frame(frame_number))
        invokeRestart("muffleWarning")
    })
    .number_of_warnings
}

