\name{clmm.control}
\alias{clmm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Set control parameters for cumulative link mixed models
}
\description{
  Set control parameters for cumulative link mixed models
}
\usage{
clmm.control(method = c("nlminb", "ucminf", "model.frame"), ..., trace = 0,
maxIter = 50, gradTol = 1e-4, maxLineIter = 50,  useMatrix = FALSE,
innerCtrl = c("warnOnly", "noWarn", "giveError"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{
    the optimizer used to maximize the marginal likelihood function.
  }
  \item{\dots}{control arguments passed on to the optimizer; see
    \code{\link[ucminf]{ucminf}} for details.
    \code{ucminf} for details.
  }
  \item{trace}{numerical, if > 0 information is printed about and during
    the outer optimization process, if < 0 information is also printed
    about the inner optimization process. Defaults to \code{0}.
  }
  \item{maxIter}{the maximum number of Newton updates of the inner
    optimization. \code{50}.
  }
  \item{gradTol}{the maximum absolute gradient of the inner
    optimization.
  }
  \item{maxLineIter}{the maximum number of step halfings allowed if
    a Newton(-Raphson) step over shoots during the inner optimization.
  }
  \item{useMatrix}{if \code{TRUE}, a general implementation of the
    Laplace approximation using the Matrix package is used, while if
    \code{FALSE} (default), a C implementation of the Laplace
    approximation valid only for models with a single random effects
    term is used when possible.
    \code{TRUE} is not valid for models fitted with quadrature methods.
  }
  \item{innerCtrl}{the use of warnings/errors if the inner optimization
    fails to converge.
  }
}
\value{
  a list of control parameters
}
\author{
  Rune Haubo B Christensen
}
\seealso{
  \code{\link{clmm}}
  }
\keyword{models}

