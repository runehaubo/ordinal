\name{clmm2.control}
\alias{clmm2.control}
\title{Set control parameters for cumulative link mixed models}
\description{
  Set control parameters for cumulative link mixed models
}
\usage{
clmm2.control(method = c("ucminf", "nlminb", "model.frame"), ...,
             trace = 0, maxIter = 50, gradTol = 1e-4,
             maxLineIter = 50,
             innerCtrl = c("warnOnly", "noWarn", "giveError"))
}
\arguments{
  \item{method}{
    the optimizer used to maximize the marginal likelihood function.
  }
  \item{\dots}{control arguments passed on to the chosen optimizer; see
    \code{\link[ucminf]{ucminf}}, \code{\link{optim}}, and
    \code{\link{nlminb}} for details.
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
  \item{innerCtrl}{the use of warnings/errors if the inner optimization
    fails to converge.
  }
}
\details{
  When the default optimizer, \code{ucminf} is used, the default values
  of that optimizers control options are changed to \code{grtol = 1e-5}
  and \code{grad = "central"}. 
}
\value{
  a list of control parameters.
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link{clmm2}}
}
\keyword{models}
