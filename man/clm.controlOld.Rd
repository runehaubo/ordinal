\name{clm2.control}
\alias{clm2.control}
\title{Set control parameters for cumulative link models}
\description{
  Set control parameters for cumulative link models
}
\usage{
clm2.control(method = c("ucminf", "Newton", "nlminb", "optim",
            "model.frame"), ..., convTol = 1e-4,
            trace = 0, maxIter = 100, gradTol = 1e-5,
            maxLineIter = 10)
}
\arguments{
  \item{method}{
    the optimizer used to maximize the likelihood
    function. \code{"Newton"} only works for models without \code{scale},
    structured thresholds and flexible link functions,
    but is considerably faster than the other
    optimizers when applicable. \code{model.frame} simply returns a list
    of model frames with the location, scale and nominal model
    frames. \code{"optim"} uses the \code{"BFGS"} method.
  }
  \item{\dots}{control arguments passed on to the chosen optimizer; see
    \code{\link[ucminf]{ucminf}}, \code{\link{optim}}, and
    \code{\link{nlminb}} for details.
  }
  \item{convTol}{convergence criterion on the size of the maximum
    absolute gradient.
  }
  \item{trace}{numerical, if > 0 information is printed about and during
    the optimization process. Defaults to \code{0}.
  }
  \item{maxIter}{the maximum number of Newton-Raphson iterations.
    Defaults to \code{100}.
  }
  \item{gradTol}{the maximum absolute gradient. This is the termination
    criterion and defaults to \code{1e-5}.
  }
  \item{maxLineIter}{the maximum number of step halfings allowed if
    a Newton(-Raphson) step over shoots. Defaults to \code{10}.
  }
}
\value{
  a list of control parameters.
}
\author{Rune Haubo B Christensen}
\seealso{
\code{\link{clm2}}
}
\keyword{models}
