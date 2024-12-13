\name{clm.control}
\alias{clm.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Set control parameters for cumulative link models}
\description{
  Set control parameters for cumulative link models
}
\usage{
clm.control(method = c("Newton", "model.frame", "design", "ucminf", "nlminb",
   "optim"), 
   sign.location = c("negative", "positive"), 
   sign.nominal = c("positive", "negative"), 
   ..., trace = 0L,
   maxIter = 100L, gradTol = 1e-06, maxLineIter = 15L, relTol = 1e-6,
   tol = sqrt(.Machine$double.eps), maxModIter = 5L,
   convergence = c("warn", "silent", "stop", "message"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{method}{\code{"Newton"} fits the model by maximum likelihood and
    \code{"model.frame"} cause \code{\link{clm}} to return the
    \code{model.frame}, \code{"design"} causes \code{\link{clm}} to
    return a list of design matrices etc. that can be used with
    \code{\link{clm.fit}}. \code{ucminf}, \code{nlminb} and \code{optim} refer 
    to general purpose optimizers.
  }
  \item{sign.location}{change sign of the location part of the model.
  }
  \item{sign.nominal}{change sign of the nominal part of the model.
  }
  \item{trace}{numerical, if \code{> 0} information is printed about and during
    the optimization process. Defaults to \code{0}.
  }
  \item{maxIter}{the maximum number of Newton-Raphson iterations.
    Defaults to \code{100}.
  }
  \item{gradTol}{the maximum absolute gradient; defaults to \code{1e-6}.
  }
  \item{maxLineIter}{the maximum number of step halfings allowed if
    a Newton(-Raphson) step over shoots. Defaults to \code{15}.
  }
  \item{relTol}{relative convergence tolerence: relative change in the
    parameter estimates between Newton iterations. Defaults to \code{1e-6}.
  }
  \item{tol}{numerical tolerence on eigenvalues to determine
    negative-definiteness of Hessian. If the Hessian of a model fit is
    negative definite, the fitting algorithm did not converge. If the
    Hessian is singular, the fitting algorithm did converge albeit not
    to a \emph{unique} optimum, so one or more parameters are not
    uniquely determined even though the log-likelihood value is.
  }
  \item{maxModIter}{the maximum allowable number of consecutive
    iterations where the Newton step needs to be modified to be a decent
    direction. Defaults to \code{5}.
  }
  \item{convergence}{action to take if the fitting algorithm did not
    converge.
  }
  \item{\dots}{control arguments parsed on to \code{\link[ucminf]{ucminf}},
    \code{\link{nlminb}} or \code{\link{optim}}.
  }
}
\value{
  a list of control parameters.
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link{clm}}
}
\keyword{models}
