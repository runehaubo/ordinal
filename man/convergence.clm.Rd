\name{convergence}
\alias{convergence}
\alias{convergence.clm}
\alias{print.convergence.clm}
\title{Check convergence of cumulative link models}
\description{
  Check the accuracy of the parameter estimates of cumulative link
  models. The number of correct decimals and number of significant
  digits is given for the maximum likelihood estimates of the parameters
  in a cumulative link model fitted with \code{\link{clm}}.
}
\usage{

convergence(object, ...)

\method{convergence}{clm}(object, digits = max(3, getOption("digits") - 3),
   tol = sqrt(.Machine$double.eps), ...)

}
\arguments{
  \item{object}{for the \code{clm} method an object of class
    \code{"clm"}, i.e., the result of a call to \code{clm}.
  }
  \item{digits}{the number of digits in the printed table.
  }
  \item{tol}{numerical tolerence to judge if the Hessian is positive
    definite from its smallest eigenvalue.
  }
  \item{...}{arguments to a from methods. Not used by the \code{clm} method.
  }
}
\value{
  Convergence information. In particular a table where the \code{Error}
  column gives the numerical error in the parameter estimates. These
  numbers express how far the parameter estimates in the fitted model
  are from the true maximum likelihood estimates for this
  model. The \code{Cor.Dec} gives the number of correct decimals with
  which the the parameters are determined and the \code{Sig.Dig} gives
  the number of significant digits with which the parameters are
  determined.

  The number denoted \code{logLik.error} is the error in the value of
  log-likelihood in the fitted model at the parameter values of that
  fit. An accurate determination of the log-likelihood is essential for
  accurate likelihood ratio tests in model comparison.
}
\details{
  The number of correct decimals is defined as...

  The number of significant digits is defined as ...

  The number of correct decimals and the number of significant digits
  are determined from the numerical errors in the parameter
  estimates. The numerical errors are determined from the Method
  Independent Error Theorem (Elden et al, 2004) and is based on the
  Newton step evaluated at convergence.

}
\references{
  Elden, L., Wittmeyer-Koch, L. and Nielsen, H. B. (2004) \emph{Introduction
    to Numerical Computation --- analysis and Matlab illustrations.}
  Studentliteratur.
}
%% \seealso{
%% }
\examples{

## Simple model:
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)
convergence(fm1)

}
\author{Rune Haubo B Christensen}
\keyword{models}
