\name{ranef}
\alias{ranef}
\alias{condVar}
\alias{ranef.clmm}
\alias{condVar.clmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract conditional modes and conditional variances from clmm objects
}
\description{
  The ranef function extracts the conditional modes of the random
  effects from a clmm object. That is, the modes of the distributions
  for the random effects given the observed data and estimated model
  parameters. In a Bayesian language they are posterior modes.

  The conditional variances are computed from the second order
  derivatives of the conditional distribution of the random
  effects. Note that these variances are computed at a fixed value of
  the model parameters and thus do not take the uncertainty of the
  latter into account.
}
\usage{
ranef(object, ...)

condVar(object, ...)

\method{ranef}{clmm}(object, condVar=FALSE, ...)

\method{condVar}{clmm}(object, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a \code{\link{clmm}} object.
  }
  \item{condVar}{
    an optional logical argument indicating of conditional variances
    should be added as attributes to the conditional modes.
  }
  \item{\dots}{
    currently not used by the \code{clmm} methods.
  }
}
\details{
  The \code{ranef} method returns a list of \code{data.frame}s; one for
  each distinct grouping factor. Each \code{data.frame} has as many rows
  as there are levels for that grouping factor and as many columns as
  there are random effects for each level. For example a model can
  contain a random intercept (one column) or a random
  intercept and a random slope (two columns) for the same grouping
  factor.

  If conditional variances are requested, they are returned in the same
  structure as the conditional modes (random effect
  estimates/predictions).
}
\value{
  The \code{ranef} method returns a list of \code{data.frame}s with the
  random effects predictions/estimates computed as conditional
  modes. If \code{condVar = TRUE} a \code{data.frame} with the
  conditional variances is stored as an attribute on each
  \code{data.frame} with conditional modes.

  The \code{condVar} method returns a list of \code{data.frame}s with
  the conditional variances. It is a convenience function that simply
  computes the conditional modes and variances, then extracts and
  returns only the latter.
}
\author{
  Rune Haubo B Christensen
}
\examples{

fm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)

## Extract random effect estimates/conditional modes:
re <- ranef(fm1, condVar=TRUE)

## Get conditional variances:
attr(re$judge, "condVar")
## Alternatively:
condVar(fm1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

