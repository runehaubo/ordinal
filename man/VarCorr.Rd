\name{VarCorr}
\alias{VarCorr}
\alias{VarCorr.clmm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Extract variance and correlation parameters
}
\description{
  The VarCorr function extracts the variance and (if present)
  correlation parameters for random effect terms in a
  cumulative link mixed model (CLMM) fitted with \code{clmm}.
}
\usage{
VarCorr(x, ...)

\method{VarCorr}{clmm}(x, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{a \code{\link{clmm}} object.
  }
  \item{\dots}{
    currently not used by the \code{clmm} method.
  }
}
\details{
  The \code{VarCorr} method returns a list of \code{data.frame}s; one for
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

  A list of matrices with variances in the diagonal and correlation
  parameters in the off-diagonal --- one matrix for each random effects term
  in the model. Standard deviations are provided as attributes to the
  matrices.

}
\author{
  Rune Haubo B Christensen
}
\examples{

fm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)
VarCorr(fm1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

