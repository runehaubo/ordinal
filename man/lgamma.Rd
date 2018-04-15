\name{lgamma}
\alias{plgamma}
\alias{dlgamma}
\alias{glgamma}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  The log-gamma distribution
%%  ~~function to do ... ~~
}
\description{
  Density, distribution function and gradient of density for the
  log-gamma distribution. 
  These are implemented in C
  for speed and care is taken that the correct results are provided for
  values of \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or just
  extremely small or large values.

  The log-gamma is a flexible location-scale distribution on the real
  line with an extra parameter, \eqn{\lambda}. For \eqn{\lambda = 0} the
  distribution equals the normal or Gaussian distribution, and for
  \eqn{\lambda} equal to 1 and -1, the Gumbel minimum and maximum
  distributions are obtained.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

plgamma(q, lambda, lower.tail = TRUE)

dlgamma(x, lambda, log = FALSE)

glgamma(x, lambda)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x,q}{
    numeric vector of quantiles.
  }
  \item{lambda}{
    numerical scalar
  }
%%   \item{location}{
%%     numeric scalar.
%%   }
%%   \item{scale}{
%%     numeric scalar.
%%   }
  \item{lower.tail}{
    logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \leq x]}{P[X <= x]} otherwise, \eqn{P[X > x]}.
  }
  \item{log}{
    logical; if \code{TRUE}, probabilities p are given as log(p).
  }
}
\details{
  If \eqn{\lambda < 0} the distribution is right skew, if
  \eqn{\lambda = 0} the distribution is symmetric (and equals the normal
  distribution), and if \eqn{\lambda > 0} the distribution is left
  skew.
  
  The log-gamma distribution function is defined as \ldots pending.
  
  The density and gradient of the density are defined as\ldots pending.
  
  These distribution functions, densities and gradients are used in the
  Newton-Raphson algorithms in fitting cumulative link models with
  \code{\link{clm2}} and cumulative link mixed models with
  \code{\link{clmm2}} using the log-gamma link. 
}
\value{
  \code{plgamma} gives the distribution function, \code{dlgamma}
  gives the density and \code{glgamma} gives the gradient of the
  density. 
}
\references{
  Genter, F. C. and Farewell, V. T. (1985) Goodness-of-link testing in
  ordinal regression models. \emph{The Canadian Journal of Statistics},
  13(1), 37-44.
}
\seealso{
  Gradients of densities are also implemented for the normal, logistic, 
  cauchy, cf. \code{\link[=gnorm]{gfun}} and the Gumbel distribution,
  cf. \code{\link[=dgumbel]{gumbel}}.
}
\author{
  Rune Haubo B Christensen
}
\examples{

## Illustrating the link to other distribution functions: 
x <- -5:5
plgamma(x, lambda = 0) == pnorm(x)
all.equal(plgamma(x, lambda = -1), pgumbel(x)) ## TRUE, but:
plgamma(x, lambda = -1) == pgumbel(x)
plgamma(x, lambda = 1) == pgumbel(x, max = FALSE)

dlgamma(x, lambda = 0) == dnorm(x)
dlgamma(x, lambda = -1) == dgumbel(x)
dlgamma(x, lambda = 1) == dgumbel(x, max = FALSE)

glgamma(x, lambda = 0) == gnorm(x)
all.equal(glgamma(x, lambda = -1), ggumbel(x)) ## TRUE, but:
glgamma(x, lambda = -1) == ggumbel(x)
all.equal(glgamma(x, lambda = 1), ggumbel(x, max = FALSE)) ## TRUE, but:
glgamma(x, lambda = 1) == ggumbel(x, max = FALSE)
## There is a loss of accuracy, but the difference is very small: 
glgamma(x, lambda = 1) - ggumbel(x, max = FALSE)

## More examples:
x <- -5:5
plgamma(x, lambda = .5)
dlgamma(x, lambda = .5)
glgamma(x, lambda = .5)

}
\keyword{distribution}

