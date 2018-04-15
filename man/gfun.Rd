\name{gfun}
\alias{gnorm}
\alias{glogis}
\alias{gcauchy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Gradients of common densities
%%  ~~function to do ... ~~
}
\description{
  Gradients of common density functions in their standard forms, i.e.,
  with zero location (mean) and unit scale. These are implemented in C
  for speed and care is taken that the correct results are provided for
  the argument being \code{NA}, \code{NaN}, \code{Inf}, \code{-Inf} or
  just extremely small or large.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{

gnorm(x)

glogis(x)

gcauchy(x)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    numeric vector of quantiles.
  }
}
\details{
  The gradients are given by:
  \itemize{
    \item{gnorm: If \eqn{f(x)} is the normal density with mean 0 and
      spread 1, then the gradient is \deqn{f'(x) = -x f(x)}
    }
    \item{glogis: If \eqn{f(x)} is the logistic density with mean 0 and
      scale 1, then the gradient is
      \deqn{f'(x) = 2 \exp(-x)^2 (1 + \exp(-x))^{-3} -
	\exp(-x)(1+\exp(-x))^{-2}}
    }
    \item{pcauchy: If
      \eqn{f(x) = [\pi(1 + x^2)^2]^{-1}}{f(x) =1 / [pi (1 + x^2)^2]}
      is the cauchy density with mean 0 and scale 1, then the gradient
      is
      \deqn{f'(x) = -2x [\pi(1 + x^2)^2]^{-1}}{f'(x) = -2x / [pi (1 +
	x^2)^2]}
    }
  }
  
  These gradients are used in the Newton-Raphson algorithms in fitting
  cumulative link models with \code{\link{clm}} and cumulative link
  mixed models with \code{\link{clmm}}.
}
\value{
  a numeric vector of gradients.
}
\seealso{
  Gradients of densities are also implemented for the extreme value
  distribtion (\code{\link[=dgumbel]{gumbel}}) and the the log-gamma distribution
  (\code{\link[=lgamma]{log-gamma}}). 
}
\author{
  Rune Haubo B Christensen
}
\examples{

x <- -5:5
gnorm(x)
glogis(x)
gcauchy(x)

}
\keyword{distribution}

