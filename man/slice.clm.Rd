\name{slice}
\alias{slice}
\alias{slice.clm}
\alias{plot.slice.clm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Slice the likelihood of a clm
}
\description{
  Slice likelihood and plot the slice. This is usefull for illustrating
  the likelihood surface around the MLE (maximum likelihood estimate)
  and provides graphics to substantiate (non-)convergence of a model
  fit. Also, the closeness of a quadratic approximation to the
  log-likelihood function can be inspected for relevant parameters. A
  slice is considerably less computationally demanding than a profile.
}
\usage{
slice(object, ...)

\method{slice}{clm}(object, parm = seq_along(par), lambda = 3,
     grid = 100, quad.approx = TRUE, ...)

\method{plot}{slice.clm}(x, parm = seq_along(x),
    type = c("quadratic", "linear"), plot.mle = TRUE,
    ask = prod(par("mfcol")) < length(parm) && dev.interactive(), ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{for the \code{clm} method an object of class
    \code{"clm"}, i.e., the result of a call to \code{clm}.
  }
  \item{x}{
    a \code{slice.clm} object, i.e., the result of
    \code{slice(clm.object)}.
  }
  \item{parm}{
    for \code{slice.clm} a numeric or character vector indexing
    parameters, for \code{plot.slice.clm} only a numeric vector is
    accepted. By default all parameters are selected.
  }
  \item{lambda}{
    the number of curvature units on each side of the MLE the slice
    should cover.
  }
  \item{grid}{
    the number of values at which to compute the log-likelihood for each
    parameter.
  }
  \item{quad.approx}{
    compute and include the quadratic approximation to the
    log-likelihood function?
  }
  \item{type}{
    \code{"quadratic"} plots the log-likelihood function which is
    approximately quadratic, and \code{"linear"} plots the
    signed square root of the log-likelihood function which is
    approximately linear.
  }
  \item{plot.mle}{
    include a vertical line at the MLE (maximum likelihood estimate)
    when \code{type = "quadratic"}? Ignored for \code{type = "linear"}.
  }
  \item{ask}{
    logical; if \code{TRUE}, the user is asked before each plot, see
    \code{\link{par}}\code{(ask=.)}.
  }
  \item{\dots}{
    further arguments to \code{plot.default} for the plot method. Not
    used in the slice method.
  }
}
%% \details{ bla
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  The \code{slice} method returns a list of \code{data.frame}s with one
  \code{data.frame} for each parameter slice. Each \code{data.frame}
  contains in the first column the values of the parameter and in the
  second column the values of the (positive) log-likelihood
  \code{"logLik"}. A third column is present if \code{quad.approx = TRUE}
  and contains the corresponding quadratic approximation to the
  log-likelihood. The original model fit is included as the attribute
  \code{"original.fit"}.

  The \code{plot} method produces a plot of the likelihood slice for
  each parameter.

}
\author{
  Rune Haubo B Christensen
}
\examples{

## fit model:
fm1 <- clm(rating ~ contact + temp, data = wine)
## slice the likelihood:
sl1 <- slice(fm1)

## three different ways to plot the slices:
par(mfrow = c(2,3))
plot(sl1)
plot(sl1, type = "quadratic", plot.mle = FALSE)
plot(sl1, type = "linear")

## Verify convergence to the optimum:
sl2 <- slice(fm1, lambda = 1e-5, quad.approx = FALSE)
plot(sl2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

