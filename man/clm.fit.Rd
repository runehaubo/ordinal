\name{clm.fit}
\alias{clm.fit}
\alias{clm.fit.default}
\alias{clm.fit.factor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit Cumulative Link Models
%%  ~~function to do ... ~~
}
\description{
  A direct fitter of cumulative link models.
}
\usage{

clm.fit(y, ...)

\method{clm.fit}{default}(y, ...)

\method{clm.fit}{factor}(y, X, S, N, weights = rep(1, nrow(X)),
     offset = rep(0, nrow(X)), S.offset = rep(0, nrow(X)),
     control = list(), start, doFit=TRUE,
     link = c("logit", "probit", "cloglog", "loglog", "cauchit"),
     threshold = c("flexible", "symmetric", "symmetric2", "equidistant"),
     ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{for the default method a list of model components. For the
    factor method the response variable; a factor, preferably and ordered
    factor.
  }
  \item{X, S, N}{optional design matrices for the regression parameters,
    scale parameters and nominal parameters respectively.
  }
  \item{weights}{optional case weights.
  }
  \item{offset}{an optional offset.
  }
  \item{S.offset}{an optional offset for the scale part of the model.
  }
  \item{control}{a list of control parameters, optionally a call to
    \code{\link{clm.control}}.
  }
  \item{start}{an optional list of starting values of the form
    \code{c(alpha, beta, zeta)} for the thresholds and nominal effects
    (\code{alpha}), regression parameters (\code{beta}) and scale
    parameters (\code{zeta}).
  }
  \item{doFit}{logical for whether the model should be fit or the model
    environment should be returned.
  }
  \item{link}{the link function.
  }
  \item{threshold}{the threshold structure, see further at
    \code{\link{clm}}.
  }
  \item{\dots}{currently not used.}
}
\details{
  This function does almost the same thing that \code{\link{clm}} does:
  it fits a cumulative link model. The main differences are that
  \code{clm.fit} does not setup design matrices from formulae and only
  does minimal post processing after parameter estimation.

  Compared to \code{\link{clm}}, \code{clm.fit} does little to warn the
  user of any problems with data or model. However, \code{clm.fit} will
  attempt to identify column rank defecient designs. Any unidentified
  parameters are indicated in the \code{aliased} component of the fit.

  \code{clm.fit.factor} is not able to check if all thresholds are
  increasing when nominal effects are specified since it needs access to
  the terms object for the nominal model. If the terms object for the
  nominal model (\code{nom.terms}) is included in \code{y}, the default
  method is able to chech if all thresholds are increasing.

  %% In contrast to \code{\link{clm}}, \code{clm.fit} allows non-positive
  %% weights.
}

\value{
  A list with the following components:
  \code{aliased, alpha, coefficients, cond.H, convergence, df.residual,
    edf, fitted.values, gradient, Hessian, logLik, maxGradient, message,
    n, niter, nobs, tJac, vcov}
  and optionally
  \code{beta, zeta}
  These components are documented in \code{\link{clm}}.
}
%% \references{ bla
%%   %% ~put references to the literature/web site here ~
%% }
\author{
  Rune Haubo B Christensen
}
%% \note{ bla
%% %%  ~~further notes~~
%% }
%%
%% %% ~Make other sections like Warning with \section{Warning }{....} ~
%%
\seealso{
  \code{\link{clm}}
}
\examples{

## A simple example:
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)
## get the model frame containing y and X:
mf1 <- update(fm1, method="design")
names(mf1)
res <- clm.fit(mf1$y, mf1$X) ## invoking the factor method
stopifnot(all.equal(coef(res), coef(fm1)))
names(res)

## Fitting with the default method:
mf1$control$method <- "Newton"
res2 <- clm.fit(mf1)
stopifnot(all.equal(coef(res2), coef(fm1)))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}
