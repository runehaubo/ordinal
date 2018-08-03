\name{confint.clm2}
\alias{confint.clm2}
\alias{confint.profile.clm2}
\alias{profile.clm2}
\alias{plot.profile.clm2}
\title{
  Confidence intervals and profile likelihoods for parameters in
  cumulative link models
}
\description{
  Computes confidence intervals from the profiled likelihood for one or
  more parameters in a fitted cumulative link model, or plots the
  profile likelihood function.
}
\usage{
\method{confint}{clm2}(object, parm, level = 0.95, whichL = seq_len(p),
        whichS = seq_len(k), lambda = TRUE, trace = 0, \dots)

\method{confint}{profile.clm2}(object, parm = seq_along(Pnames), level = 0.95, \dots)

\method{profile}{clm2}(fitted, whichL = seq_len(p), whichS = seq_len(k),
        lambda = TRUE, alpha = 0.01, maxSteps = 50, delta = LrootMax/10,
        trace = 0, stepWarn = 8, \dots)

\method{plot}{profile.clm2}(x, parm = seq_along(Pnames), level = c(0.95, 0.99),
        Log = FALSE, relative = TRUE, fig = TRUE, n = 1e3, ..., ylim = NULL)
}
\arguments{
  \item{object}{
    a fitted \code{\link{clm2}} object or a \code{profile.clm2} object.
  }
  \item{fitted}{
    a fitted \code{\link{clm2}} object.
  }
  \item{x}{a \code{profile.clm2} object.
  }
  \item{parm}{not used in \code{confint.clm2}.

    For \code{confint.profile.clm2}:
    a specification of which parameters are to be given confidence
    intervals, either a vector of numbers or a vector of names. If
    missing, all parameters are considered.

    For \code{plot.profile.clm2}:
    a specification of which parameters the profile likelihood are to be
    plotted for, either a vector of numbers or a vector of names. If
    missing, all parameters are considered.
  }
  \item{level}{
    the confidence level required.
  }
  \item{whichL}{
    a specification of which \emph{location} parameters are to be given confidence
    intervals, either a vector of numbers or a vector of names. If
    missing, all location parameters are considered.

  }
  \item{whichS}{
    a specification of which \emph{scale} parameters are to be given confidence
    intervals, either a vector of numbers or a vector of names. If
    missing, all scale parameters are considered.
  }
  \item{lambda}{
    logical. Should profile or confidence intervals be computed for the
    link function parameter? Only used when one of the flexible link
    functions are used; see the \code{link}-argument in
    \code{\link{clm2}}.
  }
  \item{trace}{
    logical.  Should profiling be traced?
  }
  \item{alpha}{Determines the range of profiling. By default the
    likelihood is profiled in the 99\% confidence interval region as
    determined by the profile likelihood.
  }
  \item{maxSteps}{the maximum number of profiling steps in each
    direction (up and down) for each parameter.
  }
  \item{delta}{the length of profiling steps. To some extent this
    parameter determines the degree of accuracy of the profile
    likelihood in that smaller values, i.e. smaller steps gives a higher
    accuracy. Note however that a spline interpolation is used when
    constructing confidence intervals so fairly long steps can provide
    high accuracy.
  }
  \item{stepWarn}{a warning is issued if the no. steps in each direction
    (up or down) for a parameter is less than \code{stepWarn} (defaults
    to 8 steps) because this indicates an unreliable profile.
  }
  \item{Log}{should the profile likelihood be plotted on the log-scale?
  }
  \item{relative}{should the relative or the absolute likelihood be
    plotted?
  }
  \item{fig}{should the profile likelihood be plotted?
  }
  \item{n}{the no. points used in the spline interpolation of the
    profile likelihood.
  }
  \item{ylim}{overrules default y-limits on the plot of the profile
    likelihood.
  }
  \item{\dots}{
    additional argument(s) for methods including \code{range} (for the
    hidden function \code{profileLambda}) that sets
    the range of values of \code{lambda} at which the likelihood should
    be profiled for this parameter.
  }

}
\value{
  \code{confint}:
  A matrix (or vector) with columns giving lower and upper confidence
  limits for each parameter. These will be labelled as (1-level)/2 and
  1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).
  The parameter names are preceded with \code{"loc."} or \code{"sca."}
  to indicate whether the confidence interval applies to a location or a
  scale parameter.

  \code{plot.profile.clm2} invisibly returns the profile object.

}
\details{
  These \code{confint} methods call
  the appropriate profile method, then finds the
  confidence intervals by interpolation of the profile traces.
  If the profile object is already available, this should be used as the
  main argument rather than the fitted model object itself.

  In \code{plot.profile.clm2}: at least one of \code{Log} and
  \code{relative} arguments have to be \code{TRUE}.

}
\author{Rune Haubo B Christensen}
\seealso{
\code{\link{profile}} and \code{\link{confint}}
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))

## More manageable data set:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))

m1 <- clm2(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")

## profile
pr1 <- profile(m1)
par(mfrow = c(2, 2))
plot(pr1)

m9 <- update(m1, link = "log-gamma")
pr9 <- profile(m9, whichL = numeric(0), whichS = numeric(0))
par(mfrow = c(1, 1))
plot(pr9)

plot(pr9, Log=TRUE, relative = TRUE)
plot(pr9, Log=TRUE, relative = TRUE, ylim = c(-4, 0))
plot(pr9, Log=TRUE, relative = FALSE)

## confint
confint(pr9)
confint(pr1)

## Extend example from polr in package MASS:
## Fit model from polr example:
if(require(MASS)) {
    fm1 <- clm2(Sat ~ Infl + Type + Cont, scale = ~ Cont, weights = Freq,
                data = housing)
    pr1 <- profile(fm1)
    confint(pr1)
    par(mfrow=c(2,2))
    plot(pr1)
}

}
\keyword{internal}
