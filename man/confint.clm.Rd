\name{confint}
\alias{confint.clm}
\alias{confint.profile.clm}
\alias{profile.clm}
\alias{plot.profile.clm}
\title{
  Confidence intervals and profile likelihoods for parameters in
  cumulative link models
}
\description{
  Computes confidence intervals from the profiled likelihood for one or
  more parameters in a cumulative link model, or plots the
  profile likelihood.
}
\usage{

\method{confint}{clm}(object, parm, level = 0.95,
        type = c("profile", "Wald"), trace = FALSE, ...)

\method{confint}{profile.clm}(object, parm = seq_len(nprofiles),
        level = 0.95, ...)

\method{profile}{clm}(fitted, which.beta = seq_len(nbeta),
        which.zeta = seq_len(nzeta), alpha = 0.001,
        max.steps = 50, nsteps = 8, trace = FALSE, step.warn = 5,
        control = list(), ...)

\method{plot}{profile.clm}(x, which.par = seq_len(nprofiles),
        level = c(0.95, 0.99), Log = FALSE, relative = TRUE, root =
        FALSE, fig = TRUE, approx = root, n = 1e3,
        ask = prod(par("mfcol")) < length(which.par) && dev.interactive(),
        ..., ylim = NULL)
}
\arguments{
  \item{object, fitted, x}{
    a fitted \code{\link{clm}} object or a \code{profile.clm} object.
  }
  \item{parm, which.par, which.beta, which.zeta}{
    a numeric or character vector indicating which regression
    coefficients should be profiled. By default all coefficients are
    profiled. Ignored for \code{confint.clm} where all parameters are
    considered.
  }
  \item{level}{
    the confidence level. For the \code{plot} method a vector of levels
    for which horizontal lines should be drawn.
  }
  \item{type}{
    the type of confidence interval.
  }
  \item{trace}{
    if \code{trace} is \code{TRUE} or positive, information about
    progress is printed.
  }
  \item{Log}{
    should the profile likelihood be plotted on the log-scale?
  }
  \item{relative}{
    should the relative or the absolute likelihood be plotted?
  }
  \item{root}{
    should the (approximately linear) likelihood root statistic be
    plotted?
  }
  \item{approx}{
    should the Gaussian or quadratic approximation to the (log)
    likelihood be included?
  }
  \item{fig}{
    should the profile likelihood be plotted?
  }
  \item{ask}{
    logical; if \code{TRUE}, the user is asked before each plot, see
    \code{\link{par}}\code{(ask=.)}.
  }
  \item{n}{
    the no. points used in the spline interpolation of the profile
    likelihood.
  }
  \item{ylim}{overrules default y-limits on the plot of the profile
    likelihood.
  }
  \item{alpha}{
    The likelihood is profiled in the 100*(1-alpha)\% confidence region
    as determined by the profile likelihood.
  }
  \item{control}{
    a list of control parameters for \code{\link{clm}}. Possibly use
    \code{\link{clm.control}} to set these.
  }
  %%\item{lambda}{
  %%  logical. Should profile or confidence intervals be computed for the
  %%  link function parameter? Only used when one of the flexible link
  %%  functions are used; see the \code{link}-argument in
  %%  \code{\link{clm}}.
  %%}
  \item{max.steps}{
    the maximum number of profiling steps in each direction for each
    parameter.
  }
  \item{nsteps}{
    the (approximate) number of steps to take in each direction of the
    profile for each parameter.
    The step length is determined accordingly assuming a quadratic
    approximation to the log-likelihood function.
    The actual number of steps will often be close to \code{nsteps}, but
    will deviate when the log-likelihood functions is irregular.
  }
  \item{step.warn}{
    a warning is issued if the number of steps in each
    direction (up or down) for a parameter is less than
    \code{step.warn}. If few steps are taken, the profile will be
    unreliable and derived confidence intervals will be inaccurate.
  }
  \item{\dots}{
    additional arguments to be parsed on to methods.
  }

}
\value{
  \code{confint}:
  A matrix with columns giving lower and upper confidence
  limits for each parameter. These will be labelled as (1-level)/2 and
  1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).

  \code{plot.profile.clm} invisibly returns the profile object, i.e., a
  list of \code{\link{data.frame}}s with an \code{lroot} component for
  the likelihood root statistic and a  matrix \code{par.vals} with
  values of the parameters.
}
\details{
  These \code{confint} methods call
  the appropriate profile method, then finds the
  confidence intervals by interpolation of the profile traces.
  If the profile object is already available, this should be used as the
  main argument rather than the fitted model object itself.
  }
\author{Rune Haubo B Christensen}
\seealso{
\code{\link{profile}} and \code{\link{confint}}
}
\examples{

## Accurate profile likelihood confidence intervals compared to the
## conventional Wald intervals:
fm1 <- clm(rating ~ temp * contact, data = wine)
confint(fm1) ## type = "profile"
confint(fm1, type = "Wald")
pr1 <- profile(fm1)
confint(pr1)

## plotting the profiles:
par(mfrow = c(2, 2))
plot(pr1, root = TRUE) ## check for linearity
par(mfrow = c(2, 2))
plot(pr1)
par(mfrow = c(2, 2))
plot(pr1, approx = TRUE)
par(mfrow = c(2, 2))
plot(pr1, Log = TRUE)
par(mfrow = c(2, 2))
plot(pr1, Log = TRUE, relative = FALSE)
## Not likely to be useful but allowed for completeness:
par(mfrow = c(2, 2))
plot(pr1, Log = FALSE, relative = FALSE)

## Example from polr in package MASS:
## Fit model from polr example:
if(require(MASS)) {
    fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq,
               data = housing)
    pr1 <- profile(fm1)
    confint(pr1)
    par(mfrow=c(2,2))
    plot(pr1)
}

}
\keyword{models}
