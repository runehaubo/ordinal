\name{profile.clmm2}
\alias{profile.clmm2}
\alias{confint.clmm2}
\alias{confint.profile.clmm2}
\alias{profile.clmm2}
\alias{plot.profile.clmm2}
\title{
  Confidence intervals and profile likelihoods for the standard
  deviation for the random term in cumulative link mixed models
}
\description{
  Computes confidence intervals from the profiled likelihood for
  the standard devation for the random term in a fitted cumulative link
  mixed model, or plots the associated profile likelihood function.
}
\usage{
\method{confint}{profile.clmm2}(object, parm = seq_along(Pnames), level = 0.95, \dots)

\method{profile}{clmm2}(fitted, alpha = 0.01, range, nSteps = 20, trace = 1, \dots)

\method{plot}{profile.clmm2}(x, parm = seq_along(Pnames), level = c(0.95, 0.99),
        Log = FALSE, relative = TRUE, fig = TRUE, n = 1e3, ..., ylim = NULL)
}
\arguments{
  \item{object}{
    a fitted \code{profile.clmm2} object.
  }
  \item{fitted}{
    a fitted \code{\link{clmm2}} object.
  }
  \item{x}{a \code{profile.clmm2} object.
  }
  \item{parm}{
    For \code{confint.profile.clmm2}:
    a specification of which parameters are to be given confidence
    intervals, either a vector of numbers or a vector of names. If
    missing, all parameters are considered.
    Currently only \code{"stDev"} or \code{1} are supported.

    For \code{plot.profile.clmm2}:
    a specification of which parameters the profile likelihood are to be
    plotted for, either a vector of numbers or a vector of names. If
    missing, all parameters are considered.
    Currently only \code{"stDev"} or \code{1} are supported.
  }
  \item{level}{
    the confidence level required. Observe that the model has to be
    profiled in the appropriate region; otherwise the limits are
    \code{NA}.
  }
  \item{trace}{
    logical. Should profiling be traced? Defaults to \code{TRUE} due to
    the time consuming nature of the computation.
  }
  \item{alpha}{Determines the range of profiling. By default the
    likelihood is profiled approximately in the 99\% confidence interval
    region as determined by the Wald approximation. This is usually
    sufficient for 95\% profile likelihood confidence limits.
  }
  \item{range}{if range is specified, this overrules the range
    computation based on \code{alpha}. \code{range} should be all
    positive and \code{stDev} is profiled in \code{range(range)}.
  }
  \item{nSteps}{the number of points at which to profile the likelihood
    function. This determines the resolution and accuracy of the profile
    likelihood function; higher values gives a higher resolution, but
    also longer computation times.
  }
  \item{Log}{should the profile likelihood be plotted on the log-scale?
  }
  \item{relative}{should the relative or the absolute likelihood be
    plotted?
  }
  \item{fig}{should the profile likelihood be plotted?
  }
  \item{n}{the no. points used in the spline interpolation of the
    profile likelihood for plotting.
  }
  \item{ylim}{overrules default y-limits on the plot of the profile
    likelihood.
  }
  \item{\dots}{
    additional argument(s), e.g. graphical parameters for the
    \code{plot} method.
  }

}
\details{
  A \code{confint.clmm2} method deliberately does not exist due to the
  time consuming nature of the computations. The user is required to
  compute the profile object first and then call \code{confint} on the
  profile object to obtain profile likelihood confidence intervals.

  In \code{plot.profile.clm2}: at least one of \code{Log} and
  \code{relative} arguments have to be \code{TRUE}.
}
\value{
  \code{confint}:
  A matrix with columns giving lower and upper confidence
  limits. These will be labelled as (1-level)/2 and
  1 - (1-level)/2 in \% (by default 2.5\% and 97.5\%).

  \code{plot.profile.clm2} invisibly returns the profile object.
}
\author{Rune Haubo B Christensen}
\seealso{
\code{\link{profile}} and \code{\link{confint}}
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))

if(require(lme4)) { ## access cbpp data
    cbpp2 <- rbind(cbpp[,-(2:3)], cbpp[,-(2:3)])
    cbpp2 <- within(cbpp2, {
        incidence <- as.factor(rep(0:1, each=nrow(cbpp)))
        freq <- with(cbpp, c(incidence, size - incidence))
    })

    ## Fit with Laplace approximation:
    fm1 <- clmm2(incidence ~ period, random = herd, weights = freq,
                 data = cbpp2, Hess = 1)

    pr.fm1 <- profile(fm1)
    confint(pr.fm1)

    par(mfrow = c(2,2))
    plot(pr.fm1)
    plot(pr.fm1, Log=TRUE, relative = TRUE)
    plot(pr.fm1, Log=TRUE, relative = FALSE)
}

}
\keyword{models}
