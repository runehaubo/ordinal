\name{ordinal-package}
\alias{ordinal-package}
\alias{ordinal}
\docType{package}
\title{
  Regression Models for Ordinal Data via Cumulative Link (Mixed) Models
}
\description{
  This package facilitates analysis of ordinal (ordered categorical
  data) via cumulative link models (CLMs) and cumulative link mixed
  models (CLMMs). Robust and efficient computational methods gives
  speedy and accurate estimation. A wide range of methods for model fits
  aids the data analysis.
}
\details{
\tabular{ll}{
Package: \tab ordinal\cr
Type: \tab Package\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package implements cumualtive link models and cumulative link
models with normally distributed random effects, denoted cumulative link
mixed (effects) models. Cumulative link models are also known as ordered
regression models, proportional odds models, proportional hazards models
for grouped survival times and ordered logit/probit/... models.

Cumulative link models are fitted with \code{\link{clm}} and the main
features are:
\itemize{
  \item{A range of standard link functions are available.}
  \item{In addition to the standard location (additive) effects, scale
    (multiplicative) effects are also allowed.}
  \item{nominal effects are allowed for any subset of the predictors ---
    these effects are also known as partial proportional odds effects
    when using the logit link.}
  \item{Restrictions can be imposed on the thresholds/cut-points, e.g.,
    symmetry or equidistance.}
  \item{A (modified) Newton-Raphson algorithm provides the maximum
    likelihood estimates of the parameters. The estimation scheme is robust,
    fast and accurate.}
  \item{Rank-deficient designs are identified and unidentified
    coefficients exposed in \code{print} and \code{summary} methods as
    with \code{\link{glm}}.}
  \item{A suite of standard methods are available including \code{anova},
    \code{add}/\code{drop}-methods, \code{step}, \code{profile},
    \code{confint}.}
  \item{A \code{slice} method facilitates illustration of
    the likelihood function and a \code{convergence} method summarizes
    the accuracy of the model estimation.}
  \item{The \code{predict} method can predict probabilities, response
    class-predictions and cumulative probabilities, and it provides
    standard errors and confidence intervals for the predictions.}
}

Cumulative link mixed models are fitted with \code{\link{clmm}} and the
main features are:
\itemize{
  \item{Any number of random effect terms can be included.}
  \item{The syntax for the model formula resembles that of \code{\link[lme4]{lmer}} from the \code{lme4} package.}
  \item{Nested random effects, crossed random effects and partially
    nested/crossed random effects are allowed.}
  \item{Estimation is via maximum likelihood using the Laplace
    approximation or adaptive Gauss-Hermite quadrature (one random
    effect).}
  \item{Vector-valued and correlated random effects such as random
    slopes (random coefficient models) are fitted with the Laplace
    approximation.}
  \item{Estimation employs sparse matrix methods from the
    \code{\link[Matrix]{Matrix}} package. }
  \item{During model fitting a Newton-Raphson algorithm updates the
    conditional modes of the random effects a large number of times. The
    likelihood function is optimized with a general purpose optimizer.}
}

A major update of the package in August 2011 introduced new and improved
implementations of \code{\link{clm}} and \code{\link{clmm}}. The old
implementations are available with \code{\link{clm2}} and
\code{\link{clmm2}}. At the time of writing there is functionality in
\code{clm2} and \code{clmm2} not yet available in \code{clm} and
\code{clmm}. This includes flexible link functions (log-gamma and
Aranda-Ordaz links) and a profile method for random effect variance
parameters in CLMMs. The new implementations are expected to take over
the old implementations at some point, hence the latter will eventually
be \code{\link[=.Deprecated]{deprecated}} and
\code{\link[=.Defunct]{defunct}}.

}
\author{
  Rune Haubo B Christensen

  Maintainer: Rune Haubo B Christensen <rune.haubo@gmail.com>
}
%% \references{
%%   ~~ Literature or other references for background information ~~
%% }
\keyword{ package }
%% \seealso{
%% ~~ Optional links to other man pages, e.g. ~~
%% %% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%% }
\examples{

## A simple cumulative link model:
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)

## A simple cumulative link mixed model:
fmm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)
summary(fmm1)

}
