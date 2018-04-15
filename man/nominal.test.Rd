\name{nominal_test}
\alias{nominal_test}
\alias{scale_test}
\alias{nominal_test.clm}
\alias{scale_test.clm}
\title{
  Likelihood ratio tests of model terms in scale and nominal formulae
}
\description{
  Add all model terms to scale and nominal formulae and perform
  likelihood ratio tests. These tests can be viewed as goodness-of-fit
  tests. With the logit link, \code{nominal_test} provides likelihood
  ratio tests of the proportional odds assumption. The \code{scale_test}
  tests can be given a similar interpretation.
}
\usage{
nominal_test(object, ...)

\method{nominal_test}{clm}(object, scope, trace=FALSE, ...)

scale_test(object, ...)

\method{scale_test}{clm}(object, scope, trace=FALSE, ...)

}
\arguments{
  \item{object}{for the \code{clm} method an object of class
    \code{"clm"}, i.e., the result of a call to \code{clm}.
  }
  \item{scope}{
    a formula or character vector specifying the terms to add to scale
    or nominal. In \code{nominal_test} terms in scope already in
    \code{nominal} are ignored. In \code{scale_test} terms in scope
    already in \code{scale} are ignored.

    In \code{nominal_test} the default is to add all terms
    from \code{formula} (location part) and \code{scale} that are not
    also in \code{nominal}.

    In \code{scale_test} the default is to add
    all terms from \code{formula} (location part) that are not also in
    \code{scale}.
  }
  \item{trace}{
    if \code{TRUE} additional information may be given on the fits as
    they are tried.
  }
  \item{\dots}{
    arguments passed to or from other methods.
  }
}
\value{
  A table of class \code{"anova"} containing columns for the change
  in degrees of freedom, AIC, the likelihood ratio statistic and a
  p-value based on the asymptotic chi-square distribtion of the
  likelihood ratio statistic under the null hypothesis.
}
\details{
  The definition of AIC is only up to an additive constant because the
  likelihood function is only defined up to an additive constant.
}
\author{Rune Haubo B Christensen}
\examples{

## Fit cumulative link model:
fm <- clm(rating ~ temp + contact, data=wine)
summary(fm)
## test partial proportional odds assumption for temp and contact:
nominal_test(fm)
## no evidence of non-proportional odds.
## test if there are signs of scale effects:
scale_test(fm)
## no evidence of scale effects.

## tests of scale and nominal effects for the housing data from MASS:
if(require(MASS)) {
    fm1 <- clm(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    scale_test(fm1)
    nominal_test(fm1)
    ## Evidence of multiplicative/scale effect of 'Cont'. This is a breach
    ## of the proportional odds assumption.
}

}
\keyword{models}
