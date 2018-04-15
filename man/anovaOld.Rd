\name{anova.clm2}
%%\alias{anova}
\alias{anova.clm2}
\alias{anova.clmm2}
\title{Likelihood ratio test of cumulative link models}
\description{
  Comparison of cumulative link models in likelihood ratio tests.
  The models may differ by terms in location, scale and nominal
  formulae, in link, threshold function and random effect structure.
}
\usage{
\method{anova}{clm2}(object, ..., test = c("Chisq", "none"))
\method{anova}{clmm2}(object, ..., test = c("Chisq", "none"))
}
\arguments{
  \item{object}{a \code{\link{clm2}} object.
  }
  \item{\dots}{one or more additional \code{\link{clm2}} objects.
  }
  \item{test}{if \code{test = "none"} the p-value for the likelihood
    ratio test is suppressed.
  }
}
\value{
  The method returns an object of class \code{Anova} (for printing) and
  \code{data.frame} with the following elements
  \item{Model}{character description of the cumulative link models being
    compared. Location, scale and nominal formulae are separated by
    "|"s in this order.
  }
  \item{Resid.df}{the residual degrees of freedom
  }
  \item{-2logLik}{twice the negative log likelihood (proportional to the
    deviance)}
  \item{Test}{indication of which models are being compared.
  }
  \item{DF}{the difference in the degrees of freedom in the models being
    compared, i.e. the degrees of freedom for the chi-squared test.
  }
  \item{LR stat.}{the likelihood ratio statistic.
  }
  \item{Pr(Chi)}{the p-value from the likelihood ratio test. Absent if
    \code{test = "none"}.
  }
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link[ordinal]{clm2}}, \code{\link[=addterm.clm2]{addterm}},
  \code{\link[ordinal:addtermOld]{dropterm}} and
  \code{\link[=anova]{anova.default}}
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))
m1 <- clm2(SURENESS ~ PROD, scale = ~PROD, data = soup,
          link = "logistic")

## anova
anova(m1, update(m1, scale = ~.-PROD))
mN1 <- clm2(SURENESS ~ 1, nominal = ~PROD, data = soup,
           link = "logistic")
anova(m1, mN1)
anova(m1, update(m1, scale = ~.-PROD), mN1)

## Fit model from polr example:
if(require(MASS)) {
    fm1 <- clm2(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    anova(fm1, update(fm1, scale =~ Cont))
}

}
\keyword{models}
