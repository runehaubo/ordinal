\name{anova.clm}
%%\alias{anova}
\alias{anova.clm}
\title{ANODE Tables and Likelihood ratio test of cumulative link models}
\description{
  Type I, II, and III analysis of deviance (ANODE) tables for 
  cumulative link models and
  comparison of cumulative link models with likelihood ratio tests.
  Models may differ by terms in location, scale and nominal
  formulae, in link, threshold function.
}
\usage{
\method{anova}{clm}(object, ..., type = c("I", "II", "III", "1", "2", "3"))
}
\arguments{
  \item{object}{a \code{\link{clm}} object.
  }
  \item{\dots}{optionally one or more additional \code{\link{clm}} objects.
  }
  \item{type}{the type of hypothesis test if \code{anova} is called with a 
  single model; ignored if more than one model is passed to the method.
  }
}
\details{
  The ANODE table returned when \code{anova} is called with a single model apply only to 
  terms in \code{formula}, that is, terms in \code{nominal} and \code{scale} are
  ignored. 
}
\value{
  An analysis of deviance table based on Wald chi-square test if called with a
  single model and a comparison of
  models with likelihood ratio tests if called with more than one model.
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link[ordinal]{clm}}
}
\examples{

## Analysis of deviance tables with Wald chi-square tests:
fm <- clm(rating ~ temp * contact, scale=~contact, data=wine)
anova(fm, type="I")
anova(fm, type="II")
anova(fm, type="III")

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
