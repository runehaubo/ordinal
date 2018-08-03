\name{addterm.clm2}
\alias{addterm.clm2}
\alias{dropterm.clm2}
\title{
Try all one-term additions to and deletions from a model
}
\description{
Try fitting all models that differ from the current model by adding or
deleting a single term from those supplied while maintaining
marginality.
}
\usage{
\method{addterm}{clm2}(object, scope, scale = 0, test = c("none", "Chisq"),
        k = 2, sorted = FALSE, trace = FALSE,
        which = c("location", "scale"), \dots)
\method{dropterm}{clm2}(object, scope, scale = 0, test = c("none", "Chisq"),
        k = 2, sorted = FALSE, trace = FALSE,
        which = c("location", "scale"), \dots)
}
\arguments{
  \item{object}{
    A \code{\link{clm2}} object.
  }
  \item{scope}{
    for \code{addterm}:
    a formula specifying a maximal model which should include the current
    one. All additional terms in the maximal model with all marginal terms
    in the original model are tried.
    For \code{dropterm}:
    a formula giving terms which might be dropped. By default, the model
    formula. Only terms that can be dropped and maintain marginality are
    actually tried.
  }
  \item{scale}{
    used in the definition of the AIC statistic for selecting the
    models. Specifying \code{scale} asserts that the dispersion is known.
  }
  \item{test}{
    should the results include a test statistic relative to the original
    model?  The Chisq test is a likelihood-ratio test.
  }
  \item{k}{
    the multiple of the number of degrees of freedom used for the penalty.
    Only \code{k=2} gives the genuine AIC: \code{k = log(n)} is sometimes referred
    to as BIC or SBC.
  }
  \item{sorted}{
    should the results be sorted on the value of AIC?
  }
  \item{trace}{
    if \code{TRUE} additional information may be given on the fits as they are tried.
  }
  \item{which}{should additions or deletions occur in location or scale
    models?
  }
  \item{\dots}{
    arguments passed to or from other methods.
}}
\value{
  A table of class \code{"anova"} containing columns for the change
  in degrees of freedom, AIC and the likelihood ratio statistic. If
  \code{test = "Chisq"} a column also contains the
  p-value from the Chisq test.
}
\details{
The definition of AIC is only up to an additive constant because the
likelihood function is only defined up to an additive constant.
}
\author{Rune Haubo B Christensen}
\seealso{
  \code{\link[ordinal]{clm2}}, \code{\link[=anova.clm2]{anova}},
  \code{\link{addterm.default}} and \code{\link{dropterm.default}}
}
\examples{

options(contrasts = c("contr.treatment", "contr.poly"))

if(require(MASS)) { ## dropterm, addterm, housing
    mB1 <- clm2(SURENESS ~ PROD + GENDER + SOUPTYPE,
                scale = ~ COLD, data = soup, link = "probit",
                Hess = FALSE)
    dropterm(mB1, test = "Chi")       # or
    dropterm(mB1, which = "location", test = "Chi")
    dropterm(mB1, which = "scale", test = "Chi")
    addterm(mB1, scope = ~.^2, test = "Chi", which = "location")
    addterm(mB1, scope = ~ . + GENDER + SOUPTYPE,
            test = "Chi", which = "scale")
    addterm(mB1, scope = ~ . + AGEGROUP + SOUPFREQ,
            test = "Chi", which = "location")

    ## Fit model from polr example:
    fm1 <- clm2(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    addterm(fm1, ~ Infl + Type + Cont, test= "Chisq", which = "scale")
    dropterm(fm1, test = "Chisq")
}

}
\keyword{internal}
