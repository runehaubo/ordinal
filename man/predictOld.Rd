\name{predict.clm2}
\alias{predict.clm2}
\alias{predict.clmm2}
\title{Predict Method for CLM fits}
\description{
  Obtains predictions from a cumulative link (mixed) model.
}
\usage{
\method{predict}{clm2}(object, newdata, ...)

%% \method{predict}{clmm}(object, newdata, ...)
}
\arguments{
  \item{object}{a fitted object of class inheriting from
    \code{clm2} including \code{clmm2} objects.}
  \item{newdata}{optionally, a data frame in which to look for variables
    with which to predict. Observe that the response variable should
    also be present.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
  This method does not duplicate the behavior of
  \code{predict.polr} in package \code{MASS} which produces a
  matrix instead of a vector of predictions. The behavior of
  \code{predict.polr} can be mimiced as shown in the examples.

  If \code{newdata} is not supplied, the fitted values are obtained. For
  \code{clmm2} fits this means predictions that are controlled for the
  observed value of the random effects. If the predictions for a
  random effect of zero, i.e. an average 'subject', are wanted, the same
  data used to fit the model should be supplied in the \code{newdata}
  argument. For \code{clm2} fits those two sets of predictions are
  identical.
}
\value{
  A vector of predicted probabilities.
}
\author{Rune Haubo B Christensen}
\seealso{
    \code{\link[ordinal]{clm2}}, \code{\link[ordinal]{clmm2}}.
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))

## More manageable data set for less voluminous printing:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))
dat26

m1 <- clm2(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")
predict(m1)

mN1 <-  clm2(sureness ~ 1, nominal = ~prod, data = dat26,
            weights = wghts)
predict(mN1)

predict(update(m1, scale = ~.-prod))


#################################
## Mimicing the behavior of predict.polr:
if(require(MASS)) {
    ## Fit model from polr example:
    fm1 <- clm2(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    predict(fm1)

    set.seed(123)
    nlev <- 3
    y <- gl(nlev, 5)
    x <- as.numeric(y) + rnorm(15)
    fm.clm <- clm2(y ~ x)
    fm.polr <- polr(y ~ x)

    ## The equivalent of predict.polr(object, type = "probs"):
    (pmat.polr <- predict(fm.polr, type = "probs"))
    ndat <- expand.grid(y = gl(nlev,1), x = x)
    (pmat.clm <- matrix(predict(fm.clm, newdata = ndat), ncol=nlev,
                        byrow = TRUE))
    all.equal(c(pmat.clm), c(pmat.polr), tol = 1e-5) # TRUE

    ## The equivalent of predict.polr(object, type = "class"):
    (class.polr <- predict(fm.polr))
    (class.clm <- factor(apply(pmat.clm, 1, which.max)))
    all.equal(class.clm, class.polr) ## TRUE
}

}
\keyword{models}
