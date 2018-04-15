\name{clm2}
\alias{clm2}
\title{Cumulative link models}
\description{
  A new improved implementation of CLMs is available in \code{\link{clm}}.

  Fits cumulative link models with an additive model for the location
  and a multiplicative model for the scale. The function allows for
  structured thresholds. A popular special case of a CLM is the
  proportional odds model. In addition to the standard link functions,
  two flexible link functions, "Arandar-Ordaz" and "log-gamma" are
  available, where an extra link function parameter provides additional
  flexibility. A subset of the predictors can be allowed to have nominal
  rather than ordinal effects. This has been termed "partial
  proportional odds" when the link is the logistic.
}
\usage{
clm2(location, scale, nominal, data, weights, start, subset,
    na.action, contrasts, Hess = TRUE, model,
    link = c("logistic", "probit", "cloglog", "loglog",
    "cauchit", "Aranda-Ordaz", "log-gamma"), lambda,
    doFit = TRUE, control,
    threshold = c("flexible", "symmetric", "equidistant"), ...)
}
\arguments{
\item{location}{
  a formula expression as for regression models, of the form
  \code{response ~ predictors}. The response should be a factor
  (preferably an ordered factor), which will be interpreted as an
  ordinal response with levels ordered as in the factor.
  The model must have an intercept: attempts to remove one will
  lead to a warning and will be ignored. An offset may be used.  See the
  documentation of \code{\link{formula}} for other details.
}
\item{scale}{
  a optional formula expression as for the location part, of the form
  \code{ ~ predictors}, i.e. with an empty left hand side.
  An offset may be used. See the
  documentation of \code{\link{formula}} for other details.
}
\item{nominal}{
  an optional formula of the form \code{ ~ predictors}, i.e. with an
  empty left hand side. The effects of the predictors in this formula are
  assumed to nominal.
}
\item{data}{
  an optional data frame in which to interpret the variables occurring
  in the formulas.
}
\item{weights}{
  optional case weights in fitting. Defaults to 1.
}
\item{start}{
  initial values for the parameters in the format
  \code{c(alpha, beta, log(zeta), lambda)}.
}
\item{subset}{
  expression saying which subset of the rows of the data should  be used
  in the fit. All observations are included by default.
}
\item{na.action}{
  a function to filter missing data. Applies to terms in all three formulae.
}
\item{contrasts}{
  a list of contrasts to be used for some or all of
  the factors appearing as variables in the model formula.
}
\item{Hess}{
  logical for whether the Hessian (the inverse of the observed
  information matrix)
  should be computed.
  Use \code{Hess = TRUE} if you intend to call \code{summary} or
  \code{vcov} on the fit and \code{Hess = FALSE} in all other instances
  to save computing time. The argument is ignored if
  \code{method = "Newton"} where the Hessian is always computed and
  returned. Defaults to \code{TRUE}.
}
\item{model}{
  logical for whether the model frames should be part of the returned
  object.
}
\item{link}{link function, i.e. the type of location-scale distribution
  assumed for the latent distribution. The \code{Aranda-Ordaz} and
  \code{log-gamma} links add additional flexibility with a link function
  parameter, \code{lambda}. The \code{Aranda-Ordaz} link
  (Aranda-Ordaz, 1983) equals the logistic
  link, when \code{lambda = 1} and approaches the \code{loglog} link when
  \code{lambda} approaches zero. The \code{log-gamma} link (Genter and
  Farewell, 1985) equals the
  \code{loglog} link when \code{lambda = 1}, the \code{probit} link when
  \code{lambda = 0} and the \code{cloglog} link when \code{lambda =
    -1}.
}
\item{lambda}{numerical scalar: the link function parameter. Used in
  combination with link \code{Aranda-Ordaz} or \code{log-gamma} and
  otherwise ignored. If lambda is specified, the model is estimated with
  lambda fixed at this value and otherwise lambda is estimated by
  ML. For \code{Aranda-Ordaz} lambda has to be positive; \code{> 1e-5}
  for numerical reasons.
}
\item{doFit}{logical for whether the model should be fit or the model
  environment should be returned.
}
\item{control}{a call to \code{\link{clm2.control}}.
}
\item{threshold}{specifies a potential structure for the thresholds
  (cut-points). \code{"flexible"} provides the standard unstructured
  thresholds, \code{"symmetric"} restricts the distance between the
  thresholds to be symmetric around the central one or two thresholds
  for odd or equal numbers or thresholds respectively, and
  \code{"equidistant"} restricts the distance between consecutive
  thresholds to the same value.
}
\item{\dots}{
  additional arguments are passed on to \code{\link{clm2.control}} and
  possibly further on to the optimizer, which can lead to surprising
  error or warning messages when mistyping arguments etc.
}

}
\details{
  There are methods for the standard model-fitting functions, including
  \code{\link{summary}}, \code{\link{vcov}},
  \code{\link[ordinal]{predict}},
  \code{\link[=anova.clm2]{anova}}, \code{\link{logLik}},
  \code{\link[=profile.clm2]{profile}},
  \code{\link[=profile.clm2]{plot.profile}},
  \code{\link[=confint.clm2]{confint}},
  \code{\link[=update.clm2]{update}},
  \code{\link[=addterm.clm2]{dropterm}},
  \code{\link[=addterm.clm2]{addterm}}, and an \code{extractAIC} method.

  The design of the implementation is inspired by an idea proposed by
  Douglas Bates in the talk "Exploiting sparsity in model matrices"
  presented at the DSC conference in Copenhagen, July 14 2009.
  Basically an environment is set up with all the information needed to
  optimize the likelihood function. Extractor functions are then used to
  get the value of likelihood at current or given parameter values and
  to extract current values of the parameters. All computations are
  performed inside the environment and relevant variables are updated
  during the fitting process. After optimizer termination relevant
  variables are extracted from the environment and the remaining are
  discarded.

  Some aspects of \code{clm2}, for instance, how starting values are
  obtained, and of the associated methods are
  inspired by \code{\link[MASS]{polr}} from package \code{MASS}.

}
\value{
  If \code{doFit = FALSE} the result is an environment
  representing the model ready to be optimized.
  If \code{doFit = TRUE} the result is an
  object of class \code{"clm2"} with the following components:

  \item{beta}{the parameter estimates of the location part.
  }
  \item{zeta}{the parameter estimates of the scale part on the log
    scale; the scale parameter estimates on the original scale are given
    by \code{exp(zeta)}.
  }
  \item{Alpha}{vector or matrix of the threshold parameters.
  }
  \item{Theta}{vector or matrix of the thresholds.
  }
  \item{xi}{vector of threshold parameters, which, given a
    threshold function (e.g. \code{"equidistant"}), and possible nominal
    effects define the class boundaries, \code{Theta}.
  }
  \item{lambda}{the value of lambda if lambda is supplied or estimated,
    otherwise missing.
  }
  \item{coefficients}{the coefficients of the intercepts
    (\code{theta}), the location  (\code{beta}), the scale
    (\code{zeta}), and the link function parameter (\code{lambda}).
  }
  \item{df.residual}{the number of residual degrees of freedoms,
    calculated using the weights.
  }
  \item{fitted.values}{vector of fitted values for each
    observation. An observation here is each of the scalar elements of
    the multinomial table and not a multinomial vector.
  }
  \item{convergence}{\code{TRUE} if the gradient
    based convergence criterion is met and \code{FALSE} otherwise.
  }
  \item{gradient}{vector of gradients for all the parameters
    at termination of the optimizer.
  }
  \item{optRes}{list with results from the optimizer. The contents of the
    list depends on the choice of optimizer.
  }
  \item{logLik}{the log likelihood of the model at
    optimizer termination.
  }
  \item{Hessian}{if the model was fitted with \code{Hess = TRUE}, this
    is the Hessian matrix of the parameters at the optimum.
  }
  \item{scale}{\code{model.frame} for the scale model.
  }
  \item{location}{\code{model.frame} for the location model.
  }
  \item{nominal}{\code{model.frame} for the nominal model.
  }
  \item{edf}{the (effective) number of degrees of freedom used by the
    model.
  }
  \item{start}{the starting values.
  }
  \item{convTol}{convergence tolerance for the maximum absolute
    gradient of the parameters at termination of the optimizer.
  }
  \item{method}{character, the optimizer.
  }
  \item{y}{the response variable.
  }
  \item{lev}{the names of the levels of the response variable.
  }
  \item{nobs}{the (effective) number of observations, calculated as the
    sum of the weights.
  }
  \item{threshold}{character, the threshold function used in the model.
  }
  \item{estimLambda}{\code{1} if lambda is estimated in one of the
    flexible link functions and \code{0} otherwise.
  }
  \item{link}{character, the link function used in the model.
  }
  \item{call}{the matched call.
  }
  \item{contrasts}{contrasts applied to terms in location and scale
    models.
  }
  \item{na.action}{the function used to filter missing data.
  }
}
\author{Rune Haubo B Christensen}
\references{
  Agresti, A. (2002) \emph{Categorical Data Analysis.} Second edition.  Wiley.

  Aranda-Ordaz, F. J. (1983) An Extension of the Proportional-Hazards
  Model for Grouped Data. \emph{Biometrics}, 39, 109-117.

  Genter, F. C. and Farewell, V. T. (1985) Goodness-of-link testing in
  ordinal regression models. \emph{The Canadian Journal of Statistics},
  13(1), 37-44.

  Christensen, R. H. B., Cleaver, G. and Brockhoff, P. B. (2011)
  Statistical and Thurstonian models for the A-not A protocol with and
  without sureness. \emph{Food Quality and Preference, 22},
  pp. 542-549.
}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))

## A tabular data set:
(tab26 <- with(soup, table("Product" = PROD, "Response" = SURENESS)))
dimnames(tab26)[[2]] <- c("Sure", "Not Sure", "Guess", "Guess", "Not Sure", "Sure")
dat26 <- expand.grid(sureness = as.factor(1:6), prod = c("Ref", "Test"))
dat26$wghts <- c(t(tab26))

m1 <- clm2(sureness ~ prod, scale = ~prod, data = dat26,
          weights = wghts, link = "logistic")

## print, summary, vcov, logLik, AIC:
m1
summary(m1)
vcov(m1)
logLik(m1)
AIC(m1)
coef(m1)
coef(summary(m1))

## link functions:
m2 <- update(m1, link = "probit")
m3 <- update(m1, link = "cloglog")
m4 <- update(m1, link = "loglog")
m5 <- update(m1, link = "cauchit", start = coef(m1))
m6 <- update(m1, link = "Aranda-Ordaz", lambda = 1)
m7 <- update(m1, link = "Aranda-Ordaz")
m8 <- update(m1, link = "log-gamma", lambda = 1)
m9 <- update(m1, link = "log-gamma")

## nominal effects:
mN1 <-  clm2(sureness ~ 1, nominal = ~ prod, data = dat26,
            weights = wghts, link = "logistic")
anova(m1, mN1)

## optimizer / method:
update(m1, scale = ~ 1, method = "Newton")
update(m1, scale = ~ 1, method = "nlminb")
update(m1, scale = ~ 1, method = "optim")
\dontshow{
update(m1, scale = ~ 1, method = "model.frame")
update(m1, location = ~.-prod, scale = ~ 1,
       nominal = ~ prod, method = "model.frame")
}

## threshold functions
mT1 <- update(m1, threshold = "symmetric")
mT2 <- update(m1, threshold = "equidistant")
anova(m1, mT1, mT2)

## Extend example from polr in package MASS:
## Fit model from polr example:
if(require(MASS)) {
    fm1 <- clm2(Sat ~ Infl + Type + Cont, weights = Freq, data = housing)
    fm1
    summary(fm1)
    ## With probit link:
    summary(update(fm1, link = "probit"))
    ## Allow scale to depend on Cont-variable
    summary(fm2 <- update(fm1, scale =~ Cont))
    anova(fm1, fm2)
    ## which seems to improve the fit
}

#################################
## It is possible to fit multinomial models (i.e. with nominal
## effects) as the following example shows:
if(require(nnet)) {
    (hous1.mu <- multinom(Sat ~ 1, weights = Freq, data = housing))
    (hous1.clm <- clm2(Sat ~ 1, weights = Freq, data = housing))

    ## It is the same likelihood:
    all.equal(logLik(hous1.mu), logLik(hous1.clm))

    ## and the same fitted values:
    fitHous.mu <-
        t(fitted(hous1.mu))[t(col(fitted(hous1.mu)) == unclass(housing$Sat))]
    all.equal(fitted(hous1.clm), fitHous.mu)

    ## The coefficients of multinom can be retrieved from the clm2-object
    ## by:
    Pi <- diff(c(0, plogis(hous1.clm$xi), 1))
    log(Pi[2:3]/Pi[1])

    ## A larger model with explanatory variables:
    (hous.mu <- multinom(Sat ~ Infl + Type + Cont, weights = Freq, data = housing))
    (hous.clm <- clm2(Sat ~ 1, nominal = ~ Infl + Type + Cont, weights = Freq,
                      data = housing))

    ## Almost the same likelihood:
    all.equal(logLik(hous.mu), logLik(hous.clm))

    ## And almost the same fitted values:
    fitHous.mu <-
        t(fitted(hous.mu))[t(col(fitted(hous.mu)) == unclass(housing$Sat))]
    all.equal(fitted(hous.clm), fitHous.mu)
    all.equal(round(fitted(hous.clm), 5), round(fitHous.mu), 5)
}

}
\keyword{models}
