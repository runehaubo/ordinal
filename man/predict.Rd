\name{predict.clm}
\alias{predict.clm}
\title{Predict Method for CLM fits}
\description{
  Obtains predictions from a cumulative link model.
}
\usage{

\method{predict}{clm}(object, newdata, se.fit = FALSE, interval = FALSE,
           level = 0.95,
           type = c("prob", "class", "cum.prob", "linear.predictor"),
           na.action = na.pass, ...)

}
\arguments{
  \item{object}{a fitted object of class inheriting from
    \code{clm}.}
  \item{newdata}{optionally, a data frame in which to look for variables
    with which to predict. Note that all predictor variables should be
    present having the same names as the variables used to fit the
    model. If the response variable is present in \code{newdata}
    predictions are obtained for the levels of the response as given by
    \code{newdata}. If the response variable is omitted from
    \code{newdata} predictions are obtained for all levels of the
    response variable for each of the rows of \code{newdata}.
  }
  \item{se.fit}{should standard errors of the predictions be provided?
    Not applicable and ignored when \code{type = "class"}.
  }
  \item{interval}{should confidence intervals for the predictions be
    provided?  Not applicable and ignored when \code{type = "class"}.
  }
  \item{level}{the confidence level.
  }
  \item{type}{the type of predictions. \code{"prob"} gives
    probabilities, \code{"class"} gives predicted response class
    membership defined as highest probability prediction,
    \code{"cum.prob"} gives cumulative probabilities (see details)
    and \code{"linear.predictor"} gives predictions on the scale of the
    linear predictor including the boundary categories.
  }
  \item{na.action}{function determining what should be done with missing
    values in \code{newdata}. The default is to predict \code{NA}.
  }
  \item{\dots}{further arguments passed to or from other methods.
  }
}
\details{

  If \code{newdata} is omitted and \code{type = "prob"} a vector of
  fitted probabilities are returned identical to the result from
  \code{fitted}.

  If \code{newdata} is supplied and the response
  variable is omitted, then predictions, standard errors and intervals
  are matrices rather than vectors with the same number of rows as
  \code{newdata} and with one column for each response class. If
  \code{type = "class"} predictions are always a vector.

  If \code{newdata} is omitted, the way missing values in the original fit are handled
  is determined by the \code{na.action} argument of that fit. If
  \code{na.action = na.omit} omitted cases will not appear in the
  residuals, whereas if \code{na.action = na.exclude}
  they will appear (in predictions, standard
  errors or interval limits), with residual value \code{NA}. See also
  \code{\link{napredict}}.

  If \code{type = "cum.prob"} or \code{type = "linear.predictor"} there
  will be two sets of predictions, standard errors and intervals; one
  for j and one for j-1 (in the usual notation) where j = 1, ..., J index
  the response classes.

  If newdata is supplied and the response variable is omitted, then
  \code{predict.clm} returns much the same thing as \code{predict.polr}
  (matrices of predictions). Similarly, if \code{type = "class"}.

  If the fit is rank-deficient, some of the columns of the design matrix
  will have been dropped. Prediction from such a fit only makes sense if
  newdata is contained in the same subspace as the original data. That
  cannot be checked accurately, so a warning is issued
  (cf. \code{\link{predict.lm}}).

}
\value{
  A list containing the following components
  \item{fit}{predictions or fitted values if \code{newdata} is not
    supplied.
  }
  \item{se.fit}{if \code{se.fit=TRUE} standard errors of the predictions
    otherwise \code{NULL}.
  }
  \item{upr, lwr}{if \code{interval=TRUE} lower and upper confidence
    limits.}

}
\author{Rune Haubo B Christensen}
\seealso{
    \code{\link[ordinal]{clm}}, \code{\link[ordinal]{clmm}}.
}
\examples{

## simple model:
fm1 <- clm(rating ~ contact + temp, data=wine)
summary(fm1)

## Fitted values with standard errors and confidence intervals:
predict(fm1, se.fit=TRUE, interval=TRUE) # type="prob"
## class predictions for the observations:
predict(fm1, type="class")

newData <- expand.grid(temp = c("cold", "warm"),
                       contact = c("no", "yes"))

## Predicted probabilities in all five response categories for each of
## the four cases in newData:
predict(fm1, newdata=newData, type="prob")
## now include standard errors and intervals:
predict(fm1, newdata=newData, se.fit=TRUE, interval=TRUE, type="prob")


}
\keyword{models}
