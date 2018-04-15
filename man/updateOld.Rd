\name{update.clm2}
\alias{update.clm2}
\alias{update.clmm2}
\title{Update method for cumulative link models}
\description{
  Update method for cumulative link models fitted with \code{clm2}.
  This makes it possible to use e.g.
  \code{update(obj, location = ~ . - var1, scale = ~ . + var2)}
}
\usage{
\method{update}{clm2}(object, formula., location, scale, nominal,...,
        evaluate = TRUE)
\method{update}{clmm2}(object, formula., location, scale, nominal,...,
        evaluate = TRUE)
}
\arguments{
  \item{object}{a \code{\link{clm2}} object.
  }
  \item{formula.}{not used---unfortunately this argument is part of the
    default method.
  }
  \item{location}{an optional new formula for the location; see
    \code{\link{update.formula}} for details.
  }
  \item{scale}{an optional new formula for the scale; see
    \code{\link{update.formula}} for details.
  }
  \item{nominal}{an optional new formula for nominal effects; see
    \code{\link{update.formula}} for details.
  }
  \item{\dots}{additional arguments to the call, or arguments with
    changed values.
  }
  \item{evaluate}{if true evaluate the new call else return the call.
  }
}
\value{
  If \code{evaluate = TRUE} the fitted object is returned,
  otherwise the updated call.
}
\author{Rune Haubo B Christensen}
\examples{
options(contrasts = c("contr.treatment", "contr.poly"))

m1 <-  clm2(SURENESS ~ PROD, scale = ~PROD, data = soup,
           link = "logistic")

m2 <- update(m1, link = "probit")
m3 <- update(m1, link = "cloglog")
m4 <- update(m1, link = "loglog")
anova(m1, update(m1, scale = ~.-PROD))
mT1 <- update(m1, threshold = "symmetric")

}
\keyword{models}
