\name{wine}
\alias{wine}
\title{
  Bitterness of wine
}
\description{
  The \code{wine} data set is adopted from Randall(1989) and from a
  factorial experiment on factors determining the bitterness of
  wine. Two treatment factors (temperature and contact) each have two
  levels. Temperature and contact between juice and skins can be
  controlled when cruching grapes during wine production. Nine judges
  each assessed wine from two bottles from each of the four treatment
  conditions, hence there are 72 observations in all.
}
\usage{
wine
}
\format{
\describe{
  \item{\code{response}}{
    scorings of wine bitterness on a 0---100 continuous scale.
  }
  \item{\code{rating}}{
    ordered factor with 5 levels; a grouped version of \code{response}.
  }
  \item{\code{temp}}{
    temperature: factor with two levels.
  }
  \item{\code{contact}}{
    factor with two levels (\code{"no"} and \code{"yes"}).
  }
  \item{\code{bottle}}{
    factor with eight levels.
  }
  \item{\code{judge}}{
    factor with nine levels.
  }
}}
\source{
  Data are adopted from Randall (1989).
}
\references{
  Randall, J (1989). The analysis of sensory data by generalised linear
  model. \emph{Biometrical journal 7}, pp. 781--793.

  Tutz, G. and W. Hennevogl (1996). Random effects in ordinal regression
  models. \emph{Computational Statistics & Data Analysis 22},
  pp. 537--557.
}
\examples{

head(wine)
str(wine)

## Variables 'rating' and 'response' are related in the following way:
(intervals <- seq(0,100, by = 20))
all(wine$rating == findInterval(wine$response, intervals)) ## ok

## A few illustrative tabulations:
## Table matching Table 5 in Randall (1989):
temp.contact.bottle <- with(wine, temp:contact:bottle)[drop=TRUE]
xtabs(response ~ temp.contact.bottle + judge, data = wine)

## Table matching Table 6 in Randall (1989):
with(wine, {
  tcb <- temp:contact:bottle
  tcb <- tcb[drop=TRUE]
  table(tcb, rating)
})
## or simply: with(wine, table(bottle, rating))

## Table matching Table 1 in Tutz & Hennevogl (1996):
tab <- xtabs(as.numeric(rating) ~ judge + temp.contact.bottle,
             data = wine)
colnames(tab) <-
  paste(rep(c("c","w"), each = 4), rep(c("n", "n", "y", "y"), 2),
        1:8, sep=".")
tab


## A simple model:
m1 <- clm(rating ~ temp * contact, data = wine)
summary(m1)

}

\keyword{datasets}
