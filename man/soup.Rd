\name{soup}
\alias{soup}
\title{
  Discrimination study of packet soup
}
\description{
  The \code{soup} data frame has 1847 rows and 13 variables. 185
  respondents participated in an A-not A discrimination test with
  sureness. Before experimentation the respondents were familiarized
  with the reference product and during experimentation, the respondents
  were asked to rate samples on an ordered scale with six categories
  given by combinations of (reference, not reference) and (sure, not
  sure, guess) from 'referene, sure' = 1 to 'not reference, sure' = 6.
  %given by the levels of the \code{SURENESS} variable.
}
\usage{
soup
}
\format{
\describe{
  \item{\code{RESP}}{
    factor with 185 levels: the respondents in the study.
  }
  \item{\code{PROD}}{
    factor with 2 levels: index reference and test products.
  }
  \item{\code{PRODID}}{
    factor with 6 levels: index reference and the five test product
    variants.
  }
  \item{\code{SURENESS}}{
    ordered factor with 6 levels: the respondents ratings of soup
    samples.
  }
  \item{\code{DAY}}{
    factor with two levels: experimentation was split over two days.
  }
  \item{\code{SOUPTYPE}}{
    factor with three levels: the type of soup regularly consumed by the
    respondent.
  }
  \item{\code{SOUPFREQ}}{
    factor with 3 levels: the frequency with which the respondent
    consumes soup.
  }
  \item{\code{COLD}}{
    factor with two levels: does the respondent have a cold?
  }
  \item{\code{EASY}}{
    factor with ten levels: How easy did the respondent find the
    discrimation test? 1 = difficult, 10 = easy.
  }
  \item{\code{GENDER}}{
    factor with two levels: gender of the respondent.
  }
  \item{\code{AGEGROUP}}{
    factor with four levels: the age of the respondent.
  }
  \item{\code{LOCATION}}{
    factor with three levels: three different locations where
    experimentation took place. 
  }
%%   \item{\code{SEQ}}{
%%     integer vector: the sequence at which experimentation took
%%     place. Numbering restarted at the second day of experimentation.
%%   }
}}
\source{
  Data are produced by Unilever Research. Permission to publish
  the data is granted. 
}
\references{
  Christensen, R. H. B., Cleaver, G. and Brockhoff, P. B.(2011)
  Statistical and Thurstonian models for the A-not A protocol with and
  without sureness. \emph{Food Quality and Preference, 22},
  pp. 542-549. 
}

\keyword{datasets}
