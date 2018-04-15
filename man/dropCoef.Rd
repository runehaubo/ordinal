\name{drop.coef}
\alias{drop.coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Ensure Full Rank Design Matrix
}
\description{
  Coefficients (columns) are dropped from a design matrix to ensure that
  it has full rank.
}
\usage{
drop.coef(X, silent = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    a design matrix, e.g., the result of \code{\link{model.matrix}}
    possibly of less than full column rank, i.e., with redundant
    parameters. Works for \code{ncol(X) >= 0} and \code{nrow(X) >= 0}.
  }
  \item{silent}{
    should a message not be issued if X is column rank deficient?
  }
}
\details{
  Redundant columns of the design matrix are identified with the
  LINPACK implementation of the \code{\link{qr}} decomposition and
  removed. The returned design matrix will have \code{qr(X)$rank}
  columns.
}
\value{
  The design matrix \code{X} without redundant columns.
}
\author{
  Rune Haubo B Christensen
}
\seealso{
  \code{\link{qr}} and \code{\link{lm}}
}
\examples{

X <- model.matrix( ~ PRODID * DAY, data = soup)
ncol(X)
newX <- drop.coef(X)
ncol(newX)

## Essentially this is being computed:
qr.X <- qr(X, tol = 1e-7, LAPACK = FALSE)
newX <- X[, qr.X$pivot[1:qr.X$rank], drop = FALSE]
## is newX of full column rank?
ncol(newX) == qr(newX)$rank
## the number of columns being dropped:
ncol(X) - ncol(newX)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{models}

