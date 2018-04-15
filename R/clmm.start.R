## This file contains:
## Functions to compute starting values for clmm()s.

clmm.start <- function(frames, link, threshold) {
  ## get starting values from clm:
  fit <- with(frames,
              clm.fit(y=y, X=X, weights=wts, offset=off, link=link,
                      threshold=threshold))

  ## initialize variance parameters to zero:
  start <- c(fit$par, rep(0, length(frames$grList)))
  return(start)
}

