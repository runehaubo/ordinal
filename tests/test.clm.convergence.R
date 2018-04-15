library(ordinal)


## Testing that errors in chol() are caught soon enough:
cy <- with(wine, which(temp == "cold" & contact == "yes"))
wine2 <- subset(wine, subset=(!1:nrow(wine) %in% cy))
wine2[c(9, 15, 46), "rating"] <- NA
fm1 <- clm(rating ~ temp, scale=~contact, nominal=~contact,
           data=wine2)
fm1 <- try(clm(rating ~ temp, scale=~contact, nominal=~contact,
               data=wine2, control=list(gradTol=1e-12)), silent=TRUE)
fm2 <- try(clm(rating ~ temp, scale=~contact, nominal=~contact,
               data=wine2, control=list(gradTol=1e-15)), silent=TRUE)
## These gave errors in version 2014.11-12.
stopifnot(!inherits(fm1, "try-error"))
stopifnot(!inherits(fm2, "try-error"))
summary(fm1)
summary(fm2)

## Error in convergence.clm() due to bad evaluation of model
## environment with update(object, doFit=FALSE):
wine3 <- wine
set.seed(1234)
wts <- runif(nrow(wine3), 0, 2)
fm3 <- clm(rating ~ temp + contact, data=wine3,
           weights=wts)
c0 <- convergence(fm3)
set.seed(1234)
fm3 <- clm(rating ~ temp + contact, data=wine3,
           weights=runif(nrow(wine3), 0, 2))
c1 <- convergence(fm3)
c0$info$logLik.Error
c1$info$logLik.Error
all.equal(c0$info$logLik.Error, c1$info$logLik.Error)
## In version 2014.11-14:
## > wine3 <- wine
## > set.seed(1234)
## > wts <- runif(nrow(wine3), 0, 2)
## > fm3 <- clm(rating ~ temp + contact, data=wine3,
## +            weights=wts)
## > c0 <- convergence(fm3)
## > set.seed(1234)
## > fm3 <- clm(rating ~ temp + contact, data=wine3,
## +            weights=runif(nrow(wine3), 0, 2))
## > c1 <- convergence(fm3)
## > c0$info$logLik.Error
## [1] "<1e-10"
## > c1$info$logLik.Error
## [1] "4.80e+00"
## > all.equal(c0$info$logLik.Error, c1$info$logLik.Error)
## [1] "1 string mismatch"
stopifnot(c0$info$logLik.Error ==
          c1$info$logLik.Error)
