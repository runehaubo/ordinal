# test.sign.R

# Test the use of sign.location and sign.nominal in clm.control():

library(ordinal)

fm1 <- clm(rating ~ temp + contact, data=wine)
fm2 <- clm(rating ~ temp + contact, data=wine,
           sign.location="positive")
# dput(names(fm1))
keep <- c("aliased", "alpha", "cond.H", 
          "contrasts", "convergence", "df.residual", "edf", 
          "fitted.values", "formula", "formulas", "gradient", 
          "info", "link", "logLik", "maxGradient", "message", "model", 
          "n", "niter", "nobs", "start", "terms", "Theta", "threshold", 
          "tJac", "xlevels", "y", "y.levels")
check <- mapply(function(x, y) isTRUE(all.equal(x, y)), fm1[keep], fm2[keep])
stopifnot(all(check))
stopifnot(isTRUE(all.equal(
  fm1$beta, - fm2$beta
)))

fm1 <- clm(rating ~ temp, nominal=~ contact, data=wine)
fm2 <- clm(rating ~ temp, nominal=~ contact, data=wine,
           sign.nominal="negative")
keep <- c("aliased", "beta", "cond.H", 
          "contrasts", "convergence", "df.residual", "edf", 
          "fitted.values", "formula", "formulas", "gradient", 
          "info", "link", "logLik", "maxGradient", "message", "model", 
          "n", "niter", "nobs", "start", "terms", "Theta", "threshold", 
          "tJac", "xlevels", "y", "y.levels")
# check <- mapply(function(x, y) isTRUE(all.equal(x, y)), fm1, fm2)
check <- mapply(function(x, y) isTRUE(all.equal(x, y)), fm1[keep], fm2[keep])
stopifnot(all(check))
stopifnot(isTRUE(all.equal(
  fm1$alpha[5:8], -fm2$alpha[5:8]
)))


fm1 <- clm(rating ~ temp, nominal=~ contact, data=wine)
fm2 <- clm(rating ~ temp, nominal=~ contact, data=wine,
           sign.nominal="negative", sign.location="positive")
keep <- c("aliased", "cond.H", 
          "contrasts", "convergence", "df.residual", "edf", 
          "fitted.values", "formula", "formulas", "gradient", 
          "info", "link", "logLik", "maxGradient", "message", "model", 
          "n", "niter", "nobs", "start", "terms", "Theta", "threshold", 
          "tJac", "xlevels", "y", "y.levels")
# check <- mapply(function(x, y) isTRUE(all.equal(x, y)), fm1, fm2)
check <- mapply(function(x, y) isTRUE(all.equal(x, y)), fm1[keep], fm2[keep])
stopifnot(all(check))
stopifnot(
  isTRUE(all.equal(fm1$alpha[5:8], -fm2$alpha[5:8])),
  isTRUE(all.equal(fm1$beta, -fm2$beta))
)

# Check predict method:
newData <- with(wine, expand.grid(temp=levels(temp), contact=levels(contact)))
(p1 <- predict(fm1, newdata=newData))
(p2 <- predict(fm2, newdata=newData))
stopifnot(isTRUE(all.equal(p1, p2)))

stopifnot(isTRUE(
  all.equal(predict(fm1, newdata=wine, se=TRUE, interval=TRUE),
            predict(fm2, newdata=wine, se=TRUE, interval=TRUE))
))

# Check profile and confint methods:
confint.default(fm1)
confint.default(fm2)

stopifnot(
  isTRUE(all.equal(confint(fm1), -confint(fm2)[, 2:1, drop=FALSE], 
                   check.attributes=FALSE))
)

fm1 <- clm(rating ~ temp + contact, data=wine)
fm2 <- clm(rating ~ temp + contact, data=wine,
           sign.location="positive")
pr1 <- profile(fm1)
pr2 <- profile(fm2)
stopifnot(
  isTRUE(all.equal(confint(fm1), - confint(fm2)[, 2:1], check.attributes=FALSE))
)

