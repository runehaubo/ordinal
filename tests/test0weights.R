library(ordinal)
options(contrasts = c("contr.treatment", "contr.poly"))
## library(devtools)
## r2path <- "/Users/rhbc/Documents/Rpackages/ordinal/pkg/ordinal"
## clean_dll(pkg = r2path)
## load_all(r2path)

## one zero weight:
data(wine, package="ordinal")
wts <- rep(1, nrow(wine))
wine$rating
wts[1] <- 0
fm1 <- clm(rating ~ contact + temp, data=wine, weights=wts)
fm1
fm1$n ## 72
fm1$nobs ## 71
confint(fm1)
plot(profile(fm1))
plot(slice(fm1), 5)
convergence(fm1)
drop1(fm1, test="Chi")
add1(fm1, scope=~.^2, test="Chi")
## clm_anova(fm1)
pred <- predict(fm1, newdata=wine) ## OK
step.fm1 <- step(fm1, trace=0)
fitted(fm1)
dim(model.matrix(fm1)$X)
dim(model.matrix(fm1, "B")$B1)
mf <- update(fm1, method="model.frame")
str(mf)
wts <- mf$wts
dim(model.matrix(fm1)$X[wts > 0, , drop=FALSE])

fm1b <- clm(rating ~ temp, scale=~contact, data=wine, weights=wts)
summary(fm1b)
pr <- profile(fm1b)
confint(pr)
plot(pr, 1)
fm1c <- clm(rating ~ temp, nominal=~contact, data=wine, weights=wts)
summary(fm1c)
pr <- profile(fm1c)
confint(pr)
plot(pr, 1)

## nominal.test(fm1)
## scale.test(fm1)

## zero out an entire response category:
wts2 <- 1 * with(wine, rating != "2")
fm2 <- clm(rating ~ contact + temp, data=wine, weights=wts2)
fm2
fm2$n ## 72
fm2$nobs ## 50
## Dimension of X and B1, B2 differ:
dim(model.matrix(fm2)$X)
dim(model.matrix(fm2, "B")$B1)
## Cannot directly evaluate predictions on the original data:
try(predict(fm2, newdata=wine), silent=TRUE)[1]
confint(fm2)
profile(fm2)
plot(slice(fm2), 5)
step.fm2 <- step(fm2, trace=0)
fitted(fm2)
## Scale and nominal effects:
fm2b <- clm(rating ~ temp, scale=~contact, data=wine, weights=wts2)
summary(fm2b)
pr <- profile(fm2b)
confint(pr)
plot(pr, 1)
fm2c <- clm(rating ~ temp, nominal=~contact, data=wine, weights=wts2)
summary(fm2c)
pr <- profile(fm2c)
confint(pr)
plot(pr, 1)
pred <- predict(fm2c, newdata=wine[!names(wine) %in% "rating"])
pred <- predict(fm2b, newdata=wine[!names(wine) %in% "rating"])

## nominal.test(fm2)
## scale.test(fm2)

## Different data sets (error):
try(anova(fm1, fm2), silent=TRUE)[1] ## OK

## Test clm.fit:
wts2 <- 1 * with(wine, rating != "2")
mf2 <- update(fm2, method="design")
fm3 <- with(mf2, clm.fit(y, X, weights=wts))

#################################
