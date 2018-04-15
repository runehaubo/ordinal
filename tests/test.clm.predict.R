library(ordinal)
## source("test.clm.predict.R")

## library(devtools)
## r2path <- "/Users/rhbc/Documents/Rpackages/ordinal/pkg/ordinal"
## clean_dll(pkg = r2path)
## load_all(r2path)

cy <- with(wine, which(temp == "cold" & contact == "yes"))
options("contrasts" = c("contr.treatment", "contr.poly"))
getOption("contrasts")

## Example model

wine1.clm <- clm(rating ~ temp*contact, subset = -cy, data = wine)
summary(wine1.clm)
names(wine1.clm)

wine.clm <- clm(rating~temp*contact, data=wine)
summary(wine.clm)
names(wine.clm)
## Make sure the same elements are present with a rank deficient model
## fit:
stopifnot(all(names(wine1.clm) == names(wine.clm)))

## With treatment contrasts:
options("contrasts" = c("contr.treatment", "contr.poly"))
wine.clm <- clm(rating~temp*contact, data=wine)
coef(summary(wine.clm))
head(model.matrix(wine.clm)$X)
wine.clm$contrasts
head(pred1 <- predict(wine.clm)$fit)

## With sum contrasts:
options("contrasts" = c("contr.sum", "contr.poly"))
wine.clm <- clm(rating~temp*contact, data=wine)
coef(summary(wine.clm))
head(model.matrix(wine.clm)$X)
wine.clm$contrasts
head(pred2 <- predict(wine.clm)$fit)

## Mixture of sum and treatment contrasts:
options("contrasts" = c("contr.treatment", "contr.poly"))
wine.clm <- clm(rating~temp*contact, data=wine,
                contrasts=list(temp="contr.sum"))
coef(summary(wine.clm))
head(model.matrix(wine.clm)$X)
wine.clm$contrasts
head(pred3 <- predict(wine.clm)$fit)

stopifnot(isTRUE(all.equal(pred1, pred2)))
stopifnot(isTRUE(all.equal(pred1, pred3)))

#################################
### Now for a rank deficient fit:
#################################

cy <- with(wine, which(temp == "cold" & contact == "yes"))
options("contrasts" = c("contr.treatment", "contr.poly"))
wine1.clm <- clm(rating ~ temp*contact, subset = -cy, data = wine)
coef(summary(wine1.clm))
attributes(model.matrix(wine1.clm)$X)$contrasts
wine1.clm$contrasts
head(pred4 <- predict(wine1.clm)$fit)

options("contrasts" = c("contr.sum", "contr.poly"))
wine1.clm <- clm(rating ~ temp*contact, subset = -cy, data = wine)
attributes(model.matrix(wine1.clm)$X)$contrasts
options("contrasts" = c("contr.treatment", "contr.poly"))
attributes(model.matrix(wine1.clm)$X)$contrasts
## Notice that the contrasts change in the attributes of the fit!!!
coef(summary(wine1.clm))
wine1.clm$contrasts
head(pred5 <- predict(wine1.clm)$fit)

head(cbind(pred4, pred5))
stopifnot(isTRUE(all.equal(pred4, pred5)))

options("contrasts" = c("contr.treatment", "contr.poly"))
wine1.clm <- clm(rating ~ temp*contact, subset = -cy, data = wine,
                 contrasts=list(temp="contr.sum"))
coef(summary(wine1.clm))
head(model.matrix(wine1.clm)$X)
attributes(model.matrix(wine1.clm)$X)$contrasts
wine1.clm$contrasts
head(pred6 <- predict(wine1.clm)$fit)

head(cbind(pred4, pred5, pred6))
stopifnot(isTRUE(all.equal(pred4, pred6)))
##################################################################

## Compare equality of fitted values for models with different contrasts:
options("contrasts" = c("contr.treatment", "contr.poly"))
fm1 <- clm(rating ~ temp + contact, data=wine)
fitted(fm1)
options("contrasts" = c("contr.sum", "contr.poly"))
fm2 <- clm(rating ~ temp + contact, data=wine)
fitted(fm2)
options("contrasts" = c("contr.treatment", "contr.poly"))
fm3 <- clm(rating ~ temp + contact, data=wine,
           contrasts=list(contact="contr.sum"))
fitted(fm3)
stopifnot(isTRUE(all.equal(fitted(fm1), fitted(fm2))))
stopifnot(isTRUE(all.equal(fitted(fm1), fitted(fm3))))

##################################################################
## Compare equality of fitted values for models with different
## contrasts in face of aliased coefficients:
options("contrasts" = c("contr.treatment", "contr.poly"))
cy <- with(wine, which(temp == "cold" & contact == "yes"))
Wine <- subset(wine, subset=!(temp == "cold" & contact == "yes"))
fm1 <- clm(rating ~ temp + contact, data=Wine)
options("contrasts" = c("contr.sum", "contr.poly"))
fm2 <- clm(rating ~ temp + contact, data=Wine)
options("contrasts" = c("contr.treatment", "contr.poly"))
fm3 <- clm(rating ~ temp + contact, data=Wine,
           contrasts=list(contact="contr.sum"))

stopifnot(isTRUE(all.equal(fitted(fm1), fitted(fm2))))
stopifnot(isTRUE(all.equal(fitted(fm1), fitted(fm3))))
stopifnot(isTRUE(all.equal(predict(fm1)$fit, predict(fm2)$fit)))
stopifnot(isTRUE(all.equal(predict(fm1)$fit, predict(fm3)$fit)))

#################################
## Does this also happen if the wine data has changed?
options("contrasts" = c("contr.treatment", "contr.poly"))
Wine <- subset(wine, subset=!(temp == "cold" & contact == "yes"))
fm1 <- clm(rating ~ temp + contact, data=Wine)
fit1 <- fitted(fm1)
pred1 <- predict(fm1)$fit
Wine <- wine
pred2 <- predict(fm1)$fit
stopifnot(isTRUE(all.equal(fit1, pred1)))
stopifnot(isTRUE(all.equal(fit1, pred2)))

## What if weights, say, is an expression?
## Notice that updating the model object changes it:
set.seed(123)
fm1 <- clm(rating ~ temp + contact, data=wine,
           weights=runif(nrow(wine), .5, 1.5))
fm2 <- update(fm1)
stopifnot(isTRUE(all.equal(fitted(fm1), predict(fm1)$fit)))
stopifnot(!isTRUE(all.equal(fitted(fm1), fitted(fm2))))

#################################
## Test equality of fits and predictions of models with:
##   'x + I(x^2)' and 'poly(x, 2)':
## December 25th 2014, RHBC.
data(wine)
set.seed(1)
x <- rnorm(nrow(wine), sd=2) + as.numeric(wine$rating)
range(x)

## Comparison of 'x + I(x^2)' and 'poly(x, 2)':
fm3 <- clm(rating ~ temp + x + I(x^2), data=wine)
fm4 <- clm(rating ~ temp + poly(x, 2), data=wine)
## Same model fits, but different parameterizations:
stopifnot(
    !isTRUE(all.equal(coef(fm3), coef(fm4), check.names=FALSE))
    )
stopifnot(isTRUE(all.equal(logLik(fm3), logLik(fm4))))
newData <- expand.grid(temp = levels(wine$temp),
                       x=seq(-1, 7, 3))
predict(fm3, newdata=newData)$fit
predict(fm4, newdata=newData)$fit
stopifnot(isTRUE(all.equal(fitted(fm3), fitted(fm4))))
stopifnot(isTRUE(
    all.equal(predict(fm3, newdata=newData)$fit,
              predict(fm4, newdata=newData)$fit)))
#################################
