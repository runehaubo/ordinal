# test.clm.flex.link.R

library(ordinal)

fm <- clm(rating ~ contact + temp, data=wine, link="log-gamma")
fm
summary(fm)
vcov(fm)
logLik(fm)
extractAIC(fm)
fm2 <- update(fm, link="probit")
anova(fm, fm2)
head(model.matrix(fm)$X)
head(model.frame(fm))
coef(fm)
coef(summary(fm))
nobs(fm)
terms(fm)
# profile(fm) # not implemented
confint(fm)

predict(fm, se=TRUE, interval = TRUE)
predict(fm, type="class")
newData <- expand.grid(temp = c("cold", "warm"),
                       contact = c("no", "yes"))

## Predicted probabilities in all five response categories for each of
## the four cases in newData:
predict(fm, newdata=newData, type="prob")
predict(fm, newdata=newData, type="class")

predict(fm, newdata=newData, type="prob", se.fit = TRUE, interval = TRUE)


## Aranda-Ordaz link:
fm <- clm(rating ~ contact + temp, data=wine, link="Aranda-Ordaz")
fm
summary(fm)
vcov(fm)
logLik(fm)
extractAIC(fm)
fm2 <- update(fm, link="logit")
anova(fm, fm2)
head(model.matrix(fm)$X)
head(model.frame(fm))
coef(fm)
coef(summary(fm))
nobs(fm)
terms(fm)
# profile(fm) # not implemented
confint(fm)

predict(fm, se=TRUE, interval = TRUE)
predict(fm, type="class")
newData <- expand.grid(temp = c("cold", "warm"),
                       contact = c("no", "yes"))

## Predicted probabilities in all five response categories for each of
## the four cases in newData:
predict(fm, newdata=newData, type="prob")
predict(fm, newdata=newData, type="class")

predict(fm, newdata=newData, type="prob", se.fit = TRUE, interval = TRUE)
