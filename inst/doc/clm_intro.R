### R code from vignette source 'clm_intro.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Initialize
###################################################

## Load common packages, functions and set settings:
## library(sensR)
library(ordinal)
library(xtable)
##
RUN <- FALSE    #redo computations and write .RData files
## Change options:
op <- options() ## To be able to reset settings
options("digits" = 7)
options(help_type = "html")
options("width" = 75)
options("SweaveHooks" = list(fig=function()
        par(mar=c(4,4,.5,0)+.5)))
options(continue=" ")



###################################################
### code chunk number 2: clm_intro.Rnw:243-253
###################################################
## data(wine)
tab <- with(wine, table(temp:contact, rating))
mat <- cbind(rep(c("cold", "warm"), each = 2),
             rep(c("no", "yes"), 2),
             tab)
colnames(mat) <- c("Temperature", "Contact",
                   paste("~~", 1:5, sep = ""))
xtab <- xtable(mat)
print(xtab, only.contents=TRUE, include.rownames=FALSE,
      sanitize.text.function = function(x) x)


###################################################
### code chunk number 3: comLogitModel
###################################################
getOption("SweaveHooks")[["fig"]]()
x <- seq(-4, 4, len = 100)
set.seed(12345)
theta <- (-1:1) * 1.5  + runif(3, min = -.25, max = .25)
plot(range(x), c(0,1), type = "n", xlab = "", ylab = "",
     axes = FALSE)
for(i in theta) lines(x, pnorm(-x, mean = i))
## lines(x, rep(0, length(x)))
## lines(x, rep(1, length(x)))
axis(1, labels = FALSE, lwd.ticks = 0)
axis(2, at = c(0, 1), las = 1)
mtext(expression(paste(bold(x)[i]^{T}, bold(beta))), side = 1, line = 1)
mtext(expression(gamma[ij] == P(Y[i] <= j)), side = 2, line = 2, las = 0)
## text(-3, .05, expression(j == 0))
text(-1.3, .2, expression(j == 1))
text(.3, .5, expression(j == 2))
text(1.2, .8, expression(j == 3))
## text(3, .95, expression(j == 4))


###################################################
### code chunk number 4: clm_intro.Rnw:353-370
###################################################
wine$rate <- as.numeric(wine$rating)

## The cumulative link model can be understood as combining the
## information from 1,...,J-1 binomial GLMs:
## Consider the estimates from logistic regression models and compare
## with estimates from the glm:

## coefficient estimates:
cf <- sapply(2:5, function(j) coef(glm(rate < j ~ contact, data=wine,
  family=binomial)))
round(cf, 3)
mean(cf[2,]) ## -1.238937
## standard errors:
ses <- sapply(2:5, function(j)
              coef(summary(glm(rate < j ~ contact, data=wine,
                               family=binomial)))[,2])
fm1 <- clm(rating ~ contact, data=wine)


###################################################
### code chunk number 5: clm_intro.Rnw:384-399
###################################################
## Table of coefficients comparing GLM and CLM
## j, GLM.intercept, GLM.contact, CLM.intercept, CLM.contact
dig <- 3
cf.se <- paste(format(round(cf[2,], dig-1), digits=dig), "(",
               format(round(ses[2,], dig-1), digits=dig), ")", sep="")
mat <- cbind(1:4, format(round(cf[1,], dig-1), digits=dig), cf.se,
             format(round(coef(fm1)[1:4], dig-1), digits=dig),
             c(paste(round(coef(fm1)[5], dig-1), "(",
                     round(coef(summary(fm1))[5,2], dig-1), ")",
                     sep=""), rep("", 3)))
## mat
xtab <- xtable(mat)
print(xtab, only.contents=TRUE, include.rownames=FALSE,
      include.colnames=FALSE)
#      sanitize.text.function = function(x) x)


###################################################
### code chunk number 6: clm_intro.Rnw:486-487
###################################################
str(args(clm))


###################################################
### code chunk number 7: clm_intro.Rnw:518-520
###################################################
fm1 <- clm(rating ~ contact + temp, data = wine)
summary(fm1)


###################################################
### code chunk number 8: clm_intro.Rnw:544-545
###################################################
drop1(fm1, test = "Chi")


###################################################
### code chunk number 9: clm_intro.Rnw:552-554
###################################################
fm0 <- clm(rating ~ 1, data = wine)
add1(fm0, scope = ~ contact + temp, test = "Chi")


###################################################
### code chunk number 10: clm_intro.Rnw:563-564
###################################################
confint(fm1, type = "Wald")


###################################################
### code chunk number 11: oddsRatios
###################################################
round(exp(fm1$beta), 1)


###################################################
### code chunk number 12: oddsRatioCI
###################################################
round(exp(confint(fm1, type = "Wald")), 1)


###################################################
### code chunk number 13: clm_intro.Rnw:731-745
###################################################
## freq <- c(6.5, 8.2, 11.3, 23.5, 15.6, 12.7, 22.2,
##           4.3, 6, 7.7, 13.2, 10.5, 16.3, 42.1)
## year <- factor(rep(c("1960", "1970"), each = 7))
## income <- c(0, 3, 5, 7, 10, 12, 15)
## income <- paste(income, c(rep("-", 6), "+"), c(income[-1], ""),
##                 sep = "")
## income <- data.frame(year=year, freq=freq,
##                      income=factor(rep(income, 2), ordered=TRUE,
##                        levels=income))
data(income)
tab <- xtabs(pct ~ year + income, income)
attr(tab, "class") <- NULL
attr(tab, "call") <- NULL
print(xtable(as.data.frame(tab)), only.contents = TRUE)


###################################################
### code chunk number 14: clm_intro.Rnw:775-778
###################################################
links <- c("logit", "probit", "cloglog", "loglog", "cauchit")
sapply(links, function(link) {
  clm(income ~ year, data=income, weights=pct, link=link)$logLik })


###################################################
### code chunk number 15: clm_intro.Rnw:932-934
###################################################
fm2 <- clm(rating ~ contact * temp, data = wine)
anova(fm1, fm2)


###################################################
### code chunk number 16: clm_intro.Rnw:1007-1011
###################################################
tab <- with(wine, table(temp:contact, rating))
## Get full log-likelihood:
pi.hat <- tab / rowSums(tab)
(ll.full <- sum(tab * ifelse(pi.hat > 0, log(pi.hat), 0))) ## -84.01558


###################################################
### code chunk number 17: clm_intro.Rnw:1014-1019
###################################################
## fit null-model:
fm0 <- clm(rating ~ 1, data = wine)
ll.null <- fm0$logLik
## The null or total deviance:
(Deviance <- -2 * (ll.null - ll.full)) ## 39.407


###################################################
### code chunk number 18: exaDevianceTable
###################################################
nr <- nrow(tab)
nc <- ncol(tab)
df.full <- (nr - 1) * (nc - 1) ## 12
pchisq(Deviance, df=df.full, lower.tail = FALSE)
## something is going on in the data...

## fit treatment model:
fm1 <- clm(rating ~ temp * contact, data = wine)
## residual deviance:
dev.resid <- -2 * (fm1$logLik - ll.full) ## 4.8012
pchisq(dev.resid, df=9, lower = FALSE)
## Test for treatments:
anova(fm0, fm1) ## 34.606
## Test of interaction:
fm2 <- clm(rating ~ temp + contact, data = wine)
anova(fm1, fm2)
## tests of main effects:
drop1(fm2, test = "Chi")


###################################################
### code chunk number 19: clm_intro.Rnw:1180-1195
###################################################
data(wine)
tab <- with(wine, {
  tcb <- temp:contact:bottle
  tcb <- tcb[drop = TRUE]
  table(tcb, rating)
})
mat <- cbind(rep(c("cold", "warm"), each = 4),
             rep(rep(c("no", "yes"), each=2), 2),
             1:8, tab)
colnames(mat) <-
  c("Temperature", "Contact", "Bottle",
    paste("~~", 1:5, sep = ""))
xtab <- xtable(mat)
print(xtab, only.contents=TRUE, include.rownames=FALSE,
      sanitize.text.function = function(x) x)


###################################################
### code chunk number 20: wineTableData (eval = FALSE)
###################################################
## data(wine)
## with(wine, {
##   tcb <- temp:contact:bottle
##   tcb <- tcb[drop=TRUE]
##   table(tcb, rating)
## })


###################################################
### code chunk number 21: latentCLM
###################################################
getOption("SweaveHooks")[["fig"]]()
## x <- seq(-4, 4, len = 100)
plot(x, dnorm(x), type = "l", xlab = "", ylab = "", axes = FALSE)
axis(1, labels = FALSE, lwd.ticks = 0)
segments(theta, 0, theta, .4, lty = 2)
segments(theta, 0, theta, dnorm(theta))
mtext(c(expression(theta[1]), expression(theta[2]),
        expression(theta[3])), side = 3, line = 0,
      at = theta)
mtext(c(expression(pi[i1]), expression(pi[i2]), expression(pi[i3]),
        expression(pi[i4])), at = c(-1.9, -.5, 1, 2.1), line = -1.8, side = 1)
mtext(expression(paste(bold(x)[i]^{T}, bold(beta))), side = 1,
      line = 1)


###################################################
### code chunk number 22: clm_intro.Rnw:1326-1331
###################################################
getOption("SweaveHooks")[["fig"]]()
x <- seq(-6, 6, len = 1e3)
plot(x, dlogis(x), type = "l", ylab = "Density", xlab = "",
     axes=FALSE)
axis(1); axis(2)
lines(x, dnorm(x, sd = pi/sqrt(3)), lty=2)


###################################################
### code chunk number 23: clm_intro.Rnw:1333-1337
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(x, plogis(x), type = "l", ylab="Distribution", xlab="",
     axes=FALSE)
axis(1); axis(2)
lines(x, pnorm(x, sd = pi/sqrt(3)), lty = 2)


###################################################
### code chunk number 24: clm_intro.Rnw:1367-1371
###################################################
fm1 <- clm(rating ~ contact + temp, data = wine, link = "logit")
fm2 <- clm(rating ~ contact + temp, data = wine, link = "probit")
structure(rbind(coef(fm1), coef(fm2)),
          dimnames=list(c("logit", "probit"), names(coef(fm1))))


###################################################
### code chunk number 25: clm_intro.Rnw:1374-1375
###################################################
coef(fm1) / (pi / sqrt(3))


###################################################
### code chunk number 26: clm_intro.Rnw:1404-1406
###################################################
coef(clm(rating ~ temp, data = wine, link = "probit"))["tempwarm"]
coef(clm(rating ~ temp + contact, data = wine, link = "probit"))["tempwarm"]


###################################################
### code chunk number 27: clm_intro.Rnw:1412-1414
###################################################
coef(lm(as.numeric(rating) ~ temp, data = wine))["tempwarm"]
coef(lm(as.numeric(rating) ~ contact + temp, data = wine))["tempwarm"]


###################################################
### code chunk number 28: clm_intro.Rnw:1499-1502
###################################################
lm1 <- lm(as.numeric(rating) ~ contact + temp, data =wine)
sd.lm1 <- summary(lm1)$sigma
coef(lm1)[-1] / sd.lm1


###################################################
### code chunk number 29: clm_intro.Rnw:1506-1508
###################################################
fm1 <- clm(rating ~ contact + temp, data = wine, link = "probit")
coef(fm1)[-(1:4)]


###################################################
### code chunk number 30: clm_intro.Rnw:1514-1515
###################################################
diff(coef(fm1)[1:4])


###################################################
### code chunk number 31: clm_intro.Rnw:1522-1523
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(table(as.numeric(wine$rating)), ylab="", las = 1)


###################################################
### code chunk number 32: clm_intro.Rnw:1667-1669
###################################################
fm1 <- clm(rating ~ temp + contact, data=wine)
summary(fm1)


###################################################
### code chunk number 33: clm_intro.Rnw:1672-1673
###################################################
diff(fm1$alpha)


###################################################
### code chunk number 34: clm_intro.Rnw:1680-1682
###################################################
fm2 <- clm(rating ~ temp + contact, data=wine, threshold="equidistant")
summary(fm2)


###################################################
### code chunk number 35: clm_intro.Rnw:1685-1687
###################################################
a <- round(fm2$alpha[1], 3)
b <- round(fm2$alpha[2], 3)


###################################################
### code chunk number 36: clm_intro.Rnw:1693-1694
###################################################
anova(fm1, fm2)


###################################################
### code chunk number 37: slice1
###################################################
fm1 <- clm(rating ~ temp + contact, data=wine)
slice.fm1 <- slice(fm1, lambda = 5)
par(mfrow = c(2, 3))
plot(slice.fm1)


###################################################
### code chunk number 38: slice11
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 1)


###################################################
### code chunk number 39: slice12
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 2)


###################################################
### code chunk number 40: slice13
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 3)


###################################################
### code chunk number 41: slice14
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 4)


###################################################
### code chunk number 42: slice15
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 5)


###################################################
### code chunk number 43: slice16
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 6)


###################################################
### code chunk number 44: slice2
###################################################
slice2.fm1 <- slice(fm1, lambda = 1e-5)
par(mfrow = c(2, 3))
plot(slice2.fm1)


###################################################
### code chunk number 45: slice21
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 1)


###################################################
### code chunk number 46: slice22
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 2)


###################################################
### code chunk number 47: slice23
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 3)


###################################################
### code chunk number 48: slice24
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 4)


###################################################
### code chunk number 49: slice25
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 5)


###################################################
### code chunk number 50: slice26
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 6)


###################################################
### code chunk number 51: clm_intro.Rnw:1962-1964
###################################################
fm1 <- clm(rating ~ temp + contact, data=wine)
confint(fm1)


###################################################
### code chunk number 52: clm_intro.Rnw:1970-1971
###################################################
confint(fm1, type = "Wald")


###################################################
### code chunk number 53: clm_intro.Rnw:1985-1986
###################################################
str(args(ordinal:::plot.profile.clm))


###################################################
### code chunk number 54: rootStatistic
###################################################
getOption("SweaveHooks")[["fig"]]()
pr1 <- profile(fm1, which.beta="tempwarm")
plot(pr1, root=TRUE)


###################################################
### code chunk number 55: plotRootStatistic
###################################################
getOption("SweaveHooks")[["fig"]]()
pr1 <- profile(fm1, which.beta="tempwarm")
plot(pr1, root=TRUE)


###################################################
### code chunk number 56: profileLikelihood
###################################################
pr1 <- profile(fm1, alpha=1e-4)
plot(pr1)


###################################################
### code chunk number 57: prof1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=1)


###################################################
### code chunk number 58: prof2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=2)


###################################################
### code chunk number 59: clm_intro.Rnw:2110-2115
###################################################
tab <- with(wine, table(contact, rating))
dat <- data.frame(freq =c(tab),
                  contact=rep(c("no", "yes"), 5),
                  rating = factor(rep(1:5, each=2), ordered=TRUE))
dat


###################################################
### code chunk number 60: clm_intro.Rnw:2118-2119 (eval = FALSE)
###################################################
## fm1 <- clm(rating ~ contact, weights=freq)


###################################################
### code chunk number 61: clm_intro.Rnw:2122-2135
###################################################
thresholds <- 1:4
cum.rate <- as.vector(sapply(thresholds, function(x) dat$rating <= x))
rating.factor <- gl(n=length(thresholds), k=nrow(dat),
                    length=nrow(dat) * length(thresholds))
thres.X <- model.matrix(~ rating.factor - 1)
colnames(thres.X) <- paste("t", thresholds, sep="")
old.X <- -model.matrix(~contact, dat)[, -1, drop=FALSE]
new.X <- kronecker(matrix(rep(1, length(thresholds)), nc = 1), old.X)
weights <- kronecker(matrix(rep(1, length(thresholds)), nc = 1), dat$freq)
new.X <- cbind(thres.X, new.X)
colnames(new.X)[-seq(length(thresholds))] <- colnames(old.X)
p.df <- cbind(cum.rate = 1*cum.rate, as.data.frame(new.X), weights)
p.df


###################################################
### code chunk number 62: clm_intro.Rnw:2143-2146
###################################################
glm1 <- glm(cum.rate ~ t1+t2 +t3 +t4 - 1 + contactyes,
            weights=weights, family=binomial, data=p.df)
summary(glm1)


###################################################
### code chunk number 63: misc (eval = FALSE)
###################################################
## 


