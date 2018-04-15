### R code from vignette source 'clm_tutorial.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Initialize
###################################################

## Load common packages, functions and set settings:
library(ordinal)
library(xtable)
##
RUN <- FALSE    #redo computations and write .RData files
## Change options:
op <- options() ## To be able to reset settings
options("digits" = 7)
options(help_type = "html")
## options("width" = 75)
options("SweaveHooks" = list(fig=function()
        par(mar=c(4,4,.5,0)+.5)))
options(continue=" ")



###################################################
### code chunk number 2: clm_tutorial.Rnw:153-156
###################################################
data(wine)
head(wine)
str(wine)


###################################################
### code chunk number 3: clm_tutorial.Rnw:177-191
###################################################
data(wine)
temp.contact.bottle <- with(wine, temp:contact:bottle)[drop=TRUE]
tab <- xtabs(as.numeric(rating) ~ temp.contact.bottle + judge,
             data=wine)
class(tab) <- "matrix"
attr(tab, "call") <- NULL
mat <- cbind(rep(c("cold", "warm"), each = 4),
             rep(rep(c("no", "yes"), each=2), 2),
             1:8, tab)
colnames(mat) <-
  c("Temperature", "Contact", "Bottle", 1:9)
xtab <- xtable(mat)
print(xtab, only.contents=TRUE, include.rownames=FALSE,
      sanitize.text.function = function(x) x)


###################################################
### code chunk number 4: clm_tutorial.Rnw:233-235
###################################################
fm1 <- clm(rating ~ temp + contact, data=wine)
fm1


###################################################
### code chunk number 5: clm_tutorial.Rnw:238-239
###################################################
summary(fm1)


###################################################
### code chunk number 6: clm_tutorial.Rnw:258-259
###################################################
clm.control()$gradTol


###################################################
### code chunk number 7: clm_tutorial.Rnw:278-279
###################################################
exp(coef(fm1)[5])


###################################################
### code chunk number 8: clm_tutorial.Rnw:287-289
###################################################
fm2 <- clm(rating ~ temp, data=wine)
anova(fm2, fm1)


###################################################
### code chunk number 9: clm_tutorial.Rnw:295-296
###################################################
drop1(fm1, test = "Chi")


###################################################
### code chunk number 10: clm_tutorial.Rnw:301-303
###################################################
fm0 <- clm(rating ~ 1, data=wine)
add1(fm0, scope = ~ temp + contact, test = "Chi")


###################################################
### code chunk number 11: clm_tutorial.Rnw:309-310
###################################################
confint(fm1)


###################################################
### code chunk number 12: clm_tutorial.Rnw:316-317
###################################################
confint(fm1, type="Wald")


###################################################
### code chunk number 13: clm_tutorial.Rnw:324-325
###################################################
fm.cll <- clm(rating ~ contact + temp, data=wine, link="cloglog")


###################################################
### code chunk number 14: clm_tutorial.Rnw:356-358
###################################################
fm.nom <- clm(rating ~ temp, nominal=~contact, data=wine)
summary(fm.nom)


###################################################
### code chunk number 15: clm_tutorial.Rnw:382-383
###################################################
anova(fm1, fm.nom)


###################################################
### code chunk number 16: clm_tutorial.Rnw:394-395
###################################################
fm.nom2 <- clm(rating ~ temp + contact, nominal=~contact, data=wine)


###################################################
### code chunk number 17: clm_tutorial.Rnw:398-399
###################################################
summary(fm.nom2)


###################################################
### code chunk number 18: clm_tutorial.Rnw:476-478
###################################################
fm.sca <- clm(rating ~ temp + contact, scale=~temp, data=wine)
summary(fm.sca)


###################################################
### code chunk number 19: clm_tutorial.Rnw:488-489
###################################################
exp(fm.sca$zeta)


###################################################
### code chunk number 20: clm_tutorial.Rnw:507-510
###################################################
fm.equi <- clm(rating ~ temp + contact, data=wine,
               threshold="equidistant")
summary(fm.equi)


###################################################
### code chunk number 21: clm_tutorial.Rnw:517-518
###################################################
c(with(fm.equi, tJac %*% alpha))


###################################################
### code chunk number 22: clm_tutorial.Rnw:523-524
###################################################
mean(diff(fm1$alpha))


###################################################
### code chunk number 23: clm_tutorial.Rnw:530-531
###################################################
anova(fm1, fm.equi)


###################################################
### code chunk number 24: clm_tutorial.Rnw:546-547
###################################################
predict(fm1, type = "class")


###################################################
### code chunk number 25: clm_tutorial.Rnw:554-557
###################################################
newData <- expand.grid(temp=levels(wine$temp),
                       contact=levels(wine$contact))
cbind(newData, predict(fm1, newdata=newData)$fit)


###################################################
### code chunk number 26: clm_tutorial.Rnw:564-565
###################################################
head(do.call("cbind", predict(fm1, se.fit=TRUE, interval=TRUE)))


###################################################
### code chunk number 27: clm_tutorial.Rnw:576-578
###################################################
fm.nom2 <- clm(rating ~ contact, nominal=~temp, data=wine)
summary(fm.nom2)


###################################################
### code chunk number 28: clm_tutorial.Rnw:596-597
###################################################
anova(fm1, fm.nom2)


###################################################
### code chunk number 29: clm_tutorial.Rnw:605-608
###################################################
data(soup)
fm.soup <- clm(SURENESS ~ PRODID * DAY, data=soup)
summary(fm.soup)


###################################################
### code chunk number 30: clm_tutorial.Rnw:613-614
###################################################
with(soup, table(DAY, PRODID))


###################################################
### code chunk number 31: clm_tutorial.Rnw:622-625
###################################################
mm <- model.matrix( ~ PRODID * DAY, data=soup)
ncol(mm)
qr(mm, LAPACK = FALSE)$rank


###################################################
### code chunk number 32: clm_tutorial.Rnw:657-658
###################################################
convergence(fm1)


###################################################
### code chunk number 33: clm_tutorial.Rnw:672-675
###################################################
slice.fm1 <- slice(fm1, lambda = 5)
par(mfrow = c(2, 3))
plot(slice.fm1)


###################################################
### code chunk number 34: slice11
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 1)


###################################################
### code chunk number 35: slice12
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 2)


###################################################
### code chunk number 36: slice13
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 3)


###################################################
### code chunk number 37: slice14
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 4)


###################################################
### code chunk number 38: slice15
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 5)


###################################################
### code chunk number 39: slice16
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice.fm1, parm = 6)


###################################################
### code chunk number 40: slice2
###################################################
slice2.fm1 <- slice(fm1, lambda = 1e-5)
par(mfrow = c(2, 3))
plot(slice2.fm1)


###################################################
### code chunk number 41: slice21
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 1)


###################################################
### code chunk number 42: slice22
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 2)


###################################################
### code chunk number 43: slice23
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 3)


###################################################
### code chunk number 44: slice24
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 4)


###################################################
### code chunk number 45: slice25
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 5)


###################################################
### code chunk number 46: slice26
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(slice2.fm1, parm = 6)


###################################################
### code chunk number 47: profileLikelihood
###################################################
pr1 <- profile(fm1, alpha=1e-4)
plot(pr1)


###################################################
### code chunk number 48: prof1
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=1)


###################################################
### code chunk number 49: prof2
###################################################
getOption("SweaveHooks")[["fig"]]()
plot(pr1, which.par=2)


###################################################
### code chunk number 50: misc (eval = FALSE)
###################################################
## 


