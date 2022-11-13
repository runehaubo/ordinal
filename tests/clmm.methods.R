library(ordinal)
data(wine)

#################################
## model.matrix method for clmm-objects:
fmm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)
mm <- model.matrix(fmm1)
stopifnot(inherits(mm, "matrix"),
          dim(mm) == c(72, 3))

#################################
## Test that ranef, condVar and VarCorr work as they are supposed to whether or
## not nlme and lme4 are loaded:

fm <- clmm(rating ~ temp + contact + (1|judge), data = wine)
fm
ranef(fm)
VarCorr(fm)
condVar(fm)
summary(fm)

library(nlme)
ranef(fm)
VarCorr(fm)
condVar(fm)
library(lme4)
ranef(fm)
VarCorr(fm)
condVar(fm)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=sleepstudy)
ranef(fm1)
VarCorr(fm1)

ranef(fm)
VarCorr(fm)
condVar(fm)
summary(fm)
