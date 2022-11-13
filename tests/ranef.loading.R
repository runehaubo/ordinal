# check that ranef and VarCorr work even after loading ordinal: 
library(lme4)
fm1 <- lmer(Reaction ~ Days + (Days | Subject), data=sleepstudy)
ranef(fm1)
VarCorr(fm1)
library(ordinal)
ranef(fm1)
VarCorr(fm1)
