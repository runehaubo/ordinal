library(ordinal)
data(wine)

fm1 <- clm(rating ~ temp, data=wine)
fmm1 <- clmm(rating ~ temp + (1|judge), data=wine)

## These now give identical printed results:
## Previously the printed model names were messed up when anova.clmm
## were called. 
anova(fm1, fmm1)
anova(fmm1, fm1)

anova(fm1, fmm1)

## Testing if 'test' and 'type' arguments are ignored properly: 
fm1 <- clm(rating ~ temp + contact, data=wine)
fm2 <- clm(rating ~ temp, data=wine)
anova(fm1, fm2, test="Chi")
anova(fm1, fm2, type="Chi")
anova(fm1, fm2)
## calling anova.clmm
anova(fmm1, fm1, test="Chi")
anova(fmm1, fm1, type="Chi")

