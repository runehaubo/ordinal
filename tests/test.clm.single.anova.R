# test.clm.single.anova.R

library(ordinal)

fm <- clm(rating ~ temp * contact, scale=~contact, data=wine)

anova(fm, type="I")
anova(fm, type="II")
anova(fm, type="III")
anova(fm, type=1)
anova(fm, type=2)
anova(fm, type=3)
anova(fm, type="1")
anova(fm, type="2")
anova(fm, type="3")
anova(fm, type="marginal")

# Nominal effects:
fm <- clm(rating ~ temp, nominal=~contact, data=wine)
anova(fm)

# Flexible links:
# fm <- clm(rating ~ temp + contact, link="log-gamma", data=wine)
# anova(fm1, type=1)
# anova(fm1, type=2)
# anova(fm1, type=3)

# Aliased coefficients:
fm1 <- clm(SURENESS ~ PRODID * DAY, data=soup)
anova(fm1, type=1)
anova(fm1, type=2)
anova(fm1, type=3)

# Aliased term (due to nominal effects):
fm <- clm(rating ~ temp * contact, nominal=~contact, data=wine)
anova(fm, type=1)
anova(fm, type=2)
anova(fm, type=3)

# model with all NA in vcov(object):
fm <- clm(rating ~ temp * contact, nominal=~contact, scale=~contact, data=wine)
# anova(fm, type=1) # error
# anova(fm, type=2) # error
# anova(fm, type=3) # error
all(is.na(vcov(fm)))
