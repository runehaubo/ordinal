# test.clm.single.anova.R

library(ordinal)

# WRE says "using if(requireNamespace("pkgname")) is preferred, if possible."
# even in tests:
assertError <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertError(expr, ...) else invisible()
assertWarning <- function(expr, ...)
  if(requireNamespace("tools")) tools::assertWarning(expr, ...) else invisible()

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
fm1 <- clm(rating ~ temp + contact, link="log-gamma", data=wine)
anova(fm1, type=1)
anova(fm1, type=2)
anova(fm1, type=3)

# Equivalence of tests irrespective of contrasts:
fm1 <- clm(SURENESS ~ PRODID * SOUPFREQ, data=soup)
# summary(fm1)
(an1 <- anova(fm1, type=3))
fm2 <- clm(SURENESS ~ PRODID * SOUPFREQ, data=soup,
           contrasts = list(SOUPFREQ = "contr.sum", PRODID = "contr.SAS"))
# summary(fm2)
anova(fm1, fm2)
(an2 <- anova(fm2, type=3))
stopifnot(
  isTRUE(all.equal(an1, an2, check.attributes=FALSE))
)


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
assertError(anova(fm, type=1)) # error
assertError(anova(fm, type=2)) # error
assertError(anova(fm, type=3)) # error
all(is.na(vcov(fm)))
