library(ordinal)
data(wine)

#################################
## Estimation with a single simple RE term:
## Laplace:
fmm1 <- clmm(rating ~ contact + temp + (1|judge), data=wine)
summary(fmm1)
## GHQ:
fmm.ghq <- clmm(rating ~ contact + temp + (1|judge), data=wine,
                nAGQ=-10)
summary(fmm.ghq)
## AGQ:
fmm.agq <- clmm(rating ~ contact + temp + (1|judge), data=wine,
                nAGQ=10)
summary(fmm.agq)
## tests:
## Notice warning about Laplace with multiple REs when nAGQ != 1:
fmm1 <- try(clmm(rating ~ contact + temp + (1|judge) + (1|bottle),
                 data=wine, nAGQ=10))
stopifnot(inherits(fmm1, "try-error"))

#################################
## Estimation with several RE terms:
data(soup, package="ordinal")
fmm <- clmm(SURENESS ~ PROD + (1|RESP) + (1|PROD:RESP), data=soup,
            threshold="equidistant")
summary(fmm)

#################################

## Estimation with implicit intercept:
fm1 <- clmm(rating ~ 1 + (1|judge), data = wine) 
fm2 <- clmm(rating ~ (1|judge), data = wine) 
fm3 <- clmm(rating ~ 0 + (1|judge), data = wine)
stopifnot(isTRUE(all.equal(coef(fm1), coef(fm2), tolerance=1e-5)),
          isTRUE(all.equal(coef(fm1), coef(fm3), tolerance=1e-5)))
