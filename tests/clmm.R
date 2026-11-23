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

## Both fitting paths must use the requested asymmetric link (issue #76):
for (link in c("cloglog", "loglog")) {
    fit <- clmm(rating ~ temp + contact + (1|judge), data=wine, link=link)
    fit.matrix <- clmm(rating ~ temp + contact + (1|judge),
                       data=wine, link=link,
                       control=clmm.control(useMatrix=TRUE))
    fit.legacy <- clmm2(rating ~ temp + contact, random=judge,
                        data=wine, link=link)
    expected <- as.numeric(logLik(fit.matrix))
    stopifnot(abs(as.numeric(logLik(fit)) - expected) < 1e-3,
              abs(as.numeric(logLik(fit.legacy)) - expected) < 1e-3,
              isTRUE(all.equal(coef(fit), coef(fit.matrix), tolerance=1e-4)))
}
