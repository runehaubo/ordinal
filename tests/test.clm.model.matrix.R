library(ordinal)
## source("test.clm.model.matrix.R")

## library(devtools)
## r2path <- "/Users/rhbc/Documents/Rpackages/ordinal/pkg/ordinal"
## clean_dll(pkg = r2path)
## load_all(r2path)

## Check that get_clmDesign works in standard setting:
fm1 <- clm(rating ~ temp, scale=~contact, nominal=~contact, data=wine)
contr <- c(fm1$contrasts, fm1$S.contrasts, fm1$nom.contrasts)
XX <- ordinal:::get_clmDesign(fm1$model, terms(fm1, "all"), contrasts=contr)
XX2 <- update(fm1, method="design")
(keep <- intersect(names(XX), names(XX2)))
(test <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                XX[keep], XX2[keep]))
stopifnot(all(test))

## Check that get_clmDesign works with singular fit and NAs:
cy <- with(wine, which(temp == "cold" & contact == "yes"))
wine2 <- subset(wine, subset=(!1:nrow(wine) %in% cy))
wine2[c(9, 15, 46), "rating"] <- NA
fm1 <- clm(rating ~ temp, scale=~contact, nominal=~contact,
           data=wine2)
contr <- c(fm1$contrasts, fm1$S.contrasts, fm1$nom.contrasts)
XX <- ordinal:::get_clmDesign(fm1$model, terms(fm1, "all"), contrasts=contr)
XX2 <- update(fm1, method="design")
(keep <- intersect(names(XX), names(XX2)))
(test <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                XX[keep], XX2[keep]))
stopifnot(all(test))

## In this situation update and get_clmRho give the same results:
wine2 <- wine
fm1 <- clm(rating ~ temp + contact, data=wine2) ## OK
rho1 <- ordinal:::get_clmRho.clm(fm1)
l1 <- as.list(rho1)
l2 <- as.list(update(fm1, doFit=FALSE))
(test <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                l1, l2))
stopifnot(all(test))
## If we modify the data (or other subset, weights, formulae, etc.)
## used in the model call, the results from update no longer correspond
## to the elements of the fitted model object. get_clmRho gets it
## right on the other hand:
wine2[10:13, "rating"] <- NA
l3 <- as.list(ordinal:::get_clmRho.clm(fm1))
l4 <- as.list(update(fm1, doFit=FALSE))
(test <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                l1, l3))
stopifnot(all(test)) ## same
(test <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                l3, l4))
stopifnot(sum(!test) == 8) ## not all the same anymore!
## In conclusion l1, l2, and l3 are identical. l4 is different.

#################################
## Test that checkContrasts give appropriate warnings:
contr <- c(temp="contr.sum", contact="contr.sum")
fm1 <- clm(rating ~ temp + contact, scale=~contact, data=wine) ## OK
fm1 <- clm(rating ~ temp + contact, scale=~contact, data=wine,
           contrasts=contr) ## OK
fm1 <- clm(rating ~ temp, scale=~contact, data=wine,
           contrasts=contr) ## OK
## These should give warnings:
fm1 <- clm(rating ~ temp, contrasts=c(contact="contr.sum"), data=wine)
fm1 <- clm(rating ~ temp, contrasts=contr, data=wine)
fm1 <- clm(rating ~ 1, scale=~contact, contrasts=c(temp="contr.sum"),
           data=wine)
fm1 <- clm(rating ~ 1, scale=~contact, contrasts=list(temp="contr.sum"),
           data=wine)

fm0 <- clm(rating ~ temp + contact, scale=~contact, data=wine)
ordinal:::checkContrasts(fm0$S.terms, fm0$contrasts)
ordinal:::checkContrasts(fm0$S.terms, fm0$S.contrasts)
ordinal:::checkContrasts(fm0$terms, fm0$contrasts)
ordinal:::checkContrasts(fm0$terms, fm0$S.contrasts)

#################################
## Check that clm and model.matrix respects contrast settings:
options("contrasts" = c("contr.treatment", "contr.poly"))
fm0 <- clm(rating ~ temp + contact, data=wine)
options("contrasts" = c("contr.sum", "contr.poly"))
fm1 <- clm(rating ~ temp + contact, data=wine)
stopifnot(all(model.matrix(fm0)$X[, 2] %in% c(0, 1)))
stopifnot(all(model.matrix(fm1)$X[, 2] %in% c(1, -1)))

#################################
## Check that model.matrix results do not depend on global contrast
## setting:
options("contrasts" = c("contr.sum", "contr.poly"))
fm0 <- clm(rating ~ temp + contact, scale=~contact, data=wine)
MM <- model.matrix(fm0)
options("contrasts" = c("contr.treatment", "contr.poly"))
MM2 <- model.matrix(fm0)
for(x in MM) print(head(x))
for(x in MM2) print(head(x))
stopifnot(all(mapply(all.equal, MM, MM2)))

#################################
## This gave a warning before getContrasts was implemented:
fm0 <- clm(rating ~ temp + contact, scale=~contact, data=wine)
MM <- model.matrix(fm0)
## > fm0 <- clm(rating ~ temp + contact, scale=~contact, data=wine)
## > MM <- model.matrix(fm0)
## Warning message:
## In model.matrix.default(res$S.terms, data = fullmf, contrasts.arg = getContrasts(res$S.terms,  :
##   variable 'temp' is absent, its contrast will be ignored
for(x in MM) print(head(x))

