library(ordinal)

if(require(MASS)) {
    fm1 <- clm(Sat ~ Infl + Type + Cont, data=housing, weights=Freq)
    scale_test(fm1)
    nominal_test(fm1)

    fm2 <- update(fm1, scale=~Cont)
    scale_test(fm2)
    nominal_test(fm2)
    fm3 <- update(fm1, nominal=~ Cont)
    fm3$Theta
    anova(fm2, fm3)
    fm3$alpha.mat
    summary(fm3)
}

#################################
### Testing nominal_test and scale_test:
fm1 <- clm(rating ~ temp * contact, data=wine)
## names(fm1)
fm2 <- clm(rating ~ temp * contact, data=wine, nominal=~contact)
(an <- anova(fm1, fm2))
(nm <- nominal_test(fm1))
stopifnot(isTRUE(all.equal(an[2, 6], nm["contact", 5])))

fm2 <- clm(rating ~ temp * contact, data=wine, scale=~contact)
(an <- anova(fm1, fm2))
(sc <- scale_test(fm1))
stopifnot(isTRUE(all.equal(an[2, 6], sc["contact", "Pr(>Chi)"])))

## Refit a model whose formula, data, weights and subset are local.
local_clm <- function(local_data, response, predictors, local_wt, local_subset) {
    fm <- ordinal::clm(reformulate(predictors, response), data=local_data,
                       weights=local_wt, subset=local_subset, model=FALSE)
    list(fit=fm, nominal=nominal_test(fm), scale=scale_test(fm))
}
ref_wt <- rep(1:3, length.out=nrow(wine))
ref_subset <- seq_len(nrow(wine)) %% 7L != 0L
ref <- clm(rating ~ temp + contact, data=wine,
           weights=ref_wt, subset=ref_subset, model=FALSE)
expected_nominal <- nominal_test(ref)
expected_scale <- scale_test(ref)
stopifnot(all(is.finite(expected_nominal[-1L, "LRT"])),
          all(is.finite(expected_scale[-1L, "LRT"])))

local_result <- local_clm(wine, "rating", c("temp", "contact"),
                          ref_wt, ref_subset)
stopifnot(isTRUE(all.equal(local_result$nominal, expected_nominal)),
          isTRUE(all.equal(local_result$scale, expected_scale)))

restored <- unserialize(serialize(local_result$fit, NULL))
stopifnot(isTRUE(all.equal(nominal_test(restored), expected_nominal)),
          isTRUE(all.equal(scale_test(restored), expected_scale)))

attr(restored, "call.envir") <- NULL
stopifnot(isTRUE(all.equal(nominal_test(restored), expected_nominal)),
          isTRUE(all.equal(scale_test(restored), expected_scale)))

## Recover the fitting environment even when the formula came from elsewhere.
fit_formula <- rating ~ temp + contact
local_prebuilt <- function(dataset, formula) {
    ordinal::clm(formula, data=dataset, weights=ref_wt,
                 subset=ref_subset, model=FALSE)
}
prebuilt <- local_prebuilt(wine, fit_formula)
prebuilt <- unserialize(serialize(prebuilt, NULL))
stopifnot(isTRUE(all.equal(nominal_test(prebuilt), expected_nominal)),
          isTRUE(all.equal(scale_test(prebuilt), expected_scale)))

## A same-sized data set with the same name must not replace the fit's data.
set.seed(5)
dataset <- transform(wine, rating=sample(rating))
fit_collision <- function(dataset) {
    ordinal::clm(reformulate(c("temp", "contact"), "rating"),
                 data=dataset, model=FALSE)
}
ref <- clm(rating ~ temp + contact, data=wine, model=FALSE)
collision <- fit_collision(wine)
stopifnot(isTRUE(all.equal(nominal_test(collision), nominal_test(ref))),
          isTRUE(all.equal(scale_test(collision), scale_test(ref))))

fm1 <- clm(rating ~ temp + contact,
           nominal=~temp + contact, data=wine)
fm1
try(nominal_test(fm1), silent=TRUE)[1] ## gives error OK
scale_test(fm1)
fm1 <- clm(rating ~ temp + contact,
           scale=~temp + contact, data=wine)
fm1
try(scale_test(fm1), silent=TRUE)[1] ## gives error OK
nominal_test(fm1)


## Using weights:
set.seed(123454321)
wt <- runif(nrow(wine))
fm1 <- clm(rating ~ temp * contact, data=wine, weigths=wt)
nominal_test(fm1)
scale_test(fm1)

## No nominal test for judge since that model is not identifiable:
fm1 <- clm(rating ~ judge + temp + contact, data=wine)
nominal_test(fm1)
scale_test(fm1)
fm1 <- clm(rating ~ judge + temp, nominal=~contact, data=wine)
nominal_test(fm1)
summary(fm1)

## A continuous variable:
set.seed(123454321)
x <- rnorm(nrow(wine), sd=1)
fm <- clm(rating ~ temp, nominal=~contact * x, data=wine)
nominal_test(fm)
scale_test(fm)
fm <- clm(rating ~ temp + x, nominal=~contact, data=wine)
nominal_test(fm)
scale_test(fm)
## poly:
fm <- clm(rating ~ temp + poly(x, 2), nominal=~contact, data=wine)
nominal_test(fm)
scale_test(fm)
## another combination:
fm1 <- clm(SURENESS ~ PRODID + DAY + SOUPTYPE + SOUPFREQ,
           scale=~PROD,
           nominal=~ DAY*GENDER, data=soup)
fm1
nominal_test(fm1)
scale_test(fm1)

#################################

