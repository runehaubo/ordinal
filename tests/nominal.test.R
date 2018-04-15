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

