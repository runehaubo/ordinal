library(ordinal)

#################################
## 1 categorical variable in nominal:
fm <- clm(rating ~ temp, nominal=~contact, data=wine)
fm$Theta
fm$alpha.mat
## Threshold effects:
fm <- clm(rating ~ temp, nominal=~contact, data=wine,
          threshold="symmetric")
fm$Theta
fm$alpha.mat
fm <- clm(rating ~ temp, nominal=~contact, data=wine,
          threshold="equidistant")
fm$Theta
fm$alpha.mat
## Singular fit is still ok (with a warning, though)
fm <- clm(rating ~ contact, nominal=~temp, data=wine)
fm$alpha.mat
fm$Theta

#################################
## 1 continuous variable:
set.seed(123)
x <- rnorm(nrow(wine), sd=1)
fm <- clm(rating ~ temp, nominal=~ x, data=wine)
fm$alpha.mat
fm$Theta
fm <- clm(rating ~ temp, nominal=~ poly(x, 2), data=wine)
fm$alpha.mat
fm$Theta

#################################
## 1 categorical + 1 continuous variable:
set.seed(123)
x <- rnorm(nrow(wine), sd=1)
fm <- clm(rating ~ temp, nominal=~contact + x, data=wine)
fm$alpha.mat
fm$Theta
fm <- clm(rating ~ temp, nominal=~contact + x, data=wine,
          threshold="symmetric")
fm$alpha.mat
fm$Theta
#################################
### NOTE: To get the by-threshold nominal effects of continuous terms
## use:
with(fm, t(apply(alpha.mat, 1, function(th) tJac %*% th)))
#################################
## Interactions:
fm <- clm(rating ~ temp, nominal=~contact:x, data=wine)
fm$alpha.mat
fm$Theta
fm <- clm(rating ~ temp, nominal=~contact+x+contact:x, data=wine)
fm$alpha.mat
fm$Theta
fm <- clm(rating ~ temp, nominal=~contact*x, data=wine)
fm$alpha.mat
fm$Theta
## polynomial terms:
fm <- clm(rating ~ temp, nominal=~contact + poly(x, 2), data=wine)
fm$alpha.mat
fm$Theta
## logical variables: (treated like numeric variables)
wine$Con <- as.character(wine$contact) == "yes"
fm <- clm(rating ~ temp, nominal=~Con, data=wine)
fm$Theta
fm$alpha.mat
wine$Con.num <- 1 * wine$Con
fm <- clm(rating ~ temp, nominal=~Con.num, data=wine)
fm$Theta
fm$alpha.mat
#################################
## Two continuous variables:
set.seed(321)
y <- rnorm(nrow(wine), sd=1)
fm1 <- clm(rating ~ temp, nominal=~y + x, data=wine)
fm1$alpha.mat
fm1$Theta
## summary(fm1)

#################################
## 1 categorical + 2 continuous variables:
fm1 <- clm(rating ~ temp, nominal=~y + contact + x, data=wine)
fm1$alpha.mat
fm1$Theta

fm1 <- clm(rating ~ temp, nominal=~contact + x + contact:x + y,
           data=wine)
summary(fm1)
fm1$Theta
fm1$alpha.mat
fm1 <- clm(rating ~ temp, nominal=~contact*x + y, data=wine)
fm1$Theta
fm1$alpha.mat
t(fm1$alpha.mat)
fm1

#################################
## ordered factors (behaves like numerical variables):
data(soup, package="ordinal")
fm2 <- clm(SURENESS ~ 1, nominal=~PRODID + DAY, data=soup)
fm2$Theta
fm2$alpha.mat
prodid <- factor(soup$PRODID, ordered=TRUE)
fm2 <- clm(SURENESS ~ 1, nominal=~prodid + DAY, data=soup)
fm2$alpha.mat
fm2$Theta
fm2 <- clm(SURENESS ~ 1, nominal=~prodid, data=soup)
fm2$alpha.mat
fm2$Theta
#################################
## Aliased Coefficients:
##
## Example where the interaction in the nominal effects is aliased (by
## design). Here the two Theta matrices coincide. The alpha.mat
## matrices are similar except one has an extra row with NAs:
soup2 <- soup
levels(soup2$DAY)
levels(soup2$GENDER)
xx <- with(soup2, DAY == "2" & GENDER == "Female")
## Model with additive nominal effects:
fm8 <- clm(SURENESS ~ PRODID, nominal= ~ DAY + GENDER, data=soup2, subset=!xx)
fm8$alpha.mat
fm8$Theta
## Model with non-additive, but aliased nominal effects:
fm9 <- clm(SURENESS ~ PRODID, nominal= ~ DAY * GENDER, data=soup2, subset=!xx)
fm9$alpha.mat
fm9$Theta

stopEqual <- function(x, y, ca=FALSE)
    stopifnot(isTRUE(all.equal(x, y, check.attributes=ca)))

stopEqual(fm8$alpha.mat, fm9$alpha.mat[1:3, ])
stopEqual(fm8$Theta, fm9$Theta)
stopEqual(logLik(fm8), logLik(fm9))

#################################
## Weights:
set.seed(12345)
wts <- runif(nrow(soup))
fm2 <- clm(SURENESS ~ 1, nominal=~SOUPTYPE + DAY, data=soup, weights=wts)
fm2$Theta

## Offset (correctly gives and error)
fm2 <- try(clm(SURENESS ~ 1, nominal=~SOUPTYPE + DAY + offset(wts),
               data=soup), silent=TRUE)
stopifnot(inherits(fm2, "try-error"))

#################################
### Other (misc) examples:
fm2 <- clm(SURENESS ~ 1, nominal=~SOUPTYPE + DAY, data=soup)
fm2$Theta
fm2
fm2 <- clm(SURENESS ~ 1, nominal=~SOUPTYPE * DAY, data=soup)
fm2$Theta
fm2
fm2$alpha.mat
fm2 <- clm(SURENESS ~ 1, nominal=~SOUPTYPE * DAY, data=soup,
           threshold="symmetric")
fm2$Theta
fm2$alpha.mat

#################################
### Check correctness of Theta matrix when intercept is removed in
### nominal formula:
### December 25th 2014, RHBC
fm1 <- clm(rating ~ temp, nominal=~contact-1, data=wine)
fm2 <- clm(rating ~ temp, nominal=~contact, data=wine)
stopifnot(isTRUE(all.equal(fm1$Theta, fm2$Theta)))
stopifnot(isTRUE(all.equal(fm1$logLik, fm2$logLik)))
wine2 <- wine
wine2$contact <- relevel(wine2$contact, "yes")
fm3 <- clm(rating ~ temp, nominal=~contact, data=wine2)
stopifnot(isTRUE(all.equal(coef(fm1, na.rm=TRUE), coef(fm3))))
#################################

