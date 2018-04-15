library(ordinal)
## library(devtools)
## r2path <- "/Users/rhbc/Documents/Rpackages/ordinal/pkg/ordinal"
## clean_dll(pkg = r2path)
## load_all(r2path)

#################################
## Appropriate evaluation of formulas:

## These fail and give appropriate error messages:
##  fm1 <- clm(rating ~ contact, scale=temp, data=wine)
##  fm1 <- clm(rating ~ contact, scale=~Temp, data=wine)
##  fm1 <- clm(rating ~ contact, scale="temp", data=wine)
##  sca <- "temp"
##  fm1 <- clm(rating ~ contact, scale=sca, data=wine)
##  sca <- as.formula(sca)
##  sca <- as.formula(temp)
##  sca <- with(wine, as.formula(temp))

## These all work as intended with no warnings or errors:
fm1 <- clm(rating ~ contact, scale="~temp", data=wine)
fm1 <- clm(rating ~ contact, scale=~temp, data=wine)
sca <- "~temp"
fm1 <- clm(rating ~ contact, scale=sca, data=wine)
sca <- as.formula("~temp")
fm1 <- clm(rating ~ contact, scale=sca, data=wine)
fm1 <- clm(rating ~ contact, scale=as.formula(~temp), data=wine)
fm1 <- clm(rating ~ contact, scale=as.formula("~temp"), data=wine)

#################################
## can evaluate if 'formula' is a character:
f <- "rating ~ contact + temp"
clm(f, data=wine)
clm(as.formula(f), data=wine)

#################################

### finding variables in the environment of the formula:
makeform <- function() {
  f1 <- as.formula(rating ~ temp + contact)
  rating <- wine$rating
  temp <- wine$temp
  contact <- wine$contact
  f1
}
## 'makeform' makes are formula object in the environment of the
## function makeform:
f1 <- makeform()
f1 # print
class(f1)
## If we give the data, we can evaluate the model:
fm1 <- clm(f1, data=wine)
## We can also evaluate the model because the data are available in
## the environment associated with the formula:
fm1 <- clm(f1)
## For instance, the 'rating' variable is not found in the Global
## environment; we have to evaluate the 'name' of 'rating' in the
## appropriate environment:
(try(rating, silent=TRUE))
eval(as.name("rating"), envir=environment(f1))
## If instead we generate the formula in the Global environment where
## the variables are not found, we cannot evaluate the model:
f2 <- as.formula(rating ~ temp + contact)
(try(fm2 <- clm(f2), silent=TRUE))
environment(f2) <- environment(f1)
fm2 <- clm(f2)


#################################
## Use of formula-objects in location, scale and nominal:
## Bug-report from Lluís Marco Almagro <lluis.marco@upc.edu>
## 5 May 2010 17:58
f <- formula(rating ~ temp)
fs <- formula( ~ contact)
m2 <- clm(f, scale = fs, data = wine)
summary(m2)

#################################
## Other ways to construct formulas:
set.seed(12345)
y <- factor(sample(1:4,20,replace=TRUE))
x <- rnorm(20)
data <- data.frame(y=y,x=x)
rm(x, y)
fit <- clm(data$y ~ data$x)
fit
fit <- clm(data[,1] ~ data[,2])
fit
## This previously failed, but now works:
fit <- clm(data$y ~ data$x, ~data$x)
fit

#################################
## Evaluation within other functions:
## date: January 18th 2012.
##
## The problem was raised by Stefan Herzog (stefan.herzog@unibas.ch)
## January 12th 2012 in trying to make clm work with glmulti.

fun.clm <- function(formula, data)
### This only works because clm via eclm.model.frame is careful to
### evaluate the 'formula' in the parent environment such it is not the
### character "formula" that is attempted evaluated.
  clm(formula, data = data)

fun2.clm <- function(formula, data, weights, subset) {
### This should be the safe way to ensure evaluation of clm in the
### right environment.
  mc <- match.call()
  mc[[1]] <- as.name("clm")
  eval.parent(mc)
}

fun.clm(rating ~ temp + contact, data=wine) ## works
fun2.clm(rating ~ temp + contact, data=wine) ## works

form1 <- "rating ~ temp + contact"
fun.clm(form1, data=wine) ## works
fun2.clm(form1, data=wine) ## works

form2 <- formula(rating ~ temp + contact)
fun.clm(form2, data=wine) ## works
fun2.clm(form2, data=wine) ## works
## Notice that clm is not able to get the name of the data (wine)
## correct when using fun.clm.

#################################
## Evaluation of long formulas: no line breaking in getFullForm:
data(soup, package="ordinal")

rhs <- paste(names(soup)[c(3, 5:12)], collapse=" + ")
Location <- as.formula(paste("SURENESS ~ ", rhs, sep=" "))
Scale <- as.formula("~ PROD")

fm5 <- clm(Location, scale=Scale, data=soup)
summary(fm5)

#################################
## Check that "."-notation works in formula:
## December 25th 2014, RHBC
data(wine)
wine2 <- wine[c("rating", "contact", "temp")]
str(wine2)
fm0 <- clm(rating ~ ., data=wine2)
fm1 <- clm(rating ~ contact + temp, data=wine2)
keep <- c("coefficients", "logLik", "info")
fun <- function(x, y) stopifnot(isTRUE(all.equal(x, y)))
mapply(fun, fm0[keep], fm1[keep])
#################################
