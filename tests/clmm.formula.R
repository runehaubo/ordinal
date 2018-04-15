library(ordinal)
data(wine)

#################################
## Appropriate evaluation of formulas:

## These all work as intended with no warnings or errors:
fm1 <- clmm(rating ~ contact + (1|judge), data=wine)
fm1
fm1 <- clmm("rating ~ contact + (1|judge)", data=wine)
fm1
fm1 <- clmm(as.formula("rating ~ contact + (1|judge)"), data=wine)
fm1 
fm1 <- clmm(as.formula(rating ~ contact + (1|judge)), data=wine) 
fm1

#################################

### finding variables in the environment of the formula:
makeform <- function() {
  f1 <- as.formula(rating ~ temp + contact + (1|judge))
  rating <- wine$rating
  temp <- wine$temp
  contact <- wine$contact
  judge <- wine$judge
  f1
}
## 'makeform' makes are formula object in the environment of the
## function makeform:
f1 <- makeform()
f1 # print
class(f1)
## If we give the data, we can evaluate the model:
fm1 <- clmm(f1, data=wine)
## We can also evaluate the model because the data are available in
## the environment associated with the formula:
fm1 <- clmm(f1)
## For instance, the 'rating' variable is not found in the Global
## environment; we have to evaluate the 'name' of 'rating' in the
## appropriate environment:
(try(rating, silent=TRUE))
eval(as.name("rating"), envir=environment(f1))
## If instead we generate the formula in the Global environment where
## the variables are not found, we cannot evaluate the model:
f2 <- as.formula(rating ~ temp + contact + (1|judge))
(try(fm2 <- clmm(f2), silent=TRUE))
environment(f2) <- environment(f1)
fm2 <- clmm(f2)

#################################
## Use of formula-objects
f <- formula(rating ~ temp + contact + (1|judge))
m2 <- clmm(f, data = wine)
summary(m2)

#################################
## Other ways to construct formulas:
set.seed(12345)
y <- factor(sample(1:4,20,replace=TRUE))
x <- rnorm(20)
b <- gl(5, 4, labels=letters[1:5])
data <- data.frame(y=y, x=x, b=b)
rm(x, y, b)
clmm(y ~ x + (1|b), data=data)
fit <- clmm(data$y ~ data$x + (1|data$b))
fit
fit <- clmm(data[, 1] ~ data[, 2] + (1|data[, 3]))
fit

#################################
## Evaluation within other functions:
## date: January 18th 2012.
## 
## The problem was raised by Stefan Herzog (stefan.herzog@unibas.ch)
## January 12th 2012 in trying to make clmm work with glmulti.

fun.clmm <- function(formula, data)
### This only works because clmm via eclmm.model.frame is careful to
### evaluate the 'formula' in the parent environment such it is not the
### character "formula" that is attempted evaluated.
  clmm(formula, data = data)

fun2.clmm <- function(formula, data, weights, subset) {
### This should be the safe way to ensure evaluation of clmm in the
### right environment.
  mc <- match.call()
  mc[[1]] <- as.name("clmm")
  eval.parent(mc)
}

fun.clmm(rating ~ temp + contact + (1|judge), data=wine) ## works
fun2.clmm(rating ~ temp + contact + (1|judge), data=wine) ## works

form1 <- "rating ~ temp + contact + (1|judge)"
fun.clmm(form1, data=wine) ## works
fun2.clmm(form1, data=wine) ## works

form2 <- formula(rating ~ temp + contact + (1|judge))
fun.clmm(form2, data=wine) ## works
fun2.clmm(form2, data=wine) ## works
## Notice that clmm is not able to get the name of the data (wine)
## correct when using fun.clmm.

#################################

##    ## Example 2: using clmm function
##    #
##    ## Now I want to consider judge as a random effect to account for
##    ## grouping structure of data 
##    mod2 <- clmm(rating ~ temp + contact + (1|judge), data=wine)
##    
##    ##Again, I started by using my own code to run all potential models:
##    ## put names of all your variables in this vector:
##    vl2 <- c("temp", "contact")
##    ## generate list of possible combinations of variables:
##    combos2 <- NULL
##    for(i in 1:length(vl2)) {
##      combos2 <- c(combos2, combn(vl2, i, simplify = F))
##    }
##    ## create formulae and run models one by one, saving them as model1,
##    ## model2 etc... 
##    for (i in 1:length(combos2)) {
##      vs2 <- paste(combos2[[i]], collapse=" + ")
##      f2 <- formula(paste("rating ~ ", vs2, "+(1|judge)", sep=""))
##      print(f2)
##      assign(paste("model", i, sep=""), clmm(f2, data=wine))
##    }
##    summary(model1) # etc
##    summary(model2) # etc
##    summary(model3) # etc
##    
##    models <- vector("list", length(combos2))
##    for(i in 1:length(combos2)) {
##      vs2 <- paste(combos2[[i]], collapse=" + ")
##      f2 <- formula(paste("rating ~ ", vs2, "+(1|judge)", sep=""))
##      print(f2)
##      models[[i]] <- clmm(f2, data=wine)
##      ## assign(paste("model", i, sep=""), clmm(f2, data=wine))
##    }
##    
##    ## Coefficients, AIC and BIC:
##    lapply(models, function(m) coef(summary(m)))
##    lapply(models, AIC)
##    lapply(models, BIC)
##    
##    ## library(MuMIn)
##    ## dd2 <- dredge(mod2) ## does not work
##    ## ?dredge
##    ## traceback()
##    ## mod2$formula
##    ## terms(as.formula(formula(mod2)))
##    ## 
##    ## library(lme4)
##    ## fmm1 <- lmer(response ~ temp + contact + (1|judge), data=wine)
##    ## fmm1
##    ## terms(as.formula(lme4:::formula(fmm1)))
##    ## terms(as.formula(formula(fmm1)))
