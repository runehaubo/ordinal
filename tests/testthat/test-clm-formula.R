context("Appropriate evaluation of formulae in clm()")

## These fail and give appropriate error messages:
test_that("standard formulae are interpreted correctly/give right error messages", {
    expect_error(
        fm1 <- clm(rating ~ contact, scale=temp, data=wine)
        , "object 'temp' not found")
    expect_error(
        fm1 <- clm(rating ~ contact, scale=~Temp, data=wine)
        , "object 'Temp' not found")
    expect_error(
        fm1 <- clm(rating ~ contact, scale="temp", data=wine)
        , "unable to interpret 'formula', 'scale' or 'nominal'")
    sca <- "temp"
    expect_error(
        fm1 <- clm(rating ~ contact, scale=sca, data=wine)
        , "unable to interpret 'formula', 'scale' or 'nominal'")
    ## sca <- as.formula(sca)
    ## sca <- as.formula(temp)
    ## sca <- with(wine, as.formula(temp))

    ## These all work as intended with no warnings or errors:
    fm1 <- clm(rating ~ contact, scale="~temp", data=wine)
    fm2 <- clm(rating ~ contact, scale=~temp, data=wine)
    sca <- "~temp"
    fm3 <- clm(rating ~ contact, scale=sca, data=wine)
    sca <- as.formula("~temp")
    fm4 <- clm(rating ~ contact, scale=sca, data=wine)
    fm5 <- clm(rating ~ contact, scale=as.formula(~temp), data=wine)
    fm6 <- clm(rating ~ contact, scale=as.formula("~temp"), data=wine)

    ## Test that they are all clm objects:
    for(txt in paste0("fm", 1:6))
        expect_is(eval(parse(text=txt)), "clm")


#################################
    ## can evaluate if 'formula' is a character:
    f <- "rating ~ contact + temp"
    expect_is(clm(f, data=wine), "clm")
    expect_is(clm(as.formula(f), data=wine), "clm")

#################################
})

test_that("variables are found in the right environments", {
    ## finding variables in the environment of the formula:
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
    expect_is(f1, "formula")
    ## If we give the data, we can evaluate the model:
    expect_is(fm1 <- clm(f1, data=wine), "clm")
    ## We can also evaluate the model because the data are available in
    ## the environment associated with the formula:
    expect_is(fm1 <- clm(f1), "clm")
    ## For instance, the 'rating' variable is not found in the Global
    ## environment; we have to evaluate the 'name' of 'rating' in the
    ## appropriate environment:
    (try(rating, silent=TRUE))
    expect_error(
        rating
        , "'rating' not found")
    expect_is(
        eval(as.name("rating"), envir=environment(f1))
        , "factor")
    ## If instead we generate the formula in the Global environment where
    ## the variables are not found, we cannot evaluate the model:
    f2 <- as.formula(rating ~ temp + contact)
    expect_error(
        fm2 <- clm(f2)
        )
    ## Setting the appropriate environment of the formula restores the
    ## ability to evaluate the model:
    environment(f2) <- environment(f1)
    expect_is(
        fm2 <- clm(f2)
        , "clm")

    #################################
    ## Use of formula-objects in location, scale and nominal:
    ## Bug-report from Lluís Marco Almagro <lluis.marco@upc.edu>
    ## 5 May 2010 17:58
    f <- formula(rating ~ temp)
    fs <- formula( ~ contact)
    expect_is(
        m2 <- clm(f, scale = fs, data = wine)
        , "clm")
})

test_that("data indexing works in formulae", {
#################################
    ## Other ways to construct formulas:
    set.seed(12345)
    y <- factor(sample(1:4,20,replace=TRUE))
    x <- rnorm(20)
    data <- data.frame(y=y,x=x)
    rm(x, y)
    expect_is(
        fit <- clm(data$y ~ data$x)
        , "clm")
    expect_is(
        fit <- clm(data[,1] ~ data[,2])
        , "clm")
    ## This previously failed, but now works:
    expect_is(
        fit <- clm(data$y ~ data$x, ~data$x)
        , "clm")
})

test_that("clm may be invoked within functions", {
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

    expect_is(
        fun.clm(rating ~ temp + contact, data=wine) ## works
        , "clm")
    expect_is(
        fun2.clm(rating ~ temp + contact, data=wine) ## works
        , "clm")

    form1 <- "rating ~ temp + contact"
    expect_is(
        fun.clm(form1, data=wine) ## works
        , "clm")
    expect_is(
        fun2.clm(form1, data=wine) ## works
        , "clm")

    form2 <- formula(rating ~ temp + contact)
    expect_is(
        fm1 <- fun.clm(form2, data=wine) ## works
        , "clm")
    expect_is(
        fm2 <- fun2.clm(form2, data=wine) ## works
        , "clm")
    ## Notice that clm is not able to get the name of the data (wine)
    ## correct when using fun.clm:
    expect_true(deparse(fm1$call$data) == "data")
    expect_true(deparse(fm2$call$data) == "wine")
})

test_that("no line breacking in long formulae", {
    #################################
    ## Evaluation of long formulas: no line breaking in getFullForm:

    rhs <- paste(names(soup)[c(3, 5:12)], collapse=" + ")
    Location <- as.formula(paste("SURENESS ~ ", rhs, sep=" "))
    Scale <- as.formula("~ PROD")
    expect_is(
        fm5 <- clm(Location, scale=Scale, data=soup)
        , "clm")
})

test_that("'.'-notation works in formula", {
    #################################
    ## Check that "."-notation works in formula:
    ## December 25th 2014, RHBC
    data(wine)
    wine2 <- wine[c("rating", "contact", "temp")]
    ## str(wine2)
    fm0 <- clm(rating ~ ., data=wine2)
    fm1 <- clm(rating ~ contact + temp, data=wine2)
    keep <- c("coefficients", "logLik", "info")
    fun <- function(x, y) stopifnot(isTRUE(all.equal(x, y)))
    mapply(fun, fm0[keep], fm1[keep])
    fun <- function(x, y) {expect_equal(x, y); invisible()}
    mapply(fun, fm0[keep], fm1[keep])
    #################################
})
