
context("Appropriate error and warning messages from clm()")

test_that("formula is specified in clm", {
    expect_error(clm(nominal=~contact, data=wine),
                 "Model needs a formula")
    expect_error(clm(scale=~contact, data=wine),
                 "Model needs a formula")
    expect_error(clm(),
                 "Model needs a formula")
})

test_that("response is not in scale or nominal", {
    ## No response in formula:
    expect_error(
        fm <- clm(~ temp + contact, data=wine)
        , "'formula' needs a response")
    ## response in scale:
    expect_error(
        fm <- clm(rating ~ temp, scale=rating ~ contact, data=wine)
        , "response not allowed in 'scale'")
    expect_error(
        fm <- clm(rating ~ temp, nominal=rating ~ contact, data=wine)
        , "response not allowed in 'nominal'")
    wine2 <- wine
    wine2$rate <- as.numeric(as.character(wine2$rating))
    expect_error(
        fm <- clm(rate ~ temp + contact, data=wine2)
        , "response in 'formula' needs to be a factor")
})

test_that("offset is allowed in formula, but not in scale and nominal",
{
    wine2 <- wine
    set.seed(1)
    wine2$off <- runif(nrow(wine))
    ## offset in formula is fine:
    expect_is(
        clm(rating ~ temp + contact + offset(off), data=wine2)
        , "clm")
    expect_is(
        clm(rating ~ offset(off), nominal=~contact, data=wine2)
        , "clm") ## no other terms in formula.
    ## offset in scale is also fine:
    expect_is(
        clm(rating ~ temp, scale=~contact + offset(off), data=wine2)
        , "clm")
    expect_is(
        clm(rating ~ contact + temp, scale=~offset(off), data=wine2)
        , "clm") ## no other terms in scale.
    ## offset as argument is not allowed:
    expect_error(
        clm(rating ~ temp + contact, offset=off, data=wine2)
        , "offset argument not allowed: specify 'offset' in formula or scale arguments instead")
    ## offset in nominal is not allowed:
    expect_error(
        clm(rating ~ temp, nominal=~contact + offset(off), data=wine2)
        , "offset not allowed in 'nominal'")
    expect_error(
        clm(rating ~ temp, nominal=~1 + offset(off), data=wine2)
        , "offset not allowed in 'nominal'")
})


test_that("Intercept is needed and assumed", {
    expect_is(
        fm <- clm(rating ~ 1, data=wine)
        , "clm")
    expect_warning(
        fm <- clm(rating ~ -1 + temp, data=wine)
        , "an intercept is needed and assumed in 'formula'")
    expect_warning(
        fm <- clm(rating ~ 0 + temp, data=wine)
        , "an intercept is needed and assumed in 'formula'")
    expect_warning(
        fm <- clm(rating ~ 0, data=wine)
        , "an intercept is needed and assumed in 'formula'")
    ## and similar with scale (+nominal)
})

## test_that("", {
##
## })
