context("Test that clm.predict gives warnings if prevars is absent")

fm1 <- clm(rating ~ temp + contact, data=wine)
newData <- expand.grid(temp=levels(wine$temp),
                       contact=levels(wine$contact))
expect_false(givesWarnings(
    predict(fm1, newdata=newData)
    ))
attr(fm1$terms, "predvars") <- NULL
expect_warning(
    predict(fm1, newdata=newData)
    , "terms object does not have a predvars attribute")

