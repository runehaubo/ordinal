context("Contrast specification")

test_that("clm gives contrast warnings when it should", {
    ## No warnings:
    ## Different combinations of terms i various formulae. Note that the
    ## contrasts apply to e.g. 'contact' in both 'formula' and 'scale':
    contr <- c(temp="contr.sum", contact="contr.sum")
    expect_false(givesWarnings(
        fm1 <- clm(rating ~ temp + contact, scale=~contact, data=wine) ## OK
        ))
    expect_false(givesWarnings(
        fm1 <- clm(rating ~ temp + contact, scale=~contact, data=wine,
                   contrasts=contr) ## OK
        ))
    expect_false(givesWarnings(
        fm1 <- clm(rating ~ temp, scale=~contact, data=wine,
                   contrasts=contr) ## OK
        ))
    expect_false(givesWarnings(
        fm1 <- clm(rating ~ temp, nominal=~contact, data=wine,
                   contrasts=contr) ## OK
        ))
    expect_false(givesWarnings(
        fm1 <- clm(rating~1, scale=~temp, nominal=~contact, data=wine,
                   contrasts=contr) ## OK
        ))

    ## These should give warnings:
    ## A warning is given if a variable is not present in any of the
    ## formulae:
    expect_warning(
        fm <- clm(rating ~ temp, contrasts=c(contact="contr.sum"), data=wine)
        , "variable 'contact' is absent: its contrasts will be ignored")
    expect_warning(
        fm <- clm(rating ~ temp, contrasts=contr, data=wine)
        , "variable 'contact' is absent: its contrasts will be ignored")
    expect_warning(
        fm <- clm(rating ~ 1, scale=~contact, contrasts=c(temp="contr.sum"),
                  data=wine)
        , "variable 'temp' is absent: its contrasts will be ignored")
    expect_warning(
        fm <- clm(rating ~ 1, scale=~contact, contrasts=list(temp="contr.sum"),
                  data=wine)
        , "variable 'temp' is absent: its contrasts will be ignored")

})

test_that("checkContrasts gives when it should", {
    ## No warnings:
    fm0 <- clm(rating ~ temp + contact, scale=~contact, data=wine)
    expect_false(
        givesWarnings(checkContrasts(fm0$S.terms, fm0$S.contrasts))
        )
    expect_false(
        givesWarnings(checkContrasts(fm0$terms, fm0$contrasts))
        )
    expect_false(
        givesWarnings(checkContrasts(fm0$terms, fm0$S.contrasts))
        )
    expect_false(
        givesWarnings(checkContrasts(fm0$terms, fm0$S.contrasts))
        )
    ## Warning:
    expect_warning(
        checkContrasts(fm0$S.terms, fm0$contrasts)
        , "variable 'temp' is absent: its contrasts will be ignored")
})


