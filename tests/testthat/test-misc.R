context("Test of general functionality")

test_that("citation reports year", {
    txt <- citation("ordinal")
    expect_true(as.logical(grep("year", txt)))
})

