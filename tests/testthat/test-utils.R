
context("testing namedList")

a <- 1
b <- 2
c <- 3
d <- list(e=2, f=factor(letters[rep(1:2, 2)]))
g <- matrix(runif(9), 3)
h <- NULL

test_that("namedList returns a named list", {

    res <- namedList(a, b, c)
    expect_equal(names(res), c("a", "b", "c"))
    expect_equivalent(res, list(a, b, c))

    res <- namedList(a, b, c, d, g)
    expect_equal(names(res), c("a", "b", "c", "d", "g"))
    expect_equivalent(res, list(a, b, c, d, g))

    res <- namedList(a, h)
    expect_equal(names(res), c("a", "h"))
    expect_equivalent(res, list(a, h))
})
