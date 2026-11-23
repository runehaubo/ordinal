test_that("minimum Gumbel reflects around its location", {
    location <- 10
    scale <- 2
    x <- c(7, 10, 11)
    reflected <- 2 * location - x

    for (lower.tail in c(TRUE, FALSE)) {
        expect_equal(pgumbel(x, location, scale, lower.tail, max = FALSE),
                     pgumbel(reflected, location, scale, !lower.tail))
    }
    for (log in c(TRUE, FALSE)) {
        expect_equal(dgumbel(x, location, scale, log, max = FALSE),
                     dgumbel(reflected, location, scale, log))
    }
})

test_that("minimum Gumbel preserves small tail probabilities", {
    location <- -5
    scale <- 2
    expect_equal(pgumbel(location - 40 * scale, location, scale, max = FALSE) /
                     pexp(exp(-40)), 1)
    expect_equal(pgumbel(location + 4 * scale, location, scale,
                         lower.tail = FALSE, max = FALSE) /
                     pexp(exp(4), lower.tail = FALSE), 1)

    x <- c(-Inf, Inf, NA_real_)
    expect_equal(pgumbel(x, location, scale, max = FALSE), c(0, 1, NA_real_))
    expect_equal(pgumbel(x, location, scale, lower.tail = FALSE, max = FALSE),
                 c(1, 0, NA_real_))
})
