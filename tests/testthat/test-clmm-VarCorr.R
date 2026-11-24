test_that("clmm labels random effects with equal group counts", {
    set.seed(831)
    dat <- expand.grid(id = factor(1:6), item = factor(1:6),
                       x = c(-1, 1), replicate = 1:5)
    latent <- 0.4 * dat$x +
        seq(-1.5, 1.5, length.out = 6)[dat$id] +
        seq(-1, 1, length.out = 6)[dat$item] * dat$x +
        rlogis(nrow(dat))
    dat$y <- cut(latent, breaks = c(-Inf, -1, 0, 1, Inf),
                 ordered_result = TRUE)

    fit <- clmm(y ~ x + (1 | id) + (0 + x | item), data = dat)
    vc <- VarCorr(fit)
    expect_identical(names(vc), c("id", "item"))
    expect_identical(colnames(vc$id), "(Intercept)")
    expect_identical(colnames(vc$item), "x")
})
