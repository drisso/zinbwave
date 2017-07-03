context("Test zinbwave function.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("zinbwave gives same result with / without model fit", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                                colData = data.frame(bio = gl(2, 3)))

    m1 <- zinbwave(se, X="~bio")

    fit <- zinbFit(se, X="~bio")
    m2 <- zinbwave(se, fitted_model=fit)

    expect_equal(m1, m2)
})
