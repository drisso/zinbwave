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
    expect_is(m1, "SingleCellExperiment")
    expect_is(m2, "SingleCellExperiment")
})

test_that("W is the same in zinbFit and zinbwave", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    fit <- zinbFit(se, K = 2)
    m1 <- zinbwave(se, fitted_model=fit)

    expect_equivalent(getW(fit), reducedDim(m1))
})

test_that("zinbwave computes residuals and normalized values", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    fit <- zinbFit(se, K = 2)
    m1 <- zinbwave(se, fitted_model = fit, residuals = TRUE,
                   normalizedValues = TRUE, imputedValues = TRUE)

    expect_true("normalizedValues" %in% names(assays(m1)))
    expect_true("residuals" %in% names(assays(m1)))
    expect_true("imputedValues" %in% names(assays(m1)))
})
