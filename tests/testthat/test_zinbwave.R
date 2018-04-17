context("Test zinbwave function.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("zinbwave gives same result with / without model fit", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                                colData = data.frame(bio = gl(2, 3)))

    expect_warning(m1 <- zinbwave(se, X="~bio"), "No assay named `counts`")

    fit <- zinbFit(se, X="~bio")
    expect_warning(m2 <- zinbwave(se, fitted_model=fit), "No assay named `counts`")

    expect_equal(m1, m2)
    expect_is(m1, "SingleCellExperiment")
    expect_is(m2, "SingleCellExperiment")
})

test_that("W is the same in zinbFit and zinbwave", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    fit <- zinbFit(se, K = 2)
    expect_warning(m1 <- zinbwave(se, fitted_model=fit), "No assay named `counts`")

    expect_equivalent(getW(fit), reducedDim(m1))
})

test_that("zinbwave computes residuals and normalized values", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    fit <- zinbFit(se, K = 2)
    expect_warning(m1 <- zinbwave(se, fitted_model = fit, residuals = TRUE,
                   normalizedValues = TRUE, imputedValues = TRUE),
                   "No assay named `counts`")

    expect_true("normalizedValues" %in% names(assays(m1)))
    expect_true("residuals" %in% names(assays(m1)))
    expect_true("imputedValues" %in% names(assays(m1)))
})

test_that("zinbwave works with slot counts", {

    cc <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    ll <- matrix(rnorm(60), nrow=10, ncol=6)

    se <- SummarizedExperiment(assays = list(counts = cc, norm = ll),
                               colData = data.frame(bio = gl(2, 3)))

    expect_silent(m1 <- zinbwave(se, K = 0))
    expect_silent(m2 <- zinbwave(se, K = 0, which_assay = "counts"))
    expect_equal(m1, m2)

    se <- SummarizedExperiment(assays = list(norm = ll, counts = cc),
                               colData = data.frame(bio = gl(2, 3)))

    expect_silent(m1 <- zinbwave(se, K = 0))
    expect_silent(m2 <- zinbwave(se, K = 0, which_assay = "counts"))
    expect_equal(m1, m2)

})


test_that("zinbwave works without slot counts", {

    cc <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    ll <- matrix(rnorm(60), nrow=10, ncol=6)

    se <- SummarizedExperiment(assays = list(assay1 = cc, assay2 = ll),
                               colData = data.frame(bio = gl(2, 3)))

    expect_warning(m1 <- zinbwave(se, K = 0))
    expect_silent(m2 <- zinbwave(se, K = 0, which_assay = "assay1"))
    expect_equal(m1, m2)

    se <- SummarizedExperiment(assays = list(assay1 = ll, assay2 = cc),
                               colData = data.frame(bio = gl(2, 3)))

    expect_error(m1 <- zinbwave(se, K = 0), "The input matrix should contain only whole numbers")
    expect_silent(m2 <- zinbwave(se, K = 0, which_assay = "assay2"))

})
