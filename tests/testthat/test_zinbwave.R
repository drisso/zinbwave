context("Test zinbwave function.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("zinbwave gives same result with / without model fit", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                                colData = data.frame(bio = gl(2, 3)))

    expect_warning(m1 <- zinbwave(se, X="~bio", K=0), "No assay named `counts`")

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

test_that("zinbwave computes observational weihts", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    fit <- zinbFit(se, K = 2)
    expect_warning(m1 <- zinbwave(se, fitted_model = fit,
                                  observationalWeights = TRUE),
                   "No assay named `counts`")

    expect_true("weights" %in% names(assays(m1)))
    expect_true(all(assay(m1, "weights") > 0))
    expect_true(all(assay(m1, "weights") <= 1))

    w <- computeObservationalWeights(fit, assay(se))
    expect_equivalent(w, assay(m1, "weights"))
})

test_that("one-dimensional W", {
    se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    expect_silent(fit <- zinbwave(se, K = 1, which_assay = 1))

    expect_equal(NCOL(reducedDim(fit)), 1)
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

test_that("zinbwave works with subset of genes", {

    cc <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    rownames(cc) <- paste0("gene", 1:10)
    wh_genes <- c(rep(TRUE, 2), rep(FALSE, 8))
    se <- SummarizedExperiment(assays = list(counts = cc),
                               colData = data.frame(bio = gl(2, 3)),
                               rowData = data.frame(wh_genes = wh_genes))

    ## check that it works with all genes
    set.seed(123)
    expect_silent(m1 <- zinbwave(se, K=1))
    set.seed(123)
    expect_silent(m2 <- zinbwave(se, K=1, which_genes = rownames(cc)))
    set.seed(123)
    expect_silent(m3 <- zinbwave(se, K=1, which_genes = 1:10))
    set.seed(123)
    expect_silent(m4 <- zinbwave(se, K=1, which_genes = rep(TRUE, 10)))

    expect_equal(m1, m2)
    expect_equal(m1, m3)
    expect_equal(m1, m4)

    ## check that it works with both a vector and a rowData column
    set.seed(155)
    expect_silent(m1 <- zinbwave(se, K=1, which_genes = wh_genes))
    set.seed(155)
    expect_silent(m2 <- zinbwave(se, K=1, which_genes = "wh_genes"))

    expect_equal(m1, m2)

    ## check with no refit
    set.seed(155)
    expect_silent(m1 <- zinbwave(se, K=1, which_genes = wh_genes,
                                 observationalWeights = FALSE))
    set.seed(155)
    expect_silent(m2 <- zinbwave(se, K=1, which_genes = wh_genes))

    expect_equal(reducedDims(m1), reducedDims(m2))

    ## check do nothing
    m0 <- zinbwave(se, K=0, observationalWeights=FALSE)
    expect_equal(as(se, "SingleCellExperiment"), m0)

    ## check with wrong genes
    expect_error(m1 <- zinbwave(se, K=1, which_genes = c("gene1", "gene11")),
                 "index out of bounds")

    expect_error(m1 <- zinbwave(se, K=1, which_genes = "my_genes"),
                 "it must be the name of a column")

})
