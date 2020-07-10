context("Test DE functions.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("glm with weights works as expected", {
    se <- SummarizedExperiment(assays = list(counts = matrix(rpois(60, lambda=5), nrow=10, ncol=6)),
                               colData = data.frame(bio = gl(2, 3)))

    m1 <- zinbwave(se, K = 0, observationalWeights = TRUE)

    expect_true("weights" %in% names(assays(m1)))
    expect_true(all(assay(m1, "weights") > 0))
    expect_true(all(assay(m1, "weights") <= 1))

    library(edgeR)

    dge <- DGEList(counts(m1))
    dge <- calcNormFactors(dge)
    x <- factor(rep(1:2, each=3))
    design <- model.matrix(~x)
    dge$weights <- assay(m1, "weights")
    dge <- estimateDisp(dge, design)
    fit <- glmFit(dge, design)

    lrt <- glmWeightedF(fit, coef = 2)
    lrt2 <- glmWeightedF(fit, contrast = matrix(c(0, 1)))
    expect_equal(lrt$table, lrt2$table)

    lrt3 <- glmWeightedF(fit, ZI=FALSE)
    lrt4 <- glmWeightedF(fit, independentFiltering = FALSE)
    lrt5 <- glmWeightedF(fit, filter = rowMeans(counts(m1)))
    expect_equal(lrt$table, lrt5$table)
})
