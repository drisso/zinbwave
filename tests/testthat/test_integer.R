context("Test zinbModel class accessors.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("Non-integer numbers are caught", {
    set.seed(987)

    se <- SummarizedExperiment(matrix(rnorm(60, mean=5), nrow=10, ncol=6),
                               colData = data.frame(bio = gl(2, 3)))

    expect_error(m1 <- zinbwave(se), "The input matrix should contain only whole numbers.")
    expect_error(m2 <- zinbwave(se), "The input matrix should contain only whole numbers.")
})
