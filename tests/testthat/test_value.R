context("Test numerical correctness of functions.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("Estimates are reasonable when data is Poisson", {
    counts <- matrix(rpois(10000, lambda=50), nrow=100, ncol=100)
    m1 <- zinbFit(counts, commondispersion = TRUE)
    expect_true(all(getPhi(m1) < 1e-4))

    m2 <- zinbFit(counts, commondispersion = FALSE)
    expect_true(all(getPhi(m2) < 1e-4))

    expect_equivalent(round(getMu(m1), 2), round(getMu(m2), 2))
    expect_true(abs(mean(getMu(m1)) - 50) < 1)
    expect_true(mean(getPi(m1)) < 1e-2)
})

test_that("Estimates are reasonable when data is Negative Binomial", {
    counts <- matrix(rnbinom(10000, mu=50, size = 10), nrow=100, ncol=100)

    m1 <- zinbFit(counts, commondispersion = TRUE)

    expect_true(abs(mean(getMu(m1)) - 50) < 1)
    expect_true(abs(mean(getTheta(m1)) - 10) < 1)
    expect_true(mean(getPi(m1)) < 1e-2)
})
