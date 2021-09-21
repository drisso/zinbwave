context("Test zinbModel class accessors.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("dimensions", {
    set.seed(123)
    x <- zinbModel()
    expect_equal(nSamples(x), NROW(getX_mu(x)))
    expect_equal(nFeatures(x), NROW(getV_mu(x)))
    expect_equal(nFactors(x), NCOL(getW(x)))
})

test_that("getters return the promised values", {
    set.seed(123)
    x <- zinbModel()
    expect_equal(getX_mu(x), x@X[, x@which_X_mu, drop=FALSE])
    expect_equal(getX_pi(x), x@X[, x@which_X_pi, drop=FALSE])
    expect_equal(getV_mu(x), x@V[, x@which_V_mu, drop=FALSE])
    expect_equal(getV_pi(x), x@V[, x@which_V_pi, drop=FALSE])
    expect_equal(getW(x), x@W)

    expect_equal(getBeta_mu(x), x@beta_mu)
    expect_equal(getBeta_pi(x), x@beta_pi)
    expect_equal(getGamma_mu(x), x@gamma_mu)
    expect_equal(getGamma_pi(x), x@gamma_pi)
    expect_equal(getAlpha_mu(x), x@alpha_mu)
    expect_equal(getAlpha_pi(x), x@alpha_pi)

    expect_equal(getZeta(x), x@zeta)
    expect_equal(getTheta(x), exp(x@zeta))
    expect_equal(getPhi(x), exp(-x@zeta))
})
