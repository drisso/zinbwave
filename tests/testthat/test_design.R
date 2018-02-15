context("Test zinbModel class and methods.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("getX et al work with/without intercept", {

    bio <- gl(2, 3)
    gc <- rnorm(10)
    m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio), V=model.matrix(~gc),
                 which_X_pi=1L, which_V_mu=1L)

    expect_equal(NCOL(getV_mu(m)), 1)
    expect_equal(NCOL(getV_mu(m, intercept=FALSE)), 0)
    expect_equal(NCOL(getV_pi(m)), 2)
    expect_equal(NCOL(getV_pi(m, intercept=FALSE)), 1)
    expect_equal(NCOL(getX_mu(m)), 2)
    expect_equal(NCOL(getX_mu(m, intercept=FALSE)), 1)
    expect_equal(NCOL(getX_pi(m)), 1)
    expect_equal(NCOL(getX_pi(m, intercept=FALSE)), 0)


    m <- zinbFit(matrix(10, 10, 6), X=model.matrix(~bio), V=model.matrix(~gc),
                 which_X_pi=2L, which_V_mu=2L)

    expect_equal(getV_mu(m, intercept=TRUE), getV_mu(m, intercept=TRUE))
    expect_equal(getX_pi(m, intercept=TRUE), getX_pi(m, intercept=TRUE))
})

test_that("zinbFit works with genewise dispersion", {
    bio <- gl(2, 3)
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = FALSE)

    m <- zinbFit(counts, X=model.matrix(~bio), verbose = TRUE)
})

test_that("zinbFit stops if one gene has only 0 counts", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    counts <- rbind(counts, rep(0, ncol(counts)))
    expect_error(zinbFit(counts), "only 0 counts")
})

test_that("zinbFit stops if one sample has only 0 counts", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    counts <- cbind(counts, rep(0, nrow(counts)))
    expect_error(zinbFit(counts), "only 0 counts")
})

test_that("zinbFit works without X and V", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m1 <- zinbFit(counts, X = matrix(0, ncol=1, nrow=ncol(counts)))
    m2 <- zinbFit(counts, V = matrix(0, ncol=1, nrow=nrow(counts)))
    m3 <- zinbFit(counts, X = matrix(0, ncol=1, nrow=ncol(counts)),
                  V = matrix(0, ncol=1, nrow=nrow(counts)))

    expect_equal(sum(as.vector(m1@beta_mu)), 0)
    expect_equal(sum(as.vector(m1@beta_pi)), 0)
    expect_equal(sum(as.vector(m2@gamma_mu)), 0)
    expect_equal(sum(as.vector(m2@gamma_pi)), 0)
    expect_equal(sum(as.vector(m3@beta_mu)), 0)
    expect_equal(sum(as.vector(m3@beta_pi)), 0)
    expect_equal(sum(as.vector(m3@gamma_mu)), 0)
    expect_equal(sum(as.vector(m3@gamma_pi)), 0)

})

test_that("zinbFit gives the same results with matrix and SE", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    se <- SummarizedExperiment(counts)

    m1 <- zinbFit(counts)
    m2 <- zinbFit(se)
    expect_equal(m1, m2)
})

test_that("zinbFit gives the same results with matrix and formula", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    bio <- gl(2, 3)
    gcc <- rnorm(10)
    se <- SummarizedExperiment(counts, colData=data.frame(Bio=bio),
                               rowData=data.frame(GCC=gcc))

    m1 <- zinbFit(se, X = model.matrix(~bio))
    m2 <- zinbFit(se, X = "~Bio")
    expect_equivalent(m1, m2)

    m3 <- zinbFit(se, V = model.matrix(~gcc))
    m4 <- zinbFit(se, V = "~GCC")
    expect_equivalent(m3, m4)

    # misstyping
    expect_error(zinbFit(se, V = "~gc"), "V must be a matrix or a formula")

    # colData / rowData missmatch
    expect_error(zinbFit(se, V = "~BIO"), "V must be a matrix or a formula")

})

test_that("zinbFit works with K>0", {
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m <- zinbFit(counts, K = 2)
    expect_equal(dim(getW(m)), c(nSamples(m), nFactors(m)))
})

test_that("zinbSim works", {
    a <- zinbModel(n=5, J=10)
    sim <- zinbSim(a)

    expect_true(all(.is_wholenumber(sim$counts)))
    expect_true(all(.is_wholenumber(sim$dataNB)))
    expect_true(all(.is_wholenumber(sim$dataDropouts)))
})

test_that("getMu and getPi have the right dimensions", {
    bio <- gl(2, 3)
    counts <- matrix(rpois(60, lambda=5), nrow=10, ncol=6)
    m <- zinbFit(counts, X=model.matrix(~bio), commondispersion = TRUE)

    expect_equal(dim(getMu(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getLogMu(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getPi(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getLogitPi(m)), c(nSamples(m), nFeatures(m)))
    expect_equal(dim(getW(m)), c(nSamples(m), nFactors(m)))
    expect_equal(length(getPhi(m)), nFeatures(m))
    expect_equal(length(getTheta(m)), nFeatures(m))
    expect_equal(length(getZeta(m)), nFeatures(m))
})

test_that("Initialization works", {

    ## no arguments specified
    zinbModel()

    ## specify W
    mat <- matrix(rnorm(10), ncol=2)
    m <- zinbModel(W = mat)
    expect_equal(nSamples(m), nrow(mat))

    ## specify X
    m <- zinbModel(X = mat)
    expect_equal(nSamples(m), nrow(mat))

    ## specify V
    m <- zinbModel(V = mat)
    expect_equal(nFeatures(m), nrow(mat))

    ## specify different X, V for pi and mu
    m <- zinbModel(X = mat, which_X_mu=1L, which_X_pi=2L,
              V = mat, which_V_mu=2L, which_V_pi=1L)
    expect_equal(nFeatures(m), nrow(mat))
    expect_equal(nSamples(m), nrow(mat))

    ## specify O_mu
    m <- zinbModel(O_mu = mat)
    expect_equal(nSamples(m), nrow(mat))
    expect_equal(nFeatures(m), ncol(mat))

    ## specify O_pi
    m <- zinbModel(O_pi = mat)
    expect_equal(nSamples(m), nrow(mat))
    expect_equal(nFeatures(m), ncol(mat))

    ## specify empty X
    m <- zinbModel(X = matrix(0, ncol=0, nrow=10))
    expect_equal(nSamples(m), 10)

    ## specify empty V
    m <- zinbModel(V = matrix(0, ncol=0, nrow=10))
    expect_equal(nFeatures(m), 10)

    ## specify empty X and V
    m <- zinbModel(X = matrix(0, ncol=0, nrow=10),
                   V = matrix(0, ncol=0, nrow=10))
    expect_equal(nSamples(m), 10)
    expect_equal(nFeatures(m), 10)

})
