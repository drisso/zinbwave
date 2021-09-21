context("Test simulation, initialization and epsilon.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("W can be estimated from random matrix (no signal)", {
    set.seed(789)

    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(rnorm(nS*K),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(rnorm(nG), ncol = 1)
    mm = zinbModel(X = X,
                   V = V,
                   W = W,
                   gamma_mu = matrix(1, ncol = nS),
                   beta_mu = matrix(5, ncol=nG))

    my_data <- zinbSim(mm)

    expect_equal(W, getW(mm))
    expect_equal(X, getX_mu(mm))
    expect_equal(X, getX_pi(mm))
    expect_equal(V, getV_mu(mm))
    expect_equal(V, getV_pi(mm))

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts),
                              rowDat = data.frame(V),
                              colDat = data.frame(X))
    sf <- zinbFit(SE,
                 V = "~V - 1",
                 K = 2)

    round(cor(cbind(W, sf@W)), 2)
})

test_that("W can be estimated from two-dimensional signal (no signal)", {
    set.seed(789)

    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(c(rnorm(50, mean=5, sd=.1), rnorm(50, mean=1, sd=.1),
                  rnorm(10, mean=5, sd=.1), rnorm(90, mean=1, sd=.1)),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(rnorm(nG), ncol = 1)
    mm = zinbModel(X = X,
                   V = V,
                   W = W,
                   gamma_mu = matrix(1, ncol = nS),
                   beta_mu = matrix(5, ncol=nG))

    my_data <- zinbSim(mm)

    expect_equal(W, getW(mm))
    expect_equal(X, getX_mu(mm))
    expect_equal(X, getX_pi(mm))
    expect_equal(V, getV_mu(mm))
    expect_equal(V, getV_pi(mm))

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts),
                               rowDat = data.frame(V),
                               colDat = data.frame(X))
    sf = zinbFit(SE,
                 V = "~V - 1",
                 K = 2)

    round(cor(cbind(W, sf@W)), 2)
})

test_that("Initialization works with large epsilon", {
    set.seed(123)

    nS <- 100 # sample size
    nG <- 50  # number genes
    K  <- 2   # number of latent factors
    W <- matrix(c(rnorm(50, mean=5, sd=.1), rnorm(50, mean=1, sd=.1),
                  rnorm(10, mean=5, sd=.1), rnorm(90, mean=1, sd=.1)),ncol = K)
    X <- matrix(1, ncol = 1, nrow = nS)
    V <- matrix(1, ncol = 1, nrow = nG)
    mm = zinbModel(X = X,
                   V = V,
                   W = W,
                   gamma_mu = matrix(1, ncol = nS),
                   beta_mu = matrix(5, ncol=nG))

    my_data <- zinbSim(mm)

    SE <- SummarizedExperiment(assays = list(counts=my_data$counts))

    expect_silent(sf <- zinbFit(SE, K = 2, epsilon = 1e4))
    expect_silent(sf <- zinbFit(SE, K = 2, epsilon = 1e5))
    expect_silent(sf <- zinbFit(SE, K = 2, epsilon = 1e12))


})
