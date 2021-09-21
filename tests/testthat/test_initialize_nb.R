context("Test optimization without zero inflation.")
set.seed(13124)

BiocParallel::register(BiocParallel::SerialParam())

test_that("Initialization works without zero inflation", {
    set.seed(789)

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

    expect_silent(sf <- zinbFit(SE, K = 2, zeroinflation = FALSE))
    expect_silent(sf <- zinbFit(SE, K = 2, zeroinflation = FALSE))
    expect_silent(sf <- zinbFit(SE, K = 2, zeroinflation = FALSE))

})
