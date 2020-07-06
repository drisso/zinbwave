nbInitialize <- function(m, Y, nb.repeat=2,  it.max = 100,
                           BPPARAM=BiocParallel::bpparam()) {

  n <- NROW(Y)
  J <- NCOL(Y)

  if(n != nSamples(m)) {
    stop("Y needs to have ", nSamples(m), " rows (genes)")
  }

  if(J != nFeatures(m)) {
    stop("Y needs to have ", nFeatures(m), " columns (samples)")
  }

  ## 1. Define P
  P <- Y > 0

  if(any(rowSums(P) == 0)) {
    stop("Sample ", which(rowSums(P) == 0)[1], " has only 0 counts!")
  }

  if(any(colSums(P) == 0)) {
    stop("Gene ", which(colSums(P) == 0)[1], " has only 0 counts!")
  }

  ## 2. Define L
  L <- matrix(NA, nrow=n, ncol=J)
  L <- log1p(Y) - m@O_mu

  ## 3. Estimate gamma_mu and beta_mu
  iter <- 0
  beta_mu <- getBeta_mu(m)
  gamma_mu <- getGamma_mu(m)

  while (iter < nb.repeat) {

    # Optimize gamma_mu (in parallel for each sample)
    if (NCOL(getV_mu(m)) == 0) {
      iter <- nb.repeat # no need to estimate gamma_mu nor to iterate
    } else {
      Xbeta_mu <- getX_mu(m) %*% beta_mu
      gamma_mu <- matrix(unlist(bplapply(seq(n), function(i) {
        solveRidgeRegression(x=getV_mu(m),
                             y=L[i,] - Xbeta_mu[i,],
                             epsilon = getEpsilon_gamma_mu(m),
                             beta = gamma_mu[,i],
                             family="gaussian")
      } , BPPARAM=BPPARAM
      )), nrow=NCOL(getV_mu(m)))
    }

    # Optimize beta_mu (in parallel for each gene)
    if (NCOL(getX_mu(m)) == 0) {
      iter <- nb.repeat # no need to estimate gamma_mu nor to iterate
    } else {
      tVgamma_mu <- t(getV_mu(m) %*% gamma_mu)
      beta_mu <- matrix(unlist(bplapply(seq(J), function(j) {
        solveRidgeRegression(x=getX_mu(m),
                             y=L[,j] - tVgamma_mu[, j],
                             epsilon = getEpsilon_beta_mu(m),
                             beta = beta_mu[,j],
                             family="gaussian")
      }, BPPARAM=BPPARAM
      )), nrow=NCOL(getX_mu(m)))
    }

    iter <- iter+1
  }

  ## 4. Estimate W and alpha (only if K>0)
  if(nFactors(m) > 0) {

    # Compute the residual D
    D <- L - getX_mu(m) %*% beta_mu - t(getV_mu(m) %*% gamma_mu)

    # Find a low-rank approximation with trace-norm regularization
    lambda <- sqrt(getEpsilon_W(m) * getEpsilon_alpha(m))[1]
    R <- svd(D, nu=nFactors(m), nv=nFactors(m))


    # Orthogonalize to get W and alpha
    W <- (getEpsilon_alpha(m) / getEpsilon_W(m))[1]^(1/4) *
      R$u %*% diag(sqrt(R$d[1:nFactors(m)]), nrow = length(R$d[1:nFactors(m)]))
    alpha_mu <- (getEpsilon_W(m)/getEpsilon_alpha(m))[1]^(1/4) *
      diag(sqrt(R$d[1:nFactors(m)]),nrow = length(R$d[1:nFactors(m)])) %*% t(R$v)
  } else {
    W <- getW(m)
    alpha_mu <- getAlpha_mu(m)
  }


  ## 5. Initialize dispersion to 1
  zeta <- rep(0, J)

  out <-  zinbModel(X = m@X, V = m@V, O_mu = m@O_mu,
                   which_X_mu = m@which_X_mu,
                   which_V_mu = m@which_V_mu,
                   W = W, beta_mu = beta_mu,
                   gamma_mu = gamma_mu,
                   alpha_mu = alpha_mu, zeta = zeta,
                   epsilon_beta_mu = m@epsilon_beta_mu,
                   epsilon_gamma_mu = m@epsilon_gamma_mu,
                   epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                   epsilon_zeta = m@epsilon_zeta)

  return(out)
}
