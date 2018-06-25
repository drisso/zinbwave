#' @describeIn zinbFit Y is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @export
#'
#' @import SummarizedExperiment
#' @importFrom stats model.matrix as.formula
#'
#' @param which_assay numeric or character. Which assay of Y to use (only if Y
#'   is a SummarizedExperiment).
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#'
#' m <- zinbFit(se, X=model.matrix(~bio, data=colData(se)))
setMethod("zinbFit", "SummarizedExperiment",
          function(Y, X, V, K, which_assay,
                   commondispersion=TRUE, verbose=FALSE,
                   nb.repeat.initialize=2, maxiter.optimize=25,
                   stop.epsilon.optimize=.0001,
                   BPPARAM=BiocParallel::bpparam(), ...) {

              if(missing(which_assay)) {
                  if("counts" %in% assayNames(Y)) {
                      dataY <- assay(Y, "counts")
                  } else {
                      dataY <- assay(Y)
                  }
              } else {
                  if(!(is.character(which_assay) | is.numeric(which_assay))) {
                      stop("assay needs to be a numeric or character specifying which assay to use")
                  } else {
                      dataY <- assay(Y, which_assay)
                  }
              }

              if(!missing(X)) {
                  if(!is.matrix(X)) {
                      tryCatch({
                          f <- as.formula(X)
                          X <- model.matrix(f, data=colData(Y))
                      },
                      error = function(e) {
                          stop("X must be a matrix or a formula with variables in colData(Y)")
                      })
                  }
              }

              if(!missing(V)) {
                  if(!is.matrix(V)) {
                      tryCatch({
                          f <- as.formula(V)
                          V <- model.matrix(f, data=rowData(Y))
                      },
                      error = function(e) {
                          stop("V must be a matrix or a formula with variables in rowData(Y)")
                      })
                  }
              }

              # Apply zinbFit on the assay of SummarizedExperiment
              res <- zinbFit(dataY, X, V, K, commondispersion,
                             verbose, nb.repeat.initialize, maxiter.optimize,
                             stop.epsilon.optimize, BPPARAM, ...)

              return(res)
})


#' @describeIn zinbFit Y is a matrix of counts (genes in rows).
#' @export
#'
#' @param X The design matrix containing sample-level covariates, one sample per
#'   row. If missing, X will contain only an intercept. If Y is a
#'   SummarizedExperiment object, X can be a formula using the variables in the
#'   colData slot of Y.
#' @param V The design matrix containing gene-level covariates, one gene
#'   per row. If missing, V will contain only an intercept. If Y is a
#'   SummarizedExperiment object, V can be a formula using the variables in the
#'   rowData slot of Y.
#' @param K integer. Number of latent factors.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE).
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @param verbose Print helpful messages.
#' @param nb.repeat.initialize Number of iterations for the initialization of
#'   beta_mu and gamma_mu.
#' @param maxiter.optimize maximum number of iterations for the optimization
#'   step (default 25).
#' @param stop.epsilon.optimize stopping criterion in the optimization step,
#'   when the relative gain in likelihood is below epsilon (default 0.0001).
#'
#' @details By default, i.e., if no arguments other than \code{Y} are passed,
#'   the model is fitted with an intercept for the regression across-samples and
#'   one intercept for the regression across genes, both for mu and for pi.
#'
#' @details This means that by default the model is fitted with \code{X_mu =
#'   X_pi = 1_n} and \code{V_mu = V_pi = 1_J}. If the user explicitly passes the
#'   design matrices, this behavior is overwritten, i.e., the user needs to
#'   explicitly include the intercept in the design matrices.
#'
#' @details If Y is a Summarized experiment, the function uses the assay named
#'   "counts", if any, or the first assay.
#'
#' @seealso \code{\link[stats]{model.matrix}}.
#'
#' @import BiocParallel
#' @examples
#' bio <- gl(2, 3)
#' m <- zinbFit(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'              X=model.matrix(~bio))
setMethod("zinbFit", "matrix",
          function(Y, X, V, K, commondispersion=TRUE, verbose=FALSE,
                   nb.repeat.initialize=2, maxiter.optimize=25,
                   stop.epsilon.optimize=.0001,
                   BPPARAM=BiocParallel::bpparam(), ...) {

    # Check that Y contains whole numbers
    if(!all(.is_wholenumber(Y))) {
        stop("The input matrix should contain only whole numbers.")
    }

    # Transpose Y: UI wants genes in rows, internals genes in columns!
    Y <- t(Y)

    # Create a ZinbModel object
    if (verbose) {message("Create model:")}
    m <- zinbModel(n=NROW(Y), J=NCOL(Y), X=X, V=V, K=K, ...)
    if (verbose) {message("ok")}

    # Initialize the parameters
    if (verbose) {message("Initialize parameters:")}
    m <- zinbInitialize(m, Y, nb.repeat=nb.repeat.initialize, BPPARAM=BPPARAM)
    if (verbose) {message("ok")}

    # Optimize parameters
    if (verbose) {message("Optimize parameters:")}
    m <- zinbOptimize(m, Y, commondispersion=commondispersion,
                      maxiter=maxiter.optimize,
                      stop.epsilon=stop.epsilon.optimize,
                      BPPARAM=BPPARAM, verbose=verbose)
    if (verbose) {message("ok")}

    validObject(m)
    m
})

#' @describeIn zinbFit Y is a sparse matrix of counts (genes in rows).
#' @export
#'
#' @details Currently, if Y is a sparseMatrix, this calls the zinbFit method on
#'   as.matrix(Y)
#'
#' @importClassesFrom Matrix dgCMatrix
setMethod("zinbFit", "dgCMatrix",
          function(Y, ...) {
              zinbFit(as.matrix(Y), ...)
})


#' Initialize the parameters of a ZINB regression model
#'
#' The initialization performs quick optimization of the parameters with several
#' simplifying assumptions compared to the true model: non-zero counts are
#' models as log-Gaussian, zeros are modeled as dropouts. The dispersion
#' parameter is not modified.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @param nb.repeat Number of iterations for the estimation of beta_mu and
#'   gamma_mu.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu,
#'   gamma_pi, W.
#' @examples
#' Y <- matrix(rpois(60, lambda=2), 6, 10)
#' bio <- gl(2, 3)
#' time <- rnorm(6)
#' gc <- rnorm(10)
#' m <- zinbModel(Y, X=model.matrix(~bio + time), V=model.matrix(~gc),
#'              which_X_pi=1L, which_V_mu=1L, K=1)
#' m <- zinbInitialize(m, Y)
#' @export
#' @importFrom glmnet glmnet
#' @importFrom softImpute softImpute
zinbInitialize <- function(m, Y, nb.repeat=2, BPPARAM=BiocParallel::bpparam()) {

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
    L[P] <- log(Y[P]) - m@O_mu[P]

    ## 3. Define Z
    Z <- 1 - P

    ## 4. Estimate gamma_mu and beta_mu
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
                solveRidgeRegression(x=getV_mu(m)[P[i,], , drop=FALSE],
                                     y=L[i,P[i,]] - Xbeta_mu[i, P[i,]],
                                     epsilon = getEpsilon_gamma_mu(m),
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
                solveRidgeRegression(x=getX_mu(m)[P[,j], , drop=FALSE],
                                     y=L[P[,j],j] - tVgamma_mu[P[,j], j],
                                     epsilon = getEpsilon_beta_mu(m),
                                     family="gaussian")
            }, BPPARAM=BPPARAM
            )), nrow=NCOL(getX_mu(m)))
        }

        iter <- iter+1
    }

    ## 5. Estimate W and alpha (only if K>0)
    if(nFactors(m) > 0) {

        # Compute the residual D (with missing values at the 0 count)
        D <- L - getX_mu(m) %*% beta_mu - t(getV_mu(m) %*% gamma_mu)

        # Find a low-rank approximation with trace-norm regularization
        R <- softImpute::softImpute(D,
                lambda=sqrt(getEpsilon_W(m) * getEpsilon_alpha(m))[1],
                rank.max=nFactors(m))

        # Orthogonalize to get W and alpha
        W <- (getEpsilon_alpha(m) / getEpsilon_W(m))[1]^(1/4) *
            R$u %*% diag(sqrt(R$d), nrow = length(R$d))
        alpha_mu <- (getEpsilon_W(m)/getEpsilon_alpha(m))[1]^(1/4) *
            diag(sqrt(R$d), nrow = length(R$d)) %*% t(R$v)
    } else {
        W <- getW(m)
        alpha_mu <- getAlpha_mu(m)
    }

    ## 6. Estimate beta_pi, gamma_pi, and alpha_pi
    iter <- 0
    beta_pi <- getBeta_pi(m)
    gamma_pi <- getGamma_pi(m)
    alpha_pi <- getAlpha_pi(m)

    while (iter < nb.repeat) {

        # Optimize gamma_pi (in parallel for each sample)
        if (NCOL(getV_pi(m)) == 0) {
            iter <- nb.repeat # no need to estimate gamma_pi nor to iterate
        } else {
            off <- getX_pi(m) %*% beta_pi + W %*% alpha_pi
            gamma_pi <- matrix(unlist(bplapply(seq(n), function(i) {
                solveRidgeRegression(x=getV_pi(m),
                                     y=Z[i,],
                                     offset=off[i,],
                                     epsilon = getEpsilon_gamma_pi(m),
                                     family="binomial")
            }, BPPARAM=BPPARAM
            )), nrow=NCOL(getV_pi(m)))
        }

        # Optimize beta_pi and alpha_pi (in parallel for each gene)
        if (NCOL(getX_pi(m)) + nFactors(m) == 0) {
            iter <- nb.repeat # no need to estimate nor to iterate
        } else {
            tVgamma_pi <- t(getV_pi(m) %*% gamma_pi)
            XW <- cbind(getX_pi(m),W)
            s <- matrix(unlist(bplapply(seq(J), function(j) {
                solveRidgeRegression(x=XW,
                                     y=Z[,j],
                                     offset = tVgamma_pi[,j],
                                     epsilon = c(getEpsilon_beta_pi(m),
                                                 getEpsilon_alpha(m)),
                                     family="binomial")
            }, BPPARAM=BPPARAM
            )), nrow=NCOL(getX_pi(m)) + nFactors(m))
            if (NCOL(getX_pi(m))>0) {
                beta_pi <- s[1:NCOL(getX_pi(m)),,drop=FALSE]
            }
            if (nFactors(m)>0) {
                alpha_pi <-
                    s[(NCOL(getX_pi(m)) + 1):(NCOL(getX_pi(m)) + nFactors(m)),,
                      drop=FALSE]
            }
        }

        iter <- iter+1
    }

    ## 7. Initialize dispersion to 1
    zeta <- rep(0, J)

    out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu, O_pi = m@O_pi,
                     which_X_mu = m@which_X_mu, which_X_pi = m@which_X_pi,
                     which_V_mu = m@which_V_mu, which_V_pi = m@which_V_pi,
                     W = W, beta_mu = beta_mu, beta_pi = beta_pi,
                     gamma_mu = gamma_mu, gamma_pi = gamma_pi,
                     alpha_mu = alpha_mu, alpha_pi = alpha_pi, zeta = zeta,
                     epsilon_beta_mu = m@epsilon_beta_mu,
                     epsilon_gamma_mu = m@epsilon_gamma_mu,
                     epsilon_beta_pi = m@epsilon_beta_pi,
                     epsilon_gamma_pi = m@epsilon_gamma_pi,
                     epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                     epsilon_zeta = m@epsilon_zeta,
                     epsilon_min_logit = m@epsilon_min_logit)

    return(out)
}

#' Optimize the parameters of a ZINB regression model
#'
#' The parameters of the model given as argument are optimized by penalized
#' maximum likelihood on the count matrix given as argument. It is recommended
#' to call zinb_initialize before this function to have good starting point for
#' optimization, since the optimization problem is not convex and can only
#' converge to a local minimum.
#' @param m The model of class ZinbModel
#' @param Y The matrix of counts.
#' @param commondispersion Whether the dispersion is the same for all features
#'   (default=TRUE)
#' @param maxiter maximum number of iterations (default 25)
#' @param stop.epsilon stopping criterion, when the relative gain in
#'   likelihood is below epsilon (default 0.0001)
#' @param verbose print information (default FALSE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters alpha_mu, alpha_pi, beta_mu, beta_pi, gamma_mu,
#'   gamma_pi, W.
#' @examples
#' Y = matrix(10, 3, 5)
#' m = zinbModel(n=NROW(Y), J=NCOL(Y))
#' m = zinbInitialize(m, Y)
#' m = zinbOptimize(m, Y)
#' @export
zinbOptimize <- function(m, Y, commondispersion=TRUE, maxiter=25,
                         stop.epsilon=.0001, verbose=FALSE,
                         BPPARAM=BiocParallel::bpparam()) {

    total.lik=rep(NA,maxiter)
    n <- nSamples(m)
    J <- nFeatures(m)

    epsilonright <- c(getEpsilon_beta_mu(m), getEpsilon_alpha(m),
                      getEpsilon_beta_pi(m), getEpsilon_alpha(m))
    nright <- c(length(getEpsilon_beta_mu(m)), length(getEpsilon_alpha(m)),
                length(getEpsilon_beta_pi(m)), length(getEpsilon_alpha(m)))
    optimright = (sum(nright)>0)

    epsilonleft <- c(getEpsilon_gamma_mu(m),
                     getEpsilon_gamma_pi(m), getEpsilon_W(m))
    nleft <- c(length(getEpsilon_gamma_mu(m)),
               length(getEpsilon_gamma_pi(m)), length(getEpsilon_W(m)))
    optimleft = (sum(nleft)>0)

    orthog <- (nFactors(m)>0)

    # extract fixed quantities from m
    X_mu <- getX_mu(m)
    V_mu <- getV_mu(m)
    X_pi <- getX_pi(m)
    V_pi <- getV_pi(m)
    O_mu <- m@O_mu
    O_pi <- m@O_pi

    # exctract paramters from m (remember to update!)
    beta_mu <- getBeta_mu(m)
    alpha_mu <- getAlpha_mu(m)
    gamma_mu <- getGamma_mu(m)
    beta_pi <- getBeta_pi(m)
    alpha_pi <- getAlpha_pi(m)
    gamma_pi <- getGamma_pi(m)
    W <- getW(m)
    zeta <- getZeta(m)

    for (iter in seq_len(maxiter)){
        if (verbose) {message("Iteration ",iter)}

        # Evaluate total penalized likelihood
        mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                      W %*% alpha_mu + O_mu)

        logitPi <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
            W %*% alpha_pi + O_pi

        theta <- exp(zeta)

        loglik <- zinb.loglik(Y, mu, rep(theta, rep(n, J)), logitPi)

        penalty <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
            sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
            sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
            sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
            sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
            sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
            sum(getEpsilon_W(m)*t(W)^2)/2 +
            getEpsilon_zeta(m)*var(zeta)/2

        total.lik[iter] <- loglik - penalty

        if (verbose) {message("penalized log-likelihood = ",
                          total.lik[iter])}

        # If the increase in likelihood is smaller than 0.5%, stop maximization
        if(iter > 1){
            if(abs((total.lik[iter]-total.lik[iter-1]) /
                   total.lik[iter-1])<stop.epsilon)
                break
            }

        # 1. Optimize dispersion
        zeta <- zinbOptimizeDispersion(J, mu, logitPi, getEpsilon_zeta(m), Y,
                                    commondispersion=commondispersion,
                                    BPPARAM=BPPARAM)

        # Evaluate total penalized likelihood
        if (verbose) {
            pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
                sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
                sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
                sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
                sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
                sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
                sum(getEpsilon_W(m)*t(W)^2)/2 +
                getEpsilon_zeta(m)*var(zeta)/2
            message("After dispersion optimization = ",
                    zinb.loglik(Y, mu, exp(zeta), logitPi) - pen)
            }

        # 2. Optimize right factors

        if (optimright) {
            ptm <- proc.time()
            estimate <- matrix(unlist(
                bplapply(seq(J), function(j) {
                    optimright_fun(beta_mu[,j], alpha_mu[,j], beta_pi[,j], alpha_pi[,j],
                                   Y[,j], X_mu, W, V_mu[j,], gamma_mu, O_mu[,j], X_pi,
                                   V_pi[j,], gamma_pi, O_pi[,j], zeta[j], n, epsilonright)
                },BPPARAM=BPPARAM)), nrow=sum(nright))

            if (verbose) {print(proc.time()-ptm)}
            ind <- 1
            if (nright[1]>0) {
                beta_mu <- estimate[ind:(ind+nright[1]-1),,drop=FALSE]
                ind <- ind+nright[1]
                }
            if (nright[2]>0) {
                alpha_mu <- estimate[ind:(ind+nright[2]-1),,drop=FALSE]
                ind <- ind+nright[2]
            }
            if (nright[3]>0) {
                beta_pi <- estimate[ind:(ind+nright[3]-1),,drop=FALSE]
                ind <- ind+nright[3]
            }
            if (nright[4]>0) {
                alpha_pi <- estimate[ind:(ind+nright[4]-1),,drop=FALSE]
            }
        }
        # Evaluate total penalized likelihood
        if (verbose) {
            itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                          W %*% alpha_mu + O_mu)

            iterlP <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
                W %*% alpha_pi + O_pi

            pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
                sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
                sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
                sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
                sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
                sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
                sum(getEpsilon_W(m)*t(W)^2)/2 +
                getEpsilon_zeta(m)*var(zeta)/2
            message("After right optimization = ",
                    zinb.loglik(Y, itermu, exp(zeta), iterlP) - pen)
            }

        # 3. Orthogonalize
        if (orthog) {
            o <- orthogonalizeTraceNorm(W, cbind(alpha_mu,
                                                       alpha_pi),
                                      m@epsilon_W, m@epsilon_alpha)
            W <- o$U
            alpha_mu <- o$V[,1:J,drop=FALSE]
            alpha_pi <- o$V[,(J+1):(2*J),drop=FALSE]
        }

        # Evaluate total penalized likelihood
        if (verbose) {
            itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                              W %*% alpha_mu + O_mu)

            iterlP <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
                W %*% alpha_pi + O_pi

            pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
                sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
                sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
                sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
                sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
                sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
                sum(getEpsilon_W(m)*t(W)^2)/2 +
                getEpsilon_zeta(m)*var(zeta)/2

            message("After orthogonalization = ",
                    zinb.loglik(Y, itermu, exp(zeta), iterlP) - pen)
        }

        # 4. Optimize left factors
        if (optimleft) {
            ptm <- proc.time()
            estimate <- matrix(unlist(
                bplapply(seq(n), function(i) {
                    optimleft_fun(gamma_mu[,i], gamma_pi[,i], W[i,], Y[i,], V_mu, alpha_mu,
                                              X_mu[i,], beta_mu, O_mu[i,], V_pi, alpha_pi, X_pi[i,],
                                              beta_pi, O_pi[i,], zeta, epsilonleft)
                }, BPPARAM=BPPARAM)), nrow=sum(nleft))

            if (verbose) {print(proc.time()-ptm)}
            ind <- 1
            if (nleft[1]>0) {
                gamma_mu <- estimate[ind:(ind+nleft[1]-1),,drop=FALSE]
                ind <- ind+nleft[1]
            }
            if (nleft[2]>0) {
                gamma_pi <- estimate[ind:(ind+nleft[2]-1),,drop=FALSE]
                ind <- ind+nleft[2]
            }
            if (nleft[3]>0) {
                W <- t(estimate[ind:(ind+nleft[3]-1),,drop=FALSE])
                ind <- ind+nleft[3]
            }
        }

        # Evaluate total penalized likelihood
        if (verbose) {
            itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                              W %*% alpha_mu + O_mu)

            iterlP <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
                W %*% alpha_pi + O_pi

            pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
                sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
                sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
                sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
                sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
                sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
                sum(getEpsilon_W(m)*t(W)^2)/2 +
                getEpsilon_zeta(m)*var(zeta)/2

            message("After left optimization = ",
                    zinb.loglik(Y, itermu, exp(zeta), iterlP) - pen)
        }

        # 5. Orthogonalize
        if (orthog) {
            o <- orthogonalizeTraceNorm(W, cbind(alpha_mu,
                                                       alpha_pi),
                                      m@epsilon_W, m@epsilon_alpha)
            W <- o$U
            alpha_mu <- o$V[,1:J,drop=FALSE]
            alpha_pi <- o$V[,(J+1):(2*J),drop=FALSE]
        }
        # Evaluate total penalized likelihood
        if (verbose) {
            itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                              W %*% alpha_mu + O_mu)

            iterlP <- X_pi %*% beta_pi + t(V_pi %*% gamma_pi) +
                W %*% alpha_pi + O_pi

            pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
                sum(getEpsilon_alpha(m) * (alpha_pi)^2)/2 +
                sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
                sum(getEpsilon_beta_pi(m) * (beta_pi)^2)/2 +
                sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
                sum(getEpsilon_gamma_pi(m)*(gamma_pi)^2)/2 +
                sum(getEpsilon_W(m)*t(W)^2)/2 +
                getEpsilon_zeta(m)*var(zeta)/2

            message("After orthogonalization = ",
                    zinb.loglik(Y, itermu, exp(zeta), iterlP) - pen)
        }

    }

    out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu, O_pi = m@O_pi,
                     which_X_mu = m@which_X_mu, which_X_pi = m@which_X_pi,
                     which_V_mu = m@which_V_mu, which_V_pi = m@which_V_pi,
                     W = W, beta_mu = beta_mu, beta_pi = beta_pi,
                     gamma_mu = gamma_mu, gamma_pi = gamma_pi,
                     alpha_mu = alpha_mu, alpha_pi = alpha_pi, zeta = zeta,
                     epsilon_beta_mu = m@epsilon_beta_mu,
                     epsilon_gamma_mu = m@epsilon_gamma_mu,
                     epsilon_beta_pi = m@epsilon_beta_pi,
                     epsilon_gamma_pi = m@epsilon_gamma_pi,
                     epsilon_W = m@epsilon_W, epsilon_alpha = m@epsilon_alpha,
                     epsilon_zeta = m@epsilon_zeta,
                     epsilon_min_logit = m@epsilon_min_logit)

    return(out)
}

#' Optimize the dispersion parameters of a ZINB regression model
#'
#' The dispersion parameters of the model are optimized by
#' penalized maximum likelihood on the count matrix given as argument.
#' @param J The number of genes.
#' @param mu the matrix containing the mean of the negative binomial.
#' @param logitPi the matrix containing the logit of the probability parameter
#' of the zero-inflation part of the model.
#' @param epsilon the regularization parameter.
#' @param Y The matrix of counts.
#' @param commondispersion Whether or not a single dispersion for all features
#'   is estimated (default TRUE)
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @return An object of class ZinbModel similar to the one given as argument
#'   with modified parameters zeta.
#' @examples
#' Y = matrix(10, 3, 5)
#' m = zinbModel(n=NROW(Y), J=NCOL(Y))
#' m = zinbInitialize(m, Y)
#' m = zinbOptimizeDispersion(NROW(Y), getMu(m), getLogitPi(m), getEpsilon_zeta(m), Y)
#' @export
zinbOptimizeDispersion <- function(J, mu, logitPi, epsilon,
                                   Y, commondispersion=TRUE,
                                   BPPARAM=BiocParallel::bpparam()) {

    # 1) Find a single dispersion parameter for all counts by 1-dimensional
    # optimization of the likelihood
    g=optimize(f=zinb.loglik.dispersion, Y=Y, mu=mu,
               logitPi=logitPi, maximum=TRUE,interval=c(-100,100))

    zeta <- rep(g$maximum,J)

    if (!commondispersion) {

        # 2) Optimize the dispersion parameter of each sample
        locfun <- function(logt) {
            s <- sum(unlist(bplapply(seq(J),function(i) {
                zinb.loglik.dispersion(logt[i],Y[,i],mu[,i],logitPi[,i])
            }, BPPARAM=BPPARAM)))
            if (J>1) {
                s <- s - epsilon*var(logt)/2
            }
            s
        }
        locgrad <- function(logt) {
            s <- unlist(bplapply(seq(J),function(i) {
                zinb.loglik.dispersion.gradient(logt[i], Y[,i], mu[,i],
                                                logitPi[,i])
            }, BPPARAM=BPPARAM ))
            if (J>1) {
                s <- s - epsilon*(logt - mean(logt))/(J-1)
            }
            s
        }
        zeta <- optim( par=zeta, fn=locfun , gr=locgrad,
                       control=list(fnscale=-1,trace=0), method="BFGS")$par
    }
    return(zeta)
}


#' Log-likelihood of the zero-inflated negative binomial model
#'
#' Given a vector of counts, this function computes the sum of the
#' log-probabilities of the counts under a zero-inflated negative binomial
#' (ZINB) model. For each count, the ZINB distribution is parametrized by three
#' parameters: the mean value and the dispersion of the negative binomial
#' distribution, and the probability of the zero component.
#'
#' @param Y the vector of counts
#' @param mu the vector of mean parameters of the negative binomial
#' @param theta the vector of dispersion parameters of the negative binomial, or
#'   a single scalar is also possible if the dispersion parameter is constant.
#'   Note that theta is sometimes called inverse dispersion parameter (and
#'   phi=1/theta is then called the dispersion parameter). We follow the
#'   convention that the variance of the NB variable with mean mu and dispersion
#'   theta is mu + mu^2/theta.
#' @param logitPi the vector of logit of the probabilities of the zero component
#' @export
#' @return the log-likelihood of the model.
#' @importFrom copula log1pexp
#' @importFrom stats dnbinom optim optimize rbinom rnbinom runif var
#' @examples
#' n <- 10
#' mu <- seq(10,50,length.out=n)
#' logitPi <- rnorm(10)
#' zeta <- rnorm(10)
#' Y <- rnbinom(n=n, size=exp(zeta), mu=mu)
#' zinb.loglik(Y, mu, exp(zeta), logitPi)
#' zinb.loglik(Y, mu, 1, logitPi)
zinb.loglik <- function(Y, mu, theta, logitPi) {

    # log-probabilities of counts under the NB model
    logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

    # contribution of zero inflation
    lognorm <- - copula::log1pexp(logitPi)

    # log-likelihood
    sum(logPnb[Y>0]) + sum(logPnb[Y==0] + copula::log1pexp(logitPi[Y==0] -
                                                logPnb[Y==0])) + sum(lognorm)
}


#' Log-likelihood of the zero-inflated negative binomial model, for a fixed
#' dispersion parameter
#'
#' Given a unique dispersion parameter and a set of counts, together with a
#' corresponding set of mean parameters and probabilities of zero components,
#' this function computes the sum of the log-probabilities of the counts under
#' the ZINB model. The dispersion parameter is provided to the function through
#' zeta = log(theta), where theta is sometimes called the inverse dispersion
#' parameter. The probabilities of the zero components are provided through
#' their logit, in order to better numerical stability.
#'
#' @param zeta a scalar, the log of the inverse dispersion parameters of the
#'   negative binomial model
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @param logitPi a vector of logit of the probabilities of the zero component
#' @export
#' @seealso \code{\link{zinb.loglik}}.
#' @return the log-likelihood of the model.
#' @examples
#' mu <- seq(10,50,length.out=10)
#' logitPi <- rnorm(10)
#' zeta <- rnorm(10)
#' Y <- rnbinom(n=10, size=exp(zeta), mu=mu)
#' zinb.loglik.dispersion(zeta, Y, mu, logitPi)
zinb.loglik.dispersion <- function(zeta, Y, mu, logitPi) {
    zinb.loglik(Y, mu, exp(zeta), logitPi)
}

#' Derivative of the log-likelihood of the zero-inflated negative binomial model
#' with respect to the log of the inverse dispersion parameter
#'
#' @param zeta the log of the inverse dispersion parameters of the negative
#'   binomial
#' @param Y a vector of counts
#' @param mu a vector of mean parameters of the negative binomial
#' @param logitPi a vector of the logit of the probability of the zero component
#' @export
#' @return the gradient of the inverse dispersion parameters.
#' @seealso \code{\link{zinb.loglik}}, \code{\link{zinb.loglik.dispersion}}.
zinb.loglik.dispersion.gradient <- function(zeta, Y, mu, logitPi) {
    theta <- exp(zeta)

    # Check zeros in the count vector
    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE,Y0))
    has1 <- !is.na(match(TRUE,Y1))

    grad <- 0
    if (has1) {
        grad <- grad + sum( theta * (digamma(Y[Y1] + theta) - digamma(theta) +
                                        zeta - log(mu[Y1] + theta) + 1 -
                                         (Y[Y1] + theta)/(mu[Y1] + theta) ) )
    }

    if (has0) {
        logPnb <- suppressWarnings(dnbinom(0, size = theta, mu = mu[Y0],
                                           log = TRUE))
        grad <- grad + sum( theta * (zeta - log(mu[Y0] + theta) + 1 -
                        theta/(mu[Y0] + theta)) / (1+exp(logitPi[Y0] - logPnb)))
                        # *exp(- copula::log1pexp( -logPnb + logitPi[Y0])))

    }

    grad
}




#' Parse ZINB regression model
#'
#' Given the parameters of a ZINB regression model, this function parses the
#' model and computes the vector of log(mu), logit(pi), and the dimensions of
#' the different components of the vector of parameters. See
#' \code{\link{zinb.loglik.regression}} for details of the ZINB regression model
#' and its parameters.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi, b) concatenated
#' @param A.mu matrix of the model (see above, default=empty)
#' @param B.mu matrix of the model (see above, default=empty)
#' @param C.mu matrix of the model (see above, default=zero)
#' @param A.pi matrix of the model (see above, default=empty)
#' @param B.pi matrix of the model (see above, default=empty)
#' @param C.pi matrix of the model (see above, default=zero)
#' @return A list with slots \code{logMu}, \code{logitPi}, \code{dim.alpha} (a
#'   vector of length 3 with the dimension of each of the vectors \code{a.mu},
#'   \code{a.pi} and \code{b} in \code{alpha}), and \code{start.alpha} (a vector
#'   of length 3 with the starting indices of the 3 vectors in \code{alpha})
#' @seealso \code{\link{zinb.loglik.regression}}
zinb.regression.parseModel <- function(alpha, A.mu, B.mu, C.mu,
                                       A.pi, B.pi, C.pi) {

    n <- nrow(A.mu)
    logMu <- C.mu
    logitPi <- C.pi
    dim.alpha <- rep(0,3)
    start.alpha <- rep(NA,3)
    i <- 0

    j <- ncol(A.mu)
    if (j>0) {
        logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
        dim.alpha[1] <- j
        start.alpha[1] <- i+1
        i <- i+j
    }

    j <- ncol(A.pi)
    if (j>0) {
        logitPi <- logitPi + A.pi %*% alpha[(i+1):(i+j)]
        dim.alpha[2] <- j
        start.alpha[2] <- i+1
        i <- i+j
    }

    j <- ncol(B.mu)
    if (j>0) {
        logMu <- logMu + B.mu %*% alpha[(i+1):(i+j)]
        logitPi <- logitPi + B.pi %*% alpha[(i+1):(i+j)]
        dim.alpha[3] <- j
        start.alpha[3] <- i+1
    }

    return(list(logMu=logMu, logitPi=logitPi, dim.alpha=dim.alpha,
                 start.alpha=start.alpha))
}




#' Penalized log-likelihood of the ZINB regression model
#'
#' This function computes the penalized log-likelihood of a ZINB regression
#' model given a vector of counts.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi, b) concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param B.mu matrix of the model (see Details, default=empty)
#' @param C.mu matrix of the model (see Details, default=zero)
#' @param A.pi matrix of the model (see Details, default=empty)
#' @param B.pi matrix of the model (see Details, default=empty)
#' @param C.pi matrix of the model (see Details, default=zero)
#' @param C.theta matrix of the model (see Details, default=zero)
#' @param epsilon regularization parameter. A vector of the same length as
#'   \code{alpha} if each coordinate of \code{alpha} has a specific
#'   regularization, or just a scalar is the regularization is the same for all
#'   coordinates of \code{alpha}. Default=\code{O}.
#' @details The regression model is parametrized as follows: \deqn{log(\mu) =
#'   A_\mu * a_\mu + B_\mu * b + C_\mu} \deqn{logit(\Pi) = A_\pi * a_\pi + B_\pi
#'   * b} \deqn{log(\theta) = C_\theta} where \eqn{\mu, \Pi, \theta} are
#'   respectively the vector of mean parameters of the NB distribution, the
#'   vector of probabilities of the zero component, and the vector of inverse
#'   dispersion parameters. Note that the \eqn{b} vector is shared between the
#'   mean of the negative binomial and the probability of zero. The
#'   log-likelihood of a vector of parameters \eqn{\alpha = (a_\mu; a_\pi; b)}
#'   is penalized by a regularization term \eqn{\epsilon ||\alpha||^2 / 2} is
#'   \eqn{\epsilon} is a scalar, or \eqn{\sum_{i}\epsilon_i \alpha_i^2 / 2} is
#'   \eqn{\epsilon} is a vector of the same size as \eqn{\alpha} to allow for
#'   differential regularization among the parameters.
#' @return the penalized log-likelihood.
#' @export
zinb.loglik.regression <- function(alpha, Y,
                                   A.mu = matrix(nrow=length(Y), ncol=0),
                                   B.mu = matrix(nrow=length(Y), ncol=0),
                                   C.mu = matrix(0, nrow=length(Y), ncol=1),
                                   A.pi = matrix(nrow=length(Y), ncol=0),
                                   B.pi = matrix(nrow=length(Y), ncol=0),
                                   C.pi = matrix(0, nrow=length(Y), ncol=1),
                                   C.theta = matrix(0, nrow=length(Y), ncol=1),
                                   epsilon=0) {

    # Parse the model
    r <- zinb.regression.parseModel(alpha=alpha,
                                    A.mu = A.mu,
                                    B.mu = B.mu,
                                    C.mu = C.mu,
                                    A.pi = A.pi,
                                    B.pi = B.pi,
                                    C.pi = C.pi)

    # Call the log likelihood function
    z <- zinb.loglik(Y, exp(r$logMu), exp(C.theta), r$logitPi)

    # Penalty
    z <- z - sum(epsilon*alpha^2)/2
    z
}


#' Gradient of the penalized log-likelihood of the ZINB regression model
#'
#' This function computes the gradient of the penalized log-likelihood of a ZINB
#' regression model given a vector of counts.
#'
#' @param alpha the vectors of parameters c(a.mu, a.pi, b) concatenated
#' @param Y the vector of counts
#' @param A.mu matrix of the model (see Details, default=empty)
#' @param B.mu matrix of the model (see Details, default=empty)
#' @param C.mu matrix of the model (see Details, default=zero)
#' @param A.pi matrix of the model (see Details, default=empty)
#' @param B.pi matrix of the model (see Details, default=empty)
#' @param C.pi matrix of the model (see Details, default=zero)
#' @param C.theta matrix of the model (see Details, default=zero)
#' @param epsilon regularization parameter. A vector of the same length as
#'   \code{alpha} if each coordinate of \code{alpha} has a specific
#'   regularization, or just a scalar is the regularization is the same for all
#'   coordinates of \code{alpha}. Default=\code{O}.
#' @details The regression model is described in
#'   \code{\link{zinb.loglik.regression}}.
#' @seealso \code{\link{zinb.loglik.regression}}
#' @return The gradient of the penalized log-likelihood.
#' @export
zinb.loglik.regression.gradient <- function(alpha, Y,
                                     A.mu = matrix(nrow=length(Y), ncol=0),
                                     B.mu = matrix(nrow=length(Y), ncol=0),
                                     C.mu = matrix(0, nrow=length(Y), ncol=1),
                                     A.pi = matrix(nrow=length(Y), ncol=0),
                                     B.pi = matrix(nrow=length(Y), ncol=0),
                                     C.pi = matrix(0, nrow=length(Y), ncol=1),
                                     C.theta = matrix(0, nrow=length(Y), ncol=1),
                                     epsilon=0) {

    # Parse the model
    r <- zinb.regression.parseModel(alpha=alpha,
                                    A.mu = A.mu,
                                    B.mu = B.mu,
                                    C.mu = C.mu,
                                    A.pi = A.pi,
                                    B.pi = B.pi,
                                    C.pi = C.pi)

    theta <- exp(C.theta)
    mu <- exp(r$logMu)
    n <- length(Y)

    # Check zeros in the count matrix
    Y0 <- Y <= 0
    Y1 <- Y > 0
    has0 <- !is.na(match(TRUE,Y0))
    has1 <- !is.na(match(TRUE,Y1))

    # Check what we need to compute,
    # depending on the variables over which we optimize
    need.wres.mu <- r$dim.alpha[1] >0 || r$dim.alpha[3] >0
    need.wres.pi <- r$dim.alpha[2] >0 || r$dim.alpha[3] >0

    # Compute some useful quantities
    muz <- 1/(1+exp(-r$logitPi))
    clogdens0 <- dnbinom(0, size = theta[Y0], mu = mu[Y0], log = TRUE)
    # dens0 <- muz[Y0] + exp(log(1 - muz[Y0]) + clogdens0)
    # More accurate: log(1-muz) is the following
    lognorm <- -r$logitPi - copula::log1pexp(-r$logitPi)

    dens0 <- muz[Y0] + exp(lognorm[Y0] + clogdens0)

    # Compute the partial derivatives we need
    ## w.r.t. mu
    if (need.wres.mu) {
        wres_mu <- numeric(length = n)
        if (has1) {
            wres_mu[Y1] <- Y[Y1] - mu[Y1] *
                (Y[Y1] + theta[Y1])/(mu[Y1] + theta[Y1])
        }
        if (has0) {
            #wres_mu[Y0] <- -exp(-log(dens0) + log(1 - muz[Y0]) + clogdens0 + r$logtheta[Y0] - log(mu[Y0] + theta[Y0]) + log(mu[Y0]))
            # more accurate:
            wres_mu[Y0] <- -exp(-log(dens0) + lognorm[Y0] + clogdens0 +
                                    C.theta[Y0] - log(mu[Y0] + theta[Y0]) +
                                    log(mu[Y0]))
        }
    }

    ## w.r.t. pi
    if (need.wres.pi) {
        wres_pi <- numeric(length = n)
        if (has1) {
            #            wres_pi[Y1] <- -1/(1 - muz[Y1]) * linkobj$mu.eta(r$logitPi)[Y1]
            # simplification for the logit link function:
            wres_pi[Y1] <- - muz[Y1]

        }
        if (has0) {
            # wres_pi[Y0] <- (linkobj$mu.eta(r$logitPi)[Y0] - exp(clogdens0) * linkobj$mu.eta(r$logitPi)[Y0])/dens0
            # simplification for the logit link function:
            wres_pi[Y0] <- (1 - exp(clogdens0)) * muz[Y0] * (1-muz[Y0]) / dens0
        }
    }

    # Make gradient
    grad <- numeric(0)

    ## w.r.t. a_mu
    if (r$dim.alpha[1] >0) {
        istart <- r$start.alpha[1]
        iend <- r$start.alpha[1]+r$dim.alpha[1]-1
        grad <- c(grad , colSums(wres_mu * A.mu) -
                      epsilon[istart:iend]*alpha[istart:iend])
    }

    ## w.r.t. a_pi
    if (r$dim.alpha[2] >0) {
        istart <- r$start.alpha[2]
        iend <- r$start.alpha[2]+r$dim.alpha[2]-1
        grad <- c(grad , colSums(wres_pi * A.pi) -
                      epsilon[istart:iend]*alpha[istart:iend])
    }

    ## w.r.t. b
    if (r$dim.alpha[3] >0) {
        istart <- r$start.alpha[3]
        iend <- r$start.alpha[3]+r$dim.alpha[3]-1
        grad <- c(grad , colSums(wres_mu * B.mu) +
                      colSums(wres_pi * B.pi) -
                      epsilon[istart:iend]*alpha[istart:iend])
    }

    grad
}

optimright_fun <- function(beta_mu, alpha_mu, beta_pi, alpha_pi,
                           Y, X_mu, W, V_mu, gamma_mu, O_mu, X_pi,
                           V_pi, gamma_pi, O_pi, zeta, n, epsilonright) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(beta_mu, alpha_mu,
                 beta_pi, alpha_pi),
           Y=Y, A.mu=cbind(X_mu, W),
           C.mu=t(V_mu %*% gamma_mu) + O_mu,
           A.pi=cbind(X_pi, W),
           C.pi=t(V_pi %*% gamma_pi) + O_pi,
           C.theta=matrix(zeta, nrow = n, ncol = 1),
           epsilon=epsilonright,
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}

optimleft_fun <- function(gamma_mu, gamma_pi, W, Y, V_mu, alpha_mu,
                          X_mu, beta_mu, O_mu, V_pi, alpha_pi, X_pi,
                          beta_pi, O_pi, zeta, epsilonleft) {
    optim( fn=zinb.loglik.regression,
           gr=zinb.loglik.regression.gradient,
           par=c(gamma_mu, gamma_pi,
                 t(W)),
           Y=t(Y),
           A.mu=V_mu,
           B.mu=t(alpha_mu),
           C.mu=t(X_mu%*%beta_mu + O_mu),
           A.pi=V_pi,
           B.pi=t(alpha_pi),
           C.pi=t(X_pi%*%beta_pi + O_pi),
           C.theta=zeta,
           epsilon=epsilonleft,
           control=list(fnscale=-1,trace=0),
           method="BFGS")$par
}

.is_wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
    abs(x - round(x)) < tol
}
