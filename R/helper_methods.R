#' Toy dataset to check the model
#'
#' @name toydata
#' @aliases toydata
#'
#' @format A matrix of integers (counts) with 96 samples (rows) and 500 genes
#'   (columns).
NULL

#' Initialize an object of class ZinbModel
#'
#' @param X matrix. The design matrix containing sample-level covariates, one
#'   sample per row.
#' @param V matrix. The design matrix containing gene-level covariates, one gene
#'   per row.
#' @param O_mu matrix. The offset matrix for mu.
#' @param O_pi matrix. The offset matrix for pi.
#' @param which_X_mu integer. Indeces of which columns of X to use in the
#'   regression of mu.
#' @param which_V_mu integer. Indeces of which columns of V to use in the
#'   regression of mu.
#' @param which_X_pi integer. Indeces of which columns of X to use in the
#'   regression of pi.
#' @param which_V_pi integer. Indeces of which columns of V to use in the
#'   regression of pi.
#' @param W matrix. The factors of sample-level latent factors.
#' @param beta_mu matrix or NULL. The coefficients of X in the regression of mu.
#' @param gamma_mu matrix or NULL. The coefficients of V in the regression of
#'   mu.
#' @param alpha_mu matrix or NULL. The coefficients of W in the regression of
#'   mu.
#' @param beta_pi matrix or NULL. The coefficients of X in the regression of pi.
#' @param gamma_pi matrix or NULL. The coefficients of V in the regression of
#'   pi.
#' @param alpha_pi matrix or NULL. The coefficients of W in the regression of
#'   pi.
#' @param zeta numeric. A vector of log of inverse dispersion parameters.
#' @param epsilon nonnegative scalar. Regularization parameter.
#' @param epsilon_beta_mu nonnegative scalar. Regularization parameter for
#'   beta_mu.
#' @param epsilon_gamma_mu nonnegative scalar. Regularization parameter for
#'   gamma_mu.
#' @param epsilon_beta_pi nonnegative scalar. Regularization parameter for
#'   beta_pi.
#' @param epsilon_gamma_pi nonnegative scalar. Regularization parameter for
#'   gamma_pi.
#' @param epsilon_W nonnegative scalar. Regularization parameter for W.
#' @param epsilon_alpha nonnegative scalar. Regularization parameter for alpha
#'   (both alpha_mu and alpha_pi).
#' @param epsilon_zeta nonnegative scalar. Regularization parameter for zeta.
#' @param epsilon_min_logit scalar. Minimum regularization parameter for
#'   parameters of the logit model, including the intercept.
#' @param n integer. Number of samples.
#' @param J integer. Number of genes.
#' @param K integer. Number of latent factors.
#'
#' @export
#'
#'
#' @details This is a wrapper around the new() function to create an
#'   instance of class \code{ZinbModel}. Rarely, the user will need to create a
#'   \code{ZinbModel} object from scratch, as tipically this is the result of
#'   \code{\link{zinbFit}}.
#'
#' @details If any of \code{X}, \code{V}, \code{W} matrices are passed,
#'   \code{n}, \code{J}, and \code{K} are inferred. Alternatively, the user can
#'   specify one or more of \code{n}, \code{J}, and \code{K}.
#'
#' @details The regularization parameters can be set by a unique parameter
#'   \code{epsilon} or specific values for the different regularization
#'   parameters can also be provided.
#'   If only \code{epsilon} is specified, the other parameters take the
#'   following values:
#'   \itemize{
#'   \item epsilon_beta = epsilon/J
#'   \item epsilon_gamma = epsilon/n
#'   \item epsilon_W = epsilon/n
#'   \item epsilon_alpha = epsilon/J
#'   \item epsilon_zeta = epsilon
#'   }
#'   We empirically found that large values of \code{epsilon} provide a more
#'   stable estimation of \code{W}.
#'
#' @details A call with no argument has the following default values: \code{n =
#'   50}, \code{J = 100}, \code{K = 0}, \code{epsilon=J}.
#'
#' @details Although it is possible to create new instances of the class by
#'   calling this function, this is not the most common way of creating
#'   \code{ZinbModel} objects. The main use of the class is within the
#'   \code{\link{zinbFit}} function.
#'
#' @return an object of class \code{\linkS4class{ZinbModel}}.
#'
#' @examples
#' a <- zinbModel()
#' nSamples(a)
#' nFeatures(a)
#' nFactors(a)
#'
zinbModel <- function(X, V, O_mu, O_pi, which_X_mu,
                      which_X_pi, which_V_mu, which_V_pi, W, beta_mu, beta_pi,
                      gamma_mu, gamma_pi, alpha_mu, alpha_pi, zeta, epsilon,
                      epsilon_beta_mu, epsilon_gamma_mu, epsilon_beta_pi,
                      epsilon_gamma_pi, epsilon_W, epsilon_alpha,
                      epsilon_zeta, epsilon_min_logit, n, J, K) {

    # Find n (default 50), J (default 100), K (default 0)
    if (missing(n)) {
        if (!missing(X)) {
            n <- NROW(X)
        } else if (!missing(gamma_mu)) {
            n <- NCOL(gamma_mu)
        } else if (!missing(gamma_pi)) {
            n <- NCOL(gamma_pi)
        } else if (!missing(W)) {
            n <- NROW(W)
        } else if (!missing(O_mu)) {
            n <- NROW(O_mu)
        } else if (!missing(O_pi)) {
            n <- NROW(O_pi)
        } else {
            n <- 50
        }
    }
    if (missing(J)) {
        if (!missing(V)) {
            J <- NROW(V)
        } else if (!missing(beta_mu)) {
            J <- NCOL(beta_mu)
        } else if (!missing(beta_pi)) {
            J <- NCOL(beta_pi)
        } else if (!missing(alpha_mu)) {
            J <- NCOL(alpha_mu)
        } else if (!missing(alpha_pi)) {
            J <- NCOL(alpha_pi)
        } else if (!missing(O_mu)) {
            J <- NCOL(O_mu)
        } else if (!missing(O_pi)) {
            J <- NCOL(O_pi)
        } else if (!missing(zeta)) {
            J <- length(zeta)
        } else {
            J <- 100
        }
    }
    if (missing(K)) {
        if (!missing(W)) {
            K <- NCOL(W)
        } else {
            K <- 0
        }
    }

    # Set the different slots for the matrices
    if(missing(X)) {
        X <- matrix(1, nrow=n, ncol=1)
    }

    if (missing(V)) {
        V <- matrix(1, nrow=J, ncol=1)
    }

    if(missing(which_X_mu)) {
        if (NCOL(X)>0) {
            which_X_mu <- seq(NCOL(X))
        } else {
            which_X_mu <- integer(0)
        }
    }

    if (missing(which_X_pi)) {
        if (NCOL(X)>0) {
            which_X_pi <- seq(NCOL(X))
        } else {
            which_X_pi <- integer(0)
        }
    }

    if (missing(which_V_mu)) {
        if (NCOL(V)>0) {
            which_V_mu <- seq(NCOL(V))
        } else {
            which_V_mu <- integer(0)
        }
    }

    if (missing(which_V_pi)) {
        if (NCOL(V)>0) {
            which_V_pi <- seq(NCOL(V))
        } else {
            which_V_pi <- integer(0)
        }
    }

    X_mu_intercept <- FALSE
    if(length(which_X_mu) > 0) {
        if (all(X[, which_X_mu[1]]==1)) {
            X_mu_intercept <- TRUE
        }
    }


    V_mu_intercept <- FALSE
    if(length(which_V_mu) > 0) {
        if (all(V[, which_V_mu[1]]==1)) {
            V_mu_intercept <- TRUE
        }
    }

    X_pi_intercept <- FALSE
    if(length(which_X_pi) > 0) {
        if (all(X[, which_X_pi[1]]==1)) {
            X_pi_intercept <- TRUE
        }
    }

    V_pi_intercept <- FALSE
    if(length(which_V_pi) > 0) {
        if (all(V[, which_V_pi[1]]==1)) {
            V_pi_intercept <- TRUE
        }
    }

    if (missing(O_mu)) {
        O_mu <- matrix(0, nrow=n, ncol=J)
    }
    if (missing(O_pi)) {
        O_pi <- matrix(0, nrow=n, ncol=J)
    }

    if (missing(W)) {
        W <- matrix(0, nrow=n , ncol=K)
    }

    if (missing(beta_mu)) {
        beta_mu <- matrix(0, nrow=length(which_X_mu), ncol=J)
    }
    if (missing(beta_pi)) {
        beta_pi <- matrix(0, nrow=length(which_X_pi), ncol=J)
    }
    if (missing(gamma_mu)) {
        gamma_mu <- matrix(0, nrow=length(which_V_mu), ncol=n)
    }
    if (missing(gamma_pi)) {
        gamma_pi <- matrix(0, nrow=length(which_V_pi), ncol=n)
    }
    if (missing(alpha_mu)) {
        alpha_mu <- matrix(0, nrow=K , ncol=J)
    }
    if (missing(alpha_pi)) {
        alpha_pi <- matrix(0, nrow=K , ncol=J)
    }
    if (missing(zeta)) {
        zeta <- numeric(J)
    }

    # Regularization parameters
    if (missing(epsilon)) {
        epsilon <- J
    }

    if (missing(epsilon_min_logit)) {
        epsilon_min_logit <- 1e-3
    }
    if (missing(epsilon_beta_mu)) {
        epsilon_beta_mu <- epsilon/J
    }
    if (missing(epsilon_gamma_mu)) {
        epsilon_gamma_mu <- epsilon/n
    }
    if (missing(epsilon_beta_pi)) {
        epsilon_beta_pi <- epsilon/J
    }
    if (missing(epsilon_gamma_pi)) {
        epsilon_gamma_pi <- epsilon/n
    }
    if (missing(epsilon_W)) {
        epsilon_W <- epsilon/n
    }
    if (missing(epsilon_alpha)) {
        epsilon_alpha <- epsilon/J
    }
    if (missing(epsilon_zeta)) {
        epsilon_zeta <- epsilon
    }

    obj <- new(Class="ZinbModel",
               X = X, V = V, O_mu = O_mu, O_pi = O_pi, which_X_mu = which_X_mu,
               which_X_pi = which_X_pi, which_V_mu = which_V_mu,
               which_V_pi = which_V_pi, X_mu_intercept = X_mu_intercept,
               X_pi_intercept = X_pi_intercept, V_mu_intercept = V_mu_intercept,
               V_pi_intercept = V_pi_intercept, W = W, beta_mu = beta_mu,
               gamma_mu = gamma_mu, alpha_mu = alpha_mu, beta_pi = beta_pi,
               gamma_pi = gamma_pi, alpha_pi = alpha_pi, zeta = zeta,
               epsilon_beta_mu = epsilon_beta_mu,
               epsilon_gamma_mu = epsilon_gamma_mu,
               epsilon_beta_pi = epsilon_beta_pi,
               epsilon_gamma_pi = epsilon_gamma_pi,
               epsilon_W = epsilon_W, epsilon_alpha = epsilon_alpha,
               epsilon_zeta = epsilon_zeta,
               epsilon_min_logit = epsilon_min_logit)

    validObject(obj) # call of the inspector
    return(obj)
}

#' @export
#' @describeIn ZinbModel show useful info on the object.
#'
#' @param object an object of class \code{ZinbModel}.
setMethod("show", "ZinbModel",
          function(object) {
              cat(paste0("Object of class ZinbModel.\n",
                         nSamples(object), " samples; ", nFeatures(object),
                         " genes.\n",
                         NCOL(getX_mu(object)),
                         " sample-level covariate(s) (mu); ",
                         NCOL(getX_pi(object)),
                         " sample-level covariate(s) (pi);\n",
                         NCOL(getV_mu(object)),
                         " gene-level covariate(s) (mu); ",
                         NCOL(getV_pi(object)),
                         " gene-level covariate(s) (pi);\n",
                         nFactors(object), " latent factor(s).\n"))
          }
)


################################################################
# Extract various informations and variables from a ZINB model #
################################################################

#' @export
#' @describeIn ZinbModel returns the number of samples.
#' @param x an object of class \code{ZinbModel}.
setMethod("nSamples", "ZinbModel",
          function(x) {
              return(NROW(x@X))
          }
)

#' @export
#' @describeIn ZinbModel returns the number of features.
setMethod("nFeatures", "ZinbModel",
          function(x) {
              return(NROW(x@V))
          }
)

#' @export
#' @describeIn ZinbModel returns the number of latent factors.
setMethod("nFactors", "ZinbModel",
          function(x) {
              return(NCOL(x@W))
          }
)


#' @export
#' @describeIn ZinbModel returns the sample-level design matrix for mu.
#' @param intercept logical. Whether to return the intercept (ignored if the
#'   design matrix has no intercept). Default \code{TRUE}
setMethod("getX_mu", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@X_mu_intercept && !intercept) {
                  which_X_mu <- object@which_X_mu[-1]
              } else {
                  which_X_mu <- object@which_X_mu
              }
              return(object@X[, which_X_mu, drop=FALSE])
          }
)

#' @export
#' @describeIn ZinbModel returns the sample-level design matrix for pi.
setMethod("getX_pi", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@X_pi_intercept && !intercept) {
                  which_X_pi <- object@which_X_pi[-1]
              } else {
                  which_X_pi <- object@which_X_pi
              }
              return(object@X[, which_X_pi, drop=FALSE])
          }
)

#' @export
#' @describeIn ZinbModel returns the gene-level design matrix for mu.
setMethod("getV_mu", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@V_mu_intercept && !intercept) {
                  which_V_mu <- object@which_V_mu[-1]
              } else {
                  which_V_mu <- object@which_V_mu
              }
              return(object@V[, which_V_mu, drop=FALSE])
          }
)

#' @export
#' @describeIn ZinbModel returns the sample-level design matrix for pi.
setMethod("getV_pi", "ZinbModel",
          function(object, intercept=TRUE) {
              if(object@V_pi_intercept && !intercept) {
                  which_V_pi <- object@which_V_pi[-1]
              } else {
                  which_V_pi <- object@which_V_pi
              }
              return(object@V[, which_V_pi, drop=FALSE])
          }
)

#' @export
#' @describeIn ZinbModel returns the logarithm of the mean of the non-zero
#'   component.
setMethod("getLogMu", "ZinbModel",
          function(object) {
              return(getX_mu(object) %*% object@beta_mu +
                         t(getV_mu(object) %*% object@gamma_mu) +
                         object@W %*% object@alpha_mu +
                         object@O_mu)
          }
)

#' @export
#' @describeIn ZinbModel returns the mean of the non-zero component.
setMethod("getMu", "ZinbModel",
    function(object) {
        return(exp(getLogMu(object)))
    }
)

#' @export
#' @describeIn ZinbModel returns the logit-probability of zero.
#' @importFrom stats binomial
setMethod("getLogitPi", "ZinbModel",
          function(object) {
              return(getX_pi(object) %*% object@beta_pi +
                         t(getV_pi(object) %*% object@gamma_pi) +
                         object@W %*% object@alpha_pi +
                         object@O_pi)
          }
)

#' @export
#' @describeIn ZinbModel returns the probability of zero.
#' @importFrom stats binomial
setMethod("getPi", "ZinbModel",
    function(object) {
        # return(stats::binomial()$linkinv(getLogitPi(object))
        # Instead of the call to stats::binomial() in the previous line, we
        # directly compute with the exp() function which remains exact for
        # smaller values of the arguments.
        return(1/(1+exp(-getLogitPi(object))))
    }
)

#' @export
#' @describeIn ZinbModel returns the log of the inverse of the dispersion
#'   parameter.
setMethod("getZeta", "ZinbModel",
          function(object) {
              return(object@zeta)
          }
)

#' @export
#' @describeIn ZinbModel returns the dispersion parameter.
setMethod("getPhi", "ZinbModel",
          function(object) {
              return(exp(-object@zeta))
          }
)

#' @export
#' @describeIn ZinbModel returns the inverse of the dispersion parameter.
setMethod("getTheta", "ZinbModel",
          function(object) {
              return(exp(object@zeta))
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{beta_mu}.
setMethod("getEpsilon_beta_mu", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_beta_mu, length(object@which_X_mu))
              if (object@X_mu_intercept) {
                  e[1] <- 0
              }
              e
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{gamma_mu}.
setMethod("getEpsilon_gamma_mu", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_gamma_mu, length(object@which_V_mu))
              if (object@V_mu_intercept) {
                  e[1] <- 0
              }
              e
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{beta_pi}.
setMethod("getEpsilon_beta_pi", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_beta_pi, length(object@which_X_pi))
              if (object@X_pi_intercept) {
                  e[1] <- object@epsilon_min_logit
              }
              e
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{gamma_pi}.
setMethod("getEpsilon_gamma_pi", "ZinbModel",
          function(object) {
              e <- rep(object@epsilon_gamma_pi, length(object@which_V_pi))
              if (object@V_pi_intercept) {
                  e[1] <- object@epsilon_min_logit
              }
              e
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{W}.
setMethod("getEpsilon_W", "ZinbModel",
          function(object) {
              rep(object@epsilon_W, nFactors(object))
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{alpha}.
setMethod("getEpsilon_alpha", "ZinbModel",
          function(object) {
              rep(object@epsilon_alpha, nFactors(object))
          }
)

#' @export
#' @describeIn ZinbModel returns the regularization parameters for
#'   \code{zeta}.
setMethod("getEpsilon_zeta", "ZinbModel",
          function(object) {
              object@epsilon_zeta
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix W of inferred sample-level
#'   covariates.
setMethod("getW", "ZinbModel",
          function(object) {
              object@W
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix beta_mu of inferred parameters.
setMethod("getBeta_mu", "ZinbModel",
          function(object) {
              object@beta_mu
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix beta_pi of inferred parameters.
setMethod("getBeta_pi", "ZinbModel",
          function(object) {
              object@beta_pi
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix gamma_mu of inferred parameters.
setMethod("getGamma_mu", "ZinbModel",
          function(object) {
              object@gamma_mu
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix gamma_pi of inferred parameters.
setMethod("getGamma_pi", "ZinbModel",
          function(object) {
              object@gamma_pi
          }
)


#' @export
#' @describeIn ZinbModel returns the matrix alpha_mu of inferred parameters.
setMethod("getAlpha_mu", "ZinbModel",
          function(object) {
              object@alpha_mu
          }
)

#' @export
#' @describeIn ZinbModel returns the matrix alpha_pi of inferred parameters.
setMethod("getAlpha_pi", "ZinbModel",
          function(object) {
              object@alpha_pi
          }
)

#' @export
#' @describeIn nParams returns the total number of parameters in the model.
setMethod("nParams", "ZinbModel",
          function(model) {
              X_pi <- getX_pi(model)
              X_mu <- getX_mu(model)
              V_pi <- getV_pi(model)
              V_mu <- getV_mu(model)
              disp <- get

              n <- nSamples(model)
              J <- nFeatures(model)
              K <- nFactors(model)
              M_pi <- NCOL(X_pi)
              M_mu <- NCOL(X_mu)
              L_pi <- NCOL(V_pi)
              L_mu <- NCOL(V_mu)
              ndisp <- length(unique(getZeta(model)))

              J * (M_mu + M_pi) + n * (L_mu + L_pi) + 2 * K * J + n * K + ndisp
          }
)

########################
# Other useful methods #
########################

#' @export
#' @describeIn zinbSim simulate from a ZINB distribution.
#'
setMethod(
    f="zinbSim",
    signature="ZinbModel",
    definition=function(object, seed) {

        if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
            runif(1)
        }

        if (missing(seed)) {
            RNGstate <- get(".Random.seed", envir = .GlobalEnv)
        } else {
            R.seed <- get(".Random.seed", envir = .GlobalEnv)
            set.seed(seed)
            RNGstate <- structure(seed, kind = as.list(RNGkind()))
            on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
        }

        mu <- getMu(object)
        pi <- getPi(object)
        theta <- getTheta(object)
        n <- nSamples(object)
        J <- nFeatures(object)

        # Simulate negative binomial with the mean matrix and dispersion
        # parameters
        i <- seq(n*J)
        datanb <- rnbinom(length(i), mu = mu[i], size = theta[ceiling(i/n)])
        data.nb <- matrix(datanb, nrow = n)

        # Simulate the binary dropout matrix. "1" means that a dropout (zero) is
        # observed instead of the value
        i <- seq(n*J)
        datado <- rbinom(length(i), size=1, prob = pi[i])
        data.dropout <- matrix(datado, nrow = n)

        # Matrix of zero-inflated counts
        counts <- data.nb * (1 - data.dropout)

        # Fraction of zeros in the matrix
        zero.fraction <- sum(counts == 0) / (n*J)

        ret <- list(counts = t(counts), dataNB = t(data.nb),
                    dataDropouts = t(data.dropout),
                    zeroFraction = zero.fraction)
        attr(ret, "seed") <- RNGstate
        ret
    }
)

#' @export
#' @describeIn loglik return the log-likelihood of the ZINB model.
setMethod(
    f="loglik",
    signature=c("ZinbModel","matrix"),
    definition=function(model, x) {
        zinb.loglik(x, getMu(model),
                    rep(getTheta(model), rep(nSamples(model),nFeatures(model))),
                    getLogitPi(model))
    }
)

#' @export
#' @describeIn zinbAIC returns the AIC of the ZINB model.
setMethod(
    f="zinbAIC",
    signature=c("ZinbModel","matrix"),
    definition=function(model, x) {
        if ((nSamples(model) != nrow(x))|(nFeatures(model) != ncol(x))) {
            stop("x and model should have the same dimensions!")
        }
        k <- nParams(model)
        ll <- loglik(model, x)
        return(2*k - 2*ll)
    }
)

#' @export
#' @describeIn zinbBIC returns the BIC of the ZINB model.
setMethod(
    f="zinbBIC",
    signature=c("ZinbModel","matrix"),
    definition=function(model, x) {
        n <- nSamples(model)
        if ((n != nrow(x))|(nFeatures(model) != ncol(x))) {
            stop("x and model should have the same dimensions!")
        }
        k <- nParams(model)
        ll <- loglik(model, x)
        return(log(n)*k - 2*ll)
    }
)

#' @export
#' @describeIn penalty return the penalization.
setMethod(
    f="penalty",
    signature="ZinbModel",
    definition=function(model) {
        sum(getEpsilon_alpha(model)*(model@alpha_mu)^2)/2 +
        sum(getEpsilon_alpha(model)*(model@alpha_pi)^2)/2 +
        sum(getEpsilon_beta_mu(model)*(model@beta_mu)^2)/2 +
        sum(getEpsilon_beta_pi(model)*(model@beta_pi)^2)/2 +
        sum(getEpsilon_gamma_mu(model)*(model@gamma_mu)^2)/2 +
        sum(getEpsilon_gamma_pi(model)*(model@gamma_pi)^2)/2 +
        sum(getEpsilon_W(model)*t(model@W)^2)/2 +
        getEpsilon_zeta(model)*var(getZeta(model))/2
    }
)

# Copied and modified from log1pexp of the copula package by
# Marius Hofert, Ivan Kojadinovic, Martin Maechler, Jun Yan,
# Johanna G. Neslehova
# This version assumes that there are no NA's.
log1pexp <- function (x, c0 = -37, c1 = 18, c2 = 33.3)
{
    r <- exp(x)
    if (any(i <- c0 < x & (i1 <- x <= c1)))
        r[i] <- log1p(r[i])
    if (any(i <- !i1 & (i2 <- x <= c2)))
        r[i] <- x[i] + 1/r[i]
    if (any(i3 <- !i2))
        r[i3] <- x[i3]
    r
}
