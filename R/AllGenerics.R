#' Generic function that returns the number of latent factors
#'
#' Given an object that describes a dataset or a model involving latent factors,
#' this function returns the number of latent factors.
#' @param x an object that describes a dataset or a model involving latent
#'   factors
#' @return the number of latent factors
#' @export
setGeneric("nFactors", function(x) standardGeneric("nFactors"))

#' Returns the matrix of mean parameters
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of mean parameters.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of mean parameters
#' @details Note that although the user interface of \code{\link{zinbFit}}
#' requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#' samples in row and genes in column). Hence the parameter matrix returned by
#' this function is of n x J dimensions.
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getMu(a)
#' @export
setGeneric("getMu", function(object) standardGeneric("getMu"))

#' Returns the matrix of logarithm of mean parameters
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of logarithm of mean parameters.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of logarithms of mean parameters
#' @details Note that although the user interface of \code{\link{zinbFit}}
#' requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#' samples in row and genes in column). Hence the parameter matrix returned by
#' this function is of n x J dimensions.
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getLogMu(a)
#' @export
setGeneric("getLogMu", function(object) standardGeneric("getLogMu"))

#' Returns the matrix of logit of probabilities of zero
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of logit of probabilities of 0.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of logit-probabilities of 0
#' @details Note that although the user interface of \code{\link{zinbFit}}
#' requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#' samples in row and genes in column). Hence the parameter matrix returned by
#' this function is of n x J dimensions.
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getLogitPi(a)
#' @export
setGeneric("getLogitPi", function(object) standardGeneric("getLogitPi"))

#' Returns the matrix of probabilities of zero
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of probabilities of 0.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the matrix of probabilities of 0
#' @details Note that although the user interface of \code{\link{zinbFit}}
#' requires a J x n matrix, internally this is stored as a n x J matrix (i.e.,
#' samples in row and genes in column). Hence the parameter matrix returned by
#' this function is of n x J dimensions.
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getPi(a)
#' @export
setGeneric("getPi", function(object) standardGeneric("getPi"))


#' Returns the vector of dispersion parameters
#'
#' Given an object that describes a matrix of zero-inflated negative binomial
#' distributions, returns the vector of dispersion parameters \code{phi}.
#' @param object an object that describes a matrix of zero-inflated.
#'   distributions.
#' @return the vector of dispersion parameters
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getPhi(a)
#' @export
setGeneric("getPhi", function(object) standardGeneric("getPhi"))

#' Returns the vector of inverse dispersion parameters
#'
#' Given an object that describes a matrix of zero-inflated negative binomial
#' distributions, returns the vector of inverse dispersion parameters
#' \code{theta}.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the vector of inverse dispersion parameters theta
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getTheta(a)
#' @export
setGeneric("getTheta", function(object) standardGeneric("getTheta"))

#' Returns the vector of log of inverse dispersion parameters
#'
#' Given an object that describes a matrix of zero-inflated negative binomial
#' distributions, returns the vector \code{zeta} of log of inverse dispersion
#' parameters
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the vector \code{zeta} of log of inverse dispersion parameters
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getZeta(a)
#' @export
setGeneric("getZeta", function(object) standardGeneric("getZeta"))

#' Returns the low-dimensional matrix of inferred sample-level covariates W
#'
#' Given an object that contains the fit of a ZINB-WaVE model, returns the
#' matrix \code{W} of low-dimensional matrix of inferred sample-level
#' covariates.
#'
#' @param object a \code{\linkS4class{ZinbModel}} object, typically the result
#'   of \code{\link{zinbFit}}.
#' @return the matrix \code{W} of inferred sample-level covariates.
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getW(a)
#' @export
setGeneric("getW", function(object) standardGeneric("getW"))

#' Simulate counts from a zero-inflated negative binomial model
#'
#' Given an object that describes zero-inflated negative binomial distribution,
#' simulate counts from the distribution.
#' @param object an object that describes a matrix of zero-inflated negative
#'   binomial.
#' @param seed an optional integer to specify how the random number generator
#'   should be initialized with a call to \code{set.seed}. If missing, the
#'   random generator state is not changed.
#' @param ... additional arguments.
#' @return A list with the following elements.
#'   \itemize{
#'   \item{counts}{the matrix with the simulated counts.}
#'   \item{dataNB}{the data simulated from the negative binomial.}
#'   \item{dataDropouts}{the data simulated from the binomial process.}
#'   \item{zeroFraction}{the fraction of zeros.}
#'   }
#' @examples
#' a <- zinbModel(n=5, J=10)
#' zinbSim(a)
#' @export
setGeneric("zinbSim",function(object, seed, ...) standardGeneric("zinbSim"))

#' Compute the log-likelihood of a model given some data
#'
#' Given a statistical model and some data, this function computes the
#' log-likelihood of the model given the data, i.e., the log-probability of the
#' data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @param ... additional arguments.
#' @return The log-likelihood of the model given the data.
#' @examples
#' m <- zinbModel(n=5, J=10)
#' x <- zinbSim(m)
#' loglik(m, x$counts)
#' @export
setGeneric("loglik", function(model, x, ...) standardGeneric("loglik"))

#' Compute the penalty of a model
#'
#' Given a statistical model with regularization parameters, compute the
#' penalty.
#' @param model an object that describes a statistical model with regularization
#'   parameters.
#' @return The penalty of the model.
#' @examples
#' m <- zinbModel(K=2)
#' penalty(m)
#' @export
setGeneric("penalty", function(model) standardGeneric("penalty"))

#' Fit a ZINB regression model
#'
#' Given an object with the data, it fits a ZINB model.
#'
#' @param Y The data (genes in rows, samples in columns).
#' @param ... Additional parameters to describe the model, see
#'   \code{\link{zinbModel}}.
#' @return An object of class \code{ZinbModel} that has been fitted by penalized
#'   maximum likelihood on the data.
#' @export
setGeneric("zinbFit", function(Y, ...) standardGeneric("zinbFit"))

#' Returns the sample-level design matrix for mu
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the sample-level design matrix for mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the sample-level design matrix for mu
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getX_mu(a)
setGeneric("getX_mu", function(object, ...) standardGeneric("getX_mu"))

#' Returns the sample-level design matrix for pi
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the sample-level design matrix for pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the sample-level design matrix for pi
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getX_pi(a)
setGeneric("getX_pi", function(object, ...) standardGeneric("getX_pi"))

#' Returns the gene-level design matrix for mu
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the gene-level design matrix for mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the gene-level design matrix for mu
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getV_mu(a)
setGeneric("getV_mu", function(object, ...) standardGeneric("getV_mu"))

#' Returns the gene-level design matrix for pi
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the gene-level design matrix for pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the gene-level design matrix for pi
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getV_pi(a)
setGeneric("getV_pi", function(object, ...) standardGeneric("getV_pi"))

#' Returns the matrix of paramters beta_mu
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with X_mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of beta_mu parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getBeta_mu(a)
setGeneric("getBeta_mu", function(object, ...) standardGeneric("getBeta_mu"))

#' Returns the matrix of paramters beta_pi
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with X_pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of beta_pi parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getBeta_pi(a)
setGeneric("getBeta_pi", function(object, ...) standardGeneric("getBeta_pi"))

#' Returns the matrix of paramters gamma_mu
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with V_mu
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of gamma_mu parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getGamma_mu(a)
setGeneric("getGamma_mu", function(object, ...) standardGeneric("getGamma_mu"))

#' Returns the matrix of paramters gamma_pi
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with V_pi
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of gamma_pi parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getGamma_pi(a)
setGeneric("getGamma_pi", function(object, ...) standardGeneric("getGamma_pi"))

#' Returns the matrix of paramters alpha_mu
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with W for the mean part (mu)
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of alpha_mu parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getAlpha_mu(a)
setGeneric("getAlpha_mu", function(object, ...) standardGeneric("getAlpha_mu"))

#' Returns the matrix of paramters alpha_pi
#'
#' Given an object that describes a matrix of zero-inflated distributions,
#' returns the matrix of parameters associated with W for the zero part (pi)
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @param ... Additional parameters.
#' @return the matrix of alpha_pi parameters
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getAlpha_pi(a)
setGeneric("getAlpha_pi", function(object, ...) standardGeneric("getAlpha_pi"))

#' Returns the vector of regularization parameter for beta_mu
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of rows in the parameter \code{beta_mu} with the regularization parameters
#' associated to each row.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{beta_mu}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_beta_mu(a)
setGeneric("getEpsilon_beta_mu",
           function(object) standardGeneric("getEpsilon_beta_mu"))

#' Returns the vector of regularization parameter for gamma_mu
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of columns in the parameter \code{gamma_mu} with the regularization
#' parameters associated to each row.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{gamma_mu}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_gamma_mu(a)
setGeneric("getEpsilon_gamma_mu",
           function(object) standardGeneric("getEpsilon_gamma_mu"))

#' Returns the vector of regularization parameter for beta_pi
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of rows in the parameter \code{beta_pi} with the regularization parameters
#' associated to each row.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{beta_pi}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_beta_pi(a)
setGeneric("getEpsilon_beta_pi",
           function(object) standardGeneric("getEpsilon_beta_pi"))

#' Returns the vector of regularization parameter for gamma_pi
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of columns in the parameter \code{gamma_pi} with the regularization
#' parameters associated to each column.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{gamma_pi}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_gamma_pi(a)
setGeneric("getEpsilon_gamma_pi",
           function(object) standardGeneric("getEpsilon_gamma_pi"))

#' Returns the vector of regularization parameter for W
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of columns in the parameter \code{W} with the regularization
#' parameters associated to each column.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{W}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_W(a)
setGeneric("getEpsilon_W", function(object) standardGeneric("getEpsilon_W"))

#' Returns the vector of regularization parameter for alpha
#'
#' Given an object describing a ZINB model, returns a vector of size the number
#' of rows in the parameter \code{alpha} with the regularization parameters
#' associated to each row. Here \code{alpha} refers to both \code{alpha_mu} and
#' \code{alpha_pi}, which have the same size and have the same regularization.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{alpha_mu} and
#'   \code{alpha_pi}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_alpha(a)
setGeneric("getEpsilon_alpha",
           function(object) standardGeneric("getEpsilon_alpha"))

#' Returns the regularization parameter for the dispersion parameter
#'
#' The regularization parameter penalizes the variance of zeta, the log of
#' the dispersion parameters across samples.
#' @param object an object that describes a matrix of zero-inflated
#'   distributions.
#' @return the regularization parameters for \code{zeta}.
#' @export
#' @examples
#' a <- zinbModel(n=5, J=10)
#' getEpsilon_zeta(a)
setGeneric("getEpsilon_zeta",
           function(object) standardGeneric("getEpsilon_zeta"))

#' Generic function that returns the number of features
#'
#' Given an object that describes a dataset or a model, it returns the number of
#' features.
#' @param x an object that describes a dataset or a model.
#' @return the number of features.
#' @export
setGeneric(
    name = "nFeatures",
    def = function(x) {
        standardGeneric("nFeatures")
    }
)

#' Generic function that returns the number of samples
#'
#' Given an object that describes a model or a dataset, it returns the number of
#' samples.
#' @param x an object that describes a dataset or a model.
#' @return the number of samples.
#' @export
setGeneric(
    name = "nSamples",
    def = function(x) {
        standardGeneric("nSamples")
    }
)

#' Generic function that returns the total number of parameters of the model
#'
#' Given an object that describes a model or a dataset, it returns total number of
#' parameters of the model.
#' @param model an object that describes a dataset or a model.
#' @return the total number of parameters of the model.
#' @export
setGeneric(
    name = "nParams",
    def = function(model) {
        standardGeneric("nParams")
    }
)

#' Perform dimensionality reduction using a ZINB regression model with
#' gene and cell-level covariates.
#'
#' Given an object with the data, it performs dimensionality reduction using
#' a ZINB regression model with gene and cell-level covariates.
#'
#' @param Y The data (genes in rows, samples in columns). Currently implemented
#'   only for \code{SummarizedExperiment}.
#' @param ... Additional parameters to describe the model, see
#'   \code{\link{zinbModel}}.
#' @return An object of class \code{SingleCellExperiment}; the dimensionality
#'   reduced matrix is stored in the \code{reducedDims} slot and optionally
#'   normalized values and residuals are added in the list of assays.
#' @export
setGeneric("zinbwave", function(Y, ...) standardGeneric("zinbwave"))

#' Compute the AIC of a model given some data
#'
#' Given a statistical model and some data, this function computes the AIC
#' of the model given the data, i.e., the AIC of the data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @return the AIC of the model.
#' @export
setGeneric(
    name = "zinbAIC",
    def = function(model, x) {
        standardGeneric("zinbAIC")
    }
)

#' Compute the BIC of a model given some data
#'
#' Given a statistical model and some data, this function computes the BIC
#' of the model given the data, i.e., the BIC of the data under the model.
#' @param model an object that describes a statistical model.
#' @param x an object that describes data.
#' @return the BIC of the model.
#' @export
setGeneric(
    name = "zinbBIC",
    def = function(model, x) {
        standardGeneric("zinbBIC")
    }
)


