#' Class ZinbModel
#'
#' Objects of this class store all the values needed to work with a
#' zero-inflated negative binomial (ZINB) model, as described in the vignette.
#' They contain all information to fit a model by penalized maximum likelihood
#' or simulate data from a model.
#'
#' @slot X matrix. The design matrix containing sample-level covariates, one
#'   sample per row.
#' @slot V matrix. The design matrix containing gene-level covariates, one gene
#'   per row.
#' @slot O_mu matrix. The offset matrix for mu.
#' @slot O_pi matrix. The offset matrix for pi.
#' @slot which_X_mu integer. Indeces of which columns of X to use in the
#'   regression of mu.
#' @slot which_V_mu integer. Indeces of which columns of V to use in the
#'   regression of mu.
#' @slot which_X_pi integer. Indeces of which columns of X to use in the
#'   regression of pi.
#' @slot which_V_pi integer. Indeces of which columns of V to use in the
#'   regression of pi.
#' @slot X_mu_intercept logical. TRUE if X_mu contains an intercept.
#' @slot X_pi_intercept logical. TRUE if X_pi contains an intercept.
#' @slot V_mu_intercept logical. TRUE if V_mu contains an intercept.
#' @slot V_pi_intercept logical. TRUE if V_pi contains an intercept.
#' @slot W matrix. The factors of sample-level latent factors.
#' @slot beta_mu matrix or NULL. The coefficients of X in the regression of mu.
#' @slot gamma_mu matrix or NULL. The coefficients of V in the regression of mu.
#' @slot alpha_mu matrix or NULL. The coefficients of W in the regression of mu.
#' @slot beta_pi matrix or NULL. The coefficients of X in the regression of pi.
#' @slot gamma_pi matrix or NULL. The coefficients of V in the regression of pi.
#' @slot alpha_pi matrix or NULL. The coefficients of W in the regression of pi.
#' @slot zeta numeric. A vector of log of inverse dispersion parameters.
#' @slot epsilon_beta_mu nonnegative scalar. Regularization parameter for
#'   beta_mu
#' @slot epsilon_gamma_mu nonnegative scalar. Regularization parameter for
#'   gamma_mu
#' @slot epsilon_beta_pi nonnegative scalar. Regularization parameter for
#'   beta_pi
#' @slot epsilon_gamma_pi nonnegative scalar. Regularization parameter for
#'   gamma_pi
#' @slot epsilon_W nonnegative scalar. Regularization parameter for W
#' @slot epsilon_alpha nonnegative scalar. Regularization parameter for alpha
#'   (both alpha_mu and alpha_pi)
#' @slot epsilon_zeta nonnegative scalar. Regularization parameter for zeta
#' @slot epsilon_min_logit scalar. Minimum regularization parameter for
#'   parameters of the logit model, including the intercept.
#'
#' @details For the full description of the model see the model vignette.
#'   Internally, the slots are checked so that the matrices are of the
#'   appropriate dimensions: in particular, \code{X}, \code{O_mu}, \code{O_pi},
#'   and \code{W} need to have \code{n} rows, \code{V} needs to have \code{J}
#'   rows, \code{zeta} must be of length \code{J}.
#' @name ZinbModel-class
#' @import methods
#' @exportClass ZinbModel
#' @aliases ZinbModel
#'
#' @return \code{nSamples} returns the number of samples; \code{nFeatures}
#' returns the number of features; \code{nFactors} returns the number of latent
#' factors.
#'
setClass(
    Class = "ZinbModel",
    slots = list(X = "matrix",
                 V = "matrix",
                 O_mu = "matrix",
                 O_pi = "matrix",
                 which_X_mu = "integer",
                 which_V_mu = "integer",
                 which_X_pi = "integer",
                 which_V_pi = "integer",
                 X_mu_intercept = "logical",
                 V_mu_intercept = "logical",
                 X_pi_intercept = "logical",
                 V_pi_intercept = "logical",
                 W = "matrix",
                 beta_mu = "matrix",
                 gamma_mu = "matrix",
                 alpha_mu = "matrix",
                 beta_pi = "matrix",
                 gamma_pi = "matrix",
                 alpha_pi = "matrix",
                 zeta = "numeric",
                 epsilon_beta_mu = "numeric",
                 epsilon_gamma_mu = "numeric",
                 epsilon_beta_pi = "numeric",
                 epsilon_gamma_pi = "numeric",
                 epsilon_W = "numeric",
                 epsilon_alpha = "numeric",
                 epsilon_zeta = "numeric",
                 epsilon_min_logit = "numeric"
                 )
)

setValidity("ZinbModel", function(object){
    n <- NROW(getX_mu(object)) # number of samples
    J <- NROW(getV_mu(object)) # number of genes
    K <- NCOL(getW(object)) # number of latent factors

    if(K > n) {
        return("Cannot have more latent factors than samples.")
    }
    if(K > J) {
        return("Cannot have more latent factors than genes.")
    }
    if(NROW(getW(object)) != n) {
        return("W must have n rows!")
    }
    if((length(object@which_X_mu)>0) && (max(object@which_X_mu) > NCOL(object@X))) {
        return("which_X_mu: subscript out of bound!")
    }
    if((length(object@which_X_pi)>0) && (max(object@which_X_pi) > NCOL(object@X))) {
        return("which_X_pi: subscript out of bound!")
    }
    if((length(object@which_V_mu)>0) && (max(object@which_V_mu) > NCOL(object@V))) {
        return("which_V_mu: subscript out of bound!")
    }
    if((length(object@which_V_pi)>0) && (max(object@which_V_pi) > NCOL(object@V))) {
        return("which_V_pi: subscript out of bound!")
    }
    if(NCOL(getX_mu(object)) != NROW(getBeta_mu(object))){
        return("beta_mu must have the same number of rows as there are columns in X_mu!")
    }
    if(NCOL(getX_pi(object)) != NROW(getBeta_pi(object))){
        return("beta_pi must have the same number of rows as there are columns in X_pi!")
    }
    if(NCOL(getV_mu(object)) != NROW(getGamma_mu(object))){
        return("gamma_mu must have the same number of rows as there are columns in V_mu!")
    }
    if(NCOL(getV_pi(object)) != NROW(getGamma_pi(object))){
        return("gamma_pi must have the same number of rows as there are columns in V_pi!")
    }
    if(NCOL(getBeta_mu(object)) != J) {
        return("beta_mu must have J columns!")
    }
    if(NCOL(getBeta_pi(object)) != J) {
        return("beta_pi must have J columns!")
    }
    if(NCOL(getGamma_mu(object)) != n) {
        return("gamma_mu must have n columns!")
    }
    if(NCOL(getGamma_pi(object)) != n) {
        return("gamma_pi must have n columns!")
    }
    if(NCOL(getAlpha_mu(object)) != J) {
        return("alpha_mu must have J columns!")
    }
    if(NCOL(getAlpha_pi(object)) != J) {
        return("alpha_pi must have J columns!")
    }
    if(NROW(getAlpha_mu(object)) != K) {
        return("alpha_mu must have K rows!")
    }
    if(NROW(getAlpha_pi(object)) != K) {
        return("alpha_pi must have K rows!")
    }
    if(NROW(object@O_mu) != n) {
        return("O_mu must have n rows!")
    }
    if(NROW(object@O_pi) != n) {
        return("O_pi must have n rows!")
    }
    if(NCOL(object@O_mu) != J) {
        return("O_mu must have J columns!")
    }
    if(NCOL(object@O_pi) != J) {
        return("O_pi must have J columns!")
    }
    if(length(getZeta(object)) != J) {
        return("zeta must have length J!")
    }
    if((length(object@epsilon_beta_mu) != 1) || (object@epsilon_beta_mu < 0)) {
        return("epsilon_beta_mu must be a nonnegative scalar !")
    }
    if((length(object@epsilon_gamma_mu) != 1) || (object@epsilon_gamma_mu < 0)) {
        return("epsilon_gamma_mu must be a nonnegative scalar !")
    }
    if((length(object@epsilon_beta_pi) != 1) || (object@epsilon_beta_pi < 0)) {
        return("epsilon_beta_pi must be a nonnegative scalar !")
    }
    if((length(object@epsilon_gamma_pi) != 1) || (object@epsilon_gamma_pi < 0)) {
        return("epsilon_gamma_pi must be a nonnegative scalar !")
    }
    if((length(object@epsilon_W) != 1) || (object@epsilon_W < 0)) {
        return("epsilon_W must be a nonnegative scalar !")
    }
    if((length(object@epsilon_alpha) != 1) || (object@epsilon_alpha < 0)) {
        return("epsilon_alpha must be a nonnegative scalar !")
    }
    if((length(object@epsilon_zeta) != 1) || (object@epsilon_zeta < 0)) {
        return("epsilon_zeta must be a nonnegative scalar !")
    }

    if((length(object@epsilon_min_logit) != 1) || (object@epsilon_min_logit < 0)) {
        return("epsilon_min_logit must be a nonnegative scalar !")
    }
    return(TRUE)
}
)
