#' Solve ridge regression or logistic regression problems
#'
#' This function solves a regression or logistic regression problem regularized
#' by a L2 or weighted L2 penalty. Contrary to \code{lm.ridge} or \code{glmnet},
#' it works for any number of predictors.
#' @param x a matrix of covariates, one sample per row, one covariate per
#'   column.
#' @param y a vector of response (continuous for regression, 0/1 binary for
#'   logistic regression)
#' @param beta an initial solution where optimization starts (null vector by
#'   default)
#' @param epsilon a scalar or vector of regularization parameters (default
#'   \code{1e-6})
#' @param family a string to choose the type of regression (default
#'   \code{family="gaussian"})
#' @param offset a vector of offsets (default null vector)
#' @return A vector solution of the regression problem
#' @details When \code{family="gaussian"}, we solve the ridge regression problem
#'   that finds the \eqn{\beta} that minimizes: \deqn{||y - x \beta||^2 +
#'   \epsilon||\beta||^2/2 .} When \code{family="binomial"} we solve the ridge
#'   logistic regression problem \deqn{min \sum_i [-y_i (x \beta)_i +
#'   log(1+exp(x\beta)_i)) ] + \epsilon||\beta||^2/2 .} When \code{epsilon} is a
#'   vector of size equal to the size of \code{beta}, then the penalty is a
#'   weighted L2 norm \eqn{\sum_i \epsilon_i \beta_i^2 / 2}.
#' @export
solveRidgeRegression <- function(x, y, beta=rep(0,NCOL(x)), epsilon=1e-6, family=c("gaussian","binomial"), offset=rep(0,NROW(x))) {

    family <-  match.arg(family)

    # loglik
    f <- if (family == "gaussian") {
        function(b) {
            eta <- x %*% b + offset
            l <- sum((y-eta)^2)/2
            l + sum(epsilon*b^2)/2
        }
    } else if (family == "binomial") {
        function(b) {
            eta <- x %*% b + offset
            l <- sum(-y*eta + clog1pexp(eta))
            l + sum(epsilon*b^2)/2
        }
    }

    # gradient of loglik
    g <- if (family == "gaussian") {
        function(b) {
            eta <- x %*% b + offset
            l <- t(x) %*% (-y + eta)
            l + epsilon*b
        }
    } else if (family == "binomial") {
        function(b) {
            eta <- x %*% b + offset
            l <- t(x) %*% (-y + 1/(1+exp(-eta)))
            l + epsilon*b
        }
    }

    # optimize
    m <- optim( fn=f , gr=g , par=beta, control=list(trace=0) , method="BFGS")
    m$par

}
