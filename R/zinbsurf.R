#' @describeIn zinbsurf Y is a
#'   \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
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
#' @param K integer. Number of latent factors. Specify \code{K = 0} if only
#'   computing observational weights.
#' @param zeroinflation Whether or not a ZINB model should be fitted. If FALSE,
#'   a negative binomial model is fitted instead.
#' @param which_assay numeric or character. Which assay of Y to use. If missing,
#'   if `assayNames(Y)` contains "counts" then that is used. Otherwise, the
#'   first assay is used.
#' @param which_genes character. Which genes to use to estimate W (see details).
#'   Ignored if \code{fitted_model} is provided.
#' @param prop_fit numeric between 0 and 1. The proportion of cells to use for
#'   the zinbwave fit.
#' @param BPPARAM object of class \code{bpparamClass} that specifies the
#'   back-end to be used for computations. See
#'   \code{\link[BiocParallel]{bpparam}} for details.
#' @param verbose Print helpful messages.
#'
#' @details This function implements an approximate strategy, in which the full
#'   \code{zinbwave} model is fit only on a random subset of the data
#'   (controlled by the \code{prop_fit} parameter). The rest of the samples are
#'   subsequently projected onto the low-rank space. This strategy is much
#'   faster and uses less memory than the full \code{\link{zinbwave}} method. It
#'   is recommended with extremely large datasets.
#'
#' @details By default \code{zinbsurf} uses all genes to estimate \code{W}.
#'   However, we recommend to use the top 1,000 most variable genes for this
#'   step. In general, a user can specify any custom set of genes to be used to
#'   estimate \code{W}, by specifying either a vector of gene names, or a single
#'   character string corresponding to a column of the \code{rowData}.
#'
#'
#' @examples
#' se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
#'                            colData = data.frame(bio = gl(2, 3)))
#' colnames(se) <- paste0("sample", 1:6)
#' m <- zinbsurf(se, X="~bio", K = 1, prop_fit = .5, which_assay = 1)
setMethod("zinbsurf", "SummarizedExperiment",
          function(Y, X, V, K, which_assay, which_genes,
                   zeroinflation = TRUE, prop_fit = .1,
                   BPPARAM=BiocParallel::bpparam(), verbose = FALSE, ...) {

              if(prop_fit <= 0 | prop_fit >=1) {
                  stop("`prop_fit` must be in (0, 1).")
              }
              if(is.null(colnames(Y))) {
                  stop("to use `zinbsurf` `Y` needs to have colnames.")
              }

              if(K == 0) {
                  stop("please specify `K > 0`.")
              }

              out <- as(Y, "SingleCellExperiment")

              if(!missing(which_genes)) {
                  if(length(which_genes) == 1) {
                      if(!which_genes %in% colnames(rowData(Y))) {
                          stop("if `which_genes` has length 1, it must be the name of a column of `rowData(Y)`.")
                      }
                      which_genes <- rowData(Y)[,which_genes]
                  }

                  Y <- Y[which_genes,]
              }

              sample_idx <- sample(seq_len(NCOL(Y)), size = floor(NCOL(Y) * prop_fit))
              out_idx <- setdiff(seq_len(NCOL(Y)), sample_idx)

              Ysub <- Y[,sample_idx]
              Yout <- Y[,out_idx]

              if(verbose) {
                  message("Fitting zinbwave model using ", length(sample_idx),
                          " cells.")
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

                  Xsub <- X[sample_idx,]
                  Xout <- X[out_idx,]

                  fit <- zinbFit(Ysub, Xsub, V, K, BPPARAM=BPPARAM,
                                 verbose = verbose,
                                 zeroinflation = zeroinflation, ...)
                  newm <- zinbModel(n=NCOL(Yout), J=NROW(Yout),
                                    X = Xout, V = V, K = K, ...)
              } else {
                  fit <- zinbFit(Ysub, X, V, K, BPPARAM=BPPARAM,
                                 verbose = verbose,
                                 zeroinflation = zeroinflation, ...)
                  newm <- zinbModel(n=NCOL(Yout), J=NROW(Yout),
                                    X = X, V = V, K = K, ...)
              }

              if(missing(which_assay)) {
                  if("counts" %in% assayNames(Yout)) {
                      dataY <- assay(Yout, "counts")
                  } else {
                      warning("No assay named `counts`, using first assay.",
                              "Use `assay` to specify a different assay.")
                      dataY <- assay(Yout)
                  }
              } else {
                  if(!(is.character(which_assay) | is.numeric(which_assay))) {
                      stop("assay needs to be a numeric or character specifying which assay to use")
                  } else {
                      dataY <- assay(Yout, which_assay)
                  }
              }

              ## realize in memory or make non-sparse
              dataY <- as.matrix(dataY)

              if(verbose) {
                  message("Projecting the remaining ", length(out_idx),
                          " cells.")
              }

              if(zeroinflation) {
                  epsilonleft <- c(getEpsilon_gamma_mu(fit),
                                   getEpsilon_gamma_pi(fit), getEpsilon_W(fit))

                  nleft <- c(length(getEpsilon_gamma_mu(fit)),
                             length(getEpsilon_gamma_pi(fit)),
                             length(getEpsilon_W(fit)))

                  newfit <- matrix(unlist(
                      bplapply(seq_along(out_idx), function(i) {
                          optimleft_fun(getGamma_mu(newm)[,i],
                                        getGamma_pi(newm)[,i],
                                        t(getW(newm)[i,]), dataY[,i],
                                        getV_mu(fit),
                                        getAlpha_mu(fit), getX_mu(newm)[i,],
                                        getBeta_mu(fit), newm@O_mu[i,],
                                        getV_pi(fit),
                                        getAlpha_pi(fit), getX_pi(newm)[i,],
                                        getBeta_pi(fit),
                                        newm@O_pi[i,], getZeta(fit), epsilonleft)
                      }, BPPARAM = BPPARAM)), nrow=sum(nleft))
              } else {
                  epsilonleft <- c(getEpsilon_gamma_mu(fit),
                                   getEpsilon_W(fit))

                  nleft <- c(length(getEpsilon_gamma_mu(fit)),
                             length(getEpsilon_W(fit)))

                  newfit <- matrix(unlist(
                      bplapply(seq_along(out_idx), function(i) {
                          optimleft_fun_nb(getGamma_mu(newm)[,i],
                                        t(getW(newm)[i,]), dataY[,i],
                                        getV_mu(fit),
                                        getAlpha_mu(fit), getX_mu(newm)[i,],
                                        getBeta_mu(fit), newm@O_mu[i,],
                                        getZeta(fit), epsilonleft)
                      }, BPPARAM = BPPARAM)), nrow=sum(nleft))

              }

              newW <- t(newfit[-(seq_len(sum(nleft[1:2]))), ,drop = FALSE])
              W <- rbind(getW(fit), newW)
              o <- orthogonalizeTraceNorm(W, cbind(getAlpha_mu(fit),
                                                   getAlpha_pi(fit)),
                                          fit@epsilon_W, fit@epsilon_alpha)
              W <- o$U
              rownames(W) <- c(colnames(Ysub), colnames(Yout))
              colnames(W) <- paste0('W', seq_len(nFactors(fit)))

              reducedDim(out, "zinbwave") <- W[colnames(Y), , drop=FALSE]
              return(out)
          }
)
