% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zinbwave.R
\name{computeDevianceResiduals}
\alias{computeDevianceResiduals}
\title{Deviance residuals of the zero-inflated negative binomial model}
\usage{
computeDevianceResiduals(model, x, ignoreW = TRUE)
}
\arguments{
\item{model}{the zinb model}

\item{x}{the matrix of counts n cells by J genes}

\item{ignoreW}{logical, if true matrix \code{W} is ignored. Default is TRUE.}
}
\value{
the matrix of deviance residuals of the model.
}
\description{
Given a matrix of counts, this function computes the
deviance residuals under a zero-inflated negative binomial
(ZINB) model.
}
\examples{
se <- SummarizedExperiment(matrix(rpois(60, lambda=5), nrow=10, ncol=6),
    colData = data.frame(bio = gl(2, 3)))
m <- zinbFit(se, X=model.matrix(~bio, data=colData(se)),
    BPPARAM=BiocParallel::SerialParam())
computeDevianceResiduals(m, t(assay(se)))
}
