#' Orthogonalize matrices to minimize trace norm of their product
#'
#' Given two matrices U and V that can be multiplied, this function finds two
#' new matrices U2 and V2 such that their product is conserved (U*V = U2*V2) and
#' such that a||U||^2 + b||V||^2 is minimized.
#' @param U left matrix
#' @param V right matrix
#' @param a weight of the norm of U (default=1)
#' @param b weight of the norm of V (default=1)
#' @return A list with the two matrices that solve the problem in the slots U
#'   and V.
#' @export
#' @examples
#' U <- matrix(rnorm(15),5,3)
#' V <- matrix(rnorm(12),3,4)
#' o <- orthogonalizeTraceNorm(U,V)
#' norm( U%*%V - o$U%*%o$V) # should be zero
#' sum(U^2)+sum(V^2)
#' sum(o$U^2)+sum(o$V^2) # should be smaller
orthogonalizeTraceNorm <- function(U, V, a=1, b=1) {

    # do QR of U
    U.qr <- qr (U)
    U.Q <- qr.Q (U.qr)
    U.R <- qr.R (U.qr)

    # do QR of t(V)
    V.qr <- qr (t(V))
    V.Q <- qr.Q (V.qr)
    V.R <- qr.R (V.qr)

    # do SVD of the U.R %*% t(V.R) matrix to have orthog %*% diag %*% orthog
    A <- svd( U.R %*% t(V.R) )

    # Scaling factor
    s <- (a/b)^{1/4}

    # orthogonalized U
    U2 <- 1/s * U.Q %*% A$u %*% sqrt(diag(A$d,length(A$d)))

    # orthogonalized V and W
    V2 <- s * t( V.Q %*% A$v %*% sqrt(diag(A$d,length(A$d))) )

    list(U=U2, V=V2)
}
