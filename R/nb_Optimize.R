nbOptimize <- function(m, Y, commondispersion=TRUE, maxiter=25,
                         stop.epsilon=.0001, verbose=FALSE,
                         BPPARAM=BiocParallel::bpparam()) {

  total.lik=rep(NA,maxiter)
  n <- nSamples(m)
  J <- nFeatures(m)

  epsilonright <- c(getEpsilon_beta_mu(m), getEpsilon_alpha(m))

  nright <- c(length(getEpsilon_beta_mu(m)), length(getEpsilon_alpha(m)))

  optimright = (sum(nright)>0)

  epsilonleft <- c(getEpsilon_gamma_mu(m), getEpsilon_W(m))
  nleft <- c(length(getEpsilon_gamma_mu(m)), length(getEpsilon_W(m)))
  optimleft = (sum(nleft)>0)

  orthog <- (nFactors(m)>0)

  # extract fixed quantities from m
  X_mu <- getX_mu(m)
  V_mu <- getV_mu(m)
  O_mu <- m@O_mu

  # exctract paramters from m (remember to update!)
  beta_mu <- getBeta_mu(m)
  alpha_mu <- getAlpha_mu(m)
  gamma_mu <- getGamma_mu(m)
  W <- getW(m)
  zeta <- getZeta(m)

  for (iter in seq_len(maxiter)){
    if (verbose) {message("Iteration ",iter)}

    # Evaluate total penalized likelihood
    mu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                W %*% alpha_mu + O_mu)

    theta <- exp(zeta)

    loglik <- nb.loglik(Y, mu, rep(theta, rep(n, J)))

    penalty <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
               sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
               sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
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
    zeta <- nbOptimizeDispersion(J, mu, getEpsilon_zeta(m), Y,
                                   commondispersion=commondispersion,
                                   BPPARAM=BPPARAM)

    # Evaluate total penalized likelihood
    if (verbose) {
      pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
             sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
             sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
             sum(getEpsilon_W(m)*t(W)^2)/2 +
             getEpsilon_zeta(m)*var(zeta)/2
      message("After dispersion optimization = ",
              nb.loglik(Y, mu, exp(zeta)) - pen)
    }

    # 2. Optimize right factors

    if (optimright) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(J), function(j) {
          optimright_fun_nb(beta_mu[,j], alpha_mu[,j], Y[,j], X_mu,
                         W, V_mu[j,], gamma_mu, O_mu[,j], zeta[j], n, epsilonright)
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

    }
    # Evaluate total penalized likelihood
    if (verbose) {
      itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                      W %*% alpha_mu + O_mu)


      pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
             sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
             sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
             sum(getEpsilon_W(m)*t(W)^2)/2 +
             getEpsilon_zeta(m)*var(zeta)/2
      message("After right optimization = ",
              nb.loglik(Y, itermu, exp(zeta)) - pen)
    }

    # 3. Orthogonalize
    if (orthog) {
      o <- orthogonalizeTraceNorm(W, alpha_mu, m@epsilon_W, m@epsilon_alpha)
      W <- o$U
      alpha_mu <- o$V
    }

    # Evaluate total penalized likelihood
    if (verbose) {
      itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                      W %*% alpha_mu + O_mu)

      pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
             sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
             sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
             sum(getEpsilon_W(m)*t(W)^2)/2 +
             getEpsilon_zeta(m)*var(zeta)/2

      message("After orthogonalization = ",
              nb.loglik(Y, itermu, exp(zeta)) - pen)
    }

    # 4. Optimize left factors
    if (optimleft) {
      ptm <- proc.time()
      estimate <- matrix(unlist(
        bplapply(seq(n), function(i) {
          optimleft_fun_nb(gamma_mu[,i], W[i,], Y[i,], V_mu, alpha_mu,
                        X_mu[i,], beta_mu, O_mu[i,], zeta, epsilonleft)
        }, BPPARAM=BPPARAM)), nrow=sum(nleft))

      if (verbose) {print(proc.time()-ptm)}
      ind <- 1
      if (nleft[1]>0) {
        gamma_mu <- estimate[ind:(ind+nleft[1]-1),,drop=FALSE]
        ind <- ind+nleft[1]
      }
      if (nleft[2]>0) {
        W <- t(estimate[ind:(ind+nleft[2]-1),,drop=FALSE])
        ind <- ind+nleft[2]
      }
    }

    # Evaluate total penalized likelihood
    if (verbose) {
      itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                      W %*% alpha_mu + O_mu)


      pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
             sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
             sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
             sum(getEpsilon_W(m)*t(W)^2)/2 +
             getEpsilon_zeta(m)*var(zeta)/2

      message("After left optimization = ",
              nb.loglik(Y, itermu, exp(zeta)) - pen)
    }

    # 5. Orthogonalize
    if (orthog) {
      o <- orthogonalizeTraceNorm(W, alpha_mu,
                                  m@epsilon_W, m@epsilon_alpha)
      W <- o$U
      alpha_mu <- o$V[,1:J,drop=FALSE]
    }
    # Evaluate total penalized likelihood
    if (verbose) {
      itermu <- exp(X_mu %*% beta_mu + t(V_mu %*% gamma_mu) +
                      W %*% alpha_mu + O_mu)


      pen <- sum(getEpsilon_alpha(m) * (alpha_mu)^2)/2 +
             sum(getEpsilon_beta_mu(m) * (beta_mu)^2)/2 +
             sum(getEpsilon_gamma_mu(m)*(gamma_mu)^2)/2 +
             sum(getEpsilon_W(m)*t(W)^2)/2 +
             getEpsilon_zeta(m)*var(zeta)/2

      message("After orthogonalization = ",
              nb.loglik(Y, itermu, exp(zeta)) - pen)
    }

  }

  out <- zinbModel(X = m@X, V = m@V, O_mu = m@O_mu,
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




#####################################################################################################

nb.loglik <- function(Y, mu, theta) {

  # log-probabilities of counts under the NB model
  logPnb <- suppressWarnings(dnbinom(Y, size = theta, mu = mu, log = TRUE))

  sum(logPnb)

}





#####################################################################################################

nb.loglik.dispersion <- function(zeta, Y, mu){

  nb.loglik(Y, mu, exp(zeta))

}



#####################################################################################################

nb.regression.parseModel <- function(alpha, A.mu, B.mu, C.mu) {

  n <- nrow(A.mu)
  logMu <- C.mu
  dim.alpha <- rep(0,2)
  start.alpha <- rep(NA,2)
  i <- 0

  j <- ncol(A.mu)
  if (j>0) {
    logMu <- logMu + A.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[1] <- j
    start.alpha[1] <- i+1
    i <- i+j
  }


  j <- ncol(B.mu)
  if (j>0) {
    logMu <- logMu + B.mu %*% alpha[(i+1):(i+j)]
    dim.alpha[2] <- j
    start.alpha[2] <- i+1
  }

  return(list(logMu=logMu, dim.alpha=dim.alpha,
              start.alpha=start.alpha))
}








#####################################################################################################

nb.loglik.regression <- function(alpha, Y,
                                 A.mu = matrix(nrow=length(Y), ncol=0),
                                 B.mu = matrix(nrow=length(Y), ncol=0),
                                 C.mu = matrix(0, nrow=length(Y), ncol=1),
                                 C.theta = matrix(0, nrow=length(Y), ncol=1),
                                 epsilon=0) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                A.mu = A.mu,
                                B.mu = B.mu,
                                C.mu = C.mu)

  # Call the log likelihood function
  z <- nb.loglik(Y, exp(r$logMu), exp(C.theta))

  # Penalty
  z <- z - sum(epsilon*alpha^2)/2
  z
}

nb.loglik.regression.gradient <- function(alpha, Y,
                                          A.mu = matrix(nrow=length(Y), ncol=0),
                                          B.mu = matrix(nrow=length(Y), ncol=0),
                                          C.mu = matrix(0, nrow=length(Y),
                                                        ncol=1),
                                          C.theta = matrix(0, nrow=length(Y),
                                                           ncol=1),
                                          epsilon=0) {

  # Parse the model
  r <- nb.regression.parseModel(alpha=alpha,
                                A.mu = A.mu,
                                B.mu = B.mu,
                                C.mu = C.mu)
  Y=as.vector(Y)
  theta <- exp(C.theta)
  mu <- exp(r$logMu)
  n <- length(Y)

  # Check what we need to compute,
  # depending on the variables over which we optimize
  need.wres.mu <- r$dim.alpha[1] >0 || r$dim.alpha[2] >0

  # Compute the partial derivatives we need
  ## w.r.t. mu
  if (need.wres.mu) {
    wres_mu <- numeric(length = n)
    wres_mu <- Y - mu *
        (Y + theta)/(mu + theta)
    wres_mu <- as.vector(wres_mu)
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


  ## w.r.t. b
  if (r$dim.alpha[2] >0) {
    istart <- r$start.alpha[2]
    iend <- r$start.alpha[2]+r$dim.alpha[2]-1
    grad <- c(grad , colSums(wres_mu * B.mu) -
                epsilon[istart:iend]*alpha[istart:iend])
  }

  grad
}

optimright_fun_nb <- function(beta_mu, alpha_mu, Y, X_mu, W,
                           V_mu, gamma_mu, O_mu, zeta, n, epsilonright) {
  optim( fn=nb.loglik.regression,
         gr=nb.loglik.regression.gradient,
         par=c(beta_mu, alpha_mu),
         Y=Y,
         A.mu=cbind(X_mu, W),
         C.mu=t(V_mu %*% gamma_mu) + O_mu,
         C.theta=matrix(zeta, nrow = n, ncol = 1),
         epsilon=epsilonright,
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}

optimleft_fun_nb <- function(gamma_mu, W, Y, V_mu, alpha_mu,
                          X_mu, beta_mu, O_mu, zeta, epsilonleft) {
  optim( fn=nb.loglik.regression,
         gr=nb.loglik.regression.gradient,
         par=c(gamma_mu, t(W)),
         Y=t(Y),
         A.mu=V_mu,
         B.mu=t(alpha_mu),
         C.mu=t(X_mu%*%beta_mu + O_mu),
         C.theta=zeta,
         epsilon=epsilonleft,
         control=list(fnscale=-1,trace=0),
         method="BFGS")$par
}

nbOptimizeDispersion <- function(J, mu, epsilon,
                                 Y, commondispersion=TRUE,
                                 BPPARAM=BiocParallel::bpparam()) {

  # 1) Find a single dispersion parameter for all counts by 1-dimensional
  # optimization of the likelihood
  g=optimize(f=nb.loglik.dispersion, Y=Y, mu=mu,
             maximum=TRUE,interval=c(-100,100))

  zeta <- rep(g$maximum,J)

  if (!commondispersion) {

    # 2) Optimize the dispersion parameter of each sample
    locfun <- function(logt) {
      s <- sum(unlist(bplapply(seq(J),function(i) {
        nb.loglik.dispersion(logt[i],Y[,i],mu[,i])
      }, BPPARAM=BPPARAM)))
      if (J>1) {
        s <- s - epsilon*var(logt)/2
      }
      s
    }
    locgrad <- function(logt) {
      s <- unlist(bplapply(seq(J),function(i) {
        nb.loglik.dispersion.gradient(logt[i], Y[,i], mu[,i])
      }, BPPARAM=BPPARAM ))
      if (J>1) {
        s <- s - epsilon*(logt - mean(logt))/(J-1)
      }
      s
    }
    zeta <- optim(par=zeta, fn=locfun , gr=locgrad,
                  control=list(fnscale=-1,trace=0), method="BFGS")$par
  }
  return(zeta)
}

nb.loglik.dispersion.gradient <- function(zeta, Y, mu) {
  theta <- exp(zeta)

  grad <- 0
  grad <- grad + sum( theta * (digamma(Y + theta) - digamma(theta) +
                                 zeta - log(mu + theta) + 1 -
                                 (Y + theta)/(mu + theta) ) )

  grad
}
