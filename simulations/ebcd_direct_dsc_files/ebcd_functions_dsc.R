
normal_means_loglik <- function(x, s, Et, Et2) {
  idx <- is.finite(s) & s > 0
  x <- x[idx]
  s <- s[idx]
  Et <- Et[idx]
  Et2 <- Et2[idx]
  
  return(-0.5 * sum(log(2 * pi * s^2) + (1 / s^2) * (Et2 - 2 * x * Et + x^2)))
}

PolarU <- function(A) {
  svdA <- svd(A)
  out <- svdA$u %*% t(svdA$v)
  return(out)
}

greedyrotation <- function(Z=NULL, Z.perp=NULL,
                           Rl=NULL, R=NULL, l=NULL){
  
  # compute an orthogonal complement
  if(is.null(Z.perp)){
    Z.perp <- qr.Q(qr(Z), complete=TRUE)[, (ncol(Z)+1):ncol(Z)]
  }
  
  if(is.null(Rl)){
    Rl <- R %*% l
  }
  
  Z.perpT.Rl <- t(Z.perp) %*% Rl
  z <- (Z.perp %*% Z.perpT.Rl)/norm(Z.perpT.Rl, type="F")
  z <- as.vector(z)
  return(z)
}

ebcd_init <- function(X = NULL,
                      S = NULL,
                      C = NULL,
                      N = NULL) {
  if (!is.null(C)) {
    A <- C
  } else if (!is.null(S)) {
    eigS <- eigen(S)
    C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
    A <- C
  } else if (!is.null(X)) {
    N <- nrow(X)
    
    if (nrow(X) > ncol(X)) {
      S <- crossprod(X) / N
      eigS <- eigen(S)
      C <- sqrt(nrow(S)) * eigS$vectors %*% diag(sqrt(pmax(eigS$values, 0))) %*% t(eigS$vectors)
      A <- C
    } else {
      A <- X
    }
  }
  
  if (is.null(N)) {
    N <- nrow(A)
  }
  
  nrowA <- nrow(A)
  
  tau <- prod(dim(A)) / sum(A^2)
  Z <- matrix(0, nrow = nrowA, ncol = 0)
  EL <- matrix(0, nrow = ncol(A), ncol = 0)
  
  
  ebcd_obj <- list(
    A = A, N = N, nrowA = nrowA,
    tau = tau, Z = Z, EL = EL
  )
  
  
  return(ebcd_obj)
}

#alternative to irlba::irlba is RSpectra::svds()
#irlba is a randomized svd method so I get different results if I don't set a seed (I think)

compute_z_init <- function(ebcd_obj, l, R){
  # check if the estimate of l is zero
  if (sum(l^2) < 10^(-16)){
    z <- matrix(rep(0, nrow(R)), ncol = 1)
  } else {
    
  Rl <- R %*% l
  
  if(ncol(ebcd_obj$Z)==0){
    z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))
  }else{
    
    Z.perp <- qr.Q(qr(ebcd_obj$Z), complete=TRUE)[, -(1:ncol(ebcd_obj$Z))]
    z <- sqrt(ebcd_obj$nrowA) * greedyrotation(Z.perp = Z.perp, Rl = Rl)
  }
  }
  return(z)
}

ebcd_nonneg_greedy <- function(ebcd_obj,
                        Kmax = 1,
                        ebnm_fn = ebnm::ebnm_point_exponential,
                        l_tol = 1e-6,
                        maxiter = 500) {
  
  if(length(ebnm_fn)==1){
    ebnm_fn <- rep(list(ebnm_fn), Kmax)
  }
  
  for (k in 1:Kmax) {
    # initialize an additional factor
    R <- ebcd_obj$A - tcrossprod(ebcd_obj$Z, ebcd_obj$EL)
    svd1 <- irlba::irlba(R, nv = 1, nu = 0)
    dv <- svd1$d * svd1$v #
    l <- dv / sqrt(ebcd_obj$nrowA)
    
    # l initialization for non-negative constraint
    l_opt1 <- pmax(l,0)
    l_opt2 <- pmax(-l, 0)
    
    z_opt1 <- compute_z_init(ebcd_obj, l_opt1, R)
    z_opt2 <- compute_z_init(ebcd_obj, l_opt2, R)
    
    error1 <- sum((R - z_opt1 %*% t(l_opt1))^2)
    error2 <- sum((R - z_opt2 %*% t(l_opt2))^2)
    
    if (error1 < error2){
      l <- l_opt1
      z <- z_opt1
    } else {
      l <- l_opt2
      z <- z_opt2
    }
    
    v <- l^2
    if(ncol(ebcd_obj$Z)!=0){
      Z.perp <- qr.Q(qr(ebcd_obj$Z), complete=TRUE)[, -(1:ncol(ebcd_obj$Z))]
    }
    
    tau <- ebcd_obj$tau
    
    for (iter in 1:maxiter){
      # shrinkage step
      x <- c(crossprod(R, z)) / ebcd_obj$nrowA
      s <- rep(sqrt(1 / (ebcd_obj$N * tau)), times=length(x))
      e <- ebnm_fn[[k]](x = x, s = s)
      
      
      l.old <- l
      l <- e$posterior$mean # issues when this is zero
      v <- e$posterior$sd^2
      
      # I added this - check if l is estimated to be zero
      # If it is zero, then break out of optimization loop and stop adding factors
      if (sqrt(sum(l^2)) < 10^(-16)){
        print('Estimate of l is zero')
        break
      }
      
      Rl <- R %*% l
      
      # (greedy) Rotation step
      if(ncol(ebcd_obj$Z)==0){
        z <- sqrt(ebcd_obj$nrowA) * Rl / sqrt(sum(Rl^2))
      }else{
        z <- sqrt(ebcd_obj$nrowA) * greedyrotation(Z.perp = Z.perp, Rl = Rl)
      }
      
      # Precision step
      tau <- prod(dim(R)) / (sum((R - z %*% t(l))^2) + ebcd_obj$nrowA * sum(v))
      
      # check convergence and exit if converged
      if(sqrt(sum((l - l.old)^2))/sqrt(sum(l.old^2)) < l_tol){
        break
      }
    }
    
    # I added this - check if l is estimated to be zero
    # if zero, then stop adding factors
    if (sqrt(sum(l^2)) < 10^(-16)){
      break
    }
    
    ebcd_obj$EL <- cbind(ebcd_obj$EL, l)
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    ebcd_obj$tau <- tau
  }
  
  ebcd_obj$ebnm_fn <- ebnm_fn
  
  return(ebcd_obj)
}

ebcd_backfit <- function(ebcd_obj,
                         tol = 1e-6,
                         maxiter = 5000) {
  Kmax <- ncol(ebcd_obj$Z)
  if (Kmax == 0){
    print('Nothing to backfit. The number of factors is zero.')
    return(ebcd_obj)
  }
  ebcd_obj$KL <- rep(0, length = Kmax)
  ebcd_obj$obj.old <- -Inf
  ebcd_obj$vec.obj <- c()
  
  for (iter in 1:maxiter) {
    # Shrinkage step
    ebcd_obj$V <- matrix(0, nrow = nrow(ebcd_obj$EL), ncol = ncol(ebcd_obj$EL))
    
    for (k in 1:Kmax) {
      ebnm_fn <- ebcd_obj$ebnm_fn[[k]]
      x <- c(crossprod(ebcd_obj$A, ebcd_obj$Z[, k])) / ebcd_obj$nrowA
      s <- rep(sqrt(1 / (ebcd_obj$N * ebcd_obj$tau)), times=length(x))
      e <- ebnm_fn(x = x, s = s)
      
      ebcd_obj$EL[, k] <- e$posterior$mean
      ebcd_obj$V[, k] <- e$posterior$sd^2
      ebcd_obj$KL[k] <- e$log_likelihood +
        - normal_means_loglik(x, s, ebcd_obj$EL[, k], ebcd_obj$EL[, k]^2 + ebcd_obj$V[, k])
    }
    
    # Rotation step
    ebcd_obj$Z <- sqrt(ebcd_obj$nrowA) * PolarU(ebcd_obj$A %*% ebcd_obj$EL)
    
    # Precision step
    ebcd_obj$tau <- prod(dim(ebcd_obj$A)) / (sum((ebcd_obj$A - ebcd_obj$Z %*% t(ebcd_obj$EL))^2) + ebcd_obj$nrowA * sum(ebcd_obj$V))
    
    # check convergence
    ebcd_obj$obj <- -ebcd_obj$N * ncol(ebcd_obj$A) / 2 * log(2 * pi / ebcd_obj$tau) +
      -(ebcd_obj$N * ebcd_obj$tau / 2) * (
        sum(ebcd_obj$A^2) / ebcd_obj$nrowA - 2 * sum(diag(t(ebcd_obj$A) %*% ebcd_obj$Z %*% t(ebcd_obj$EL))) / ebcd_obj$nrowA + sum(ebcd_obj$EL^2) + sum(ebcd_obj$V)
      ) +
      +sum(ebcd_obj$KL)
    
    ebcd_obj$vec.obj <- c(ebcd_obj$vec.obj, ebcd_obj$obj)
    if (iter >= 10 & ebcd_obj$obj - ebcd_obj$obj.old < tol) break
    ebcd_obj$obj.old <- ebcd_obj$obj
  }
  
  return(ebcd_obj)
}

ebcd_nonneg <- function(X = NULL,
                 S = NULL,
                 C = NULL,
                 N = NULL,
                 Kmax = 5,
                 ebnm_fn = ebnm::ebnm_point_exponential,
                 l_tol_greedy = 1e-6,
                 maxiter_greedy = 500,
                 tol_backfit = 1e-6,
                 maxiter_backfit = 5000) {
  ebcd_obj <- ebcd_init(X = X, S = S, C = C, N = N) |>
    ebcd_nonneg_greedy(
      Kmax = Kmax,
      ebnm_fn = ebnm_fn,
      l_tol = l_tol_greedy,
      maxiter = maxiter_greedy
    ) |>
    ebcd_backfit(
      tol = tol_backfit,
      maxiter = maxiter_backfit
    )
  
  return(ebcd_obj)
}

transform_ebcd_Z <- function(Y, ebcd_obj){
  Y.svd <- svd(Y)
  Y.UV <- Y.svd$u %*% t(Y.svd$v)
  Z_transformed <- Y.UV %*% ebcd_obj$Z
  return(Z_transformed)
}
