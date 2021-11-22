## some example functions for analysis
# all examples have 0-1 intervention effects
# functions at the top, examples at the bottom
# each main function returns a list with:
# [1] full covariance matrix
# [2] vector u
# [3] starting inverse submatrix for starting indexes
# [4] value of objective function

#1. simple time series design
gen_tmat <- function(n_obs,
                     cov_pars = c(0.1, 0.25),
                     start_idx) {
  tpos <- seq(-1, 1, length.out = n_obs)
  sig <- matrix(NA, nrow = n_obs, ncol = n_obs)
  for (i in 1:(n_obs - 1)) {
    for (j in (i + 1):n_obs) {
      sig[i, j] <-
        cov_pars[1] * exp(-abs(tpos[i] - tpos[j]) / cov_pars[2])
      sig[j, i] <- sig[i, j]
    }
  }
  diag(sig) <- cov_pars[1]
  X <-
    matrix(c(rep(1, n_obs), I(tpos > 0 & tpos < 0.5) * 1), ncol = 2)
  #X <- cbind(diag(n_obs),matrix(I(tpos>0)*1,ncol=1))
  
  C <- matrix(c(rep(0, ncol(X) - 1), 1), ncol = 1)
  M <- t(X) %*% solve(sig) %*% X
  cM <- t(C) %*% solve(M)
  u <- cM %*% t(X)
  
  # A <- solve(sig[start_idx,start_idx])
  # val <- matrix(u[start_idx],nrow=1)%*%A%*%matrix(u[start_idx],ncol=1)
  return(list(sig, u, cbind(X, tpos)))
}


#2. cluster stepped designs
x_mat <- function(t, J) {
  X <- kronecker(rep(1, J), diag(t))
  XJ <- kronecker(diag(J), rep(1, t))
  X <- cbind(X, XJ)
  int <- c()
  for (i in 1:J) {
    int <- c(int, rep(0, t - (i - 1)), rep(1, i - 1))
  }
  X <- cbind(X, int)
  return(X[, 2:ncol(X)])
}

#function to make Sigma
cov_mat_i <- function(M, icc, R, t, J) {
  mat1 <- matrix(0, nrow = t, ncol = t)
  for (i in 1:t) {
    for (j in 1:t) {
      mat1[i, j] <- icc * R ^ (abs(i - j))
    }
  }
  ones <- matrix(1, nrow = M, ncol = M)
  sig <- kronecker(mat1, ones)
  #sig <- sig + diag(nrow(sig))
  diag(sig) <- 1
  sig <- kronecker(diag(J), sig)
  return(sig)
}

gen_clmat <- function(x_pars = c(6, 7, 10),
                      cov_pars = c(0.25, 0.5),
                      start_idx = 1:100) {
  X <- x_mat(x_pars[1], x_pars[2])
  X <- X[rep(1:nrow(X), each = x_pars[3]),]
  sig <-
    cov_mat_i(x_pars[3], cov_pars[1], cov_pars[2], x_pars[1], x_pars[2])
  
  C <- matrix(c(rep(0, ncol(X) - 1), 1), ncol = 1)
  M <- t(X) %*% solve(sig) %*% X
  cM <- t(C) %*% solve(M)
  u <- cM %*% t(X)
  
  A <- solve(sig[start_idx, start_idx])
  val <-
    matrix(u[start_idx], nrow = 1) %*% A %*% matrix(u[start_idx], ncol = 1)
  
  return(list(sig, u, A, val))
}

#3. geospatial example
gen_stmat <- function(n_dim,
                      cov_pars = c(0.1, 0.25),
                      start_idx) {
  #ndim is the size in one dimension, so total obs is n_dim^2
  xpos <- seq(-1, 1, length.out = 2 * n_dim + 1)
  xpos <- xpos[(1:n_dim) * 2]
  
  pos <- expand.grid(x = xpos, y = xpos)
  sig <- matrix(NA, nrow = nrow(pos), ncol = nrow(pos))
  cat("\nBuilding covariance matrix:\n")
  for (i in 1:(nrow(pos) - 1)) {
    for (j in (i + 1):nrow(pos)) {
      d <- sqrt((pos[i, 1] - pos[j, 1]) ^ 2 + (pos[i, 2] - pos[j, 2]) ^ 2)
      sig[i, j] <- cov_pars[1] * exp(-abs(d) / cov_pars[2])
      sig[j, i] <- sig[i, j]
    }
    cat("\rRow: ", i, " of ", nrow(pos))
  }
  diag(sig) <- cov_pars[1]
  #X <- matrix(c(rep(1,n_dim^2),ifelse(sqrt(pos$x^2 + pos$y^2)<0.5,1,0)),ncol=2)
  X <-
    matrix(c(rep(1, n_dim ^ 2), sqrt(pos$x ^ 2 + pos$y ^ 2)), ncol =
             2)
  
  C <- matrix(c(rep(0, ncol(X) - 1), 1), ncol = 1)
  M <- t(X) %*% solve(sig) %*% X
  cM <- t(C) %*% solve(M)
  u <- cM %*% t(X)
  
  # A <- solve(sig[start_idx,start_idx])
  # val <- matrix(u[start_idx],nrow=1)%*%A%*%matrix(u[start_idx],ncol=1)
  return(list(sig, u, pos, X))
}

#3. geospatial example
gen_stmat_non_linear <- function(n_dim = 20,
                                 cov_pars = c(0.1, 0.25),
                                 theta = c(1.0, 0.5, 0.2),
                                 start_idx = NULL) {
  #ndim is the size in one dimension, so total obs is n_dim^2
  xpos <- seq(-1, 1, length.out = 2 * n_dim + 1)
  xpos <- xpos[(1:n_dim) * 2]
  
  pos <- expand.grid(x = xpos, y = xpos)
  sig <- matrix(NA, nrow = nrow(pos), ncol = nrow(pos))
  cat("\nBuilding covariance matrix:\n")
  for (i in 1:(nrow(pos) - 1)) {
    for (j in (i + 1):nrow(pos)) {
      d <- sqrt((pos[i, 1] - pos[j, 1]) ^ 2 + (pos[i, 2] - pos[j, 2]) ^ 2)
      sig[i, j] <- cov_pars[1] * exp(-abs(d) / cov_pars[2])
      sig[j, i] <- sig[i, j]
    }
    cat("\rRow: ", i, " of ", nrow(pos))
  }
  diag(sig) <- cov_pars[1]
  #X <- matrix(c(rep(1,n_dim^2),ifelse(sqrt(pos$x^2 + pos$y^2)<0.5,1,0)),ncol=2)
  
  d <- sqrt(pos$x ^ 2 + pos$y ^ 2)
  beta0 = theta[0]
  beta1 = theta[1]
  beta2 = theta[2]
  X <-
    matrix(c(rep(1, n_dim ^ 2), exp(-d * beta2), -d * beta1 * exp(-d * beta2)),
           ncol = 3)
  
  #C <- matrix(c(rep(0,ncol(X)-1),1),ncol=1)
  C <- c(0, 1, 1)
  M <- t(X) %*% solve(sig) %*% X
  cM <- t(C) %*% solve(M)
  u <- cM %*% t(X)
  
  # A <- solve(sig[start_idx,start_idx])
  # val <- matrix(u[start_idx],nrow=1)%*%A%*%matrix(u[start_idx],ncol=1)
  return(list(sig, u, pos, X))
}

## GENERATING EXAMPLES
# some examples are below
# can rewrite the code if needed
# with the u and sigma outputs of these function you can also reorder the matrices
# based on starting values, function below will reorder and return new starting submatrix

# u and sigma from other functions
# new order is ordering of indexes
# m is desired sample size (so it takes 1:m for new submatrix)
reorder_mat <- function(u,
                        sigma,
                        new_order,
                        m) {
  u <- u[new_order]
  sigma <- sigma[new_order, new_order]
  A <- solve(sigma[1:m, 1:m])
  return(list(sigma, u, A))
}

# time series
gen_tmat(n_obs = 10, start_idx = 1:5)

# cluster stepped design
# option 1: 420 observations
gen_clmat()

#option 2: 100 observations
gen_clmat(x_pars = c(4, 5, 5))

# geospatial example
gen_stmat(10, start_idx = 1:50)

#example of reordering
out <- gen_tmat(n_obs = 10, start_idx = 1:5)
out2 <- reorder_mat(out[[2]],
                    out[[1]],
                    sample(1:10, 10, replace = F),
                    5)