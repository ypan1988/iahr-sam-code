require(Rcpp)
require(RcppArmadillo)

sourceCpp("rem_add_funcs.cpp")

# R versions of functions

choose_swap <- function(idx_in,A,sig,u){
  idx_out <- 1:nrow(sig)
  idx_out <- idx_out[!idx_out%in%idx_in]
  val_out <- sapply(1:length(idx_in),function(i)remove_one(A,i-1,u[idx_in]))
  rm1A <- remove_one_mat(A,which.max(val_out)-1)
  idx_in <- idx_in[-which.max(val_out)]
  val_in <- sapply(idx_out,function(i)add_one(rm1A,sig[i,i],sig[idx_in,i],u[c(idx_in,i)]))
  swap_idx <- idx_out[which.max(val_in)]
  newA <- add_one_mat(rm1A,sig[swap_idx,swap_idx],sig[idx_in,swap_idx])
  idx_in <- c(idx_in,swap_idx)
  val <- obj_fun(newA,u[idx_in])
  return(list(val,idx_in,newA))
}

choose_swap_robust <- function(idx_in,A_list,sig_list,u_list,weights){
  idx_out <- 1:nrow(sig_list[[1]])
  idx_out <- idx_out[!idx_out%in%idx_in]

  val_out_mat <- matrix(NA,nrow=length(idx_in),ncol=length(A_list))
  val_in_mat <- matrix(NA,nrow=length(idx_out),ncol=length(A_list))

  for(idx in 1:length(A_list)){
    val_out_mat[,idx] <- sapply(1:length(idx_in),function(i)
      remove_one(A_list[[idx]],i-1,u_list[[idx]][idx_in]))
  }
  val_out <- as.numeric(val_out_mat %*% weights)

  rm1A <- list()
  for(idx in 1:length(A_list)){
    rm1A[[idx]] <- remove_one_mat(A_list[[idx]],which.max(val_out)-1)
  }

  idx_in <- idx_in[-which.max(val_out)]

  for(idx in 1:length(A_list)){
    val_in_mat[,idx] <- sapply(idx_out,function(i)add_one(rm1A[[idx]],
                                                          sig_list[[idx]][i,i],
                                                          sig_list[[idx]][idx_in,i],
                                                          u_list[[idx]][c(idx_in,i)]))
  }
  val_in <- as.numeric(val_in_mat %*% weights)
  swap_idx <- idx_out[which.max(val_in)]

  newA <- list()
  for(idx in 1:length(A_list)){
    newA[[idx]] <- add_one_mat(rm1A[[idx]],sig_list[[idx]][swap_idx,swap_idx],
                               sig_list[[idx]][idx_in,swap_idx])
  }
  idx_in <- c(idx_in,swap_idx)
  val <- val_in[which.max(val_in)]
  return(list(val,idx_in,newA))
}

grad <- function(idx_in,A,sig,u,tol=1e-9, trace = TRUE){
  new_val <- obj_fun(A,u[idx_in])
  diff <- 1
  i <- 0
  while(diff > tol){
    val <- new_val
    i <- i + 1
    out <- choose_swap(idx_in,A,sig,u)
    new_val <- out[[1]]
    diff <- new_val - val
    if(diff>0){
      A <- out[[3]]
      idx_in <- out[[2]]
    }
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }
  }
  return(idx_in)
}

grad_robust <- function(idx_in,
                        C_list,
                        X_list,
                        sig_list,
                        w=NULL,
                        tol=1e-9,
                        trace = TRUE){

  if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
  if(sum(w)!=1)w <- w/sum(w)
  if(!is(w,"matrix"))w <- matrix(w,ncol=1)
  if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
  if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
  if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
  if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")

  A_list <- list()
  u_list <- list()
  for(i in 1:length(sig_list)){
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M <- t(X_list[[i]]) %*% solve(sig_list[[i]]) %*% X_list[[i]]
    cM <- t(C_list[[i]]) %*% solve(M)
    u_list[[i]] <- cM %*% t(X_list[[i]])
  }

  new_val_vec <- matrix(sapply(1:length(A_list),function(i)obj_fun(A_list[[i]], u_list[[i]][idx_in])),nrow=1)

  new_val <- as.numeric(new_val_vec %*% w)
  diff <- 1
  i <- 0
  while(diff > tol){
    val <- new_val
    i <- i + 1
    out <- choose_swap_robust(idx_in,A_list,sig_list,u_list, w)
    new_val <- out[[1]]
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }
    diff <- new_val - val
    if(diff>0){
      A_list <- out[[3]]
      idx_in <- out[[2]]
    }

  }

  #return variance
  if(trace){
    var_vals <- c()
    for(i in 1:length(X_list)){
      var_vals[i] <- optim_fun(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in])
    }
    cat("\nVariance for individual model(s):\n")
    print(var_vals)
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(var_vals*c(w)))
    }
  }



  return(idx_in)
}

optim_fun <- function(C,X,S){
  if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
  M <- t(X) %*% solve(S) %*% X
  val <- t(C) %*% solve(M) %*% C
  return(val)
}

optim_fun2 <- function(C,X,S){
  if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
  M <- t(X) %*% solve(S) %*% X
  val <- diag(solve(M))[c(C) != 0]
  return(val)
}


sourceCpp("gd_search.cpp")
grad_robust2 <- function(idx_in, C_list, X_list, sig_list, w=NULL, tol=1e-9,
                        trace = TRUE){
  if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
  if(sum(w)!=1)w <- w/sum(w)
  if(!is(w,"matrix"))w <- matrix(w,ncol=1)
  if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
  if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
  if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
  if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")

  A_list <- list()
  u_list <- list()
  for(i in 1:length(sig_list)){
    A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
    M <- t(X_list[[i]]) %*% solve(sig_list[[i]]) %*% X_list[[i]]
    cM <- t(C_list[[i]]) %*% solve(M)
    u_list[[i]] <- cM %*% t(X_list[[i]])
  }

  idx_in <- GradRobust(length(sig_list), idx_in-1, do.call(rbind, A_list),
            do.call(rbind, sig_list), do.call(cbind, u_list), w, tol, trace)
  idx_in <- idx_in + 1

  #return variance
  if(trace){
    var_vals <- c()
    for(i in 1:length(X_list)){
      var_vals[i] <- optim_fun(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in])
    }
    cat("\nVariance for individual model(s):\n")
    print(var_vals)
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(var_vals*c(w)))
    }
  }
  return(idx_in)
}

# For a given m find the optimal power vector
max_var <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)

  if (length(idx_in) != nrow(X_list[[1]]))
  idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  v0
}


# For a given m find the optimal power vector
max_power <- function(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE){

  # randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)
  #idx_in <- (1:m)*round(nrow(X_list[[1]])/m,0)

  if (length(idx_in) != nrow(X_list[[1]]))
  idx_in <- grad_robust2(idx_in, C_list, X_list, sig_list, w, 1e-9, trace)

  v0 <- c()
  for(i in 1:length(X_list)){
    v0 <- c(v0,optim_fun2(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]))
  }

  pow <- pnorm(sqrt(theta[unlist(C_list)!=0]/sqrt(v0)) - qnorm(1-alpha/2))

  pow
}

sample_size <- function(theta, alpha, pwr_target, m, C_list, X_list, sig_list, w) {
  iter <- 0
  pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w)
  while (!all(pwr_new - pwr_target > 0) & m < nrow(X_list[[1]])) {
    iter <- iter + 1
    m <- m + 1
    pwr_new <- max_power(theta, alpha, m, C_list, X_list, sig_list, w, trace = FALSE)
    cat("\nm = ", m)
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }
  return(m)
}

sample_size2 <- function(theta, alpha, pwr_target, C_list, X_list, sig_list, w) {

  cat("\nTarget power = ", pwr_target)

  lo <- max(unlist(lapply(C_list,function(i)length(unlist(i)))))*3
  hi <- nrow(X_list[[1]])
  pwr_new_lo <-NULL
  while(is.null(pwr_new_lo)){
    cat("\nlo = ", lo)
    pwr_new_lo <- tryCatch(
      max_power(theta, alpha, lo, C_list, X_list, sig_list, w, trace = FALSE),
      error=function(i)NULL)
    lo <- lo+10
  }
  pwr_new_hi <- max_power(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\nmin power = ", min(pwr_new_lo))
  cat("\nmax power = ", min(pwr_new_hi))

  if (min(pwr_new_hi) < pwr_target | min(pwr_new_lo) > pwr_target)
    stop("\ntarget power is not in range of ", min(pwr_new_lo) , " and ", min(pwr_new_hi))

  v_hi    <- max_var(theta, alpha, hi, C_list, X_list, sig_list, w, trace = FALSE)
  v_target<- (theta[unlist(C_list)!=0]/(qnorm(pwr_target) + qnorm(1-alpha/2))^2)^2
  guess   <- round(max(v_hi / v_target * hi))
  pwr_new_guess <- max_power(theta, alpha, guess, C_list, X_list, sig_list, w, trace = FALSE)

  cat("\ninitial guess = ", guess, " with power = ", min(pwr_new_guess))

  if (min(pwr_new_guess) < pwr_target) lo <- guess
  if (min(pwr_new_guess) > pwr_target) hi <- guess

  while (lo <= hi) {
    mid <- lo + round((hi - lo) / 2)
    cat("\nlo = ", lo)
    cat("  hi = ", hi)
    cat(" mid = ", mid)
    pwr_new <- max_power(theta, alpha, mid, C_list, X_list, sig_list, w, trace = FALSE)
    if (pwr_target < min(pwr_new)) hi = mid - 1
    if (pwr_target > min(pwr_new)) lo = mid + 1
    cat("\ntarget: ", pwr_target)
    cat("  minpwr: ", min(pwr_new))
  }

  return(mid)
}
