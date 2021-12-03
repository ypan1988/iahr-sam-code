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
  
  X_id <- 1:nrow(X_list[[1]])
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
    #check to see whether to remove any covariates
    new_val <- out[[1]]
    diff <- new_val - val
    
    if(diff>0){
      idx_in <- out[[2]]
      for(j in 1:length(A_list)){
        #see if any columns don't have any observations in
        rss <- colSums(X_list[[j]][idx_in,])
        if(any(rss==0)){
          if (trace) {
            cat("\nIter: ",i)
            cat(" ",diff)
            cat("  X not full rank...")
          }
          #which columns are empty
          idx_rm <- which(rss==0)
          #which rows to keep
          if(length(idx_rm)==1){
            nrm_rm <- which(X_list[[j]][,idx_rm]==0)
          } else {
            nrm_rm <- which(apply(X_list[[j]][,idx_rm],1,function(x)all(x==0)))
          }
          new_idx <- rep(NA,nrow(X_list[[j]]))
          #remove rows and columns corresponding to point to remove
          X_list[[j]] <- X_list[[j]][nrm_rm,-idx_rm]
          #fix IDs to new X
          new_idx[nrm_rm] <- 1:nrow(X_list[[j]])
          C_list[[j]] <- matrix(C_list[[j]][-idx_rm,],ncol=1)
          if(trace)cat(paste0("Removing column ",idx_rm))
          
          idx_in <- new_idx[idx_in]
          X_id <- X_id[nrm_rm]
          
          #replace all the matrices
          sig_list[[j]] <- sig_list[[j]][nrm_rm,nrm_rm]
          A_list[[j]] <- solve(sig_list[[j]][idx_in,idx_in])
          M <- t(X_list[[j]]) %*% solve(sig_list[[j]]) %*% X_list[[j]]
          cM <- t(C_list[[j]]) %*% solve(M)
          u_list[[j]] <- cM %*% t(X_list[[j]])
        } else {
          
          if (trace) {
            cat("\nIter: ",i)
            cat(" ",diff)
          }
          A_list[[j]] <- out[[3]][[j]]
        }
      }
    }
    
   
    
  }
  
  #return variance
  if(trace){
    var_vals <- c()
    for(i in 1:length(X_list)){
      var_vals[i] <- tryCatch(optim_fun(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]),error=function(i){NA})
    }
    cat("\nVariance for individual model(s):\n")
    print(var_vals)
    if(length(A_list)>1){
      cat("\n Weighted average variance:\n")
      print(sum(var_vals*c(w)))
    }
  }
  
  
  return(X_id[idx_in])
}

# grad_robust <- function(idx_in,
#                         C_list,
#                         X_list,
#                         sig_list,
#                         w=NULL,
#                         tol=1e-9,
#                         trace = TRUE){
#   
#   if(is.null(w))w <- rep(1/length(sig_list),length(sig_list))
#   if(sum(w)!=1)w <- w/sum(w)
#   if(!is(w,"matrix"))w <- matrix(w,ncol=1)
#   if(!all(unlist(lapply(C_list,function(x)is(x,"matrix")))))stop("All C_list must be matrices")
#   if(!all(unlist(lapply(sig_list,function(x)is(x,"matrix")))))stop("All sig_list must be matrices")
#   if(!all(unlist(lapply(X_list,function(x)is(x,"matrix")))))stop("All X_list must be matrices")
#   if((length(C_list)!=length(X_list))|length(X_list)!=length(sig_list))stop("Lists must be same length")
#   
#   A_list <- list()
#   u_list <- list()
#   for(i in 1:length(sig_list)){
#     A_list[[i]] <- solve(sig_list[[i]][idx_in,idx_in])
#     M <- t(X_list[[i]]) %*% solve(sig_list[[i]]) %*% X_list[[i]]
#     cM <- t(C_list[[i]]) %*% solve(M)
#     u_list[[i]] <- cM %*% t(X_list[[i]])
#   }
#   
#   new_val_vec <- matrix(sapply(1:length(A_list),function(i)obj_fun(A_list[[i]], u_list[[i]][idx_in])),nrow=1)
#   
#   new_val <- as.numeric(new_val_vec %*% w)
#   diff <- 1
#   i <- 0
#   while(diff > tol){
#     val <- new_val
#     i <- i + 1
#     out <- choose_swap_robust(idx_in,A_list,sig_list,u_list, w)
#     print(out[[1]])
#     #check to see whether to remove any covariates
#     for(j in 1:length(sig_list)){
#       #rss <- rowSums(t(X_list[[j]][out[[2]],])%*%out[[3]][[j]]%*%X_list[[j]][out[[2]],])
#       rss <- colSums(X_list[[j]][out[[2]],])
#       if(any(rss==0)){
#         
#         if (trace) {
#           cat("\nIter: ",i)
#           cat(paste0("  Removing column ",which(rss==0)))
#         }
#         
#         new_idx <- rep(NA,nrow(X_list[[j]]))
#         nrm_idx <- which(X_list[[j]][,which(rss==0)]==0)
#         # print(length(nrm_idx))
#         # 
#         # print(length(1:nrow(X_list[[j]])))
#         
#         X_list[[j]] <- X_list[[j]][X_list[[j]][,which(rss==0)]==0,]
#         new_idx[nrm_idx] <- 1:nrow(X_list[[j]])
#         X_list[[j]] <- X_list[[j]][,-which(rss==0)]
#         C_list[[j]] <- matrix(C_list[[j]][-which(rss==0),],ncol=1)
#         sig_list[[j]] <- sig_list[[j]][nrm_idx,nrm_idx]
#         M <- t(X_list[[j]]) %*% solve(sig_list[[j]]) %*% X_list[[j]]
#         cM <- t(C_list[[j]]) %*% solve(M)
#         u_list[[j]] <- cM %*% t(X_list[[j]])
#         idx_in <- new_idx[out[[2]]]
#         A_list[[j]] <- solve(sig_list[[j]][idx_in,idx_in])
#         
#       } else {
#         new_val <- out[[1]]
#         diff <- new_val - val
#         if (trace) {
#           cat("\nIter: ",i)
#           cat(" ",diff)
#         }
#         if(diff>0){
#           A_list <- out[[3]]
#           idx_in <- out[[2]]
#         }
#       }
#     }
#     
#     
#     
#     
#   }
#   
#   #return variance
#   if(trace){
#     var_vals <- c()
#     for(i in 1:length(X_list)){
#       var_vals[i] <- tryCatch(optim_fun(C_list[[i]],X_list[[i]][idx_in,],sig_list[[i]][idx_in,idx_in]),error=function(i){NA})
#     }
#     cat("\nVariance for individual model(s):\n")
#     print(var_vals)
#     if(length(A_list)>1){
#       cat("\n Weighted average variance:\n")
#       print(sum(var_vals*c(w)))
#     }
#   }
#   
#   
#   
#   return(idx_in)
# }


optim_fun <- function(C,X,S){
  if(!any(is(C,"matrix"),is(X,"matrix"),is(S,"matrix")))stop("C, X, S must be matrices")
  M <- t(X) %*% solve(S) %*% X
  # print(M)
  # print(eigen(M))
  # print(solve(M))
  val <- t(C) %*% solve(M) %*% C
  return(val)
}
