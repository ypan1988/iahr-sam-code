require(Matrix)
require(ggplot2)
require(Rcpp)
require(RcppArmadillo)

setwd("C:/Users/pany/Documents/Projects/iahr-sam-code")
sourceCpp("gd_search.cpp")
source("gd_search.R")
source("test_examples.R")

########################################
##### Example with Cluster & T data ####
########################################

#test data
dat <- gen_clmat()

#starting vector
d <- sample(c(rep(1,100),rep(0,320)),420)
idx_in <- which(d==1)
sig <- dat[[1]]
u <- dat[[2]]
A <- solve(sig[idx_in,idx_in])

d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)

sum(d-1 == d2)
length(intersect(d-1, d2))

# DO NOT RUN
# require(microbenchmark)
# mres1 <- microbenchmark(grad(idx_in,A,sig,u,tol=1e-10,F),
#                         Grad(idx_in-1,A,sig,u,tol=1e-10,F),
#                         check=NULL,times=100)

## plot results for d
Xq <- x_mat(6,7)
Xq <- Xq[rep(1:nrow(Xq),each=10),]
Xq <- cbind(Xq,0)
Xq[sort(d),ncol(Xq)] <- 1
Xq <- cbind(Xq,rep(1:7,each=60))
Xq <- cbind(Xq,rep(rep(1:6,each=10),7))
Xq <- as.data.frame(Xq)
Xqa <- aggregate(Xq[,c('V14','int')],list(Xq$V15,Xq$V16),sum)

ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=V14))+
  geom_tile()+
  geom_label(aes(label=int/10))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:6)+
  scale_y_continuous(breaks = 1:7)

## plot results for d2
Xq <- x_mat(6,7)
Xq <- Xq[rep(1:nrow(Xq),each=10),]
Xq <- cbind(Xq,0)
Xq[sort(d2),ncol(Xq)] <- 1
Xq <- cbind(Xq,rep(1:7,each=60))
Xq <- cbind(Xq,rep(rep(1:6,each=10),7))
Xq <- as.data.frame(Xq)
Xqa <- aggregate(Xq[,c('V14','int')],list(Xq$V15,Xq$V16),sum)

ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=V14))+
  geom_tile()+
  geom_label(aes(label=int/10))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:6)+
  scale_y_continuous(breaks = 1:7)

########################################
#####       Another example        #####
########################################

#test data
dat <- gen_stmat(25)

#random starting vector
d <- sample(c(rep(1,200),rep(0,425)),625)
idx_in <- which(d==1)

sig <- dat[[1]]
s1 <- solve(sig)
u <- dat[[2]]
A <- solve(sig[idx_in,idx_in])
X <- dat[[4]]
pos <- dat[[3]]

d <- grad(idx_in,A,sig,u,tol=1e-10)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10)

pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- X[,2]


ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

########################################
##### Example for the non-linear   #####
########################################

obj_fun_val <- function(d){
  solve(u[d]%*%solve(sig[d,d])%*%matrix(u[d],ncol=1))
}

#test data for non linear
dat <- gen_stmat_non_linear(cov_pars = c(0.5,1),theta = c(1,0.5,3))

#starting vector
d <- sample(c(rep(1,100),rep(0,300)),400)
idx_in <- which(d==1)
sig <- dat[[1]]
s1 <- solve(sig)
u <- dat[[2]]
A <- solve(sig[idx_in,idx_in])
X <- dat[[4]]
pos <- dat[[3]]

d <- grad(idx_in,A,sig,u,tol=1e-20,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-20,T)

length(intersect(d-1, d2))

pos$isin <- 0
pos[d,'isin'] <- 1
pos$int <- dat[[5]]

ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- X[,2]

ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")


########################################
##### Example for the robust       #####
########################################

num_dat <- 5
wts <- rep(1/num_dat, num_dat)
dat <- lapply(1:num_dat, function(i) gen_stmat_non_linear(cov_pars = c(0.5,1),theta = c(1,0.5,3)))
d <- sample(c(rep(1,100),rep(0,300)),400)
idx_in <- which(d==1)
sig_list <- lapply(1:num_dat, function(i) dat[[i]][[1]])
u_list <- lapply(1:num_dat, function(i) dat[[i]][[2]])
A_list <- lapply(1:num_dat, function(i) solve(sig_list[[i]][idx_in,idx_in]))
X_list <- lapply(1:num_dat, function(i) dat[[i]][[4]])
C_list <- lapply(1:num_dat, function(i) dat[[i]][[6]])

grad_robust_old <- function(idx_in,A_list,sig_list,u_list,weights,tol=1e-9, trace = TRUE){
  new_val_vec <- sapply(1:length(A_list),function(i)obj_fun(A_list[[i]], u_list[[i]][idx_in]))
  new_val <- as.numeric(new_val_vec %*% weights)
  diff <- 1
  i <- 0
  while(diff > tol){
    val <- new_val
    i <- i + 1
    out <- choose_swap_robust(idx_in,A_list,sig_list,u_list, weights)
    new_val <- out[[1]]
    diff <- new_val - val
    if(diff>0){
      A_list <- out[[3]]
      idx_in <- out[[2]]
    }
    if (trace) {
      cat("\nIter: ",i)
      cat(" ",diff)
    }
  }
  return(idx_in)
}

system.time(grad_robust_old(idx_in,A_list, sig_list,u_list,wts,tol=1e-20,T))
system.time(GradRobust(num_dat, idx_in-1,do.call(rbind, A_list),
                       do.call(rbind, sig_list),do.call(cbind,u_list),wts,tol=1e-20,T)
)

d = grad_robust(idx_in,C_list,X_list,sig_list,wts,tol=1e-20,T)
d2= grad_robust2(idx_in,C_list,X_list,sig_list,wts,tol=1e-20,T)

sum(d == d2)
