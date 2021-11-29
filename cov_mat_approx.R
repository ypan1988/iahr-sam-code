# function to produce covariance matrix 
# for GLMM
# inputs:
# Xb = f(X;b) the fixed effects part of the mean function evaluated at a value of beta = b
# Z the design matrix for the random effects components
# D the covariance matrix of the random effects
# model is the type of model

logit <- function(x){
  exp(x)/(1+exp(x))
}

gen_cov_mat <- function(Xb,
                        Z,
                        D,
                        family,
                        var_par=NULL){
  # assume random effects value is at zero
  f <- family
  if(!f[1]%in%c("poisson","binomial","gaussian","gamma"))stop("family must be one of Poisson, Binomial, Gaussian, Gamma")
  
  if(f[1]=="poisson"){
    if(f[2]=="log"){
      W <- diag(1/(exp(Xb)))
    }
    if(f[2]=="identity"){
      W <- diag(exp(Xb))
    }
  }
  
  if(f[1]=="binomial"){
    if(f[2]=="logit"){
      W <- diag(1/(logit(Xb)*(1-logit(Xb))*length(Xb)))
    }
    if(f[2]=="log"){
      W <- diag((1-logit(Xb))/(logit(Xb)*length(Xb)))
    }
    if(f[2]=="identity"){
      W <- diag((logit(Xb)*(1-logit(Xb)))/length(Xb))
    }
    if(f[2]=="probit"){
      W <- diag((pnorm(Xb)*(1-pnorm(Xb)))/(dnorm(Xb)*length(Xb)))
    }
  }
  
  if(f[1]=="gaussian"){
    if(f[2]=="identity"){
      if(is.null(var_par))stop("For gaussian(link='identity') provide var_par")
      W <- var_par*diag(length(Xb))
    }
    if(f[2]=="log"){
      if(is.null(var_par))stop("For gaussian(link='log') provide var_par")
      W <- diag(var_par/exp(Xb))
    }
  }
  
  if(f[1]=="gamma"){
    if(f[2]=="inverse"){
      if(is.null(var_par))stop("For gamma(link='inverse') provide var_par")
      W <- var_par*diag(length(Xb))
    }
  }
  
  if(is(D,"numeric")){
    S <- W + D * Z %*% t(Z)
  } else {
    S <- W + Z %*% D %*% t(Z)
  }
  
  
  return(S)
  
}


################################
# Functions for generating random effects design matrices
# df should be a data frame with columns corresponding to the position of each observation in each dimension
# dims should be a list indicating which columns of df go with which funciton (e.g. list(c(1,2),c(3))) indicates 2D first position (spatial)
# and 1D second poisiton (time)
# funs should be a list of lists of functions indicating the covariance function for each element of the list in dims. current options are
# "exponential" and "indicator" (more to be added), sytax should be,e.g. 
# list(list("exponential",pars=c(1,2)),list("indicator",pars=c(1,2)))
gen_re_mat <- function(df,
                       dims,
                       funs){
  if(!is(dims,"list"))stop("dims should be list") 
  if(!is(funs,"list"))stop("funs should be list") 
  if(length(dims)!=length(funs))stop("dims and funs should be same length")
  
  nD <- length(dims)
  rownames(df) <- NULL
  df_nodup <- df[!duplicated(df),]
  zdim2 <- nrow(df_nodup)
  Z <- matrix(0,nrow=nrow(df),ncol=zdim2)
  for(i in 1:zdim2){
    mat <- unlist(sapply(1:nrow(df),function(j)isTRUE(all.equal(df[j,],df_nodup[i,],check.attributes=FALSE,use.names=FALSE))))
    Z[mat,i] <- 1
  }
  
  Dlist <- list()
  for(i in 1:nD){
    Dlist[[i]] <- as.matrix(dist(df_nodup[,dims[[i]]],upper = TRUE, diag=TRUE))
  }
  
  D <- matrix(1, nrow=zdim2,ncol=zdim2)  
  
  for(i in 1:nD){
    D <- D*do.call(funs[[i]][[1]],list(Dlist[[i]],funs[[i]][[2]]))
  }
  
  return(list(Z,D))
  
}

exponential <- function(x, pars){
  pars[1]*exp(-x*pars[2])
}

indicator <- function(x, pars){
  pars[1]*I(x==0)
}


###
# the X matrix for the cluster example

x_mat <- function(t,J){
  X <- kronecker(rep(1,J),diag(t))
  XJ <- kronecker(diag(J),rep(1,t))
  X <- cbind(X,XJ)
  int <- c()
  for(i in 1:J){
    int <- c(int,rep(0,t-(i-1)),rep(1,i-1))
  }
  X <- cbind(X, int)
  return(X[,2:ncol(X)])
}




###############################################
# examples
# 
# the example with six time periods and seven clusters with 10 people per clusters
df <- expand.grid(t=1:6,J = 1:7)
df <- df[rep(1:nrow(df),each=10),]

X <- x_mat(6,7)
X <- X[rep(1:nrow(X),each=10),]
#some random parameters for the example, assume cluster means are zero
beta <- matrix(c(rnorm(5),rep(0,7),0.5),ncol=1)
Xb <- X%*%beta

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exponential",c(1,0.5)),
                     list("indicator",c(0.5))
                   ))

# linear model example
S <- gen_cov_mat(Xb,
                 Z = Sout[[1]],
                 D = Sout[[2]],
                 family=gaussian(link="identity"),
                 var_par = 1)

# binomial-logit model example
S <- gen_cov_mat(Xb,
                 Z = Sout[[1]],
                 D = Sout[[2]],
                 family=binomial(link="logit"))



#########################################################
#an example of geospatial design, no temporal dimension
n_dim <- 10
xpos <- seq(-1,1,length.out = 2*n_dim+1)
xpos <- xpos[(1:n_dim)*2]
pos <- expand.grid(x=xpos,y=xpos)

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,0.5))
                   ))

# poisson model example
X <- matrix(c(rep(1,nrow(pos)),rnorm(nrow(pos))),ncol=2)
beta <- matrix(c(-0.5,0.5),ncol=1)
Xb <- X%*%beta

S <- gen_cov_mat(Xb,
                 Z = Sout[[1]],
                 D = Sout[[2]],
                 family=poisson())
