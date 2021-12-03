### function to generate cov and X mats based on model formula

create_data <- function(formula,
                        family,
                        data,
                        theta,
                        Z,
                        D,
                        C,
                        var_par=NULL,
                        verbose=TRUE){
  
  # parse the formula
  if(length(formula)==3)stop("formula should not have dependent variable.")
  if(!all(all.vars(formula)%in%colnames(data)))stop(paste0("Not all variables named in data :",
                                                           paste0(all.vars(formula)[!all.vars(formula)%in%colnames(data)],collapse=" ")))
  mf1 <- formula[[2]]
  
  f <- as.formula(formula)
  vars1 <- all.vars(formula)
  
  # extract list of functions and list of variables for each function
  f1 <- list()
  v1 <- list()
  if(mf1[[1]]=="+"){
    sym <- "+"
    iter <- 0
    while(sym == "+"){
      iter <- iter+1
      sym <- as.character(mf1[[c(rep(2,iter-1),1)]])
      
      if(sym == "+"){
        idx <- c(rep(2,iter-1),3,1)
      } else {
        idx <- c(rep(2,iter-1),1)
      }
      
      f1[[iter]] <- as.character(mf1[[idx]])
      v1[[iter]] <- all.vars(mf1[[idx[1:(length(idx)-1)]]])
    }
  } else {
    f1[[1]] <- c(as.character(mf1[[1]]))
    v1[[1]] <- all.vars(mf1)
  }
  
  mod <- gen_model_string(f1,v1)
  if(verbose)cat("MEAN FUNCTION:\n")
  if(verbose)print(mod[[1]])
  if(length(theta)!=(length(f1)+1))stop("Length theta not equal to number of lists of parameters")
  
  
  #check functions are in our 
  #print(unlist(f1))
  if(!all(unlist(f1)%in%c("Exp","Log","Linear","Algebraic")))stop("Functions not recognised")
  
  #build the matrix!
  X <- matrix(1,nrow=nrow(data),ncol=1)
  m0 <- rep(theta[[1]][[1]],nrow(data))
  
  #cycle through the functions and add them in order and add them in columns
  for(i in 1:length(f1)){
    idx <- (length(f1)-i+1)
    
    if(is(theta[[i+1]],"list"))theta[[i+1]] <- unlist(theta[[i+1]])
    
    if(f1[[idx]] %in% c("Exp","Log")){
      m1 <- rep(0,nrow(data))
      for(j in 1:length(v1[[idx]])){
        m1 <- m1 + theta[[i+1]][[j+1]] * data[,v1[[idx]][j]]
      }
      if(f1[[idx]] == "Exp"){
        X <- cbind(X,matrix(exp(m1),ncol=1))
        for(j in 1:length(v1[[idx]])){
          X <- cbind(X,matrix(theta[[i+1]][[1]]*data[,v1[[idx]][j]]*exp(m1),ncol=1))
        }
        
        m0 <- m0 + theta[[i+1]][[1]]*exp(m1)
      }
      if(f1[[idx]] == "Log"){
        X <- cbind(X,matrix(log(m1),ncol=1))
        for(j in 1:length(v1[[idx]])){
          X <- cbind(X,matrix(theta[[i+1]][[j+1]]*data[,v1[[idx]][j]]/m1,ncol=1))
        }
        m0 <- m0 + theta[[i+1]][[1]]*log(m1)
      }
    }
    
    if(f1[[idx]] %in% c("Linear")){
      m1 <- rep(0,nrow(data))
      for(j in 1:length(v1[[idx]])){
        m1 <- m1 + theta[[i+1]][[j]] * data[,v1[[idx]][[j]]]
      }
      for(j in 1:length(v1[[idx]])){
        X <- cbind(X,matrix(data[,v1[[idx]][j]],ncol=1))
      }
      m0 <- m0 + m1
    }
    
    if(f1[[idx]] %in% c("Algebraic")){
      m1 <- rep(0,nrow(data))
      for(j in 1:length(v1[[idx]])){
        m1 <- m1 + theta[[i+1]][[(j*2 - 1)]] * data[,v1[[idx]][j]] ^ theta[[i+1]][[(j*2)]]
      }
      for(j in 1:length(v1[[idx]])){
        X <- cbind(X,matrix(-(data[,v1[[idx]][j]]^theta[[i+1]][[(j*2)]])/(m1^2),ncol=1))
        X <- cbind(X,matrix(-(theta[[i+1]][[(j*2 - 1)]]*data[,v1[[idx]][j]]^theta[[i+1]][[(j*2)]] * log(data[,v1[[idx]][j]]))/(m1^2),ncol=1))
      }
      m0 <- m0 + 1/(1+m1)
    }
    
    
  }
  
  S <- gen_cov_mat(m0,
                   Z = Z,
                   D = D,
                   family=family,
                   var_par = var_par)
  
  return(list(C,X,S))
  
}

gen_model_string <- function(f1,v1){
  count <- 1
  str <- "theta[[1]][1]"
  #par_vars <- c()
  for(i in 1:length(f1)){
    count <- 1
    idx <- (length(f1)-i+1)
    if(f1[[idx]] %in% c("Exp","Log","Linear")){
      #count <- count+1
      
      if(f1[[idx]] %in% c("Exp","Log")){
        str <- paste0(str," + theta[[",i+1,"]][",count,"]*")
        count <- count +1
      } else {
        str <- paste0(str," + ")
      }
      
      str1 <- ""
      for(j in 1:length(v1[[idx]])){
        str1 <- paste0(str1,"theta[[",i+1,"]][",count,"]*",v1[[idx]][j])
        count <- count +1
        #par_vars <- c(par_vars,count)
        if(j != length(v1[[idx]])){
          str1 <- paste0(str1," + ")
        }
      }
      if(f1[[idx]]=="Exp"){
        str <- paste0(str,"Exp(",str1,")")
      }
      if(f1[[idx]]=="Log"){
        str <- paste0(str,"Log(",str1,")")
      }
      if(f1[[idx]]=="Linear"){
        str <- paste0(str,str1)
      }
    }
    if(f1[[idx]] %in% c("Algebraic")){
      str1 <- ""
      for(j in 1:length(v1[[idx]])){
        #count <- count +1
        str1 <- paste0(str1,"theta[[",i+1,"]][",count,"]*",v1[idx][j],"^theta[[",i+1,"]][",count+1,"]")
        #par_vars <- c(par_vars,count)
        count <- count +1
        if(j != length(v1[[idx]])){
          str1 <- paste0(str1," + ")
        }
      }
      
      str <- paste0(str," + 1/(1+",str1,")")
    }
  }
  return(list(str,count))
}


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
  pars[1]*exp(-x/pars[2])
}

indicator <- function(x, pars){
  pars[1]*I(x==0)
}

exp_power <- function(x,pars){
  pars[1]^(x)
}
