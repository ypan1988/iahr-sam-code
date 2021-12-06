#power 

power <- function(theta,alpha,pwr,m,C_list,X_list,sig_list){
  #randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)
  
  #OPTION1: RUN grad robust once, output the power for given sample size
  # grad robsut returns variance v0
  # calculate power as;
  pow[i] <- pnorm(sqrt(theta[1]/sqrt(v0[i])) - qnorm(1-alpha/2))
  
  
  #OPTION 2: Output the sample size (and design), for a given power
  #do in loop for all theta where C[i]!=0
  v0 <- (theta/((qnorm(pow) + qnorm(1-alpha/2))^2))^2 
  
  #aim
  all(v1<v0)
  
  # call grad-robust
  
  # A. check variance at m=m, v1. then
  # if v1 -v0 < -tol, then then set m = m-1
  # if v1 -v0 > +tol then set m = m+1
  # run grad robust again,
  
  # other wise if abs(v1-v0)<tol, stop
  #
  #
  
  
}