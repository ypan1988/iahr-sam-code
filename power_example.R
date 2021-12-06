#power 

power <- function(theta,alpha,pwr,m,C_list,X_list,sig_list){
  #randomly generate starting position
  d <- sample(c(rep(1,m),rep(0,nrow(X_list[[1]])-m)),nrow(X_list[[1]]))
  idx_in <- which(d==1)
  
  
  #do in loop for all theta where C[i]!=0
  v0[i] <- ((qnorm(power[i]) +qnorm(1-alpha/2))/theta[i])^2 
  
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