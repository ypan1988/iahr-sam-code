####################################################################################
# TEST EXAMPLES ####################################################################
####################################################################################
require(Matrix)
require(ggplot2)
require(Rcpp)
require(RcppArmadillo)

setwd("C:/Users/pany/Documents/Projects/iahr-sam-code")
sourceCpp("gd_search.cpp")
source("gd_search.R")
source("create_data.R")

#starting vector
d <- sample(c(rep(1,100),rep(0,n_dim^2-100)),n_dim^2)
idx_in <- which(d==1)

n_dim <- 20
xpos <- seq(-1,1,length.out = 2*n_dim+1)
xpos <- xpos[(1:n_dim)*2]
pos <- expand.grid(x=xpos,y=xpos)

data <- data.frame(x=-sqrt(pos$x^2 + pos$y^2),y=rnorm(n_dim^2),z=runif(n_dim^2))

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,0.1))
                   ))
# optimal design for exponential models
dat1 <- create_data(~Exp(x)+Linear(y,z),
            family=binomial(link="logit"),
            data=data,
            theta=list(
              list(1),
              list(0.5,5),
              list(1,1)
            ),
            Z=Sout[[1]],
            D=Sout[[2]],
            C=matrix(c(0,0.5,0.5,0,0)),
            var_par = 1)



sig <- dat1[[2]]
u <- dat1[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- 0.5*exp(data$x*2)


p_exp <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

#optimal design for linear model
dat2 <- create_data(~Linear(x,y,z),
                   family=binomial(link="logit"),
                   data=data,
                   theta=list(
                     list(1),
                     list(0.5,1,1)
                   ),
                   Z=Sout[[1]],
                   D=Sout[[2]],
                   C=matrix(c(0,1,0,0)),
                   var_par = 1)


# robust for both previous designs
sig <- dat2[[2]]
u <- dat2[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- 0.5*data$x

p_linear <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

######
# test robust
# compare exp and linear models

A_list <- list(solve(dat1[[2]][idx_in,idx_in]),solve(dat2[[2]][idx_in,idx_in]))
sig_list <- list(dat1[[2]],dat2[[2]])
u_list <- list(dat1[[1]],dat2[[1]])


d <- grad_robust(idx_in,
                 A_list=A_list,
                 sig_list=sig_list,
                 u_list=u_list,
                 tol=1e-10,
                 trace=T)


pos$isin <- 0
pos[d,'isin'] <- 1
pos$int <- data$x#0.5*exp(data$x*2)


p_robust <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

##################################
# graphical comparison

ggpubr::ggarrange(p_exp,p_linear,p_robust)


####################################################################################
# TEST EXAMPLES 2 ####################################################################
####################################################################################
#starting vector
# d <- sample(c(rep(1,100),rep(0,n_dim^2-100)),n_dim^2)
# idx_in <- which(d==1)
# 
# n_dim <- 20
# xpos <- seq(-1,1,length.out = 2*n_dim+1)
# xpos <- xpos[(1:n_dim)*2]
# pos <- expand.grid(x=xpos,y=xpos)
# 
# data <- data.frame(x=-sqrt(pos$x^2 + pos$y^2),y=rnorm(n_dim^2),z=runif(n_dim^2))

# Sout <- gen_re_mat(pos,
#                    dims = list(c(1,2)),
#                    funs = list(
#                      list("exponential",c(1,0.1))
#                    ))

# optimal design for log_binomial
dat1 <- create_data(~Linear(x,y,z),
                    family=binomial(link="log"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,1,1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,1,0,0)),
                    var_par = 1)



sig <- dat1[[2]]
u <- dat1[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- 0.5*data$x


p_logb <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

#optimal design for logit-binomial model
dat2 <- create_data(~Linear(x,y,z),
                    family=binomial(link="logit"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,1,1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,1,0,0)),
                    var_par = 1)


# robust for both previous designs
sig <- dat2[[2]]
u <- dat2[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- 0.5*data$x

p_logitb <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

#optimal design for identity-binomial model
dat3 <- create_data(~Linear(x,y,z),
                    family=binomial(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,1,1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,1,0,0)),
                    var_par = 1)


# robust for all previous designs
sig <- dat3[[2]]
u <- dat3[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d3 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d3+1,'isin'] <- 1
pos$int <- 0.5*data$x

p_idenb <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

######
# test robust
# compare exp and linear models

A_list <- list(solve(dat1[[2]][idx_in,idx_in]),
               solve(dat2[[2]][idx_in,idx_in]),
               solve(dat3[[2]][idx_in,idx_in]))
sig_list <- list(dat1[[2]],
                 dat2[[2]],
                 dat3[[2]])
u_list <- list(dat1[[1]],
               dat2[[1]],
               dat3[[1]])


d <- grad_robust(idx_in,
                 A_list=A_list,
                 sig_list=sig_list,
                 u_list=u_list,
                 tol=1e-10,
                 trace=T)


pos$isin <- 0
pos[d,'isin'] <- 1
pos$int <- 0.5*data$x


p_robust2 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")

##################################
# graphical comparison

ggpubr::ggarrange(p_logb,p_logitb,p_idenb,p_robust2)
