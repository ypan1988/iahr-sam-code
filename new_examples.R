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

dat <- create_data(~Exp(x)+Linear(y,z),
            family=binomial(link="log"),
            data=data,
            theta=list(
              list(1),
              list(0.5,2),
              list(1,1)
            ),
            Z=Sout[[1]],
            D=Sout[[2]],
            C=matrix(c(0,1,1,0,0)),
            var_par = 1)


#starting vector
d <- sample(c(rep(1,100),rep(0,n_dim^2-100)),n_dim^2)
idx_in <- which(d==1)
sig <- dat[[2]]
u <- dat[[1]]
A <- solve(sig[idx_in,idx_in])

#d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)


pos$isin <- 0
pos[d2+1,'isin'] <- 1
pos$int <- 0.5*exp(data$x*2)


ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention area")+
  labs(x="x",y="y")
