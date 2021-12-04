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
n_dim <- 30

xpos <- seq(-1,1,length.out = 2*n_dim+1)
xpos <- xpos[(1:n_dim)*2]
pos <- expand.grid(x=xpos,y=xpos)
xdist <- sqrt(pos$x^2 + pos$y^2)
pos <- pos[xdist<=1,]

data <- data.frame(intc = 1, x=-sqrt(pos$x^2 + pos$y^2))

d <- sample(c(rep(1,100),rep(0,nrow(data)-100)),nrow(data))
idx_in <- which(d==1)

## compare functional form

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,0.2))
                   ))
# optimal design for exponential models
dat1 <- create_data(~Linear(intc) + Exp(x),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,5)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.5,0.5)),
                    var_par = 1)

d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]),
                  tol=1e-10,
                  trace=T)

pos$isin <- 0
pos[d1,'isin'] <- 1
pos$int <- 0.5*exp(data$x*5)

pexp1 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5 exp(-d/0.2)$)'))

# optimal design for exponential models
dat2 <- create_data(~Linear(intc) + Exp(x),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.5,0.5)),
                    var_par = 1)

d2 <- grad_robust(idx_in,
                  C=list(dat2[[1]]),
                  X_list=list(dat2[[2]]),
                  sig_list = list(dat2[[3]]))

pos$isin <- 0
pos[d2,'isin'] <- 1
pos$int <- 0.5*exp(data$x*1)

pexp2 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5 exp(-d/1)$)'))

#linear
dat3 <- create_data(~Linear(intc) + Linear(x),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(0.5),
                      list(-1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,1)),
                    var_par = 1)

d3 <- grad_robust(idx_in,
                  C=list(dat3[[1]]),
                  X_list=list(dat3[[2]]),
                  sig_list = list(dat3[[3]]))

pos$isin <- 0
pos[d3,'isin'] <- 1
pos$int <- 0.5+data$x

plin1 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5-d$)'))


#
test_cubic <- function(a,b,c){
  x <- seq(0,1,length.out=200)
  y <- a*x + b*x^2 + c*x^3
  qplot(x=x,y=y)
}

data$x2 <- data$x^2 
data$x3 <- data$x^3

dat4 <- create_data(~Linear(intc) + Linear(x,x2,x3),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(0.5),
                      list(0.5,1.5,-1.25)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.33,0.33,0.33)),
                    var_par = 1)

d4 <- grad_robust(idx_in,
                  C=list(dat3[[1]]),
                  X_list=list(dat3[[2]]),
                  sig_list = list(dat3[[3]]))

pos$isin <- 0
pos[d4,'isin'] <- 1
pos$int <- 0.5-0.5*data$x-1.5*data$x2 +1.25*data$x3

pcub1 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5-d-3d^2+2.5d^3$)'))



dr <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]],dat3[[1]],dat4[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]],dat3[[2]],dat4[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]],dat3[[3]],dat4[[3]]))

pos$isin <- 0
pos[dr,'isin'] <- 1
pos$int <- data$x

prob <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle("Robust")


ggpubr::ggarrange(pexp1,pexp2,plin1,pcub1,prob,nrow=3,ncol=2)

##############################################################
# EXAMPLE 2

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,0.05))
                   ))
# optimal design for exponential models
dat1 <- create_data(~Linear(intc) + Exp(x),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,4)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.5,0.5)),
                    var_par = 1)

d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]))

pos$isin <- 0
pos[d1,'isin'] <- 1
pos$int <- 0.5*exp(data$x*3)

pexp1 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5 exp(-d/0.25)$ and $\beta = 0.05$)'))

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,2))
                   ))
# optimal design for exponential models
dat2 <- create_data(~Linear(intc) + Exp(x),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,4)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.5,0.5)),
                    var_par = 1)

d2 <- grad_robust(idx_in,
                  C=list(dat2[[1]]),
                  X_list=list(dat2[[2]]),
                  sig_list = list(dat2[[3]]))

pos$isin <- 0
pos[d2,'isin'] <- 1
pos$int <- 0.5*exp(data$x*3)

pexp2 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'($0.5 exp(-d/0.25)$ and $\beta = 2$)'))

d3 <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]]))

pos$isin <- 0
pos[d3,'isin'] <- 1
pos$int <- 0.5*exp(data$x*3)

prob2 <- ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle("Robust")

ggpubr::ggarrange(pexp1,pexp2,prob2,nrow=2,ncol=2)

#######################################################
# two intervention locations, two time periods, one intervention only in the second period


#starting vector
n_dim <- 30

xpos <- seq(-1,1,length.out = 2*n_dim+1)
xpos <- xpos[(1:n_dim)*2]
pos <- expand.grid(x=xpos,y=xpos,t=1:2)


data <- data.frame(intc = 1, x1=-sqrt((pos$x-0.5)^2 + (pos$y-0.5)^2), x2=-sqrt((pos$x+0.5)^2 + (pos$y+0.5)^2), t=pos$t)
data$d <- 0
data[data$t==1,'d'] <- data[data$t==1,'x1']
data[data$t==2,'d'] <- apply(data[data$t==2,c('x1','x2')],1,max)

ss <- 100
d <- sample(c(rep(1,ss),rep(0,nrow(data)-ss)),nrow(data))
idx_in <- which(d==1)

## compare functional form
# as large try in parallel?

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2),3),
                   funs = list(
                     list("exponential",c(1,0.2)),
                     list("exp_power",c(0.2))
                   ),
                   parallel = TRUE)
# optimal design for exponential models
dat1 <- create_data(~Linear(intc) + Exp(d),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1),
                      list(0.5,4)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,0.5,0.5)),
                    var_par = 1)

d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]))

pos$isin <- 0
pos[d1,'isin'] <- 1
pos$int <- 0.5*exp(data$d*4)

ggplot(data=pos,aes(x=x,y=y,fill=int))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  facet_wrap(~t)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_c(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'(Two period, m=100. Intervention effect: $0.5 exp(-d/0.25)$. Cov function $0.8^{|t_i - t_j|}exp(-d_{ij}/2)$)'))

# dichotomous effect
data$int <- I(abs(data$d) < 0.3)*1

dat2 <- create_data(~Linear(intc,int),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(1,1)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(0,1)),
                    var_par = 1)

d2 <- grad_robust(idx_in,
                  C=list(dat2[[1]]),
                  X_list=list(dat2[[2]]),
                  sig_list = list(dat2[[3]]))

pos$isin <- 0
pos[d2,'isin'] <- 1
pos$int <- data$int

ggplot(data=pos,aes(x=x,y=y,fill=factor(int)))+
  geom_tile(alpha=0.8)+
  geom_point(data=pos[pos$isin==1,],aes(x=x,y=y))+
  facet_wrap(~t)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_blank())+
  scale_fill_viridis_d(name="Intervention effect")+
  labs(x="x",y="y")+
  ggtitle(TeX(r'(Two period, m=100. Dichotomous intervention effect. Cov function $0.2^{|t_i - t_j|}exp(-d_{ij}/2)$)'))
