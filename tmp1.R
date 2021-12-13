require(latex2exp)

x_mat <- function(t,J,nper=1){
  X <- kronecker(rep(1,J),diag(t))
  XJ <- kronecker(diag(J),rep(1,t))
  int <- c()
  for(i in 1:J){
    int <- c(int,rep(c(rep(0,t-(i-1)),rep(1,i-1)),nper))
  }
  if(nper>1){
    XJ <- kronecker(diag(nper),XJ)
    X <- kronecker(rep(1,nper),X)
  }
  X <- cbind(X,XJ)
  
  X <- cbind(X, int)
  return(X)
}

# the example with six time periods and seven clusters with 10 people per clusters
nT <- 6
nJ <- nT+1
M <- 10

df <- expand.grid(t=1:nT,J = 1:nJ)
# df2 <- df
# df2$J <- df2$J+nJ
# df <- rbind(df,df2)
df <- df[rep(1:nrow(df),each=M),]

X <- x_mat(nT,nJ)
# X2 <- X
X <- X[rep(1:nrow(X),each=M),]
#some random parameters for the example, assume cluster means are zero
beta <- c(rnorm(nT),rnorm(nJ-3),1.75)
#data <- as.data.frame(X[,c(1:(nT),(nT+1):(nT+3),(nT+5):ncol(X))])
# X[,1] <- 1
# data <- as.data.frame(X[,c(1,13)])

data <- as.data.frame(X)
colnames(data) <- c(paste0("t",1:(nT)),paste0("cl",c(1:nJ)),"int")
data$cons <- 1

#for the power example
data <- data[!(data$cl1==1|data$cl7==1),]
df <- df[df$J%in%c(2:6),]

ss <- 75
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.8)),
                     list("indicator",c(0.01))
                   ),
                   parallel = TRUE)

dat1 <- create_data(formula(paste0("~Linear(cons, ",paste0("t",1:(nT-1),collapse=", "),",",paste0("cl",c(2:5),collapse=", "),", int)")),
                    family=gaussian(),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-3),1)),
                    var_par = 1)

sample_size2(theta = beta, 
             alpha = 0.05,
             pwr_target = 0.8,
             C_list=list(dat1[[1]]),
             X_list=list(dat1[[2]]),
             sig_list=list(dat1[[3]]),w=NULL)

d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]),
                  tol=1e-10,
                  trace=T)

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d1),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
#print(head(Xq))
#colnames(Xq)[(ncol(Xq)-3):ncol(Xq)] <- c("int","incl","cl","t")
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
#print(Xqa)
pcl <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)


# wrapper to make it easier to run lots of examples

optim_cl_des <- function(family,
                         Sout,
                         idx_in,
                         data,
                         nT=10,nJ=11,M=20){
  # optimal design for log_binomial
  
  
  
  return(pcl)
}

# EXAMPLES 1: COMPARISON OF FUNCTION FORM OF GLM


pcl1 <- optim_cl_des(gaussian(link="identity"),
                     Sout = Sout,
                     idx_in = idx_in,
                     data=data,
                     nT=nT,nJ=nJ,M=M)
pcl1 <- pcl1+ggtitle(TeX(r'($\rho =.05$, $R=0.1$)'))

