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

# the example with six time periods and seven clusters with 10 people per clusters
nT <- 6
nJ <- nT+1
M <- 30

df <- expand.grid(t=1:nT,J = 1:nJ)
df <- df[rep(1:nrow(df),each=M),]

X <- x_mat(nT,nJ)
X <- X[rep(1:nrow(X),each=M),]
#some random parameters for the example, assume cluster means are zero
beta <- c(rep(0,nT-1),rep(0,nJ-1),0.5)
data <- as.data.frame(X[,c(1:(nT-1),(nT+1):ncol(X))])
# X[,1] <- 1
# data <- as.data.frame(X[,c(1,13)])
colnames(data) <- c(paste0("t",1:(nT-1)),paste0("cl",2:nJ),"int")





# wrapper to make it easier to run lots of examples

optim_cl_des <- function(family,
                         Sout,
                         idx_in,
                         data,
                         nT=10,nJ=11,M=20){
  # optimal design for log_binomial
  
  dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:(nT-1),collapse=", "),",",paste0("cl",2:nJ,collapse=", "),", int)")),
                      family=family,
                      data=data,
                      theta=list(
                        list(1),
                        list(beta)
                      ),
                      Z=Sout[[1]],
                      D=Sout[[2]],
                      C=matrix(c(rep(0,(ncol(data))),1)),
                      var_par = 1)
  
  
  
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
  
  return(pcl)
}

# EXAMPLES 1: COMPARISON OF FUNCTION FORM OF GLM
ss <- 200
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exponential",c(1,0.5)),
                     list("indicator",c(0.5))
                   ))

pcl_logit <- optim_cl_des(binomial(link = "logit"),
                          Sout = Sout,
                          idx_in = idx_in,
                          data=data)

pcl_pois <- optim_cl_des(poisson(link="log"),
                          Sout = Sout,
                          idx_in = idx_in,
                         data=data)

pcl_lin <- optim_cl_des(gaussian(link="identity"),
                          Sout = Sout,
                          idx_in = idx_in,
                          data=data)


ggpubr::ggarrange(pcl_lin,pcl_logit,pcl_pois)


# EXAMPLES 1: COMPARISON OF FUNCTION FORM OF GLM
Sout2 <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.2)),
                     list("indicator",c(0.05))
                   ))


pcl_logit2 <- optim_cl_des(binomial(link = "log"),
                          Sout = Sout2,
                          idx_in = idx_in,
                          data=data,
                          nT=nT,nJ=nJ,M=M)

pcl_pois2 <- optim_cl_des(poisson(link="log"),
                         Sout = Sout2,
                         idx_in = idx_in,
                         data=data,
                         nT=nT,nJ=nJ,M=M)

pcl_lin2 <- optim_cl_des(gaussian(link="identity"),
                        Sout = Sout2,
                        idx_in = idx_in,
                        data=data,
                        nT=nT,nJ=nJ,M=M)


ggpubr::ggarrange(pcl_lin2,pcl_logit2,pcl_pois2)


