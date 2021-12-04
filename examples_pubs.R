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
M <- 30

df <- expand.grid(t=1:nT,J = 1:nJ)
# df2 <- df
# df2$J <- df2$J+nJ
# df <- rbind(df,df2)
df <- df[rep(1:nrow(df),each=M),]

X <- x_mat(nT,nJ)
# X2 <- X
X <- X[rep(1:nrow(X),each=M),]
#some random parameters for the example, assume cluster means are zero
beta <- c(rep(0,nT),rep(0,nJ),0.5)
#data <- as.data.frame(X[,c(1:(nT),(nT+1):(nT+3),(nT+5):ncol(X))])
# X[,1] <- 1
# data <- as.data.frame(X[,c(1,13)])

data <- as.data.frame(X)
colnames(data) <- c(paste0("t",1:(nT)),paste0("cl",c(1:nJ)),"int")





# wrapper to make it easier to run lots of examples

optim_cl_des <- function(family,
                         Sout,
                         idx_in,
                         data,
                         nT=10,nJ=11,M=20){
  # optimal design for log_binomial
  
  dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                      family=family,
                      data=data,
                      theta=list(
                        list(beta)
                      ),
                      Z=Sout[[1]],
                      D=Sout[[2]],
                      C=matrix(c(rep(0,nT+nJ-1),1)),
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
ss <- 150
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.1)),
                     list("indicator",c(0.05))
                   ))

pcl1 <- optim_cl_des(gaussian(link="identity"),
                          Sout = Sout,
                          idx_in = idx_in,
                          data=data,
                          nT=nT,nJ=nJ,M=M)
pcl1 <- pcl1+ggtitle(TeX(r'($\rho =.05$, $R=0.1$)'))

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.1)),
                     list("indicator",c(0.25))
                   ))

pcl2 <- optim_cl_des(gaussian(link="identity"),
                     Sout = Sout,
                     idx_in = idx_in,
                     data=data,
                     nT=nT,nJ=nJ,M=M)
pcl2 <- pcl2+ggtitle(TeX(r'($\rho =0.25$, $R=0.1$)'))


Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.05))
                   ))

pcl3 <- optim_cl_des(gaussian(link="identity"),
                     Sout = Sout,
                     idx_in = idx_in,
                     data=data,
                     nT=nT,nJ=nJ,M=M)
pcl3 <- pcl3+ggtitle(TeX(r'($\rho =0.05$, $R=1$)'))

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.25))
                   ))

pcl4 <- optim_cl_des(gaussian(link="identity"),
                     Sout = Sout,
                     idx_in = idx_in,
                     data=data,
                     nT=nT,nJ=nJ,M=M)
pcl4 <- pcl3+ggtitle(TeX(r'($\rho =0.25$, $R=1$)'))

#robust
Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.1)),
                     list("indicator",c(0.05))
                   ))
dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.1)),
                     list("indicator",c(0.25))
                   ))
dat2 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.05))
                   ))
dat3 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.25))
                   ))
dat4 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]],dat3[[1]],dat4[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]],dat3[[2]],dat4[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]],dat3[[3]],dat4[[3]]))

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d1),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
pclr <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  ggtitle("Robust")

ggpubr::ggarrange(pcl1,pcl2,pcl3,pcl4,pclr,ncol=2,nrow = 3)


# EXAMPLES 1: COMPARISON OF FUNCTION FORM OF GLM
Sout2 <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.2)),
                     list("indicator",c(0.05))
                   ))


pcl_logit2 <- optim_cl_des(binomial(link = "logit"),
                          Sout = Sout2,
                          idx_in = idx_in,
                          data=data,
                          nT=nT,nJ=nJ,M=M)
pcl_logit2 <- pcl_logit2 + ggtitle("Binomial-logit")

pcl_pois2 <- optim_cl_des(binomial(link="log"),
                         Sout = Sout2,
                         idx_in = idx_in,
                         data=data,
                         nT=nT,nJ=nJ,M=M)
pcl_pois2 <- pcl_pois2 + ggtitle("Binomial-log")

pcl_bin2 <- optim_cl_des(binomial(link="identity"),
                          Sout = Sout2,
                          idx_in = idx_in,
                          data=data,
                          nT=nT,nJ=nJ,M=M)
pcl_bin2 <- pcl_bin2 + ggtitle("Binomial-identity")

pcl_lin2 <- optim_cl_des(gaussian(link="identity"),
                        Sout = Sout2,
                        idx_in = idx_in,
                        data=data,
                        nT=nT,nJ=nJ,M=M)
pcl_lin2 <- pcl_lin2 + ggtitle("Gaussian-identity")

dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=binomial(link = "logit"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
dat2 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=binomial(link = "log"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
dat3 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=binomial(link = "identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
dat4 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:3,5:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)
d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]],dat3[[1]],dat4[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]],dat3[[2]],dat4[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]],dat3[[3]],dat4[[3]]))

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d1),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
pclr2 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  ggtitle("Robust")

ggpubr::ggarrange(pcl_logit2,pcl_pois2,pcl_bin2,pcl_lin2,pclr2,ncol=2,nrow = 3)

####################################################
# another example allowing multiple clusters in the same period
# the example with six time periods and seven clusters with 10 people per clusters
nT <- 6
nJ <- nT+1
M <- 20

df <- expand.grid(t=1:nT,J = 1:(nJ*2))
df <- df[rep(1:nrow(df),each=M),]

X <- x_mat(nT,nJ,nper=2)
X <- X[rep(1:nrow(X),each=M),]
beta <- c(rep(0,nT),rep(0,nJ*2),0.5)

data <- as.data.frame(X)
colnames(data) <- c(paste0("t",1:(nT)),paste0("cl",c(1:(nJ*2))),"int")

Sout <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(0.2)),
                     list("indicator",c(0.05))
                   ))
Sout2 <- gen_re_mat(df,
                   dims = list(1,2),
                   funs = list(
                     list("exp_power",c(1)),
                     list("indicator",c(0.05))
                   ))

dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(2:(nJ*2)),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ*2-1),1)),
                    var_par = 1)

dat2 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(2:(nJ*2)),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ*2-1),1)),
                    var_par = 1)

ss <- 100
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)


d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]))


d2 <- grad_robust(idx_in,
                  C=list(dat2[[1]]),
                  X_list=list(dat2[[2]]),
                  sig_list = list(dat2[[3]]),
                  tol=1e-15)

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d1),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:(nJ*2),each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ*2))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)

pcl1 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:(nJ*2))+
  ggtitle(TeX(r'(Up to 2 clusters per sequence, $\rho =0.05$, $R=0.2$)'))

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d2),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:(nJ*2),each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ*2))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)

pcl2 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:(nJ*2))+
  ggtitle(TeX(r'(Up to 2 clusters per sequence, $\rho =0.05$, $R=1$)'))

d3 <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]]))

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d3),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:(nJ*2),each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ*2))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)

pcl3 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:(nJ*2))+
  ggtitle("Robust")

ggpubr::ggarrange(pcl1,pcl2,pcl3)


#################################################################
# EXAMPLES 3: MEAN VERSUS TIME VARYING
#

# This example does not work very well for time period specific interactions-
# i think not enough degrees of freedom





ss <- 150
d <- sample(c(rep(1,ss),rep(0,nrow(df)-ss)),nrow(df))
idx_in <- which(d==1)
#check importance of starting position - some starting positions are not full rank!
# even distribution in each cluster-period
d <- rep(c(rep(1,round(ss*M/nrow(df),0)),rep(0,M-round(ss*M/nrow(df),0))),nT*nJ)


data$dur <- 0
#data$t <- 0
for(t in 1:nT){
  #data[data[,paste0("t",t)]==1,'t'] <- t-1
  for(i in 1:nJ){
    
    if(i == 4){
      data[rowSums(data[,grepl("cl",colnames(data))])==0&data[,paste0("t",t)]==1,'dur'] <- 
        data[rowSums(data[,grepl("cl",colnames(data))])==0&data[,paste0("t",t)]==1,'int']*
        (t - (7-i))
    } else {
      data[data[,paste0("cl",i)]==1&data[,paste0("t",t)]==1,'dur'] <- data[data[,paste0("cl",i)]==1&data[,paste0("t",t)]==1,'int']*
        (t - (7-i))
    }
    
  }
}

for(i in 1:nT){
  data <- cbind(data,I(data[,'dur']==(i-1)&data[,'int']==1)*1)
  colnames(data)[ncol(data)] <- paste0("dur",i)
}
#data$intc <- 1
data$dur <- data$dur/nT

Sout <- gen_re_mat(df,
                    dims = list(1,2),
                    funs = list(
                      list("exp_power",c(0.2)),
                      list("indicator",c(0.05))
                    ))

dat1 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(2:nJ),collapse=", "),", int)")),
                    family=gaussian(link="identity"),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout2[[1]],
                    D=Sout2[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)

dat2 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(2:nJ),collapse=", "),", dur)")),
                    family=gaussian(),
                    data=data,
                    theta=list(
                      list(beta)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nT+nJ-1),1)),
                    var_par = 1)

beta2 <- c(rep(0,nT),rep(0,nJ-1),rep(0.5,nT))

dat3 <- create_data(formula(paste0("~Linear(",paste0("t",1:nT,collapse=", "),",",paste0("cl",c(1:nJ),collapse=", "),", ",
                                   paste0("dur",1:(nT),collapse=", "),")")),
                    family=gaussian(),
                    data=data,
                    theta=list(
                      list(beta2)
                    ),
                    Z=Sout[[1]],
                    D=Sout[[2]],
                    C=matrix(c(rep(0,nJ),rep(1,nT))),
                    var_par = 1)

d1 <- grad_robust(idx_in,
                  C=list(dat1[[1]]),
                  X_list=list(dat1[[2]]),
                  sig_list = list(dat1[[3]]))


d2 <- grad_robust(idx_in,
                  C=list(dat2[[1]]),
                  X_list=list(dat2[[2]]),
                  sig_list = list(dat2[[3]]),
                  tol=1e-15)

X_id <- 1:nrow(dat2[[2]])
# make suggested removals
nc <- 2
dat2[[1]] <- matrix(dat2[[1]][-nc,],ncol=1)
dat2[[3]] <- dat2[[3]][dat2[[2]][,nc]==0,dat2[[2]][,nc]==0]
X_id <- X_id[dat2[[2]][,nc]==0]
dat2[[2]] <- dat2[[2]][dat2[[2]][,nc]==0,-nc]
d <- sample(c(rep(1,ss),rep(0,nrow(dat2[[2]])-ss)),nrow(dat2[[2]]))
idx_in <- which(d==1)
####################
d2 <- X_id[d2]




d3 <- grad_robust(idx_in,
                  C=list(dat3[[1]]),
                  X_list=list(dat3[[2]]),
                  sig_list = list(dat3[[3]]))

X_id <- 1:nrow(dat3[[2]])
# make suggested removals
nc <- 16
dat3[[1]] <- matrix(dat3[[1]][-nc,],ncol=1)
dat3[[3]] <- dat3[[3]][dat3[[2]][,nc]==0,dat3[[2]][,nc]==0]
X_id <- X_id[dat3[[2]][,nc]==0]
dat3[[2]] <- dat3[[2]][dat3[[2]][,nc]==0,-nc]
d <- sample(c(rep(1,ss),rep(0,nrow(dat3[[2]])-ss)),nrow(dat3[[2]]))
idx_in <- which(d==1)
#########################
d3 <- X_id[d3]

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d2),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)

ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  ggtitle("Time invariant")

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d2),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
pm2 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  ggtitle("Time varying")

d3 <- grad_robust(idx_in,
                  C=list(dat1[[1]],dat2[[1]]),
                  X_list=list(dat1[[2]],dat2[[2]]),
                  sig_list = list(dat1[[3]],dat2[[3]]),
                  tol=1e-10,
                  trace=T)

## plot results for d
Xq <- data
Xq <- cbind(Xq,incl=0)
Xq[sort(d3),ncol(Xq)] <- 1
Xq <- cbind(Xq,cl=rep(1:nJ,each=nT*M))
Xq <- cbind(Xq,t=rep(rep(1:nT,each=M),nJ))
Xqa <- aggregate(Xq[,c('incl','int')],list(Xq[,'cl'],Xq[,'t']),sum)
pm3 <- ggplot(data=Xqa,aes(x=Group.2,y=Group.1,fill=incl))+
  geom_tile()+
  geom_label(aes(label=int/M))+
  scale_fill_viridis_c(name="N per block")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  labs(x="T",y="Cluster")+
  scale_x_continuous(breaks = 1:nT)+
  scale_y_continuous(breaks = 1:nJ)+
  ggtitle("Robust")

ggpubr::ggarrange(pm1,pm2,pm3,ncol=2,nrow = 2)




