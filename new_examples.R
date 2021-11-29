####################################################################################
# TEST EXAMPLES ####################################################################
####################################################################################

n_dim <- 10
xpos <- seq(-1,1,length.out = 2*n_dim+1)
xpos <- xpos[(1:n_dim)*2]
pos <- expand.grid(x=xpos,y=xpos)

data <- data.frame(x=sqrt(pos$x^2 + pos$y^2),y=rnorm(100),z=runif(100))

Sout <- gen_re_mat(pos,
                   dims = list(c(1,2)),
                   funs = list(
                     list("exponential",c(1,0.5))
                   ))

dat <- create_data(~Exp(x)+Linear(y,z),
            family=binomial(link="logit"),
            data=data,
            theta=list(
              list(1),
              list(0.5,-0.5),
              list(1,1)
            ),
            Z=Sout[[1]],
            D=Sout[[2]],
            C=matrix(c(0,1,1,0,0)))


#starting vector
d <- sample(c(rep(1,20),rep(0,100)),420)
idx_in <- which(d==1)
sig <- dat[[1]]
u <- dat[[2]]
A <- solve(sig[idx_in,idx_in])

d <- grad(idx_in,A,sig,u,tol=1e-10,T)
d2 <- Grad(idx_in-1,A,sig,u,tol=1e-10,T)
