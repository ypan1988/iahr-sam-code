library(glmmr)

#####
##### Markov Chain Monte Carlo Maximum Likelihood
#####

df <- nelder(~(cl(10)*t(5)) > ind(10))
df$int <- 0
df[df$cl > 5, 'int'] <- 1
des <- Design$new(
  covariance = list(
    data = df,
    formula = ~ (1|gr(cl)*ar1(t)),
    parameters = c(0.25,0.8)),
  mean.function = list(
    formula = ~ factor(t) + int - 1,
    data=df,
    parameters = c(rep(0,5),0.6),
    family = binomial())
  )
ysim <- des$sim_data()

# fits the models using MCEM 
fit1 <- des$MCML(y = ysim)
fit1


#fits the models using MCNR
fit2 <- des$MCML(y = ysim, method="mcnr")
fit2

#adds a simulated likelihood step after the MCEM algorithm
fit3 <- des$MCML(y = des$sim_data(), options = list(sim_lik_step=TRUE))
fit3


#####
##### Simulation-Based Analysis
#####

df <- nelder(~(cl(10)*t(5)) > ind(10))
df$int <- 0
df[df$cl > 5, 'int'] <- 1

mf1 <- MeanFunction$new(
  formula = ~ factor(t) + int - 1,
  data=df,
  parameters = c(rep(0,5),0.6),
  family = gaussian()
)
cov1 <- Covariance$new(
  data = df,
  formula = ~ (1|gr(cl)) + (1|gr(cl*t)),
  parameters = c(0.25,0.1)
)
des <- Design$new(
  covariance = cov1,
  mean.function = mf1,
  var_par = 1
)
# analysis using MCML mcem algorithm
test1 <- des$analysis(type="sim",
                     iter=100,
                     par=6,
                     parallel = FALSE,
                     verbose = FALSE,
                     method = "mcem",
                     m = 100)
#an analysis using the permutation test option and MCNR
test2 <- des$analysis(type="sim",
                     iter=100,
                     se.method="perm",
                     par=6,
                     parallel = FALSE,
                     verbose = FALSE,
                     options = list(
                       perm_type="unw", perm_iter=100,
                       perm_parallel=FALSE,perm_ci_steps=1000),
                     method = "mcnr",
                     m = 100)
#returning previously saved sim data
test3 <- des$analysis(type="sim_data")

#to test model misspecification we can simulate from a different model
cov2 <- Covariance$new(
  data = df,
  formula = ~ (1|gr(cl)*ar1(t)),
  parameters = c(0.25,0.8)
)
des2 <- Design$new(
  covariance = cov2,
  mean.function = mf1,
  var_par = 1
)
test4 <- des$analysis(type="sim",
                     iter=100,
                     sim_design = des2,
                     par=6,
                     parallel = FALSE,
                     verbose = FALSE,
                     method = "mcem",
                     m = 100)