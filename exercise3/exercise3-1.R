set.seed(1)
library(dplyr)

## load data: diabetes information for women of pima indian heritage (http://archive.ics.uci.edu/ml/datasets/Pima+Indians+Diabetes)
data <- read.csv('pima.csv')
summary(data)

#rename 1st column name
names(data)[1] <- "times_pregnant"

#naming x & y
y <- as.matrix(data$class_variable)
x <- data[, -9] %>%
  mutate(intercept = rep(1, n))
x <- as.matrix(x)
n <- nrow(y)
p <- ncol(x)

#GLM
GLM_est <- glm(y ~ x, family = binomial(link = probit))
GLM_est

########## Writing Gibbs Sampler##########

library(MASS)
library(truncnorm)

## prior hyperparameters befor update
k_init <- diag(rep(0.1, p)) #precision matrix
beta_init <- matrix(0.1, p) #prior guess on beta
a_init <- 1 # prior sample size for the 
b_init <- 1 # prior sum of square error for the error variance

N1 <- sum(y) # number of success
N0 <- n - N1 # number of failures

## Update parameters
n_iter <- 5000
beta <- matrix(NA, n_iter, p) #p = # of column
tau <- matrix(NA, n_iter)
beta[1, ] <- rep(0, p)
tau[1] <- 1
z <- rep(0, n)

## actual gibbs sampler
for (i in 2:n_iter){
  mu_z <- x %*% beta[i-1, ] #update mean of z based on beta ##Learn from Natalia's work about the truncnorm library. Thank you!
  z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b =0)
  z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
  k_new <- k_init + crossprod(x)
  beta_new <- solve(k_new) %*% (crossprod(x, z) + k_init %*% beta_init) # Get the betas for update
  beta[i, ] <- mvrnorm(1, beta_new, solve(tau[i-1] * k_new))
  a_new <- a_init + (n + 1) / 2 # get a new tau
  s <- t(beta_init) %*% k_init %*% beta_init + crossprod(z)
  r <- t(beta_init) %*% (k_init + crossprod(x)) %*% beta_init
  b_new <- as.numeric(b_init + 1 / 2*(s - r))
  tau[i] <- rgamma(1, a_new, b_new)
  }

burnin <- 1000
beta_post <- colMeans(betap[-(1:burnin), ])
y_hat <- ifelse(x %*% beta_post > 0, 1, 0)
accuracy.bayes <- sum(y_hat) / sum(y) * 100

########## Using MCMCpack for estimation #######
library(MCMCpack)

MCMC_logit <- MCMClogit(class_variable ~ times_pregnant + plasma_glucose + dia_blood_pressure + triceps_skinfold + serum_insulin + bmi + pedigree + age, data = data, mcmc=5000)
summary(MCMC_logit)
plot(MCMC_logit)


