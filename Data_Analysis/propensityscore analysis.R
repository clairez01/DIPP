# ------------------------------------------------------------------------------
# Propensity scoring approach
# ------------------------------------------------------------------------------

library(nimble)
library(tidyverse)
library(posterior)
library(haven)
library(ggplot2)
library(patchwork)
library(dplyr)
library(MatchIt)
library(lmtest)
library(sandwich)

# Load in CSL data
setwd("C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/data set")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")

# Log transform outcome variable
current$LOGCHG <- log(current$AVAL + 1) - log(current$BASE + 1)
historical$LOGCHG <- log(historical$AVAL + 1) - log(historical$BASE + 1)

# Remove unneeded variables
current1 <- subset(current, select = -c(CHG, AVAL, PCHG, DOSE, TRT01AN))
historical1 <- subset(historical, select = -c(CHG, AVAL, PCHG, DOSE))
historical1 <- historical1 %>% filter(!is.na(ageN))

# Combine data sets
current1$current <- 1
historical1$current <- 0
df <- rbind(current1, historical1)

# Match historical subjects to current subjects
match_obj <- matchit(current ~ BASE + ageN + sexN + BMI + haetypeN + trtgroup,
                     data = df, method = "nearest", distance ="glm", link = "logit",
                     ratio = 1,
                     caliper = .1,
                     replace = FALSE)
summary(match_obj)

# Plot balance between matched subjects
plot(match_obj, type = "jitter", interactive = FALSE)
plot(summary(match_obj), abs = FALSE)

# Obtain matched historical subjects
matched_hist <- match.data(match_obj) %>% filter(current == 0)

# Assign data to variables
y0 <- matched_hist$LOGCHG
y  <- current1$LOGCHG
X0 <- cbind(1, matched_hist$trtgroup, matched_hist$BASE)
X  <- cbind(1, current1$trtgroup, current1$BASE)


# --------------------------------------------------
# data analysis step using power prior

p <- ncol(X)
n <- nrow(X)
n0 <- nrow(X0)

# Hyperparameters
beta_init_prec = diag(0.1, p)
beta_init_mean = rep(0, p)
tau_init_shape = .1
tau_init_rate = .1
xi_shape1 = 2
xi_shape2 = 2
a0_shape1 = 2
a0_shape2 = 2
a0 = 1

# Initial Values
tau_init = 1
beta_init = rep(0, p)


pp.a0fixed <- nimbleCode({
  ## Compute prior mean / precision of beta | tau
  crossprod_X0[1:p,1:p]  <- t(X0[1:n0, 1:p]) %*% X0[1:n0, 1:p]
  crossprod_X0y0[1:p] <- t(X0[1:n0, 1:p]) %*% y0[1:n0]
  omegatilde0[1:p,1:p]     <- a0*crossprod_X0[1:p,1:p] + beta_init_prec[1:p,1:p]
  beta_init_prec_mean[1:p] <- beta_init_prec[1:p,1:p] %*% beta_init_mean[1:p]
  beta_mean[1:p] <- solve( omegatilde0[1:p,1:p],  
                           a0*crossprod_X0y0[1:p] + beta_init_prec_mean[1:p] )
  beta_prec[1:p,1:p] <- tau * omegatilde0[1:p,1:p]
  
  # ## Compute prior shape and rate parameters on precision
  tau_shape <- 0.5 * (tau_init_shape + a0 * n0)
  tau_rate  <- 0.5 * (tau_init_rate
                      + a0 * sum( (y0[1:n0])^2 )
                      + inprod(beta_init_mean[1:p], (beta_init_prec[1:p,1:p] %*% beta_init_mean[1:p]))
                      - inprod(beta_mean[1:p], omegatilde0[1:p,1:p] %*% beta_mean[1:p])
  )
  ## Joint prior for regression coefficients and precision | c0
  beta[1:p] ~ dmnorm(beta_mean[1:p], prec = beta_prec[1:p,1:p])
  tau       ~ dgamma(tau_shape, rate = tau_rate)
  
  ## Current data likelihood
  y_mean[1:n]         <- X[1:n, 1:p] %*% beta[1:p]
  y_prec_mtx[1:n,1:n] <- tau * eye[1:n,1:n] 
  y[1:n] ~ dmnorm(y_mean[1:n], prec = y_prec_mtx[1:n,1:n])
})


# compile nimble code
const4 <- list(n = n, n0 = n0, p = p)
data4  <- list(y = y, X = X, y0 = y0, X0 = X0
               , beta_init_mean = beta_init_mean
               , beta_init_prec = beta_init_prec
               , tau_init_shape = tau_init_shape
               , tau_init_rate = tau_init_rate
               , eye = diag(n)
               , a0 = a0
)
init4  <- list(beta = beta_init, tau = tau_init)
model4    <- nimbleModel(pp.a0fixed, constants = const4, data = data4, inits = init4)
cmodel4   <- compileNimble(model4)
mcmcConf4 <- configureMCMC(cmodel4, monitors = c('beta', 'tau'))
mcmc4     <- buildMCMC(mcmcConf4)
cmcmc.pp.a0fixed <- compileNimble(mcmc4, project = cmodel4)


niter  = 200000
burnin = 2000
thin   = 2

set.seed(8090)
smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
print(summary.pp.a0fixed)

# Compute DIC
logL <- function(y, X, beta1, beta2, beta3, tau) {
  beta_est <- c(beta1, beta2, beta3)
  n <- length(y)
  log_likelihood <- (n/2) * log(tau/(2 * pi)) - sum((tau/2) * (y - X %*% beta_est)^2)
  return(log_likelihood)
}

beta1_est <- mean(smpl.pp.a0fixed[,1])
beta2_est <- mean(smpl.pp.a0fixed[,2])
beta3_est <- mean(smpl.pp.a0fixed[,3])
tau_est   <- mean(smpl.pp.a0fixed[,4])

log_likelihoods <- sapply(1:niter, function(i) {
  logL(y, X, smpl.pp.a0fixed[i, 1], smpl.pp.a0fixed[i, 2], smpl.pp.a0fixed[i, 3], smpl.pp.a0fixed[i, 4])
})

pDIC <- 2*(logL(y, X, beta1_est, beta2_est, beta3_est, tau_est) - mean(log_likelihoods))
DIC  <- -2*logL(y, X, beta1_est, beta2_est, beta3_est, tau_est) + 2*pDIC
