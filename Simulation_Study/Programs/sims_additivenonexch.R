# Load libraries
library(nimble)
library(tidyverse)
library(posterior)
library(haven)
library(hdbayes)
library(cmdstanr)

## Directory to store files
save.dir <- '/work/users/c/l/clairez1/Paper1sims3/ALL_model'

# Load in grid of simulation parameters
grid_add <- readRDS("/work/users/c/l/clairez1/Paper1sims3/grid_additive.rds")
grid <- subset(grid_add, n0 == 20 | n0 == 50)

## get task ID
id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if ( is.na(id ) )
  id <- 1

## get this job's sim parameters
grid.id      <- grid[id, ]
seed.id      <- as.integer(grid.id$seed)
hist.seed.id <- as.integer(grid.id$hist.seed)
n.id         <- as.integer(grid.id$n)
n0.id        <- as.integer(grid.id$n0)
q.id         <- as.numeric(grid.id$q)
prob.unexch.id <- as.numeric(grid.id$prob.unexch)

## Obtain file name based on id
filename <- file.path(save.dir, paste0('id_', id, '_', 'n_', n.id, '_n0_', n0.id, 
                                       '_q_', q.id, '_probunexch_', prob.unexch.id, '.rds'))

set.seed(seed.id)

# ------------------------------------------------------------------------------
# Load in data

current <- readRDS("/work/users/c/l/clairez1/Paper1sims3/simdata_curr.rds")
historical <- readRDS("/work/users/c/l/clairez1/Paper1sims3/simdata_hist_add.rds")

curr_i <- current[1:n.id,]
hist <- subset(historical, n0 == n0.id & q == q.id & prob.unexch == prob.unexch.id)

df_i <- as.data.frame(cbind(curr_i$y, curr_i$trtgroup, curr_i$BASE))
df0  <- as.data.frame(cbind(hist$y, hist$trtgroup, hist$BASE))

y  <- curr_i$y
X  <- cbind(1, curr_i$trtgroup, curr_i$BASE)
y0 <- hist$y
X0 <- cbind(1, hist$trtgroup, hist$BASE)
p  <- ncol(X)

# Hyperparameters
beta_init_prec = diag(0.1, p)
beta_init_mean = rep(0, p)
tau_init_shape = 0.1
tau_init_rate  = 0.1
xi_shape1      = 2
xi_shape2      = 2
a0_shape1      = 2
a0_shape2      = 2
a0             = 0.5

# Initial Values
xi_init        = 0.5
c0_init_prob   = 0.5
tau_init       = 1
beta_init      = rep(0, p)
a0_init        = 0.5

# ------------------------------------------------------------------------------
# DIPP and NDIPP models
# ------------------------------------------------------------------------------

dipp <- nimbleCode({
  ## Prior on exchangeability probability
  xi ~ dbeta(xi_shape1, xi_shape2)
  ## Conditional prior on latent binary variable | exchangeability probability
  for ( i in 1:n0 ) {
    c0[i] ~ dbern(xi)
    c0X0[i, 1:p] <- c0[i] * X0[i, 1:p]
    c0y0[i]      <- c0[i] * y0[i]
  }
  n0exch <- sum(c0[1:n0])
  ## Compute prior mean / precision of beta | c0, tau
  crossprod_c0X0[1:p,1:p]  <- t(c0X0[1:n0, 1:p]) %*% c0X0[1:n0, 1:p]
  crossprod_c0X0_c0y0[1:p] <- t(c0X0[1:n0, 1:p]) %*% c0y0[1:n0]
  omegatilde0[1:p,1:p]     <- a0*crossprod_c0X0[1:p,1:p] + beta_init_prec[1:p,1:p]
  beta_init_prec_mean[1:p] <- beta_init_prec[1:p,1:p] %*% beta_init_mean[1:p]
  beta_mean[1:p] <- solve( omegatilde0[1:p,1:p],  
                           a0*crossprod_c0X0_c0y0[1:p] + beta_init_prec_mean[1:p] )
  beta_prec[1:p,1:p] <- tau * omegatilde0[1:p,1:p]
  
  # ## Compute prior shape and rate parameters on precision
  tau_shape <- 0.5 * (tau_init_shape + a0 * n0exch)
  tau_rate  <- 0.5 * (tau_init_rate
                      + a0 * sum( (c0y0[1:n0])^2 )
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

ndipp <- nimbleCode({
  a0 ~ dbeta(a0_shape1, a0_shape2)
  
  ## Prior on exchangeability probability
  xi ~ dbeta(xi_shape1, xi_shape2)
  ## Conditional prior on latent binary variable | exchangeability probability
  for ( i in 1:n0 ) {
    c0[i] ~ dbern(xi)
    c0X0[i, 1:p] <- c0[i] * X0[i, 1:p]
    c0y0[i]      <- c0[i] * y0[i]
  }
  n0exch <- sum(c0[1:n0])
  ## Compute prior mean / precision of beta | c0, tau
  crossprod_c0X0[1:p,1:p]  <- t(c0X0[1:n0, 1:p]) %*% c0X0[1:n0, 1:p]
  crossprod_c0X0_c0y0[1:p] <- t(c0X0[1:n0, 1:p]) %*% c0y0[1:n0]
  omegatilde0[1:p,1:p]     <- a0*crossprod_c0X0[1:p,1:p] + beta_init_prec[1:p,1:p]
  beta_init_prec_mean[1:p] <- beta_init_prec[1:p,1:p] %*% beta_init_mean[1:p]
  beta_mean[1:p] <- solve( omegatilde0[1:p,1:p],  
                           a0*crossprod_c0X0_c0y0[1:p] + beta_init_prec_mean[1:p] )
  beta_prec[1:p,1:p] <- tau * omegatilde0[1:p,1:p]
  
  # ## Compute prior shape and rate parameters on precision
  tau_shape <- 0.5 * (tau_init_shape + a0 * n0exch)
  tau_rate  <- 0.5 * (tau_init_rate
                      + a0 * sum( (c0y0[1:n0])^2 )
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

## -----------------------------------------------------------------------------
## Custom samplers for DIPP and NDIPP
## -----------------------------------------------------------------------------

sampler_dipp_conjugate <- nimbleFunction(
  name = 'sampler_dipp_conjugate'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y0             <- model$y0
    X0             <- model$X0
    a0             <- model$a0
    y0y0           <- (t(y0) %*% y0)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X0)
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    y              <- model$y
    X              <- model$X
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    
    c0             <- model[['c0']]
    C0             <- diag(c0)
    c0X0           <- C0 %*% X0
    c0y0           <- (c0 * y0)
    X0C0X0         <- (t(c0X0) %*% c0X0)
    X0C0y0         <- (t(c0X0) %*% c0y0)[, 1]
    y0c0y0         <- (t(c0y0) %*% c0y0)[1,1]
    
    ## Compute posterior parameters
    prec_n      <- XtX + prec0 + a0*X0C0X0
    Uprec_n     <- chol(prec_n)
    Uprec_n_inv <- backsolve(Uprec_n, eye_p)    ## efficient inverse of upper triangular matrix
    cov_n       <- Uprec_n_inv %*% t(Uprec_n_inv)
    mu_n        <- (cov_n %*% ( Xty + prec0_mu0 + a0*X0C0y0 ))[, 1]
    shape_n     <- 0.5 * (n + a0*sum(c0) + shape0)
    rate_n      <- 0.5 * (rate0 + a0*y0c0y0 + yty + mu0t_prec0_mu0 - (t(mu_n) %*% prec_n %*% mu_n)[1,1] )
    
    ## Sample tau
    tau_new <- rgamma(1, shape_n, rate=rate_n)
    
    ## Sample beta | tau
    beta_new <- rmnorm_chol(1, mean = mu_n, cholesky = sqrt(tau_new) * Uprec_n, prec_param = TRUE)
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['beta']] <<- beta_new
    model[['tau']] <<- tau_new
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

sampler_ndipp_conjugate <- nimbleFunction(
  name = 'sampler_ndipp_conjugate'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y0             <- model$y0
    X0             <- model$X0
    y0y0           <- (t(y0) %*% y0)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X0)
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    y              <- model$y
    X              <- model$X
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    a0             <- model[['a0']]
    c0             <- model[['c0']]
    C0             <- diag(c0)
    c0X0           <- C0 %*% X0
    c0y0           <- (c0 * y0)
    X0C0X0         <- (t(c0X0) %*% c0X0)
    X0C0y0         <- (t(c0X0) %*% c0y0)[, 1]
    y0c0y0         <- (t(c0y0) %*% c0y0)[1,1]
    
    ## Compute posterior parameters
    prec_n      <- XtX + prec0 + a0*X0C0X0
    Uprec_n     <- chol(prec_n)
    Uprec_n_inv <- backsolve(Uprec_n, eye_p)    ## efficient inverse of upper triangular matrix
    cov_n       <- Uprec_n_inv %*% t(Uprec_n_inv)
    mu_n        <- (cov_n %*% ( Xty + prec0_mu0 + a0*X0C0y0 ))[, 1]
    shape_n     <- 0.5 * (n + a0*sum(c0) + shape0)
    rate_n      <- 0.5 * (rate0 + a0*y0c0y0 + yty + mu0t_prec0_mu0 - (t(mu_n) %*% prec_n %*% mu_n)[1,1] )
    
    ## Sample tau
    tau_new <- rgamma(1, shape_n, rate=rate_n)
    
    ## Sample beta | tau
    beta_new <- rmnorm_chol(1, mean = mu_n, cholesky = sqrt(tau_new) * Uprec_n, prec_param = TRUE)
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['beta']] <<- beta_new
    model[['tau']] <<- tau_new
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

## -----------------------------------------------------------------------------
## Compile DIPP and NDIPP
## -----------------------------------------------------------------------------

dipp.const <- list(n  = n.id, n0 = n0.id, p  = p)

dipp.data <- list(y = y,  X = X,  y0 = y0,  X0 = X0,
  beta_init_mean = beta_init_mean,
  beta_init_prec = beta_init_prec,
  tau_init_shape = tau_init_shape,
  tau_init_rate  = tau_init_rate,
  xi_shape1 = xi_shape1,
  xi_shape2 = xi_shape2,
  a0 = a0,
  eye = diag(n.id)
)

dipp.init <- list(  c0 = rbinom(n0.id, 1, c0_init_prob),  xi = xi_init)

dipp.model    <- nimbleModel(dipp, constants = dipp.const, data = dipp.data, inits = dipp.init)
dipp.cmodel   <- compileNimble(dipp.model)
dipp.mcmcConf <- configureMCMC(dipp.model, monitors = c('beta', 'xi', 'tau', 'c0', 'xi'))
dipp.mcmcConf$removeSampler(c('beta', 'tau'))
dipp.mcmcConf$addSampler(c('beta', 'tau'), 'sampler_dipp_conjugate')
dipp.mcmc  <- buildMCMC(dipp.mcmcConf)
dipp.cmcmc <- compileNimble(dipp.mcmc, project = dipp.cmodel)


ndipp.const <- list(n  = n.id, n0 = n0.id, p  = p)

ndipp.data <- list(y = y,  X = X,  y0 = y0,  X0 = X0,
  beta_init_mean = beta_init_mean,
  beta_init_prec = beta_init_prec,
  tau_init_shape = tau_init_shape,
  tau_init_rate  = tau_init_rate,
  xi_shape1 = xi_shape1,
  xi_shape2 = xi_shape2,
  a0_shape1 = a0_shape1,
  a0_shape2 = a0_shape2,
  eye = diag(n.id)
)

ndipp.init <- list(c0 = rbinom(n0.id, 1, c0_init_prob),  xi = xi_init,  a0 = a0_init)

ndipp.model    <- nimbleModel(ndipp, constants = ndipp.const, data = ndipp.data, inits = ndipp.init)
ndipp.cmodel   <- compileNimble(ndipp.model)
ndipp.mcmcConf <- configureMCMC(ndipp.model, monitors = c('beta', 'xi', 'tau', 'c0', 'n0exch', 'a0'))
ndipp.mcmcConf$removeSampler(c('beta', 'tau'))
ndipp.mcmcConf$addSampler(c('beta', 'tau'), 'sampler_ndipp_conjugate')
ndipp.mcmc  <- buildMCMC(ndipp.mcmcConf)
ndipp.cmcmc <- compileNimble(ndipp.mcmc, project = ndipp.cmodel)

# ------------------------------------------------------------------------------
# POWER PRIOR MODELS
# ------------------------------------------------------------------------------

pp <- nimbleCode({
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


npp <- nimbleCode({
  a0 ~ dbeta(a0_shape1, a0_shape2)
  
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

## -----------------------------------------------------------------------------
## Custom samplers for power prior models
## -----------------------------------------------------------------------------

sampler_pp_conjugate <- nimbleFunction(
  name = 'sampler_pp_conjugate'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    
    y0             <- model$y0
    X0             <- model$X0
    a0             <- model$a0
    y0y0           <- (t(y0) %*% y0)[1,1]
    X0y0           <- (t(X0) %*% y0)[, 1]
    X0X0           <- t(X0) %*% X0
    n0             <- length(y0)
    p              <- ncol(X0)
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    
    y              <- model$y
    X              <- model$X
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    n              <- length(y)
    
    ## Compute posterior parameters
    prec_n      <- XtX + prec0 + a0 * X0X0
    Uprec_n     <- chol(prec_n)
    Uprec_n_inv <- backsolve(Uprec_n, eye_p)    ## efficient inverse of upper triangular matrix
    cov_n       <- Uprec_n_inv %*% t(Uprec_n_inv)
    mu_n        <- (cov_n %*% ( Xty + prec0_mu0 + a0*X0y0 ))[, 1]
    shape_n     <- 0.5 * (n + a0*n0 + shape0)
    rate_n      <- 0.5 * (rate0 + a0*y0y0 + yty + mu0t_prec0_mu0 - (t(mu_n) %*% prec_n %*% mu_n)[1,1] )
    
    ## Sample tau
    tau_new <- rgamma(1, shape_n, rate=rate_n)
    
    ## Sample beta | tau
    beta_new <- rmnorm_chol(1, mean = mu_n, cholesky = sqrt(tau_new) * Uprec_n, prec_param = TRUE)
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['beta']] <<- beta_new
    model[['tau']] <<- tau_new
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

sampler_npp_conjugate <- nimbleFunction(
  name = 'sampler_npp_conjugate'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y0             <- model$y0
    X0             <- model$X0
    X0X0           <- t(X0) %*% X0
    X0y0           <- (t(X0) %*% y0)[, 1]
    y0y0           <- (t(y0) %*% y0)[1,1]
    n0             <- length(y0)
    p              <- ncol(X0)
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    a0 <- model[['a0']]
    y              <- model$y
    X              <- model$X
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    n              <- length(y)
    ## Compute posterior parameters
    prec_n      <- XtX + prec0 + a0*X0X0
    Uprec_n     <- chol(prec_n)
    Uprec_n_inv <- backsolve(Uprec_n, eye_p)    ## efficient inverse of upper triangular matrix
    cov_n       <- Uprec_n_inv %*% t(Uprec_n_inv)
    mu_n        <- (cov_n %*% ( Xty + prec0_mu0 + a0*X0y0 ))[, 1]
    shape_n     <- 0.5 * (n + a0*n0 + shape0)
    rate_n      <- 0.5 * (rate0 + a0*y0y0 + yty + mu0t_prec0_mu0 - (t(mu_n) %*% prec_n %*% mu_n)[1,1] )
    
    ## Sample tau
    tau_new <- rgamma(1, shape = shape_n, rate = rate_n)
    
    ## Sample beta | tau
    beta_new <- rmnorm_chol(1, mean = mu_n, cholesky = sqrt(tau_new) * Uprec_n, prec_param = TRUE)
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['beta']] <<- beta_new
    model[['tau']] <<- tau_new
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

## -----------------------------------------------------------------------------
## Compile power prior and NPP samplers
## -----------------------------------------------------------------------------

pp.const <- list(n = n.id, n0 = n0.id, p = p)
pp.data  <- list(y = y, X = X, y0 = y0, X0 = X0
              , beta_init_mean = beta_init_mean
              , beta_init_prec = beta_init_prec
              , tau_init_shape = tau_init_shape
              , tau_init_rate = tau_init_rate
              , eye = diag(n.id)
              , a0 = a0
)
pp.init <-  list(beta = beta_init, tau = tau_init)

pp.model    <- nimbleModel(pp, constants = pp.const, data = pp.data, inits = pp.init)
pp.cmodel   <- compileNimble(pp.model)
pp.mcmcConf <- configureMCMC(pp.model, monitors = c('beta', 'tau'))
pp.mcmcConf$removeSampler(c('beta', 'tau'))
pp.mcmcConf$addSampler(c('beta', 'tau'), 'sampler_pp_conjugate')
pp.mcmc     <- buildMCMC(pp.mcmcConf)
pp.cmcmc    <- compileNimble(pp.mcmc, project = pp.cmodel)


npp.const <- list(n  = n.id,  n0 = n0.id,  p  = p)

npp.data <- list(  y = y,  X = X,
  y0 = y0,  X0 = X0,
  beta_init_mean = beta_init_mean,
  beta_init_prec = beta_init_prec,
  tau_init_shape = tau_init_shape,
  tau_init_rate = tau_init_rate,
  eye = diag(n.id),
  a0_shape1 = a0_shape1,
  a0_shape2 = a0_shape2
)

npp.init <- list(beta = beta_init,  tau  = tau_init,  a0   = a0_init)

npp.model    <- nimbleModel(npp, constants = npp.const, data = npp.data, inits = npp.init)
npp.cmodel   <- compileNimble(npp.model)
npp.mcmcConf <- configureMCMC(npp.cmodel, monitors = c('beta', 'tau', 'a0'))
npp.mcmcConf$removeSampler(c('beta', 'tau'))
npp.mcmcConf$addSampler(c('beta', 'tau'), 'sampler_npp_conjugate')
npp.mcmc  <- buildMCMC(npp.mcmcConf)
npp.cmcmc <- compileNimble(npp.mcmc, project = npp.cmodel)

# ------------------------------------------------------------------------------
# Begin sim code

# Simulation parameters
Nsims  = 500
niter  = 10000
burnin = 2000
thin   = 2

# Create empty data frame to store simulation results
results.all <- data.frame()

start.time <- Sys.time()
for (i in 1:Nsims) {
  
  block  <- (id - 1) %% (10000/Nsims) + 1
  simID  <- (block - 1) * Nsims + i
  end    <- simID * n.id
  start  <- end - n.id + 1
  
  curr_i <- current[start:end,]
  df_i   <- as.data.frame(cbind(curr_i$y, curr_i$trtgroup, curr_i$BASE))
  y      <- curr_i$y
  X      <- cbind(1, curr_i$trtgroup, curr_i$BASE)
  
  # Replace data in nimble model
  dipp.cmodel$y  <- y
  dipp.cmodel$X  <- X
  ndipp.cmodel$y <- y
  ndipp.cmodel$X <- X
  pp.cmodel$y    <- y
  pp.cmodel$X    <- X
  npp.cmodel$y   <- y
  npp.cmodel$X   <- X
  
  # DIPP
  dipp.smpl    <- runMCMC(dipp.cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  dipp.summary <- summarize_draws(dipp.smpl, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
  dipp.results <- data.frame(dipp.summary %>% filter(grepl("beta", variable) | variable %in% c("xi", "tau", "n0exch")) )
  dipp.results$method <- "dipp"
  
  # NDIPP
  ndipp.smpl    <- runMCMC(ndipp.cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  ndipp.summary <- summarize_draws(ndipp.smpl, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
  ndipp.results <- data.frame(ndipp.summary %>% filter(grepl("beta", variable) | variable %in% c("xi", "tau", "n0exch", "a0")) )
  ndipp.results$method <- "ndipp"
  
  # PP
  pp.smpl    <- runMCMC(pp.cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  pp.summary <- summarize_draws(pp.smpl, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
  pp.results <- data.frame(pp.summary %>% filter(grepl("beta", variable) | variable %in% c("tau")) )
  pp.results$method <- "pp"
  
  # NPP
  npp.smpl    <- runMCMC(npp.cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  npp.summary <- summarize_draws(npp.smpl, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
  npp.results <- data.frame(npp.summary %>% filter(grepl("beta", variable) | variable %in% c("tau", "a0")) )
  npp.results$method <- "npp"
  
  
  # LEAP
  sink(tempfile())
  fit <- glm.leap(formula = V1 ~ V2 + V3, family = gaussian(link = "identity"),
                  data.list = list(df_i, df0), K = 2,
                  iter_warmup = burnin, iter_sampling = niter, chains = 1)
  sink()
  
  leap.trt <- data.frame(variable = "trteffect", mean = mean(fit$V2), sd = sd(fit$V2),
                             X2.5. = quantile(fit$V2, 0.025), X97.5. = quantile(fit$V2, 0.975)  )
  leap.prob <- data.frame(variable = "exch.param", mean = mean(fit$`probs[1]`), sd = sd(fit$`probs[1]`),
                             X2.5. = quantile(fit$`probs[1]`, 0.025), X97.5. = quantile(fit$`probs[1]`, 0.975)  )
  leap.results <- rbind(leap.trt, leap.prob)
  leap.results$method <- "leap"
  
  results.all <- rbind(results.all, dipp.results, ndipp.results, pp.results, npp.results, leap.results)
  cat("Iteration:", i, "\n")
}

end.time <- Sys.time()
end.time - start.time

# End sim code
# ------------------------------------------

## SAVE THE RESULTS

lst <- list(
  'simscen' = grid.id
  , 'id'      = id
  , 'simres'  = results.all
)

saveRDS(lst, filename)