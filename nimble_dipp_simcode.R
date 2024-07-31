
# Load libraries
library(nimble)
library(tidyverse)
library(posterior)
library(haven)

## Directory to store files
save.dir <- '/work/users/c/l/clairez1/Paper1sims2/DIPP'

# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper1sims2/grid.rds")

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

# Load in CSL data
setwd("/nas/longleaf/home/clairez1/Paper1 files/CSL Datasets")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")

# Log transform outcome variable
current$LOGCHG <- log(current$AVAL + 1) - log(current$BASE + 1)
historical$LOGCHG <- log(historical$AVAL + 1) - log(historical$BASE + 1)

# Analysis of actual data
y  <- current$LOGCHG
X  <- cbind(1, current$trtgroup, current$BASE)
lm <- summary(lm(y ~ 0 + X)) # regression coefficients
lmbeta <- lm$coefficients[,1]
tau <- 1/(lm$sigma)^2
p  <- ncol(X)

# Parameters
beta <- cbind(lmbeta, (lmbeta * q.id)) # 1st vector = exch values, 2nd vector = unexch values
sd   <- lm$sigma

# Generate current data by sampling covariates with replacement
sampled_curr <- current[sample(nrow(current), n.id, replace = TRUE), ]
X  <- cbind(1, sampled_curr$trtgroup, sampled_curr$BASE) 
y  <- rnorm(n.id, X %*% beta[,1], sd)

# Initialize variables to store sampled data corresponding to the indices
X0_exch   <- NULL
y0_exch   <- NULL
X0_unexch <- NULL
y0_unexch <- NULL

min.dist.exch   <- Inf
min.dist.unexch <- Inf

set.seed(hist.seed.id)
for (i in 1:10000) {
  n.unexch <- round(n0.id * prob.unexch.id)
  n.exch <- n0.id - n.unexch
  
  if (n.exch > 0) {
    sampled.hist.exch <- historical[sample(nrow(historical), n.exch, replace = TRUE), ]
    X0.exch <- cbind(1, sampled.hist.exch$trtgroup, sampled.hist.exch$BASE)
    y0.exch <- rnorm(n.exch, X0.exch %*% beta[, 1], sd)
    lm.exch <- summary(lm(y0.exch ~ 0 + X0.exch))$coefficients[, 1]
    coeff.dist.exch <- (lm.exch[1] - (beta[1, 1] + beta[2, 1]))^2 + (lm.exch[2] - beta[3, 1])^2
  } else { 
    coeff.dist.exch <- Inf
  }
  
  if (n.unexch > 0) {
    sampled.hist.unexch <- historical[sample(nrow(historical), n.unexch, replace = TRUE), ]
    X0.unexch <- cbind(1, sampled.hist.unexch$trtgroup, sampled.hist.unexch$BASE)
    y0.unexch <- rnorm(n.unexch, X0.unexch %*% beta[, 2], sd)
    lm.unexch <- summary(lm(y0.unexch ~ 0 + X0.unexch))$coefficients[, 1]
    coeff.dist.unexch <- (lm.unexch[1] - (beta[1, 2] + beta[2, 2]))^2 + (lm.unexch[2] - beta[3, 2])^2
  } else { 
    coeff.dist.unexch <- Inf
  }
  
  # Update dataset if MLEs are closer
  if (!is.na(coeff.dist.exch) && coeff.dist.exch < min.dist.exch) {
    min.dist.exch <- coeff.dist.exch
    X0_exch <- X0.exch
    y0_exch <- y0.exch
  } 
  if (!is.na(coeff.dist.unexch) && coeff.dist.unexch < min.dist.unexch) {
    min.dist.unexch <- coeff.dist.unexch
    X0_unexch <- X0.unexch
    y0_unexch <- y0.unexch
  } 
  
  X0 <- rbind(X0_exch, X0_unexch)
  y0 <- c(y0_exch, y0_unexch)
}

# Hyperparameters
beta_init_prec = diag(0.1, p)
beta_init_mean = rep(0, p)
tau_init_shape = .1
tau_init_rate = .1
xi_shape1 = 2
xi_shape2 = 2
a0_shape1 = 2
a0_shape2 = 2
a0 = 0.5

# Initial Values
xi_init = 0.5
c0_init_prob = 0.5
tau_init = 1
beta_init = rep(0, p)
a0_init = 0.5

# ------------------------------------------------------------------------------
# NIMBLE CODE FOR MODELS
# ------------------------------------------------------------------------------

dipp.a0fixed <- nimbleCode({
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
## Custom samplers
## -----------------------------------------------------------------------------

sampler_dipp_a0fixed_conjugate <- nimbleFunction(
  name = 'sampler_dipp_a0fixed_conjugate'  ## give name for sampler
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

## -----------------------------------------------------------------------------
## Compile DIPP sampler with fixed a0
## -----------------------------------------------------------------------------

const <- list(n = n.id, n0 = n0.id, p = p)
data  <- list(
  y = y, X = X, y0 = y0, X0 = X0
  , beta_init_mean = beta_init_mean
  , beta_init_prec = beta_init_prec
  , tau_init_shape = tau_init_shape
  , tau_init_rate = tau_init_rate
  , xi_shape1 = xi_shape1
  , xi_shape2 = xi_shape2
  , a0 = a0
  , eye = diag(n.id)
)
init <- list(c0 = rbinom(n0.id, 1, c0_init_prob), xi = xi_init)

model    <- nimbleModel(dipp.a0fixed, constants = const, data = data, inits = init)
cmodel   <- compileNimble(model)
mcmcConf <- configureMCMC(model, monitors = c('beta', 'xi', 'tau', 'c0', 'n0exch'))
mcmcConf$removeSampler(c('beta', 'tau'))
mcmcConf$addSampler(c('beta', 'tau'), 'sampler_dipp_a0fixed_conjugate')
mcmc     <- buildMCMC(mcmcConf)
cmcmc    <- compileNimble(mcmc, project = cmodel)


#-------------------------------------------------------------------------------

# Begin sim code

# Simulation parameters
Nsims  = 2000
niter  = 10000
burnin = 2000
thin   = 2

# Create empty data frame to store simulation results
results.all <- data.frame()

start <- Sys.time()
for (i in 1:Nsims) {
  
  # Generate new current data
  sampled_curr <- current[sample(nrow(current), n.id, replace = TRUE), ]
  X  <- cbind(1, sampled_curr$trtgroup, sampled_curr$BASE) 
  y  <- rnorm(n.id, X %*% beta[,1], sd)
  
  # Replace data in nimble model
  cmodel$y   <- y
  cmodel$X   <- X
  
  smpl.dipp    <- runMCMC(cmcmc, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
  summary.dipp <- summarize_draws(smpl.dipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
  results.dipp <- data.frame(summary.dipp %>% filter(grepl("beta", variable) | variable %in% c("xi", "tau", "n0exch")))
  results.dipp$method <- "dipp"
  
  results <- rbind(results.dipp)
  results.all <- rbind(results.all, results)
  cat("Iteration:", i, "\n")
}

end <- Sys.time()
end - start

# End sim code
# ------------------------------------------

## SAVE THE RESULTS

lst <- list(
  'simscen' = grid.id
  , 'id'      = id
  , 'simres'  = results.all
)

saveRDS(lst, filename)



