# ------------------------------------------------------------------------------
# Title: Data Analysis of COMPACT Studies
# Date : April 2024
# ------------------------------------------------------------------------------

# load libraries
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
library(posterior)
library(bayesplot)

# Load in CSL data
setwd("C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/data set")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")
#historical <- historical[is.na(historical$ageN) == FALSE,]

# Log transform outcome variable
current$LOGCHG <- log(current$AVAL + 1) - log(current$BASE + 1)
historical$LOGCHG <- log(historical$AVAL + 1) - log(historical$BASE + 1)

# Assign data to y, X, y0, and x0
y  <- current$LOGCHG
X  <- cbind(1, current$trtgroup, current$BASE)
#X  <- cbind(1, current$trtgroup, current$BASE, current$ageN, current$sexN, current$BMI)
y0 <- historical$LOGCHG
X0 <- cbind(1, historical$trtgroup, historical$BASE)
#X0 <- cbind(1, historical$trtgroup, historical$BASE, historical$ageN, historical$sexN, historical$BMI)

# Analysis of actual data
lm <- summary(lm(y ~ 0 + X)) # regression coefficients
lmbeta <- lm$coefficients[,1]
tau <- 1/(lm$sigma)^2
p  <- ncol(X)
n  <- length(y)
n0 <- length(y0)

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
xi = 0.5

# ------------------------------------------------------------------------------
# Create Summary Plots

# Change theme
theme_set(theme_bw())

curr <- as.data.frame(cbind(current$LOGCHG, current$trtgroup, current$BASE, 1))
hist <- as.data.frame(cbind(historical$LOGCHG, historical$trtgroup, historical$BASE, 0))
ps_hist <- as.data.frame(cbind(matched_hist$LOGCHG, matched_hist$trtgroup, matched_hist$BASE, 2))

all.df <- rbind(curr, hist, ps_hist)
all.df$V4 <- factor(all.df$V4)
all.active <- all.df[all.df$V2 == 1,]
all.placebo <- all.df[all.df$V2 == 0,]

# Density plots of outcomes for historical and current data
dens.active <- ggplot(all.active, aes(x = V1, group = V4, fill = V4)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = "Data Source", 
                    values = c("blue", "red", "lightblue"), 
                    labels = c("Historical (COMPACT Phase 2)", "Current (COMPACT Phase 3)", "Historical (after PS matching)")) +
  labs(x = "Outcome", y = "Density", title = "Active Arms (HAEGARDA)") +
  xlim(-2, 3) + ylim(0, 1.5)

dens.placebo <- ggplot(all.placebo, aes(x = V1, group = V4, fill = V4)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(name = "Data Source", 
                    values = c("red"), 
                    labels = c("Current (Phase 3)")) +
  labs(x = "Outcome", y = "Density", title = "Placebo Arm") +
  xlim(-2, 3) + ylim(0, 1.5)

dens.active + dens.placebo

# ------------------------------------------------------------------------------
# Custom samplers

sampler_dipp_betatau_a0ran_conjugate <- nimbleFunction(
  name = 'lm_conjugate_sampler'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y              <- model$y
    X              <- model$X
    y0             <- model$y0
    X0             <- model$X0
    #a0             <- model$a0
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    y0y0           <- (t(y0) %*% y0)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X)
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

sampler_dipp_betatau_a0fix_conjugate <- nimbleFunction(
  name = 'lm_conjugate_sampler'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y              <- model$y
    X              <- model$X
    y0             <- model$y0
    X0             <- model$X0
    a0             <- model$a0
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    y0y0           <- (t(y0) %*% y0)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X)
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

sampler_c0_a0ran_conjugate <- nimbleFunction(
  name = 'c0_conjugate_sampler'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y              <- model$y
    X              <- model$X
    y0             <- model$y0
    X0             <- model$X0
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    xi_shape1      <- model$xi_shape1
    xi_shape2      <- model$xi_shape2
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X0)
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    
    a0     <- model[['a0']]
    c0     <- model[['c0']]
    beta   <- model[['beta']]
    tau    <- model[['tau']]
    c0_new <- c0
    
    ## Sample xi
    xi_shape1_new <- sum(c0) + xi_shape1
    xi_shape2_new <- n0 - sum(c0) + xi_shape2
    xi_new <- rbeta(1, xi_shape1_new, xi_shape2_new)
    
    for (i in 1:n0) {
      L0i <- (tau / (2*3.1415926))^(1/2) * exp((-tau/2) * (y0[i] - (X0[i, ] %*% beta))^2)
      
      # Compute component of Z(c0 = 1)
      c01 <- c0
      c01[i] <- 1
      C01 <- diag(c01)
      
      c0y0_1   <- c01 * y0
      c0X0_1   <- C01 %*% X0
      X0C0X0_1 <- t(c0X0_1) %*% c0X0_1
      X0c0y0_1 <- (t(c0X0_1) %*% c0y0_1)[, 1]
      y0c0y0_1 <- (t(c0y0_1) %*% c0y0_1)[1,1]
      
      prec_n_1     <- a0 * X0C0X0_1 + prec0
      Uprec_n_1    <- chol(prec_n_1)
      Uprec_n_inv1 <- backsolve(Uprec_n_1, eye_p)
      cov_n_1      <- Uprec_n_inv1 %*% t(Uprec_n_inv1)
      mu_n_1       <- cov_n_1 %*% (a0 * X0c0y0_1 + prec0_mu0)
      
      term_1      <- 0.5*(rate0 + a0*y0c0y0_1 + mu0t_prec0_mu0 - (t(mu_n_1) %*% prec_n_1 %*% mu_n_1)[1,1])
      
      # Compute component of Z(c0 = 0)
      c00 <- c0
      c00[i] <- 0
      C00 <- diag(c00)
      
      c0y0_0   <- c00 * y0
      c0X0_0   <- C00 %*% X0
      X0C0X0_0 <- t(c0X0_0) %*% c0X0_0
      X0c0y0_0 <- (t(c0X0_0) %*% c0y0_0)[, 1]
      y0c0y0_0 <- (t(c0y0_0) %*% c0y0_0)[1,1]
      
      prec_n_0     <- a0 * X0C0X0_0 + prec0
      Uprec_n_0    <- chol(prec_n_0)
      Uprec_n_inv0 <- backsolve(Uprec_n_0, eye_p)
      cov_n_0      <- Uprec_n_inv0 %*% t(Uprec_n_inv0)
      mu_n_0       <- cov_n_0 %*% (a0 * X0c0y0_0 + prec0_mu0)
      
      term_0      <- 0.5*(rate0 + a0*y0c0y0_0 + mu0t_prec0_mu0 - (t(mu_n_0) %*% prec_n_0 %*% mu_n_0)[1,1])
      
      Z_ratio_term1 <- 1 + (a0 * (t(X0[i, ]) %*% cov_n_0 %*% X0[i, ]))[1,1]
      Z_ratio_term2 <- gamma(0.5*(a0 * sum(c01) + shape0)) / gamma(0.5*(a0 * sum(c00) + shape0))
      Z_ratio_term3 <- exp(0.5 * (a0 * sum(c00) + shape0) * log(term_0) - 0.5 * (a0 * sum(c01) + shape0) * log(term_1))
      
      Z_ratio  <- Z_ratio_term1^(-0.5) * Z_ratio_term2 * Z_ratio_term3 * (2*3.1415926)^(-a0/2)
      xi_tilde <- ((1 + Z_ratio * ((1 - xi_new) / xi_new) * L0i^(-a0) )[1,1])^(-1)
      
      ## Sample c0[i]
      c0_new[i] <- rbinom(1, 1, xi_tilde)
    }
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['c0']] <<- c0_new
    model[['xi']] <<- xi_new
    model[['n0exch']] <<- sum(c0_new)
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

sampler_c0_a0fix_conjugate <- nimbleFunction(
  name = 'c0_conjugate_sampler'  ## give name for sampler
  , contains = sampler_BASE  ## doesn't change
  , setup = function(model, mvSaved, target, control) {  ## arguments always the same
    ## Place only things that are fixed here
    y              <- model$y
    X              <- model$X
    y0             <- model$y0
    X0             <- model$X0
    a0             <- model$a0
    mu0            <- model$beta_init_mean
    prec0          <- model$beta_init_prec
    shape0         <- model$tau_init_shape
    rate0          <- model$tau_init_rate
    xi_shape1      <- model$xi_shape1
    xi_shape2      <- model$xi_shape2
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    n              <- length(y)
    n0             <- length(y0)
    p              <- ncol(X0)
    prec0_mu0      <- (prec0 %*% mu0)[, 1]
    mu0t_prec0_mu0 <- (t(mu0) %*% prec0_mu0)[1,1]
    eye_p          <- diag(1, p)
  }
  ## This is the function that runs when the sampler is called
  , run = function() {
    c0     <- model[['c0']]
    beta   <- model[['beta']]
    tau    <- model[['tau']]
    c0_new <- c0
    
    ## Sample xi
    xi_shape1_new <- sum(c0) + xi_shape1
    xi_shape2_new <- n0 - sum(c0) + xi_shape2
    xi_new <- rbeta(1, xi_shape1_new, xi_shape2_new)
    
    for (i in 1:n0) {
      L0i <- (tau / (2*3.1415926))^(1/2) * exp((-tau/2) * (y0[i] - (X0[i, ] %*% beta))^2)
      
      # Compute component of Z(c0 = 1)
      c01 <- c0
      c01[i] <- 1
      C01 <- diag(c01)
      
      c0y0_1   <- c01 * y0
      c0X0_1   <- C01 %*% X0
      X0C0X0_1 <- t(c0X0_1) %*% c0X0_1
      X0c0y0_1 <- (t(c0X0_1) %*% c0y0_1)[, 1]
      y0c0y0_1 <- (t(c0y0_1) %*% c0y0_1)[1,1]
      
      prec_n_1     <- a0 * X0C0X0_1 + prec0
      Uprec_n_1    <- chol(prec_n_1)
      Uprec_n_inv1 <- backsolve(Uprec_n_1, eye_p)
      cov_n_1      <- Uprec_n_inv1 %*% t(Uprec_n_inv1)
      mu_n_1       <- cov_n_1 %*% (a0 * X0c0y0_1 + prec0_mu0)
      
      term_1      <- 0.5*(rate0 + a0*y0c0y0_1 + mu0t_prec0_mu0 - (t(mu_n_1) %*% prec_n_1 %*% mu_n_1)[1,1])
      
      # Compute component of Z(c0 = 0)
      c00 <- c0
      c00[i] <- 0
      C00 <- diag(c00)
      
      c0y0_0   <- c00 * y0
      c0X0_0   <- C00 %*% X0
      X0C0X0_0 <- t(c0X0_0) %*% c0X0_0
      X0c0y0_0 <- (t(c0X0_0) %*% c0y0_0)[, 1]
      y0c0y0_0 <- (t(c0y0_0) %*% c0y0_0)[1,1]
      
      prec_n_0     <- a0 * X0C0X0_0 + prec0
      Uprec_n_0    <- chol(prec_n_0)
      Uprec_n_inv0 <- backsolve(Uprec_n_0, eye_p)
      cov_n_0      <- Uprec_n_inv0 %*% t(Uprec_n_inv0)
      mu_n_0       <- cov_n_0 %*% (a0 * X0c0y0_0 + prec0_mu0)
      
      term_0      <- 0.5*(rate0 + a0*y0c0y0_0 + mu0t_prec0_mu0 - (t(mu_n_0) %*% prec_n_0 %*% mu_n_0)[1,1])
      
      Z_ratio_term1 <- 1 + (a0 * (t(X0[i, ]) %*% cov_n_0 %*% X0[i, ]))[1,1]
      Z_ratio_term2 <- gamma(0.5*(a0 * sum(c01) + shape0)) / gamma(0.5*(a0 * sum(c00) + shape0))
      Z_ratio_term3 <- exp(0.5 * (a0 * sum(c00) + shape0) * log(term_0) - 0.5 * (a0 * sum(c01) + shape0) * log(term_1))
      
      Z_ratio  <- Z_ratio_term1^(-0.5) * Z_ratio_term2 * Z_ratio_term3 * (2*3.1415926)^(-a0/2)
      xi_tilde <- ((1 + Z_ratio * ((1 - xi_new) / xi_new) * L0i^(-a0) )[1,1])^(-1)
      
      ## Sample c0[i]
      c0_new[i] <- rbinom(1, 1, xi_tilde)
    }
    
    ## Store sampled value and recalculate (note: <<- is necessary to store samples)
    model[['c0']] <<- c0_new
    model[['xi']] <<- xi_new
    model[['n0exch']] <<- sum(c0_new)
    model$calculate(target)  ## always needed at end
    nimCopy(from = model, to = mvSaved, row = 1, nodes = target, logProb = TRUE)  ## always needed at end
  }
  ## Leave this as is
  , methods = list(
    reset = function() {}
  )
)

lm_conjugate_sampler <- nimbleFunction(
  name = 'lm_conjugate_sampler'  ## give name for sampler
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
    y              <- model$y
    X              <- model$X
    XtX            <- t(X) %*% X
    Xty            <- (t(X) %*% y)[, 1]
    yty            <- (t(y) %*% y)[1,1]
    n              <- length(y)
    a0 <- model[['a0']]
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

# ------------------------------------------------------------------------------
# Nimble code for DIPP and power prior

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


dipp.xifixed <- nimbleCode({
  
  ## Conditional prior on latent binary variable | exchangeability probability
  a0 ~ dbeta(a0_shape1, a0_shape2)
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


dipp.fixed <- nimbleCode({
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
## Compile all samplers
## -----------------------------------------------------------------------------

# NDIPP
const <- list(n = n, n0 = n0, p = p)
data  <- list(
  y = y, X = X, y0 = y0, X0 = X0
  , beta_init_mean = beta_init_mean
  , beta_init_prec = beta_init_prec
  , tau_init_shape = tau_init_shape
  , tau_init_rate = tau_init_rate
  , xi_shape1 = xi_shape1
  , xi_shape2 = xi_shape2
  , eye = diag(n)
  , a0_shape1 = a0_shape1
  , a0_shape2 = a0_shape2
)
init <- list(c0 = rbinom(n0, 1, c0_init_prob), xi = xi_init, a0 = a0_init)
model    <- nimbleModel(ndipp, constants = const, data = data, inits = init)
cmodel   <- compileNimble(model)
mcmcConf <- configureMCMC(model, monitors = c('beta', 'xi', 'tau', 'c0', 'n0exch', 'a0'))
mcmcConf$removeSampler(c('beta', 'tau'))
mcmcConf$addSampler(c('beta', 'tau'), 'sampler_dipp_betatau_a0ran_conjugate')
mcmc     <- buildMCMC(mcmcConf)
cmcmc.ndipp <- compileNimble(mcmc, project = cmodel)

# DIPP, random a0, fixed xi
const1 <- list(n = n, n0 = n0, p = p)
data1  <- list(
  y = y, X = X, y0 = y0, X0 = X0
  , beta_init_mean = beta_init_mean
  , beta_init_prec = beta_init_prec
  , tau_init_shape = tau_init_shape
  , tau_init_rate = tau_init_rate
  , eye = diag(n)
  , a0_shape1 = a0_shape1
  , a0_shape2 = a0_shape2
  , xi = xi
)
init1 <- list(c0 = rbinom(n0, 1, c0_init_prob), a0 = a0_init)
model1    <- nimbleModel(dipp.xifixed, constants = const1, data = data1, inits = init1)
cmodel1   <- compileNimble(model1)
mcmcConf1 <- configureMCMC(model1, monitors = c('beta', 'tau', 'c0', 'n0exch', 'a0'))
mcmcConf1$removeSampler(c('beta', 'tau'))
mcmcConf1$addSampler(c('beta', 'tau'), 'sampler_dipp_betatau_a0ran_conjugate')
mcmc1     <- buildMCMC(mcmcConf1)
cmcmc.dipp.xifixed <- compileNimble(mcmc1, project = cmodel1)

# DIPP
const2 <- list(n = n, n0 = n0, p = p)
data2  <- list(
  y = y, X = X, y0 = y0, X0 = X0
  , beta_init_mean = beta_init_mean
  , beta_init_prec = beta_init_prec
  , tau_init_shape = tau_init_shape
  , tau_init_rate = tau_init_rate
  , xi_shape1 = xi_shape1
  , xi_shape2 = xi_shape2
  , eye = diag(n)
  , a0 = a0
)
init2 <- list(c0 = rbinom(n0, 1, c0_init_prob), xi = xi_init)
model2    <- nimbleModel(dipp, constants = const2, data = data2, inits = init2)
cmodel2   <- compileNimble(model2)
mcmcConf2 <- configureMCMC(model2, monitors = c('beta', 'xi', 'tau', 'c0', 'n0exch'))
mcmcConf2$removeSampler(c('beta', 'tau'))
mcmcConf2$addSampler(c('beta', 'tau'), 'sampler_dipp_betatau_a0fix_conjugate')
mcmc2     <- buildMCMC(mcmcConf2)
cmcmc.dipp <- compileNimble(mcmc2, project = cmodel2)

# normalized power prior
const3 <- list(n = n, n0 = n0, p = p)
data3  <- list(y = y, X = X, y0 = y0, X0 = X0
               , beta_init_mean = beta_init_mean
               , beta_init_prec = beta_init_prec
               , tau_init_shape = tau_init_shape
               , tau_init_rate = tau_init_rate
               , eye = diag(n)
               , a0_shape1 = a0_shape1
               , a0_shape2 = a0_shape2
)
init3  <- list(beta = beta_init, tau = tau_init, a0 = a0_init)
model3    <- nimbleModel(npp, constants = const3, data = data3, inits = init3)
cmodel3   <- compileNimble(model3)
mcmcConf3 <- configureMCMC(cmodel3, monitors = c('beta', 'tau', 'a0'))
mcmc3     <- buildMCMC(mcmcConf3)
cmcmc.pp.a0random <- compileNimble(mcmc3, project = cmodel3)

# power prior
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
model4    <- nimbleModel(pp, constants = const4, data = data4, inits = init4)
cmodel4   <- compileNimble(model4)
mcmcConf4 <- configureMCMC(cmodel4, monitors = c('beta', 'tau'))
mcmc4     <- buildMCMC(mcmcConf4)
cmcmc.pp.a0fixed <- compileNimble(mcmc4, project = cmodel4)

# DIPP, fixed a0 and xi
const6 <- list(n = n, n0 = n0, p = p)
data6  <- list(y = y, X = X, y0 = y0, X0 = X0
               , beta_init_mean = beta_init_mean
               , beta_init_prec = beta_init_prec
               , tau_init_shape = tau_init_shape
               , tau_init_rate = tau_init_rate
               , eye = diag(n)
               , a0 = a0
               , xi = xi
)
init6  <- list(beta = beta_init, tau = tau_init)
model6    <- nimbleModel(dipp.fixed, constants = const6, data = data6, inits = init6)
cmodel6   <- compileNimble(model6)
mcmcConf6 <- configureMCMC(cmodel6, monitors = c('beta', 'tau', 'c0', 'n0exch'))
mcmcConf6$removeSampler(c('beta', 'tau'))
mcmcConf6$addSampler(c('beta', 'tau'), 'sampler_dipp_betatau_a0fix_conjugate')
mcmc6     <- buildMCMC(mcmcConf6)
cmcmc.dipp.fixed <- compileNimble(mcmc6, project = cmodel6)


# ------------------------------------------------------------------------------
# Function to compute DIC

logL <- function(y, X, beta, tau) {
  n <- length(y)
  mu <- as.vector(X %*% beta)
  
  (n / 2) * log(tau / (2 * pi)) - sum((tau / 2) * (y - mu)^2)
}

DIC <- function(y, X, data) {
  # extract beta columns automatically
  beta_idx <- grep("^beta\\[", colnames(data))
  beta_mat <- as.matrix(data[, beta_idx, drop = FALSE])
  tau      <- data[, "tau"]
  
  niter <- nrow(beta_mat)
  
  # posterior means
  beta_est <- colMeans(beta_mat)
  tau_est  <- mean(tau)
  
  # log-likelihood at each iteration
  log_likelihoods <- vapply(
    seq_len(niter),
    function(i) logL(y, X, beta_mat[i, ], tau[i]),
    numeric(1)
  )
  
  # effective number of parameters
  pDIC <- 2 * (
    logL(y, X, beta_est, tau_est) -
      mean(log_likelihoods)
  )
  
  # DIC
  DIC <- -2 * logL(y, X, beta_est, tau_est) + 2 * pDIC
  
  return(DIC)
}

# ------------------------------------------------------------------------------
# Get posterior samples from each prior method

niter  = 200000
burnin = 2000
thin   = 2

# a0 = 0.3, xi ~ Beta(2,2)
set.seed(8090)
cmodel2$a0 <- 0.3
smpl.dipp <- runMCMC(cmcmc.dipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp <- summarize_draws(smpl.dipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp)

# a0 = 0.5, xi ~ Beta(2,2)
set.seed(8090)
cmodel2$a0 <- 0.5
smpl.dipp <- runMCMC(cmcmc.dipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp <- summarize_draws(smpl.dipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp)

# a0 = 0.8, xi ~ Beta(2,2)
set.seed(8090)
cmodel2$a0 <- 0.8
smpl.dipp <- runMCMC(cmcmc.dipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp <- summarize_draws(smpl.dipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp)

# a0 = 1, xi ~ Beta(2,2)
set.seed(8090)
cmodel2$a0 <- 1
smpl.dipp <- runMCMC(cmcmc.dipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp <- summarize_draws(smpl.dipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp)

# a0 ~ Beta(2,2), xi ~ Beta(2,2)
set.seed(8090)
smpl.ndipp <- runMCMC(cmcmc.ndipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.ndipp <- summarize_draws(smpl.ndipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.ndipp)

# a0 ~ Beta(5,2), xi ~ Beta(2,2)
set.seed(8090)
cmodel$a0_shape1 <- 5
smpl.ndipp <- runMCMC(cmcmc.ndipp, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.ndipp <- summarize_draws(smpl.ndipp, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.ndipp)

# a0 ~ Beta(2,2), xi = 0.5
set.seed(8090)
smpl.dipp.xifixed <- runMCMC(cmcmc.dipp.xifixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp.xifixed <- summarize_draws(smpl.dipp.xifixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp.xifixed)

# a0 = 1, xi = 0.5
set.seed(8090)
cmodel6$a0 <- 1
cmodel6$xi <- 0.5
smpl.dipp.fixed <- runMCMC(cmcmc.dipp.fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp.fixed <- summarize_draws(smpl.dipp.fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp.fixed)

# a0 = 1, xi = 0.8
set.seed(8090)
cmodel6$a0 <- 1
cmodel6$xi <- 0.8
smpl.dipp.fixed2 <- runMCMC(cmcmc.dipp.fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp.fixed2 <- summarize_draws(smpl.dipp.fixed2, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.dipp.fixed2)

# power prior, a0 = 0.3
set.seed(8090)
cmodel4$a0 <- 0.3
smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.pp.a0fixed)

# power prior, a0 = 0.5
set.seed(8090)
cmodel4$a0 <- 0.5
smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.pp.a0fixed)

# power prior, a0 = 0.8
set.seed(8090)
cmodel4$a0 <- 0.8
smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.pp.a0fixed)

# power prior, a0 = 1
set.seed(8090)
cmodel4$a0 <- 1
smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.pp.a0fixed)

# power prior, a0 ~ beta(2,2)
set.seed(8090)
smpl.pp.a0random <- runMCMC(cmcmc.pp.a0random, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0random <- summarize_draws(smpl.pp.a0random, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.pp.a0random)


# -------------------------------------------------------------------------------------
# samples with pre-processed historical data

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

# Obtain matched historical subjects
matched_hist <- match.data(match_obj) %>% filter(current == 0)

# Assign data to variables
y0 <- matched_hist$LOGCHG
X0 <- cbind(1, matched_hist$trtgroup, matched_hist$BASE)
n0 <- length(y0)



# take samples to create posterior density plots

set.seed(8090)
smpl.dipp.a0fixed <- runMCMC(cmcmc.dipp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp.a0fixed <- summarize_draws(smpl.dipp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 

smpl.dipp.a0random <- runMCMC(cmcmc.dipp.a0random, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.dipp.a0random <- summarize_draws(smpl.dipp.a0random, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 

smpl.pp.a0fixed <- runMCMC(cmcmc.pp.a0fixed, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0fixed <- summarize_draws(smpl.pp.a0fixed, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 

smpl.pp.a0random <- runMCMC(cmcmc.pp.a0random, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.pp.a0random <- summarize_draws(smpl.pp.a0random, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 


## Plot posterior density
data.dipp.a0fixed <- data.frame(smpl.dipp.a0fixed)
data.dipp.a0random <- data.frame(smpl.dipp.a0random)
data.dipp.xifixed <- data.frame(smpl.dipp.xifixed)
data.pp.a0fixed <- data.frame(smpl.pp.a0fixed)
data.pp.a0random <- data.frame(smpl.pp.a0random)
data.ref <- data.frame(smpl.ref)

all_data <- bind_rows(
  data.dipp.a0fixed %>% mutate(Type = "dipp.a0fixed"),
  data.dipp.a0random %>% mutate(Type = "dipp.a0random"),
  data.pp.a0fixed %>% mutate(Type = "pp.a0fixed"),
  data.pp.a0random %>% mutate(Type = "pp.a0random"),
  data.ref %>% mutate(Type = "ref")
)

b1 <- ggplot(all_data, aes(x = beta.1., color = Type, linetype = Type)) + 
  geom_density(size = 0.6) +
  labs(title = expression(paste("Intercept (", italic(β[1]), ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0fixed" = "black", "dipp.a0random" = "grey50", 
                                "pp.a0fixed" = "black", "pp.a0random" = "grey50",
                                "ref" = "red"),
                     labels = c("dipp.a0fixed" = expression(paste("DIPP, ", a[0], " fixed")), 
                                "dipp.a0random" = expression(paste("DIPP, ", a[0], " random")), 
                                "pp.a0fixed" = expression(paste("Power prior, ", a[0], " fixed")), 
                                "pp.a0random" = expression(paste("Power prior, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "dashed", "solid")))) +
  scale_linetype_manual(values = c("dipp.a0fixed" = "solid", "dipp.a0random" = "solid", 
                                   "pp.a0fixed" = "dashed", "pp.a0random" = "dashed", "ref" = "solid"))

b2 <- ggplot(all_data, aes(x = beta.2., color = Type, linetype = Type)) + 
  geom_density(size = 0.6) +
  labs(title = expression(paste("Treatment Effect (", italic(β[2]), ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0fixed" = "black", "dipp.a0random" = "grey50", 
                                "pp.a0fixed" = "black", "pp.a0random" = "grey50",
                                "ref" = "red"),
                     labels = c("dipp.a0fixed" = expression(paste("DIPP, ", a[0], " fixed")), 
                                "dipp.a0random" = expression(paste("DIPP, ", a[0], " random")), 
                                "pp.a0fixed" = expression(paste("Power prior, ", a[0], " fixed")), 
                                "pp.a0random" = expression(paste("Power prior, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "dashed", "solid")))) +
  scale_linetype_manual(values = c("dipp.a0fixed" = "solid", "dipp.a0random" = "solid", 
                                   "pp.a0fixed" = "dashed", "pp.a0random" = "dashed", "ref" = "solid"))

b3 <- ggplot(all_data, aes(x = beta.3., color = Type, linetype = Type)) + 
  geom_density(size = 0.6) +
  labs(title = expression(paste("Baseline C1-INH (", italic(β[3]), ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0fixed" = "black", "dipp.a0random" = "grey50", 
                                "pp.a0fixed" = "black", "pp.a0random" = "grey50",
                                "ref" = "red"),
                     labels = c("dipp.a0fixed" = expression(paste("DIPP, ", a[0], " fixed")), 
                                "dipp.a0random" = expression(paste("DIPP, ", a[0], " random")), 
                                "pp.a0fixed" = expression(paste("Power prior, ", a[0], " fixed")), 
                                "pp.a0random" = expression(paste("Power prior, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "dashed", "solid")))) +
  scale_linetype_manual(values = c("dipp.a0fixed" = "solid", "dipp.a0random" = "solid", 
                                   "pp.a0fixed" = "dashed", "pp.a0random" = "dashed", "ref" = "solid"))


tau_plot <- ggplot(all_data, aes(x = tau, color = Type, linetype = Type)) + 
  geom_density(size = 0.6) + 
  labs(title = expression(paste("Precision (", tau, ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0fixed" = "black", "dipp.a0random" = "grey50", 
                                "pp.a0fixed" = "black", "pp.a0random" = "grey50",
                                "ref" = "red"),
                     labels = c("dipp.a0fixed" = expression(paste("DIPP, ", a[0], " fixed")), 
                                "dipp.a0random" = expression(paste("DIPP, ", a[0], " random")), 
                                "pp.a0fixed" = expression(paste("Power prior, ", a[0], " fixed")), 
                                "pp.a0random" = expression(paste("Power prior, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "dashed", "solid")))) +
  scale_linetype_manual(values = c("dipp.a0fixed" = "solid", "dipp.a0random" = "solid", 
                                   "pp.a0fixed" = "dashed", "pp.a0random" = "dashed", "ref" = "solid"))


a0 <- ggplot(all_data, aes(x = a0, color = Type, linetype = Type)) + 
  geom_density(size = 0.6) + 
  labs(title = expression(paste("Global Discounting Parameter (", italic(a[0]), ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0random" = "grey50", "pp.a0random" = "grey50"),
                     labels = c("dipp.a0random" = expression(paste("DIPP, ", a[0], " random")), 
                                "pp.a0random" = expression(paste("Power prior, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "dashed")))) +
  scale_linetype_manual(values = c("dipp.a0random" = "solid", "pp.a0random" = "dashed"))


xi <- ggplot(all_data, aes(x = xi, color = Type, linetype = Type)) + 
  geom_density() + 
  labs(title = expression(paste("Latent Parameter for Individual Discounting Parameter (", italic(ξ), ")")),
       y = "Density") +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(name = "Prior", 
                     values = c("dipp.a0fixed" = "black", "dipp.a0random" = "grey50"),
                     labels = c("dipp.a0fixed" = expression(paste("DIPP, ", a[0], " fixed")), 
                                "dipp.a0random" = expression(paste("DIPP, ", a[0], " random")),
                                "ref" = "Reference Prior"),
                     guide = guide_legend(override.aes = list(linetype = c("solid", "solid")))) +
  scale_linetype_manual(values = c("dipp.a0fixed" = "solid", "dipp.a0random" = "solid"))


b1 + b2 + b3 + tau_plot + a0 + xi
 





# ------------------------------------------------------------------------------
# Analysis with no historical data included, ie. a0 = 0

ref.nimble <- nimbleCode({
  ## Compute prior mean / precision of beta | tau
  omegatilde0[1:p,1:p]  <- beta_init_prec[1:p,1:p]
  beta_init_prec_mean[1:p] <- beta_init_prec[1:p,1:p] %*% beta_init_mean[1:p]
  beta_mean[1:p] <- solve( omegatilde0[1:p,1:p], beta_init_prec_mean[1:p] )
  beta_prec[1:p,1:p] <- tau * omegatilde0[1:p,1:p]
  
  # ## Compute prior shape and rate parameters on precision
  tau_shape <- 0.5 * (tau_init_shape)
  tau_rate  <- 0.5 * (tau_init_rate
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

const5 <- list(n = n, p = p)
data5  <- list(y = y, X = X
               , beta_init_mean = beta_init_mean
               , beta_init_prec = beta_init_prec
               , tau_init_shape = tau_init_shape
               , tau_init_rate = tau_init_rate
               , eye = diag(n)
)

init5  <- list(beta = beta_init, tau = tau_init)

model5    <- nimbleModel(ref.nimble, constants = const5, data = data5, inits = init5)
cmodel5   <- compileNimble(model5)
mcmcConf5 <- configureMCMC(cmodel5, monitors = c('beta', 'tau'))
mcmc5     <- buildMCMC(mcmcConf5)
cmcmc5    <- compileNimble(mcmc5, project = cmodel5)

set.seed(8090)
model5$tau_init_shape <- 2
model5$tau_init_rate  <- 2
smpl.ref    <- runMCMC(cmcmc5, niter = burnin + thin * niter, nburnin = burnin, thin = thin)
summary.ref <- summarize_draws(smpl.ref, "mean", "sd", ~quantile(.x, probs = c(0.025, 0.975))) 
DIC(y, X, smpl.ref)

beta1_est <- mean(smpl.ref[,1])
beta2_est <- mean(smpl.ref[,2])
beta3_est <- mean(smpl.ref[,3])
tau_est   <- mean(smpl.ref[,4])

log_likelihoods <- sapply(1:nrow(smpl.ref), function(i) {
  logL(y, X, smpl.ref[i, 1], smpl.ref[i, 2], smpl.ref[i, 3], smpl.ref[i, 4])
})

pDIC <- 2*(logL(y, X, beta1_est, beta2_est, beta3_est, tau_est) - mean(log_likelihoods))
DIC  <- -2*logL(y, X, beta1_est, beta2_est, beta3_est, tau_est) + 2*pDIC

DIC
