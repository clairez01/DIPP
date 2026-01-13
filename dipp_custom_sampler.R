

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









