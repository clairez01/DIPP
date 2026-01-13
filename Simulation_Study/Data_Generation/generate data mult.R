
save.dir <- '/work/users/c/l/clairez1/Paper1sims3'


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

grid   <- readRDS("/work/users/c/l/clairez1/Paper1_sims_revisions/grid.rds")
grid_s <- grid[,2:4]
grid_unique <- grid_s[ !duplicated(grid_s[c("n0", "q", "prob.unexch")]), ]
grid_sorted <- grid_unique[ order(grid_unique$n0, grid_unique$prob.unexch, grid_unique$q), ]

set.seed(44040)

# Parameters
beta <- lmbeta
sd   <- lm$sigma

# ===========================================================
# Generate 10,000 CURRENT datasets stacked into one data frame
# ===========================================================
Nsims <- 10000
n     <- 50
sampled_curr <- current[sample(nrow(current), n*Nsims, replace = TRUE), ]
X  <- cbind(1, sampled_curr$trtgroup, sampled_curr$BASE) 
y  <- rnorm(n*Nsims, X %*% beta, sd)

i <- rep(1:Nsims, each = n)

master_current <- data.frame(
  i        = i,
  y        = y,
  trtgroup = sampled_curr$trtgroup,
  BASE     = sampled_curr$BASE
)

saveRDS(master_current, file = file.path(save.dir, "simdata_curr_mult.rds"))

# -----------------------------------------------------------
# Function to generate historical data for ONE grid scenario
# -----------------------------------------------------------
generate_hist_for_scenario <- function(n0.id, q.id, prob.unexch.id,
                                       historical, lmbeta, sd) {
  
  beta <- cbind(lmbeta, (lmbeta * q.id))
  
  # initialize
  X0_exch <- NULL
  y0_exch <- NULL
  X0_unexch <- NULL
  y0_unexch <- NULL
  X0_exch_best <- NULL
  y0_exch_best <- NULL
  X0_unexch_best <- NULL
  y0_unexch_best <- NULL
  min.dist.exch <- Inf
  min.dist.unexch <- Inf
  
  n.unexch <- round(n0.id * prob.unexch.id)
  n.exch   <- n0.id - n.unexch
  
  for (i in 1:10000) {
    
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
  
  master_historical <- data.frame(
    y0          = y0,
    trtgroup    = X0[,2],
    BASE        = X0[,3],
    n0          = n0.id,
    q           = q.id,
    prob.unexch = prob.unexch.id
  )
}


# -----------------------------------------------------------
# Apply it to your grid
# -----------------------------------------------------------
set.seed(404040)
results_list <- lapply(seq_len(nrow(grid_sorted)), function(i) {
  row <- grid_sorted[i, ]
  generate_hist_for_scenario(
    n0.id          = row$n0,
    q.id           = row$q,
    prob.unexch.id = row$prob.unexch,
    historical     = historical,
    lmbeta         = lmbeta,
    sd             = lm$sigma
  )
})

hist_df <- do.call(rbind, results_list)

saveRDS(hist_df, file = file.path(save.dir, "simdata_hist_mult.rds"))
