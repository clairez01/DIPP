# Load libraries
library(nimble)
library(tidyverse)
library(posterior)
library(haven)
library(hdbayes)
library(cmdstanr)

## Directory to store files
save.dir <- '/work/users/c/l/clairez1/Paper1sims3/LEAP_model'

# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper1_sims_revisions/grid.rds")

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
historical <- readRDS("/work/users/c/l/clairez1/Paper1sims3/simdata_hist_mult.rds")

curr_i <- current[1:n.id,]
hist <- subset(historical, n0 == n0.id & q == q.id & prob.unexch == prob.unexch.id)

df_i <- as.data.frame(cbind(curr_i$y, curr_i$trtgroup, curr_i$BASE))
df0  <- as.data.frame(cbind(hist$y, hist$trtgroup, hist$BASE))

# Begin sim code

# Simulation parameters
Nsims  = 2000
niter  = 10000
burnin = 2000
thin   = 2

# Create empty data frame to store simulation results
results.all <- data.frame()

start.time <- Sys.time()
for (i in 1:Nsims) {
  
  block  <- (id - 1) %% 5 + 1
  simID  <- (block - 1) * 2000 + i
  end    <- simID * n.id
  start  <- end - n.id + 1
  curr_i <- current[start:end,]
  df_i   <- as.data.frame(cbind(curr_i$y, curr_i$trtgroup, curr_i$BASE))
  
  sink(tempfile())
  fit <- glm.leap(
    formula = V1 ~ V2 + V3,
    family = gaussian(link = "identity"),
    data.list = list(df_i, df0),
    K = 2,
    iter_warmup = burnin,
    iter_sampling = niter,
    chains = 1
  )
  sink()
  
  leap.results <- data.frame(variable = "trteffect", mean = mean(fit$V2), sd = sd(fit$V2),
                              X2.5. = quantile(fit$V2, 0.025), X97.5. = quantile(fit$V2, 0.975)  )
  
  results.all <- rbind(results.all, leap.results)
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