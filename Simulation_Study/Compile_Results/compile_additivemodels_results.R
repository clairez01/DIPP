# Load in libraries
library(haven)
library(tidyverse)
library(mvtnorm)

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
p  <- ncol(X)

# "true" parameters
beta <- lm$coefficients[,1]
sd   <- lm$sigma
tau  <- 1/(sd)^2



# Set the working directory
setwd("/work/users/c/l/clairez1/Paper1sims3/ALL_model")
# Load in grid of simulation parameters
grid_add <- readRDS("/work/users/c/l/clairez1/Paper1sims3/grid_additive.rds")
grid <- subset(grid_add, n0 == 20 | n0 == 50)


# List all .rds files in the directory
file_list <- list.files(pattern = "\\.rds$")
file_list <- file_list[file_list != "grid.rds"]

numbers <- as.numeric(gsub("id_(\\d+)_.*", "\\1", file_list))

# Sort file list based on the extracted numbers
sorted_file_list <- file_list[order(numbers, file_list)]

# Print the sorted file list
print(sorted_file_list)

# Define a function to combine datasets for a given range of IDs

combine_datasets <- function(start_id, end_id) {
  for (id in start_id:end_id) {
    grid.id <- grid[id,]
    n <- grid.id$n
    n0 <- grid.id$n0
    q <- grid.id$q
    pr <- grid.id$prob.unexch
    
    # Check if the file exists for the current ID
    file <- paste0('id_', id, '_', 'n_', n, '_n0_', n0, '_q_', q, '_probunexch_', pr, '.rds')
    if (file.exists(file)) {
      # Load the data from the .rds file
      data <- readRDS(file)$simres
      # Extract ID 
      id <- readRDS(file)$id
      # Extract scenario
      scen <- readRDS(file)$simscen
      # Append dataframe to the list
      dfs[[id]] <- data
    }
  }
  
  if (length(dfs) == 0) {
    return(NULL)
  } else {
    
    # Combine the datasets
    results.all <- do.call(rbind, dfs)
    
    # Perform computations and summaries
    dipp  <- results.all %>% filter(variable == "beta[2]" & method == "dipp")
    ndipp <- results.all %>% filter(variable == "beta[2]" & method == "ndipp")
    pp    <- results.all %>% filter(variable == "beta[2]" & method == "pp")
    npp   <- results.all %>% filter(variable == "beta[2]" & method == "npp")
    leap  <- results.all %>% filter(variable == "trteffect" & method == "leap")
    
    # dipp
    bias     <- mean(dipp$mean - beta[2])
    MSE      <- mean( (dipp$mean - beta[2])^2 )
    covprob  <- sum(dipp$X2.5. <= beta[2] & dipp$X97.5. >= beta[2]) / nrow(dipp)
    width    <- mean(dipp$X97.5. - dipp$X2.5.)
    dipp.summary  <- data.frame(cbind(mean(dipp$mean), sd(dipp$mean), bias, MSE, covprob, width))
    dipp.summary$method <- "dipp"
    
    # ndipp
    bias    <- mean(ndipp$mean - beta[2])
    MSE     <- mean((ndipp$mean - beta[2])^2)
    covprob <- sum(ndipp$X2.5. <= beta[2] & ndipp$X97.5. >= beta[2]) / nrow(ndipp)
    width   <- mean(ndipp$X97.5. - ndipp$X2.5.)
    ndipp.summary <- data.frame(cbind(mean(ndipp$mean), sd(ndipp$mean), bias, MSE, covprob, width))
    ndipp.summary$method <- "ndipp"
    
    # power prior
    bias    <- mean(pp$mean - beta[2])
    MSE     <- mean((pp$mean - beta[2])^2)
    covprob <- sum(pp$X2.5. <= beta[2] & pp$X97.5. >= beta[2]) / nrow(pp)
    width   <- mean(pp$X97.5. - pp$X2.5.)
    pp.summary <- data.frame(cbind(mean(pp$mean), sd(pp$mean), bias, MSE, covprob, width))
    pp.summary$method <- "pp"
    
    # npp
    bias    <- mean(npp$mean - beta[2])
    MSE     <- mean((npp$mean - beta[2])^2)
    covprob <- sum(npp$X2.5. <= beta[2] & npp$X97.5. >= beta[2]) / nrow(npp)
    width   <- mean(npp$X97.5. - npp$X2.5.)
    npp.summary <- data.frame(cbind(mean(npp$mean), sd(npp$mean), bias, MSE, covprob, width))
    npp.summary$method <- "npp"
    
    # LEAP 
    bias    <- mean(leap$mean - beta[2])
    MSE     <- mean((leap$mean - beta[2])^2)
    covprob <- sum(leap$X2.5. <= beta[2] & leap$X97.5. >= beta[2]) / nrow(leap)
    width   <- mean(leap$X97.5. - leap$X2.5.)
    leap.summary <- data.frame(cbind(mean(leap$mean), sd(leap$mean), bias, MSE, covprob, width))
    leap.summary$method <- "leap"
    
    sim.summary <- rbind(dipp.summary, ndipp.summary, pp.summary, npp.summary, leap.summary)
    sim.summary$n0 <- scen$n0
    sim.summary$q  <- scen$q
    sim.summary$prob.unexch <- scen$prob.unexch
    sim.summary$param <- c("trteffect")
    colnames(sim.summary) <- c("mean", "sd", "bias", "MSE", "Coverage Prob.", "CI Width", "method", "n0", "q", "prob.unexch", "parameter")
    
    # Store summarized sims
    sim.summary <- rbind(sim.summary, sim.summary)
    return(sim.summary)
  }
  
  # Clear the dfs list for the next group
  dfs <- list()
}


# Initialize an empty list to store dataframes
dfs <- list()
results.all <- list()
sim.summary <- list()

# Loop through each group of 5 IDs
compiled.results <- list()
for (j in 1:34) {
  compiled.results <- rbind(compiled.results, combine_datasets(20*j - 19, 20*j))
}


mean.20 <- mean(compiled.results$bias[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$bias[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$bias[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$bias[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$MSE[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$MSE[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$MSE[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$MSE[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$`Coverage Prob.`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$`Coverage Prob.`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$`CI Width`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$`CI Width`[compiled.results$parameter == "trteffect" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50

setwd("/nas/longleaf/home/clairez1/Paper1_sims_rev2/compiled_sim_results")
saveRDS(compiled.results, file = "sim_summary_add_model.rds")
