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
setwd("/work/users/c/l/clairez1/Paper1sims2/PSIPP")
# Load in grid of simulation parameters
grid <- readRDS("/work/users/c/l/clairez1/Paper1sims2/grid.rds")

# List all .rds files in the directory
file_list <- list.files(pattern = "\\.rds$")

numbers <- as.numeric(gsub("id_(\\d+)_.*", "\\1", file_list))

# Sort file list based on the extracted numbers
sorted_file_list <- file_list[order(numbers, file_list)]

# Print the sorted file list
print(sorted_file_list)


# Initialize an empty list to store dataframes
dfs <- list()
results.all <- list()
sim.summary <- list()

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
    # Power prior
    beta1.pp  <- results.all %>% filter(variable == "beta[1]" & method == "psipp")
    beta2.pp  <- results.all %>% filter(variable == "beta[2]" & method == "psipp")
    beta3.pp  <- results.all %>% filter(variable == "beta[3]" & method == "psipp")
    tau.pp    <- results.all %>% filter(variable == "tau" & method == "psipp")
    
    # Compute bias
    beta1.bias  <- mean(beta1.pp$mean - beta[1])
    beta2.bias  <- mean(beta2.pp$mean - beta[2])
    beta3.bias  <- mean(beta3.pp$mean - beta[3])
    
    # Compute MSE
    beta1.MSE  <- mean( (beta1.pp$mean - beta[1])^2 )
    beta2.MSE  <- mean( (beta2.pp$mean - beta[2])^2 )
    beta3.MSE  <- mean( (beta3.pp$mean - beta[3])^2 )
    
    # Compute coverage probability
    beta1.covprob  <- sum(beta1.pp$X2.5. <= beta[1] & beta1.pp$X97.5. >= beta[1]) / nrow(beta1.pp)
    beta2.covprob  <- sum(beta2.pp$X2.5. <= beta[2] & beta2.pp$X97.5. >= beta[2]) / nrow(beta2.pp)
    beta3.covprob  <- sum(beta3.pp$X2.5. <= beta[3] & beta3.pp$X97.5. >= beta[3]) / nrow(beta3.pp)
    
    # Compute 95% credible interval width
    beta1.width  <- mean(beta1.pp$X97.5. - beta1.pp$X2.5.)
    beta2.width  <- mean(beta2.pp$X97.5. - beta2.pp$X2.5.)
    beta3.width  <- mean(beta3.pp$X97.5. - beta3.pp$X2.5.)
    
    beta1.summary  <- cbind(mean(beta1.pp$mean), sd(beta1.pp$mean), beta1.bias, beta1.MSE, beta1.covprob, beta1.width)
    beta2.summary  <- cbind(mean(beta2.pp$mean), sd(beta2.pp$mean), beta2.bias, beta2.MSE, beta2.covprob, beta2.width)
    beta3.summary  <- cbind(mean(beta3.pp$mean), sd(beta3.pp$mean), beta3.bias, beta3.MSE, beta3.covprob, beta3.width)
    tau.summary    <- cbind(mean(tau.pp$mean), sd(tau.pp$mean), NA, NA, NA, NA)
    
    sim.summary.pp <- data.frame(rbind(beta1.summary, beta2.summary, beta3.summary, tau.summary))
    sim.summary.pp$method <- "pp"
    sim.summary.pp$n0 <- scen$n0
    sim.summary.pp$q  <- scen$q
    sim.summary.pp$prob.unexch <- scen$prob.unexch
    sim.summary.pp$param <- c("beta1", "beta2", "beta3", "tau")
    colnames(sim.summary.pp) <- c("mean", "sd", "bias", "MSE", "Coverage Prob.", "CI Width", "method", "n0", "q", "prob.unexch", "parameter")
    
    # Store summarized sims
    sim.summary <- rbind(sim.summary, sim.summary.pp)
    return(sim.summary)
  }
  
  # Clear the dfs list for the next group
  dfs <- list()
}

# Loop through each group of 5 IDs
compiled.results <- list()
for (j in 1:120) {
  compiled.results <- rbind(compiled.results, combine_datasets(5*j - 4, 5*j))
}

mean.20 <- mean(compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.30 <- mean(compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00])
mean.40 <- mean(compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00] <- mean.30
compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00] <- mean.40
compiled.results$bias[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.30 <- mean(compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00])
mean.40 <- mean(compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00] <- mean.30
compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00] <- mean.40
compiled.results$MSE[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.30 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00])
mean.40 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00] <- mean.30
compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00] <- mean.40
compiled.results$`Coverage Prob.`[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


mean.20 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00])
mean.30 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00])
mean.40 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00])
mean.50 <- mean(compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00])

compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 20 & compiled.results$prob.unexch == 0.00] <- mean.20
compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 30 & compiled.results$prob.unexch == 0.00] <- mean.30
compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 40 & compiled.results$prob.unexch == 0.00] <- mean.40
compiled.results$`CI Width`[compiled.results$parameter == "beta2" & compiled.results$n0 == 50 & compiled.results$prob.unexch == 0.00] <- mean.50


setwd("/nas/longleaf/home/clairez1/R/Paper1sims2_files")
saveRDS(compiled.results, file = "sim_summary_psipp.rds")
