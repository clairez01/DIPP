library(nimble)
library(tidyverse)
library(posterior)
library(haven)
library(hdbayes)
library(cmdstanr)
library(patchwork)
library(dplyr)
library(MatchIt)
library(lmtest)
library(sandwich)


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
b  <- scale(c(current$BASE, historical$BASE))
lm <- summary(lm(y ~ 0 + X)) # regression coefficients
lmbeta <- lm$coefficients[,1]
tau <- 1/(lm$sigma)^2
p  <- ncol(X)

# ------------------------------------------------------------------------------
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
# ------------------------------------------------------------------------------

df <- as.data.frame(cbind(current$LOGCHG, current$trtgroup, current$BASE))
df0 <- as.data.frame(cbind(historical$LOGCHG, historical$trtgroup, historical$BASE))

df <- as.data.frame(cbind(current1$LOGCHG, current1$trtgroup, current1$BASE))
df0 <- as.data.frame(cbind(matched_hist$LOGCHG, matched_hist$trtgroup, matched_hist$BASE))


fit <- glm.leap(
  formula = V1 ~ V2 + V3,
  family = gaussian(link = "identity"),
  data.list = list(df, df0),
  K = 2,
  iter_warmup = 2000,
  iter_sampling = 100000,
  chains = 1
)

fit %>% summarize_draws() %>% print(n = 200)

data.frame(
  variable = "trteffect",
  mean = mean(fit$V2),
  sd = sd(fit$V2),
  X2.5. = quantile(fit$V2, 0.025),
  X97.5. = quantile(fit$V2, 0.975)
)

data.frame(
  variable = "prob",
  mean = mean(fit$`probs[1]`),
  sd = sd(fit$`probs[1]`),
  X2.5. = quantile(fit$`probs[1]`, 0.025),
  X97.5. = quantile(fit$`probs[1]`, 0.975)
)

lm(current$CHG ~ current$trtgroup + current$BASE)

# Compute DIC
compute_DIC_gaussian <- function(fit_draws, formula, data) {
  
  # Convert draws into a data.frame
  if (is.list(fit_draws) && !is.data.frame(fit_draws)) {
    draws_df <- as.data.frame(fit_draws)
  } else {
    draws_df <- as.data.frame(fit_draws)
  }
  
  # model matrix and response
  mf <- model.frame(formula, data)
  y  <- model.response(mf)
  X  <- model.matrix(formula, mf)
  n  <- length(y)
  S  <- nrow(draws_df)
  
  # coefficient names from model matrix
  coef_names <- colnames(X)
  beta_draws <- as.matrix(draws_df[, coef_names, drop = FALSE])
  
  # --- find sigma (std dev) ---
  sigma_draws <- NULL
  
  # 1) Try standard names first
  if (is.null(sigma_draws)) {
    for (nm in c("sigma","sd","sigma_y")) {
      if (nm %in% names(draws_df)) {
        sigma_draws <- draws_df[[nm]]
        break
      }
    }
  }
  
  # 2) Try precision (tau)
  if (is.null(sigma_draws)) {
    for (nm in c("tau","tau_y","precision")) {
      if (nm %in% names(draws_df)) {
        sigma_draws <- 1 / sqrt(draws_df[[nm]])
        break
      }
    }
  }
  
  # 3) Try *your* format: dispersion[1] = σ²
  if (is.null(sigma_draws)) {
    disp_name <- "dispersion[1]"
    if (disp_name %in% names(draws_df)) {
      sigma_draws <- sqrt(draws_df[[disp_name]])
    }
  }
  
  if (is.null(sigma_draws)) {
    stop("Could not find residual scale; add support for the appropriate column name.")
  }
  
  # compute deviance for each draw
  D_vec <- numeric(S)
  for (s in seq_len(S)) {
    beta_s  <- beta_draws[s, ]
    mu_s    <- as.numeric(X %*% beta_s)
    sigma_s <- sigma_draws[s]
    rss     <- sum((y - mu_s)^2)
    D_vec[s] <- n * log(2*pi*sigma_s^2) + rss / sigma_s^2
  }
  
  D_bar <- mean(D_vec)
  
  # posterior mean parameters
  beta_bar  <- colMeans(beta_draws)
  sigma_bar <- mean(sigma_draws)
  mu_bar    <- as.numeric(X %*% beta_bar)
  rss_bar   <- sum((y - mu_bar)^2)
  
  D_hat <- n * log(2*pi*sigma_bar^2) + rss_bar / sigma_bar^2
  
  DIC <- 2 * D_bar - D_hat
  pD  <- D_bar - D_hat
  
  list(DIC = DIC, pD = pD, D_bar = D_bar, D_hat = D_hat)
}

dic_res <- compute_DIC_gaussian(fit, formula = V1 ~ V2 + V3, data = df)

# if fit is a data.frame or matrix of draws:
# dic_res <- compute_DIC_gaussian(as.data.frame(fit_draws_matrix), formula = V1 ~ V2 + V4, data = df)

dic_res$DIC
dic_res$pD
