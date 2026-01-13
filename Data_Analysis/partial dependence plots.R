# load libraries
library(nimble)
library(tidyverse)
library(posterior)
library(haven)
library(ggplot2)
library(patchwork)
library(dplyr)

# Load in CSL data
setwd("C:/Users/clair/OneDrive/Documents/UNC Chapel Hill/PhD Biostatistics/Dissertation/data set")
current <- read_sas("csl830_current.sas7bdat")
historical <- read_sas("csl830_external.sas7bdat")

# Log transform outcome variable
current$LOGCHG <- log(current$AVAL + 1) - log(current$BASE + 1)
historical$LOGCHG <- log(historical$AVAL + 1) - log(historical$BASE + 1)

# Assign data to y, X, y0, and x0
y  <- current$LOGCHG
X  <- cbind(1, current$trtgroup, current$BASE)
y0 <- historical$LOGCHG
X0 <- cbind(1, historical$trtgroup, historical$BASE)

c0_means <- c(0.614, 0.643, 0.617, 0.600, 0.611, 0.662, 0.628, 0.624, 0.620,
              0.650, 0.563, 0.602, 0.601, 0.616, 0.662, 0.584, 0.607, 0.589)

# Determine threshold for top half
threshold <- sort(c0_means, decreasing = TRUE)[9]

# Assign 1 to top half, 0 to bottom half
historical$exch1 <- ifelse(c0_means >= threshold, 1, 0)

df <- historical
# Create a group factor with proper labels
df <- df %>%
  mutate(group1 = factor(case_when(
    exch1 == 0 ~ "Less exchangeable group",
    exch1 == 1 ~ "More exchangeable group"
  )))

# Add a line for all data
df_all <- df %>% mutate(group1 = factor("All data"))
df_plot <- bind_rows(df, df_all)

# Smoothed plot
dipp.plot <- ggplot(df_plot, aes(x = BASE, y = LOGCHG, color = group1)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  labs( x = "Baseline C1-INH Functional Activity",
        y = "Outcome",
        color = "Group"
  ) +
  ggtitle(expression("DIPP with " * a[0] * " = 0.3 and " * xi * "~ Beta(2, 2)")) +
  theme_minimal() +
  scale_color_manual(values = c("All data" = "black",
                                "Less exchangeable group" = "red",
                                "More exchangeable group" = "blue"))




# NDIPP c0 means
c0_means2 <- c(0.508, 0.580, 0.509, 0.482, 0.509, 0.607, 0.551, 0.519, 0.512,
               0.568, 0.436, 0.473, 0.483, 0.510, 0.636, 0.451, 0.480, 0.457)

# Determine threshold for top half
threshold2 <- sort(c0_means2, decreasing = TRUE)[9]

# Assign 1 to top half, 0 to bottom half
historical$exch2 <- ifelse(c0_means2 >= threshold2, 1, 0)

# Create group factor for plotting
df2 <- historical %>%
  mutate(group2 = factor(case_when(
    exch2 == 0 ~ "Less exchangeable group",
    exch2 == 1 ~ "More exchangeable group"
  )))

# Add a line for all data
df2_all <- df2 %>% mutate(group2 = factor("All data"))
df2_plot <- bind_rows(df2, df2_all)

# Smoothed plot
ndipp.plot <- ggplot(df2_plot, aes(x = BASE, y = LOGCHG, color = group2)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess", se = FALSE, size = 1) +
  labs(
    x = "Baseline C1-INH Functional Activity",
    y = "Outcome",
    color = "Group"
  ) +
  ggtitle(expression("NDIPP with " * a[0] * "~ Beta(2, 2) " * xi * "= 0.5")) +
  theme_minimal() +
  scale_color_manual(values = c("All data" = "black",
                                "Less exchangeable group" = "red",
                                "More exchangeable group" = "blue"))


dipp.plot + ndipp.plot
