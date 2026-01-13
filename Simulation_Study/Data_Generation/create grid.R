save.dir <- '/work/users/c/l/clairez1/Paper1sims3'

grid_add <- expand.grid(
  n0 = c(20, 30, 40, 50),
  q  = c(-1, -0.5, 0.5, 1),
  prob.unexch = c(0.25, 0.5, 0.75, 1)
)
grid_extra <- expand.grid(
  n0 = c(20, 30, 40, 50),
  q  = 0,
  prob.unexch = 0
)

grid <- rbind(grid_extra, grid_add)

# Seed for generating historical data
grid$hist.seed <- sample(seq_len(1000 * nrow(grid)), nrow(grid), replace = FALSE)

# Create 20 rows of each scenario to divide up sims
grid <- grid[rep(1:nrow(grid), each = 5),]

# Seed for each divided-up scenario
grid$seed <- sample(seq_len(1000 * nrow(grid)), nrow(grid), replace = FALSE)

saveRDS(grid, file = file.path(save.dir, "grid_additive.rds"))
