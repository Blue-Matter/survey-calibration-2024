
# Run simulation
library(tidyverse)
library(snowfall)
sfInit(parallel = TRUE, cpus = 10)

# Simulation design ----
sim_design <- expand.grid(
  sigma_O = 2,
  sigma_E = 0.3,
  kappa = c(5, 25), # Low kappa = high range = low correlation at short distances
  nyears = 10,
  log_rho = c(0, 0.2, 1),
  st = "rw",
  drift = c(0, 0.1),
  phi = c(0.2, 4)
) %>%
  mutate(range = sqrt(8)/kappa,
         OM = paste("OM", 1:nrow(.)))
#readr::write_csv(sim_design, "simulation/sim_design.csv")

source("99-functions-simulation.R")

# Operating model ----
cell_width <- 0.02
nd_y <- sfLapply(1:nrow(sim_design), sim_OM, sim_design = sim_design, seed = 2, cell_width = cell_width)
saveRDS(nd_y, file = "simulation/operating_model.rds")
#nd_y <- readRDS(file = "simulation/operating_model.rds")

for (i in 1:nrow(sim_design)) {

  message("OM ", i)
  nsim <- 100
  om <- dplyr::select(nd_y[[i]], X, Y, year, pred, cell)

  # Simulate data ----
  sfExport(list = c("simulate_data", "sim_design", "i", "om"))
  sfClusterSetupRNGstream(seed = 1)
  dat_list <- sfLapply(1:nsim, function(...) {
    simulate_data(om, phi = sim_design$phi[i], nyears = sim_design$nyears[i])
  })
  saveRDS(dat_list, paste0("simulation/dat_list_", i, ".rds"))
  #dat_list <- readRDS(paste0("simulation/dat_list_", i, ".rds"))

  # Fit models ---
  sfExport(list = c("fit_model_index", "om", "i", "cell_width"))
  index <- sfLapply(dat_list, function(j) {
    fit_model_index(j, grid_pred = dplyr::select(om, X, Y, year), cell_area = cell_width^2)
  })
  saveRDS(index, paste0("simulation/results_", i, ".rds"))
}

sfStop()

