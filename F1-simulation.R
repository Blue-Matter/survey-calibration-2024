
# Run simulation
library(tidyverse)
library(snowfall)
library(rlecuyer)

sfInit(parallel = TRUE, cpus = 10)

# Simulation design ----
sim_design <- expand.grid(
  sigma_O = 3,
  sigma_E = 1,
  kappa = c(25, 50), # Low kappa = high range = low correlation at short distances
  nyears = 10,
  log_rho = c(0, 0.2, 1),
  st = "rw",
  drift = 0,
  cluster = FALSE,
  phi = c(0.2, 4),
  ar_rho = c(1, 0.5)
) %>%
  mutate(range = sqrt(8)/kappa,
         OM = paste("OM", 1:nrow(.)))
readr::write_csv(sim_design, "simulation/sim_design.csv")

source("99-functions-simulation.R")

# Operating model ----
cell_width <- 0.02
nd_y <- sfLapply(1:nrow(sim_design), sim_OM, sim_design = sim_design, seed = 2, cell_width = cell_width)
saveRDS(nd_y, file = "simulation/operating_model.rds")
nd_y <- readRDS(file = "simulation/operating_model.rds")

for (i in 1:nrow(sim_design)) {

  message("OM ", i)
  nsim <- 100
  om <- dplyr::select(nd_y[[i]], X, Y, year, mu, cell)

  # Simulate data ----
  sfExport(list = c("simulate_data", "sim_design", "i", "om"))
  sfClusterSetupRNGstream(seed = 1)
  dat_list <- sfLapply(1:nsim, function(...) {
    simulate_data(om, phi = sim_design$phi[i], cluster = sim_design$cluster[i])
  })
  saveRDS(dat_list, paste0("simulation/dat_list_", i, ".rds"))
  dat_list <- readRDS(paste0("simulation/dat_list_", i, ".rds"))

  # Run retrospective ----
  sfExport(list = c("fit_model_index", "om", "i", "cell_width"))
  tictoc::tic()
  index_ret <- do_ret_sim(dat_list)
  tictoc::toc()
  saveRDS(index_ret, paste0("simulation/results_ret_", i, ".rds"))

  # Retrospective w/o gear effect
  tictoc::tic()
  index_ret_norho <- do_ret_sim(dat_list, FALSE)
  tictoc::toc()
  saveRDS(index_ret_norho, paste0("simulation/results_ret_", i, "_norho.rds"))

  # Run gap ----
  tictoc::tic()
  index_gap <- do_gap_sim(dat_list)
  tictoc::toc()
  saveRDS(index_gap, paste0("simulation/results_gap_", i, ".rds"))

  # Run gap w/o gear effect ----
  tictoc::tic()
  index_gap_norho <- do_gap_sim(dat_list, FALSE)
  tictoc::toc()
  saveRDS(index_gap_norho, paste0("simulation/results_gap_", i, "_norho.rds"))

}

sfStop()

