
# Run simulation
library(tidyverse)
library(snowfall)
library(rlecuyer)

sfInit(parallel = TRUE, cpus = 15)

# Simulation design ----
relative_range <- range_dist <- c(0.4, 0.55, 0.8) # Relative range equal to range because sqrt(Area) = 1
kappa <- sqrt(8)/range_dist

sim_design <- expand.grid(
  sigma_O = 3,
  sigma_E = 1,
  kappa = c(5, 25), # Low kappa = high range = low correlation at short distances
  nyears = 10,
  log_rho = c(0, 0.2, 1),
  drift = 0,
  cluster = FALSE,
  phi = c(0.2, 4),
  ar_rho = c(1, 0.5)
) %>%
  mutate(range = sqrt(8)/kappa,
         OM = paste("OM", 1:nrow(.)))
readr::write_csv(sim_design, file = "simulation_simfield/sim_design.csv")
sim_design <- read_csv("simulation_simfield/sim_design.csv")

source("99-functions-simulation.R")

# Operating model ----
nsim <- 100

cell_width <- 0.02
scenario <- 1:nrow(sim_design)
intercept <- log(50)

if (FALSE) {

  sfExport(list = "sim_OM")
  for (i in scenario) {
    om_i <- sfClusterApplyLB(1:nsim, function(x, i, sim_design, cell_width, intercept) {
      sim_OM(i, sim_design = sim_design, seed = x * 365, cell_width = cell_width, intercept = intercept)
    }, i = i, sim_design = sim_design, cell_width = cell_width, intercept = intercept)
    message("OM ", i)
    saveRDS(om_i, file = paste0("simulation_simfield/operating_model_", i, ".rds"))
  }

  # Filter data ----
  sfExport(list = c("simulate_data"))
  for (i in scenario) {
    message("OM ", i)
    om_i <- readRDS(file = paste0("simulation_simfield/operating_model_", i, ".rds"))

    sfClusterSetupRNGstream(seed = 1)
    dat_list <- sfClusterApplyLB(om_i, function(om) {
      list(
        nd_y = simulate_data(om[["nd_y"]], cluster = FALSE, resample = FALSE),
        mesh = om[["mesh"]]
      )
    })
    saveRDS(dat_list, paste0("simulation_simfield/dat_list_", i, ".rds"))
  }
}

# Fit models ----
sfExport(list = "fit_model_index")
for (i in seq(1, length(scenario))) {
  message("OM ", i)
  dat_list <- readRDS(paste0("simulation_simfield/dat_list_", i, ".rds"))

  # Starting values of estimation model from operating model (not used)
  start_rho <- list(
    b_j = c(intercept, sim_design$log_rho[i]),
    ln_kappa = matrix(log(sim_design$kappa[i]), 2, 1),
    ln_tau_O = log(sim_design$sigma_O[i]),
    ln_tau_E = log(sim_design$sigma_E[i]),
    ln_phi = log(sim_design$phi[i])
  )
  start_norho <- start_rho
  start_norho$b_j <- intercept

  # Run retrospective ----
  message("Part 1")
  tictoc::tic()
  index_ret <- do_ret_sim(dat_list, cell_width = cell_width)
  tictoc::toc()
  saveRDS(index_ret, paste0("simulation_simfield/results_ret_", i, ".rds"))

  # Retrospective w/o gear effect
  message("Part 2")
  tictoc::tic()
  index_ret_norho <- do_ret_sim(dat_list, FALSE, cell_width = cell_width)
  tictoc::toc()
  saveRDS(index_ret_norho, paste0("simulation_simfield/results_ret_", i, "_norho.rds"))

  # Run gap ----
  message("Part 3")
  tictoc::tic()
  index_gap <- do_gap_sim(dat_list, cell_width = cell_width)
  tictoc::toc()
  saveRDS(index_gap, paste0("simulation_simfield/results_gap_", i, ".rds"))

  # Run gap w/o gear effect ----
  message("Part 4")
  tictoc::tic()
  index_gap_norho <- do_gap_sim(dat_list, FALSE, cell_width = cell_width)
  tictoc::toc()
  saveRDS(index_gap_norho, paste0("simulation_simfield/results_gap_", i, "_norho.rds"))

}

sfStop()

