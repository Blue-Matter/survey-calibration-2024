
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
  drift = 0,
  cluster = c(FALSE, TRUE),
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
nd_y <- readRDS(file = "simulation/operating_model.rds")

for (i in 1:nrow(sim_design)) {

  message("OM ", i)
  nsim <- 100
  om <- dplyr::select(nd_y[[i]], X, Y, year, pred, cell)

  # Simulate data ----
  #sfExport(list = c("simulate_data", "sim_design", "i", "om"))
  #sfClusterSetupRNGstream(seed = 1)
  #dat_list <- sfLapply(1:nsim, function(...) {
  #  simulate_data(om, phi = sim_design$phi[i], nyears = sim_design$nyears[i], cluster = sim_design$cluster[i])
  #})
  #saveRDS(dat_list, paste0("simulation/dat_list_", i, ".rds"))
  dat_list <- readRDS(paste0("simulation/dat_list_", i, ".rds"))


  # Run retrospective ----
  sfExport(list = c("fit_model_index", "om", "i", "cell_width"))
  tictoc::tic()
  index_ret <- sfLapply(dat_list, function(dat) {
    require(dplyr)

    yret <- 6:10

    out <- lapply(yret, function(y) {
      dat_y <- dplyr::filter(dat, year <= y)
      grid_pred <- dplyr::select(om, X, Y, year) %>%
        dplyr::filter(year <= y)
      try(fit_model_index(dat_y, grid_pred = grid_pred, cell_area = cell_width^2, get_se = TRUE), silent = TRUE)
    })

    index <- lapply(1:length(yret), function(y) {
      if (!is.character(out[[y]])) {
        dplyr::mutate(out[[y]]$index, yret = yret[y])
      } else {
        data.frame()
      }
    }) %>%
      bind_rows()

    sd_report <- lapply(1:length(yret), function(y) {

      if (!is.character(out[[y]])) {
        out[[y]]$sd_report %>%
          as.data.frame() %>%
          dplyr::mutate(yret = yret[y], par_name = rownames(.))
      } else {
        data.frame()
      }

    }) %>%
      bind_rows()

    return(list(index = index, sd_report = sd_report))
  })
  tictoc::toc()
  saveRDS(index_ret, paste0("simulation/results_ret_", i, ".rds"))

  tictoc::tic()
  index_gap <- sfLapply(dat_list, function(dat) {
    require(dplyr)

    ygap <- 6:9

    out <- lapply(1:length(ygap), function(y) {
      ymissing <- seq(ygap[1], ygap[y])
      dat_y <- dplyr::filter(dat, !year %in% ymissing)
      grid_pred <- dplyr::select(om, X, Y, year)
      try(fit_model_index(dat_y, grid_pred = grid_pred, cell_area = cell_width^2, get_se = TRUE), silent = TRUE)
    })

    index <- lapply(1:length(ygap), function(y) {
      ymissing <- seq(ygap[1], ygap[y])
      if (!is.character(out[[y]])) {
        dplyr::mutate(out[[y]]$index, impute = year %in% ymissing, ngap = length(ymissing))
      } else {
        data.frame()
      }
    }) %>%
      bind_rows()

    sd_report <- lapply(1:length(ygap), function(y) {
      ymissing <- seq(ygap[1], ygap[y])
      ngap <- length(ymissing)
      if (!is.character(out[[y]])) {
        out[[y]]$sd_report %>%
          as.data.frame() %>%
          dplyr::mutate(ngap = ngap, par_name = rownames(.))
      } else {
        data.frame()
      }

    }) %>%
      bind_rows()

    return(list(index = index, sd_report = sd_report))
  })
  tictoc::toc()
  saveRDS(index_gap, paste0("simulation/results_gap_", i, ".rds"))

}

sfStop()

