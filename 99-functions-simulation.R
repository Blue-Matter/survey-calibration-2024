
sim_OM <- function(i, sim_design, seed = 10392,
                   y_new_vessel = ceiling(0.5 * nyears),
                   cell_width = 0.02,
                   intercept = 0) {
  set.seed(seed)
  require(tidyverse)
  require(sdmTMB)

  sigma_O <- sim_design$sigma_O[i]
  sigma_E <- sim_design$sigma_E[i]
  range <- sim_design$range[i]
  nyears <- sim_design$nyears[i]
  log_rho <- sim_design$log_rho[i]

  drift <- sim_design$drift[i]
  stopifnot(drift == 0) # for (y in 2:nyears) epsilon_st[, y, ] <- drift + epsilon_st[, y-1, ] + delta_st[, y, ]
  phi <- sim_design$phi[i]

  # Operating model ----
  nd <- expand.grid(X = seq(0, 1, cell_width), Y = seq(0, 1, cell_width)) %>%
    mutate(cell = 1:nrow(.))

  nd2 <- lapply(1:nyears, function(y) {
    mutate(nd, year = y)
  }) %>%
    bind_rows() %>%
    mutate(log_rho = ifelse(year > y_new_vessel, .env$log_rho, 0))

  mesh <- sdmTMB::make_mesh(nd2, xy_cols = c("X", "Y"), type = "cutoff", cutoff = 0.04)

  nd_y <- sdmTMB_simulate(
    ~ 1 + log_rho + cell,
    data = nd2,
    mesh = mesh,
    family = nbinom2(),
    time = "year",
    range = range,
    rho = sim_design$ar_rho[i],
    B = c(intercept, 1, 0),
    sigma_O = sigma_O,
    sigma_E = sigma_E,
    phi = phi,
    seed = seed
  ) %>%
    mutate(eta_pop = intercept +omega_s + epsilon_st,
           mu_pop = exp(eta_pop))

  return(list(nd_y = nd_y, mesh = mesh))
}


simulate_data <- function(om, phi = 1, #nyears = 10,
                          y_new_vessel = ceiling(0.5 * nyears),
                          cluster = FALSE, seed, n_annual = 100,
                          resample = TRUE,
                          ...) {
  if (!missing(seed)) set.seed(seed)
  require(tidyverse)

  nyears <- max(om$year)

  dat <- lapply(1:nyears, function(y) {
    if (cluster) {
      cell_y1 <- sample(unique(om$cell), 0.5 * n_annual)
      cell_y2 <- pmin(cell_y1 + 1, max(om$cell))
      cell_y <- c(cell_y1, cell_y2)
    } else {
      cell_y <- sample(unique(om$cell), n_annual)
    }
    dplyr::filter(om, cell %in% cell_y, year == y)
  }) %>%
    bind_rows() %>%
    mutate(new_vessel = year > y_new_vessel)

  if (resample) {
    require(MASS)
    dat$observed <- MASS::rnegbin(nrow(dat), mu = dat$mu, theta = phi)
  }

  return(dat)
}

fit_model_index <- function(dat, mesh_om = NULL,
                            formula = observed ~ 1 + new_vessel,
                            grid_pred, cell_area = 1, start = list(), get_se = FALSE, ...) {

  require(tidyverse)
  require(sdmTMB)

  y_data <- unique(dat$year)
  y_all <- seq(min(dat$year), max(dat$year))
  y_extra <- setdiff(y_all, y_data)
  extra_time <- if (length(y_extra)) y_extra else NULL

  if (is.null(mesh_om)) {
    mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), type = "cutoff", cutoff = 0.04, mesh = NULL)
  } else {
    mesh <- sdmTMB::make_mesh(dat, xy_cols = c("X", "Y"), mesh = mesh_om$mesh)
  }

  fit <- sdmTMB::sdmTMB(
    formula,
    dat,
    mesh,
    time = "year",
    extra_time = extra_time,
    spatiotemporal = "rw",
    family = nbinom2(),
    control = sdmTMB::sdmTMBcontrol(start = start)
    #...
  )

  if (packageVersion('sdmTMB') >= "1.2.0") {
    index <- sdmTMB::get_index(fit, newdata = grid_pred, bias_correct = TRUE, area = cell_area)
  } else {
    pred <- predict(
      fit,
      newdata = grid_pred,
      return_tmb_object = get_se
    )

    if (get_se) {
      index <- get_index(pred, bias_correct = TRUE, area = cell_area)
    } else {
      index <- pred %>%
        summarise(est = sum(exp(est) * cell_area), .by = year)
    }
  }

  sd_report <- summary(fit$sd_report, select = c("fixed", "report"))

  return(list(index = index, sd_report = sd_report))

}

do_ret_sim <- function(dat_list, est_gear = TRUE, cell_width, start = list(), use_om_mesh = TRUE, do_parallel = TRUE) {

  if (do_parallel) {
    fn <- snowfall::sfClusterApplyLB
  } else {
    fn <- base::lapply
  }

  index_ret <- fn(dat_list, function(dat_i, est_gear, cell_width, start) {
    require(dplyr)

    if (est_gear) {
      formula <- observed ~ 1 + new_vessel
    } else {
      formula <- observed ~ 1
    }

    yret <- 6:10
    dat <- dat_i[["nd_y"]]
    mesh_om <- dat_i[["mesh"]]

    out <- lapply(yret, function(y) {
      dat_y <- dplyr::filter(dat, year <= y)
      grid_pred <- lapply(seq_len(y), function(yy) {
        dplyr::distinct(as.data.frame(mesh_om$loc_xy)) %>% mutate(year = yy)
      }) %>%
        bind_rows() %>%
        mutate(new_vessel = FALSE)
      try(
        fit_model_index(
          dat_y, mesh_om = NULL,
          formula, grid_pred = grid_pred, cell_area = cell_width^2,
          start = start, get_se = TRUE
        ),
        silent = TRUE
      )
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
  }, est_gear = est_gear, cell_width = cell_width, start = start)

  return(index_ret)
}

do_gap_sim <- function(dat_list, est_gear = TRUE, cell_width, start = list(), do_parallel = TRUE) {
  if (do_parallel) {
    fn <- snowfall::sfClusterApplyLB
  } else {
    fn <- base::lapply
  }

  index_gap <- fn(dat_list,function(dat_i, est_gear, cell_width, start) {
    require(dplyr)

    if (est_gear) {
      formula <- observed ~ 1 + new_vessel
    } else {
      formula <- observed ~ 1
    }

    ygap <- 6:9
    dat <- dat_i[["nd_y"]]
    mesh_om <- dat_i[["mesh"]]

    out <- lapply(1:length(ygap), function(y) {
      ymissing <- seq(ygap[1], ygap[y])
      dat_y <- dplyr::filter(dat, !year %in% ymissing)
      grid_pred <- lapply(seq_len(max(dat$year)), function(yy) {
        dplyr::distinct(as.data.frame(mesh_om$loc_xy)) %>% mutate(year = yy)
      }) %>%
        bind_rows() %>%
        mutate(new_vessel = FALSE)
      try(
        fit_model_index(
          dat_y, mesh_om = NULL,
          formula, grid_pred = grid_pred, cell_area = cell_width^2,
          start = start, get_se = TRUE
        ),
        silent = TRUE
      )
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
  }, est_gear = est_gear, cell_width = cell_width, start = start)

  return(index_gap)
}
