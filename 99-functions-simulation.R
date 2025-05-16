
sim_OM <- function(i, sim_design, seed = 10392,
                   y_new_vessel = ceiling(0.5 * nyears),
                   cell_width = 0.02) {
  set.seed(seed)
  require(tidyverse)
  require(SimSurvey)

  sigma_O <- sim_design$sigma_O[i]
  sigma_E <- sim_design$sigma_E[i]
  range <- sim_design$range[i]
  nyears <- sim_design$nyears[i]
  log_rho <- sim_design$log_rho[i]
  st <- match.arg(as.character(sim_design$st[i]), choices = c("rw", "iid"))
  drift <- sim_design$drift[i]

  # Operating model ----
  nd <- expand.grid(X = seq(0, 1, cell_width), Y = seq(0, 1, cell_width)) %>%
    mutate(cell = 1:nrow(.))

  omega_matern_fn <- SimSurvey::sim_ays_covar(
    sd = sigma_O,
    lambda = 1, # nu
    range = range,
    phi_age = 0,
    phi_year = 0
  )
  epsilon_matern_fn <- SimSurvey::sim_ays_covar(
    sd = sigma_E,
    lambda = 1, # nu
    range = range,
    phi_age = 0,
    phi_year = 0
  )

  omega_s <- omega_matern_fn(nd[, c("X", "Y")], ages = 1, years = 1, cell = seq(1, nrow(nd)))
  delta_st <- epsilon_matern_fn(nd[, c("X", "Y")], ages = 1, years = seq(1, nyears), cell = seq(1, nrow(nd)))

  if (st == "rw") {
    epsilon_st <- array(NA, dim(delta_st), dimnames = list(age = 1, year = 1:nyears, cell = 1:nrow(nd)))
    epsilon_st[, 1, ] <- delta_st[, 1, ]
    for (y in 2:nyears) epsilon_st[, y, ] <- drift + epsilon_st[, y-1, ] + delta_st[, y, ]
  } else if (st == "iid") {
    epsilon_st <- delta_st
  }

  nd$omega_s <- omega_s[1, 1, ]
  nd_y <- right_join(nd, reshape2::melt(epsilon_st[1, , ])) %>%
    rename(epsilon_st = value) %>%
    mutate(log_rho = ifelse(year > y_new_vessel, .env$log_rho, 0)) %>%
    mutate(log_pred = omega_s + epsilon_st + log_rho,
           pred = exp(log_pred),
           log_pred_pop = omega_s + epsilon_st,
           pred_pop = exp(log_pred_pop))

  return(nd_y)
}


simulate_data <- function(om, phi = 1, nyears = 10,
                          y_new_vessel = ceiling(0.5 * nyears),
                          cluster = FALSE, seed, ...) {
  if (!missing(seed)) set.seed(seed)
  require(tidyverse)
  require(MASS)

  n_annual <- 100

  dat <- lapply(1:nyears, function(y) {
    if (cluster) {
      cell_y1 <- sample(unique(om$cell), ifelse(cluster, 0.5, 1) * n_annual)
      cell_y2 <- pmin(cell_y1 + 1, max(om$cell))
      cell_y <- c(cell_y1, cell_y2)
    } else {
      cell_y <- sample(unique(om$cell), ifelse(cluster, 0.5, 1) * n_annual)
    }
    dplyr::filter(om, cell %in% cell_y, year == y)
  }) %>%
    bind_rows() %>%
    mutate(new_vessel = year > y_new_vessel)

  dat$obs <- MASS::rnegbin(nrow(dat), mu = dat$pred, theta = phi)

  return(dat)
}

fit_model_index <- function(dat, grid_pred, cell_area = 1, get_se = FALSE, ...) {

  require(tidyverse)
  require(sdmTMB)

  y_data <- unique(dat$year)
  y_all <- seq(min(dat$year), max(dat$year))
  y_extra <- setdiff(y_all, y_data)
  extra_time <- if (length(y_extra)) y_extra else NULL

  mesh <- sdmTMB::make_mesh(dat, c("X", "Y"), n_knots = 100)
  fit <- sdmTMB(
    obs ~ 1 + new_vessel,
    dat,
    mesh,
    time = "year",
    extra_time = extra_time,
    spatiotemporal = "rw",
    family = nbinom2(),
    ...
  )

  newdata <- grid_pred %>%
    mutate(new_vessel = FALSE)

  pred <- predict(
    fit,
    newdata = newdata,
    return_tmb_object = get_se
  )

  if (get_se) {
    index <- get_index(pred, bias_correct = TRUE, area = cell_area)
  } else {
    index <- pred %>%
      summarise(est = sum(exp(est) * cell_area), .by = year)
  }
  sd_report <- summary(fit$sd_report, select = c("fixed", "report"))

  return(list(index = index, sd_report = sd_report))

}
