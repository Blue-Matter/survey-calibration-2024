
fit_sdmTMB <- function(dat) {

  require(sdmTMB)

  all_years <- seq(min(dat$year), max(dat$year))
  extra_time <- setdiff(all_years, unique(dat$year))
  if (!length(extra_time)) extra_time <- NULL

  mesh <- sdmTMB::make_mesh(
    data = dat,
    xy_cols = c("X", "Y"),
    type = "cutoff",
    cutoff = 20
  )

  fit <- sdmTMB::sdmTMB(
    formula = number ~ poly(log_depth, 2, raw = TRUE) + data.series,
    data = dat,
    mesh = mesh,
    time = "year",
    family = nbinom2(),
    extra_time = extra_time,
    spatial = "on",
    spatiotemporal = "rw",
    offset = log(dat$swept.area),
    anisotropy = TRUE
  )

  return(fit)
}


do_retrospective <- function(yret, dat, sp, cpus = length(yret)) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(yret, function(y, dat) {
    dat_fit <- dplyr::filter(dat, year <= y)
    fit <- fit_sdmTMB(dat_fit)
    saveRDS(fit, file = paste0("validation/models/", sp, "_ret_fit_", y, ".rds"))
    return(fit)
  }, dat = dat)

  return(fit_list)
}



do_retrospective_gap <- function(gap = 0, yret = max(dat$year), dat, sp, year_new_vessel = 1995, cpus = length(gap)) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(gap, function(ygap, dat, year_new_vessel) {
    dat_fit <- dplyr::filter(dat, .data$year < year_new_vessel | .data$year >= year_new_vessel + ygap)
    fit <- fit_sdmTMB(dat_fit)
    saveRDS(fit, file = paste0("validation/models/", sp, "_gap_fit_", ygap, ".rds"))
    return(fit)
  }, dat = dat, year_new_vessel = year_new_vessel)

  return(fit_list)
}


make_index <- function(fit, grid, year_new_vessel = 1995) {

  require(dplyr)

  all_years <- unique(fit$data$year)
  y <- seq(min(all_years), max(all_years))

  newdata <- sdmTMB::replicate_df(grid, time_name = "year", time_values = y)

  p <- predict(fit, newdata = newdata, return_tmb_object = TRUE)
  index <- sdmTMB::get_index(p, bias_correct = TRUE, area = newdata$area_km2)
  index$year_ret <- max(index$year)

  if (is.null(fit$extra_time)) {
    index$impute <- FALSE
  } else {
    index$impute <- index$year %in% fit$extra_time
  }
  index$gap <- sum(index$impute)

  return(index)
}


make_index_ret <- function(fit_ret, grid, cpus = length(fit_ret)) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  on.exit(snowfall::sfStop())

  index_list <- sfLapply(fit_ret, make_index, grid = grid)

  do.call(rbind, index_list)
}

get_vessel_effect <- function(fit_list) {

  vessel_effect <- lapply(1:length(fit_list), function(i) {
    summary(fit_list[[i]]$sd_report, select = "fixed") %>%
      as.data.frame() %>%
      mutate(var = rownames(.)) %>%
      mutate(sequence = i)
  }) %>%
    bind_rows() %>%
    filter(var == "b_j.3") %>%
    mutate(lwr = Estimate - 1.96 * `Std. Error`, upr = Estimate + 1.96 * `Std. Error`)

  return(vessel_effect)
}

