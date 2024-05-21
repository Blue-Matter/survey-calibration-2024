

expand_number_len <- function(dat_lf, dat_rho, length_bin) {
  lf_cols <- names(dat_lf)[grep("lf", names(dat_lf))]

  d <- dat_lf %>%
    select(year, vessel, trip, set, NAFOdiv, strat, strat.area, long.start, lat.start, set.depth.mean, data.series, dist.towed, number, nmeasured,
           starts_with("lf")) %>%
    reshape2::melt(measure.vars = lf_cols) %>%
    mutate(
      value.expand = ifelse(number > 0, number/nmeasured * value, 0),
      len = strsplit(as.character(.data$variable), "lf") %>% sapply(getElement, 2) %>% as.numeric() %>% round()
    )

  if (missing(dat_rho)) {
    d$rho <- 1
    d$value.calibrate <- d$value.expand * d$rho
  } else {
    d <- left_join(
      d,
      select(dat_rho, bin, rho),
      by = c("len" = "bin")
    ) %>%
      mutate(value.calibrate = value.expand * rho)
  }

  if (missing(length_bin)) {
    d$bin <- d$len
  } else {
    d <- d %>%
      mutate(bin = findInterval(len, length_bin))
  }

  dout <- d %>%
    summarise(number.size = sum(value.expand),
              number.calibrate = sum(value.calibrate),
              .by = c(year, vessel, trip, set, NAFOdiv, strat, strat.area, long.start, lat.start,
                      set.depth.mean, data.series, dist.towed, number, nmeasured, bin)) %>%
    mutate(rho.eff = ifelse(number.size > 0, number.calibrate/number.size, NA))

  if (missing(length_bin)) {
    dout$length <- dout$bin
  } else {
    dout <- dout %>%
      mutate(length = length_bin[.data$bin])
  }

  return(dout)
}

make_design_index <- function(dat) {

  index_design <- dat %>%
    filter(!is.na(strat.area)) %>%
    summarise(
      index_strat = mean(number),
      var_strat = var(number),
      n = n(),
      nsamp = unique(strat.area)/mean(swept.area),
      .by = c(year, strat, strat.area)
    ) %>%
    mutate(var_strat = ifelse(is.na(var_strat), 0, var_strat)) %>%
    mutate(area_total = sum(strat.area), .by = year) %>%
    summarise(
      est = sum(index_strat * nsamp),
      var = sum(nsamp * (nsamp - n)/n * var_strat)/sum(nsamp)/sum(nsamp), # See SimSurvey appendix
      .by = year
    ) %>%
    mutate(
      se = sqrt(var),
      cv = se/est,
      lwr = est - 1.96 * se,
      upr = est + 1.96 * se
    )
  return(index_design)
}

fit_sdmTMB <- function(dat, formula = number ~ poly(log_depth, 2, raw = TRUE) + data.series,
                       offset = log(dat$swept.area)) {

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
    formula = formula,
    data = dat,
    mesh = mesh,
    time = "year",
    family = nbinom2(),
    extra_time = extra_time,
    spatial = "on",
    spatiotemporal = "rw",
    offset = offset,
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

