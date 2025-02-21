

expand_number_len <- function(dat_lf, dat_rho, length_bin) {
  lf_cols <- names(dat_lf)[grep("lf", names(dat_lf))]

  # Numbers at length are not integers because they were standardized to 15 minute tows
  # Simply round to nearest whole number
  d <- dat_lf %>%
    select(year, vessel, trip, set, NAFOdiv, strat, strat.area, long.start, lat.start, set.depth.mean, data.series, dist.towed, number, nmeasured,
           starts_with("lf")) %>%
    reshape2::melt(measure.vars = lf_cols) %>%
    mutate(
      nmeasured = sum(value), .by = c(year, trip, set),
    ) %>%
    mutate(
      value.expand = ifelse(number > 0, number/nmeasured * value, 0),
      len = strsplit(as.character(.data$variable), "lf") %>% sapply(getElement, 2) %>%
        as.numeric() %>% round(),
    ) %>%
    mutate(nmeasured.expand = sum(value.expand), .by = c(year, trip, set))

  if (missing(dat_rho)) {
    d$rho <- 1
    d$value.expand.rho <- d$rho * d$value.expand
  } else {
    d <- left_join(
      d,
      select(dat_rho, bin, rho),
      by = c("len" = "bin")
    ) %>%
      mutate(value.expand.rho = rho * value.expand)
  }

  if (missing(length_bin)) {
    d$bin <- d$len
  } else {
    d <- d %>%
      mutate(bin = findInterval(len, length_bin))
  }
  dout <- d %>%
    summarise(value.expand = sum(value.expand),
              value.expand.rho = sum(value.expand.rho),
              value = sum(value),
              .by = c(year, vessel, trip, set, NAFOdiv, strat, strat.area, long.start, lat.start,
                      set.depth.mean, data.series, dist.towed, number, nmeasured, nmeasured.expand, bin)) %>%
    mutate(rho.eff = ifelse(value.expand > 0, value.expand.rho/value.expand, 1))


  if (missing(length_bin)) {
    dout$length <- dout$bin
  } else {
    dout <- dout %>%
      mutate(length = length_bin[.data$bin])
  }

  return(dout)
}

make_design_index <- function(dat, simsurvey = FALSE) {

  set_df <- dat %>%
    filter(!is.na(strat.area)) %>%
    summarise(
      number = sum(number),
      .by = c(year, strat, strat.area, vessel, trip, set, swept.area)
    )

  if (simsurvey) {


    require(data.table)

    # Copied from SimSurvey::strat_means
    strat_means <- function(data = NULL, metric = NULL, strat_groups = NULL, survey_groups = NULL, confidence = 95) {
      require(data.table)
      Nh <- strat_area <- tow_area <- Wh <- total <- sumYh <- nh <- gh <- meanYh <- varYh <- meanYst_lcl <- meanYst <- varYst <- df <- meanYst_ucl <- sumYst <- N <- sumYst_lcl <- sumYst_ucl <- NULL
      lc <- (100 - confidence)/200
      uc <- (100 - confidence)/200 + (confidence/100)
      d <- copy(data)
      d <- d[, c(strat_groups, metric), with = FALSE]
      setnames(d, names(d), c(strat_groups, "metric"))
      setkeyv(d, strat_groups)
      strat_tab <- d[, list(sumYh = sum(metric), meanYh = mean(metric),
                            varYh = stats::var(metric), nh = .N), by = strat_groups]
      strat_tab$varYh <- ifelse(strat_tab$nh == 1, 0, strat_tab$varYh)
      strat_tab[, `:=`(Nh, strat_area/tow_area)]
      strat_tab[, `:=`(Wh, Nh/sum(Nh)), by = survey_groups]
      strat_tab[, `:=`(total, Nh * sumYh/nh)]
      strat_tab[, `:=`(gh, Nh * (Nh - nh)/nh)]
      survey_tab <- strat_tab[, list(n = sum(nh), N = sum(Nh),
                                     meanYst = sum(Wh * meanYh), varYst = (1/((sum(Nh))^2)) *
                                       sum(gh * varYh), df = ((sum(gh * varYh))^2)/(sum((gh^2 *
                                                                                           varYh^2)/(nh - 1)))), by = survey_groups]

      # Use qnorm() instead of qt() because df is NA (1 set in many strata)
      survey_tab[, `:=`(meanYst_lcl, (meanYst - (sqrt(varYst)) *
                                        abs(stats::qnorm(lc))))]
      survey_tab[, `:=`(meanYst_ucl, (meanYst + (sqrt(varYst)) *
                                        abs(stats::qnorm(lc))))]
      survey_tab[, `:=`(sumYst, N * meanYst)]
      survey_tab[, `:=`(sumYst_lcl, (sumYst - abs(stats::qnorm(lc)) * N * sqrt(varYst)))]
      survey_tab[, `:=`(sumYst_ucl, (sumYst + abs(stats::qnorm(lc)) * N * sqrt(varYst)))]
      survey_tab[sapply(survey_tab, is.nan)] <- NA

      survey_tab <- survey_tab[, c(survey_groups, "n", "N", "df",
                                   "varYst", "meanYst", "meanYst_lcl", "meanYst_ucl", "sumYst",
                                   "sumYst_lcl", "sumYst_ucl"), with = FALSE]
      survey_tab$varYst <- sqrt(survey_tab$varYst)
      setnames(survey_tab, names(survey_tab), c(survey_groups,
                                                "sets", "sampling_units", "df", "sd", "mean", "mean_lcl",
                                                "mean_ucl", "total", "total_lcl", "total_ucl"))
      survey_tab
    }

    do_simsurvey <- set_df %>%
      rename(strat_area = strat.area) %>%
      mutate(tow_area = swept.area) %>%
      #mutate(tow_area = swept.area/1e3) %>%
      #select(year, set, strat, strat_area, bin, number.calibrate, tow_area) %>%
      data.table::as.data.table() %>%
      strat_means(
        metric = "number",
        strat_groups = c("year", "strat", "strat_area", "tow_area"),
        survey_groups = "year"
      )

    index_design <- do_simsurvey %>%
      as.data.frame() %>%
      rename(se = sd, est = total, lwr = total_lcl, upr = total_ucl) %>%
      mutate(cv = se/est)

  } else {

    index_design <- set_df %>%
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

  }
  return(index_design)
}

fit_sdmTMB <- function(dat, formula = number ~ poly(log_depth, 2, raw = TRUE) + data.series,
                       offset = log(dat$swept.area), family = nbinom2(),
                       spatiotemporal = "rw",
                       ...) {

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
    family = family,
    extra_time = extra_time,
    spatiotemporal = spatiotemporal,
    offset = offset,
    anisotropy = TRUE,
    ...
  )

  return(fit)
}


do_retrospective <- function(yret, dat, sp, cpus = length(yret), ...) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(yret, function(y, dat, ...) {
    dat_fit <- dplyr::filter(dat, year <= y)
    fit <- fit_sdmTMB(dat_fit, ...)
    saveRDS(fit, file = paste0("validation/models/", sp, "_ret_fit_", y, ".rds"))
    return(fit)
  }, dat = dat, ...)

  return(fit_list)
}

do_null_test <- function(ydummy, dat, sp, cpus = length(ydummy), ...) {
  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(ydummy, function(y, dat, ...) {

    if (!is.na(y)) {
      dat_fit <- dplyr::mutate(dat, dummy = year >= y)
      fit <- fit_sdmTMB(
        dat_fit,
        number ~ poly(log_depth, 2, raw = TRUE) + dummy,
        ...
      )
    } else {
      dat_fit <- dat
      fit <- fit_sdmTMB(
        dat_fit,
        number ~ poly(log_depth, 2, raw = TRUE),
        ...
      )
    }
    saveRDS(fit, file = paste0("validation/models/", sp, "_null_fit_", y, ".rds"))
    return(fit)
  }, dat = dat, ...)

  return(fit_list)
}

do_retrospective_gap <- function(gap = 0, yret = max(dat$year), dat, sp, year_new_vessel = 1995, cpus = length(gap), ...) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(gap, function(ygap, dat, year_new_vessel) {
    dat_fit <- dplyr::filter(dat, .data$year < year_new_vessel | .data$year >= year_new_vessel + ygap)
    fit <- fit_sdmTMB(dat_fit, ...)
    saveRDS(fit, file = paste0("validation/models/", sp, "_gap_fit_", ygap, ".rds"))
    return(fit)
  }, dat = dat, year_new_vessel = year_new_vessel, ...)

  return(fit_list)
}


make_index <- function(fit, grid, year_new_vessel = 1995, pred_only = FALSE) {

  require(dplyr)

  all_years <- unique(fit$data$year)
  y <- seq(min(all_years), max(all_years))

  newdata <- sdmTMB::replicate_df(grid, time_name = "year", time_values = y)

  p <- predict(fit, newdata = newdata, return_tmb_object = TRUE)
  if (pred_only) return(p)
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


plot_rho_size <- function(m, bin = seq(1, 200, 1)) {

  est <- data.frame(bin = bin, offset_ratio = 1, id = factor(NA))

  if (inherits(m, "glmmTMB")) {

    pred <- predict(m, est, re.form = NULL, se.fit = TRUE, cov.fit = TRUE)
    est$log_rho <- pred$fit
    est$se <- pred$se.fit

  } else if (inherits(m, "sdmTMB")) {

    est <- predict(m, est, re_form = NA, re_form_iid = ~ 0, se_fit = TRUE) %>%
      rename(log_rho = est, se = est_se)

  } else {
    stop("No predict function written for object m")
  }

  g <- est %>%
    mutate(rho = exp(log_rho)) %>%
    ggplot(aes(bin)) +
    geom_line(aes(y = log_rho)) +
    geom_line(aes(y = log_rho + 1.96 * se), linetype = 2) +
    geom_line(aes(y = log_rho - 1.96 * se), linetype = 2) +
    labs(x = "Length", y = expression(log(rho)))
  g
}

do_fit_size <- function(dat, formula = number.int ~ poly(log_depth, 2, raw = TRUE),
                        length_bin = unique(dat$length), sp, cpus = length(length_bin),
                        type = c("est", "calibrated", "ret", "gap", "norho"), family = nbinom2(), ...) {

  type <- match.arg(type)

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  sfLibrary(sdmTMB)
  sfExport(list = "fit_sdmTMB")
  on.exit(snowfall::sfStop())

  fit_list <- sfLapply(length_bin, function(i, dat, sp, formula, type, family, ...) {
    dat_fit <- dplyr::filter(dat, .data$length == i)

    fit <- fit_sdmTMB(
      dat_fit,
      formula = formula,
      offset = dat_fit$offset,
      family = family,
      ...
    )
    saveRDS(fit, paste0("validation_size/models/", sp, "_fit_", type, "_", i, ".rds"))
    return(fit)
  }, dat = dat, sp = sp, formula = formula, type = type, family = family, ...)

  return(fit_list)
}

make_index_size <- function(fit_list, length_bin, grid, cpus = length(length_bin)) {

  require(snowfall)

  snowfall::sfInit(TRUE, cpus = cpus)
  on.exit(snowfall::sfStop())

  index_list <- sfLapply(fit_list, make_index, grid = grid)

  ind <- Map(function(x, y) dplyr::mutate(x, length = y), x = index_list, y = length_bin) %>%
    dplyr::bind_rows()

  return(ind)
}
