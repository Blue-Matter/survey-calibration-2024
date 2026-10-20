
library(tidyverse)
#library(NAFOdown)

#ggplot2::theme_set(theme_nafo())
theme_set(theme_bw())

.ggsave <- function(filename, ...) {
  filename <- paste0("validation/figures/cod_", filename)
  ggplot2::ggsave(filename, ...)
}

source("99-functions.R")

sp <- c("cod", "ghl")

for (s in sp) {

  if (s == "cod") {

    .ggsave <- function(filename, ...) {
      filename <- paste0("validation/figures/cod_", filename)
      ggplot2::ggsave(filename, ...)
    }

    dat <- readr::read_csv("validation/Cod_438_EngelCampelenSurveyData.csv") %>%
      mutate(npue = number/swept.area) %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L")) %>%
      filter(season == "fall") %>%
      filter(npue < 1e4) %>%
      mutate(long.start = -1 * long.start) %>%
      mutate(log_depth = log(set.depth.mean)) %>%
      sdmTMB::add_utm_columns(ll_names = c("long.start", "lat.start")) %>%
      mutate(doorspread = 1e3 * swept.area/dist.towed) %>%
      mutate(data.series = factor(data.series, levels = c("Engel", "Campelen"))) %>%
      mutate(calibration = ifelse(data.series == "Engel", exp(2.072), 1))

  } else if (s == "ghl") {

    .ggsave <- function(filename, ...) {
      filename <- paste0("validation/figures/ghl_", filename)
      ggplot2::ggsave(filename, ...)
    }

    dat <- readr::read_csv("validation/GHL_892_EngelCampelenSurveyData.csv") %>%
      mutate(npue = number/swept.area) %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L")) %>%
      filter(season == "fall") %>%
      filter(npue < 5e4) %>%
      mutate(long.start = -1 * long.start) %>%
      mutate(log_depth = log(set.depth.mean)) %>%
      sdmTMB::add_utm_columns(ll_names = c("long.start", "lat.start")) %>%
      mutate(data.series = factor(data.series, levels = c("Engel", "Campelen")))%>%
      mutate(calibration = ifelse(data.series == "Engel", exp(0.8913), 1))
  }

  # Proportion of zeros per year
  dat %>%
    summarize(pzero = mean(npue == 0), .by = year)

  # Vessel + gear summary
  g <- dat %>% summarise(n = n(), .by = c(year, data.series, vessel, NAFOdiv)) %>%
    ggplot(aes(year, n, linetype = data.series)) +
    geom_point() +
    geom_line() +
    facet_grid(vars(vessel), vars(NAFOdiv)) +
    expand_limits(y = 0)

  #### Spatial figure of set data ----
  nfld <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
    sf::st_crop(xmin = -59, xmax = -46, ymin = 45, ymax = 56)
  g <- dat %>%
    filter(npue > 0) %>%
    ggplot() +
    geom_sf(data = nfld) +
    coord_sf() +
    geom_point(aes(long.start, lat.start, colour = log(npue)), size = 0.5) +
    scale_colour_viridis_c() +
    #facet_wrap(vars(year)) +
    labs(x = "Longitude", y = "Latitude", colour = "log(NPUE)",
         subtitle = expression("Area: 391422"~km^2)) +
    theme(legend.position = "bottom", panel.spacing = unit(0, "in"))
  .ggsave("npue_map.png", g, height = 5, width = 4)

  g <- ggplot(dat, aes(npue, fill = data.series)) +
    facet_wrap(vars(year)) +
    geom_histogram()

  g <- ggplot(dat, aes(npue, fill = data.series)) +
    #facet_grid(vars(year), vars(NAFOdiv)) +
    geom_histogram()

  g <- ggplot(dat, aes(set.depth.mean, npue, colour = data.series)) +
    geom_point() +
    facet_wrap(vars(year)) +
    #coord_cartesian(xlim = c(0, 800)) +
    geom_vline(xintercept = 1000, linetype = 3) +
    labs(x = "Depth (m)", y = "NPUE", colour = "Gear")

  #### Data set to fit full model ----
  dat_fit <- filter(dat, year >= 1985, year <= 2020)

  grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
    mutate(X = 1e-3 * X, Y = 1e-3 * Y,
           depth = -1 * depth) %>%
    filter(.data$depth > 0) %>%
    filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
    mutate(log_depth = log(depth),
           data.series = "Campelen")


  ### Estimate calibration factor
  fit <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE) + data.series)
  saveRDS(fit, paste0("validation/models/", s, "_fit_full.rds"))
  #fit <- readRDS(paste0("validation/models/", s, "_fit_full.rds"))

  index <- make_index(fit, grid)
  readr::write_csv(index, paste0("validation/index/", s, "_index.csv"))
  #index <- readr::read_csv(paste0("validation/index/", s, "_index.csv"))

  ### Provide calibration factor from CF
  fit_rho <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE),
                        offset = log(dat_fit$swept.area/dat_fit$calibration))
  saveRDS(fit_rho, paste0("validation/models/", s, "_fit_calibration.rds"))
  #fit_rho <- readRDS(paste0("validation/models/", s, "_fit_calibration.rds"))

  index_rho <- make_index(fit_rho, grid)
  readr::write_csv(index_rho, paste0("validation/index/", s, "_index_calibration.csv"))
  #index_rho <- readr::read_csv(paste0("validation/index/", s, "_index_calibration.csv"))

  ### No calibration factor
  fit_norho <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE))
  saveRDS(fit_norho, paste0("validation/models/", s, "_fit_norho.rds"))
  #fit_norho <- readRDS(paste0("validation/models/", s, "_fit_norho.rds"))

  index_norho <- make_index(fit_norho, grid)
  readr::write_csv(index_norho, paste0("validation/index/", s, "_index_norho.csv"))
  #index_norho <- readr::read_csv(paste0("validation/index/", s, "_index_norho.csv"))


  ### Plot and compare indices
  index_compare <- rbind(
    index_rho %>% mutate(Model = "Comparative rho"),
    index %>% mutate(Model = "Estimate rho"),
    index_norho %>% mutate(Model = "No rho")
  )
  g <- ggplot(index_compare, aes(year, est, linetype = Model, shape = Model)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    labs(x = "Year", y = "Campelen Index") +
    scale_shape_manual(values = c(16, 1, 4)) +
    scale_linetype_manual(values = c(1, 2, 3)) +
    theme(legend.position = "bottom")
  .ggsave("index_no_comparative.png", g, height = 4, width = 6)

  g <- ggplot(index_compare, aes(year, log(est))) +
    geom_point(shape = 1) +
    geom_line(linewidth = 0.1) +
    geom_linerange(aes(ymin = log(lwr), max = log(upr))) +
    #expand_limits(y = 0) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    facet_wrap(vars(Model), ncol = 1) +
    labs(x = "Year", y = "Campelen log(Index)") +
    scale_shape_manual(values = c(16, 1, 4)) +
    scale_linetype_manual(values = c(1, 2, 3)) +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(colour = NA, fill = NA))
  .ggsave("index_no_comparative_ci.png", g, height = 7, width = 4)

  #### Do retrospective ----
  yret <- 1995:2000
  fit_ret <- do_retrospective(yret, dat = dat_fit, sp = s)
  #fit_ret <- lapply(yret, function(y) {
  #  readRDS(paste0("validation/models/", s, "_ret_fit_", y, ".rds"))
  #})

  index_ret <- make_index_ret(fit_ret, grid)
  readr::write_csv(index_ret, paste0("validation/index/", s, "_index_ret.csv"))
  #index_ret <- readr::read_csv(paste0("validation/index/", s, "_index_ret.csv"))

  vessel_effect <- get_vessel_effect(fit_ret) %>%
    mutate(year_ret = .env$yret[.data$sequence])

  if (s == "cod") {
    log_rho_comparative <- data.frame(
      x = c(-Inf, Inf),
      y1 = 2.072 + 1.96 * 0.08256,
      y2 = 2.072 - 1.96 * 0.08256
    )
  } else {
    log_rho_comparative <- data.frame(
      x = c(-Inf, Inf),
      y1 = 0.8913 + 1.96 * 0.0616,
      y2 = 0.8913 - 1.96 * 0.0616
    )
  }

  g <- ggplot(vessel_effect, aes(year_ret, Estimate)) +
    geom_ribbon(data = log_rho_comparative, fill = "grey80", aes(x, ymin = y2, ymax = y1), inherit.aes = FALSE) +
    geom_hline(yintercept = log_rho_comparative$y1[1], linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$y2[1], linetype = 3) +
    geom_hline(yintercept = ifelse(s == "cod", 2.072, 0.8913), linetype = 2) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    labs(x = "Terminal year", y = "Campelen gear effect")
  .ggsave("retrospective_vessel_effect.png", g, height = 3, width = 4)

  g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
    geom_point() +
    geom_line() +
    expand_limits(y = 0) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    coord_cartesian(xlim = c(1985, 2000)) +
    labs(x = "Year", y = "Campelen Index", colour = "Terminal year")
  .ggsave("index_retrospective.png", g, height = 3, width = 6)

  g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    expand_limits(y = 0) +
    facet_wrap(vars(year_ret)) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
    guides(colour = "none") +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_retrospective_ci.png", g, height = 3.5, width = 6)

  ### Fit with gap
  fit_gap <- do_retrospective_gap(gap = 1:5, dat = dat_fit, sp = s)
  fit_gap <- lapply(1:5, function(i) {
    readRDS(paste0("validation/models/", s, "_gap_fit_", i, ".rds"))
  })

  index_gap <- make_index_ret(fit_gap, grid)
  readr::write_csv(index_gap, paste0("validation/index/", s, "_index_gap.csv"))
  #index_gap <- readr::read_csv(paste0("validation/index/", s, "_index_gap.csv"))

  g <- ggplot(index_gap, aes(year, est)) +
    geom_point(aes(shape = impute, colour = factor(gap))) +
    geom_line(linetype = 2, aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = gap), aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = gap), aes(colour = factor(gap))) +
    expand_limits(y = 0) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
    guides(shape = "none") +
    labs(x = "Year", y = "Campelen Index", colour = "Survey gap")
  .ggsave("index_gap.png", g, height = 3, width = 6)

  g <- ggplot(index_gap, aes(year, log(est))) +
    geom_point(aes(shape = impute, colour = factor(gap))) +
    geom_line(linetype = 2, aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = gap), aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = gap), aes(colour = factor(gap))) +
    #expand_limits(y = 0) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    coord_cartesian(xlim = c(1985, 2000)) +
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
    guides(shape = "none") +
    labs(x = "Year", y = "Campelen log(Index)", colour = "Survey gap")
  .ggsave("index_gap2.png", g, height = 3, width = 6)

  g <- index_gap %>%
    mutate(gap = paste(gap, "year gap")) %>%
    ggplot(aes(year, est, colour = factor(gap), shape = impute, group = factor(gap))) +
    geom_point() +
    geom_line(linewidth = 0.5) +
    geom_linerange(aes(ymin = lwr, max = upr), linewidth = 0.25) +
    expand_limits(y = 0) +
    facet_wrap(vars(gap)) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
    guides(colour = "none", shape = "none") +
    scale_shape_manual(values = c(16, 1)) +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_gap_ci.png", g, height = 3.5, width = 6)


  vessel_effect <- get_vessel_effect(fit_gap)

  g <- ggplot(vessel_effect, aes(sequence, Estimate)) +
    geom_ribbon(data = log_rho_comparative, fill = "grey80", aes(x, ymin = y2, ymax = y1), inherit.aes = FALSE) +
    geom_hline(yintercept = log_rho_comparative$y1[1], linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$y2[1], linetype = 3) +
    geom_hline(yintercept = ifelse(s == "cod", 2.072, 0.8913), linetype = 2) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    labs(x = "Survey gap", y = "Campelen gear effect")
  .ggsave("gap_vessel_effect.png", g, height = 3, width = 4)

  #### Design-based index ----

  # Compare strat area with values in prediction grid
  grid_area <- grid %>%
    summarise(area_km2 = sum(area_km2), .by = c(strat, div))

  survey_area <- dat_fit %>%
    summarise(strat.area = mean(strat.area, na.rm = TRUE), .by = c(strat, NAFOdiv)) %>%
    rename(div = NAFOdiv) %>%
    left_join(grid_area) %>%
    filter(!is.na(strat), !is.na(area_km2)) %>%
    mutate(ratio = area_km2/strat.area)

  g <- rbind(
    dat_fit %>%
      make_design_index() %>%
      mutate(Calibrated = FALSE),
    dat_fit %>%
      mutate(number = number * calibration) %>%
      make_design_index() %>%
      mutate(Calibrated = TRUE)
  ) %>%
    ggplot(aes(year, est, linetype = Calibrated, shape = Calibrated)) +
    #geom_linerange(aes(ymin = lwr, ymax = upr)) +
    geom_point() +
    geom_line() +
    expand_limits(y = 0) +
    geom_vline(xintercept = 1994.5, linetype = 2) +
    labs(x = "Year", y = "Design index") +
    scale_linetype_manual(values = c(3, 1)) +
    scale_shape_manual(values = c(1, 16)) +
    theme(legend.position = "bottom")
  .ggsave("index_design.png", g, height = 4.5, width = 6)



  #### NULL test ----
  grid_null <- mutate(grid, dummy = TRUE)

  ydummy <- c(2005:2010, NA)
  fit_null <- do_null_test(ydummy, dat = dat_fit, sp = s)
  #fit_null <- lapply(ydummy, function(y) {
  #  readRDS(paste0("validation/models/", s, "_null_fit_", y, ".rds"))
  #})

  index_null <- make_index_ret(fit_null, grid_null)
  readr::write_csv(index_null, paste0("validation/index/", s, "_index_null.csv"))
  index_null <- readr::read_csv(paste0("validation/index/", s, "_index_null.csv")) %>%
    mutate(ydummy = rep(.env$ydummy, each = length(1995:2020)))

  g <- index_null %>%
    filter(!is.na(ydummy)) %>%
    ggplot(aes(year, est)) +
    #geom_line(data = index_null %>% filter(is.na(ydummy)), linewidth = 1) +
    geom_point(aes(colour = factor(ydummy))) +
    geom_line(aes(colour = factor(ydummy))) +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Campelen Index", colour = "Dummy year change") +
    theme(panel.spacing = unit(0, "in"))
  .ggsave("index_null.png", g, height = 3.5, width = 6)

  g <- index_null %>%
    filter(!is.na(ydummy)) %>%
    mutate(ydummylabel = paste(ydummy, "dummy change")) %>%
    ggplot(aes(year, est, colour = ydummylabel, shape = year < ydummy)) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = lwr, ymax = upr)) +
    facet_wrap(vars(ydummylabel)) +
    expand_limits(y = 0) +
    guides(colour = "none", shape = "none") +
    scale_shape_manual(values = c(16, 1)) +
    labs(x = "Year", y = "Campelen Index", colour = "Dummy year change") +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_null_ci.png", g, height = 3.5, width = 6)

  vessel_effect <- get_vessel_effect(fit_null[1:6]) %>%
    mutate(ydummy = .env$ydummy[1:6])

  g <- ggplot(vessel_effect, aes(ydummy, Estimate)) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Dummy year change", y = "Campelen gear effect")
  .ggsave("null_vessel_effect.png", g, height = 3, width = 4)

}



