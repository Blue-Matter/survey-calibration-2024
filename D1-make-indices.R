
#2009

library(tidyverse)
library(survdat)
library(sf)
library(sdmTMB)

source('99-functions.R')

strata <- sf::st_read(dsn = system.file("extdata", "strata.shp", package = "survdat"), quiet = TRUE)
area <- sf::st_read(dsn = system.file("extdata", "EPU.shp", package = "survdat"), quiet = TRUE)

#plot(strata)
#splot(area)


# Average swept area for Albatross
#readr::read_csv("validation_nefsc/NEFSC_spring_03282025.csv") %>%
#  filter(SV_VESSEL == "AL") %>%
#  filter(!is.na(AREA_SWEPT_WINGS_MEAN_KM2)) %>%
#  pull(AREA_SWEPT_WINGS_MEAN_KM2) %>%
#  mean(na.rm = TRUE)

# Average swept area for HB
#readr::read_csv("validation_nefsc/NEFSC_spring_03282025.csv") %>%
#  filter(SV_VESSEL == "HB") %>%
#  filter(!is.na(AREA_SWEPT_WINGS_MEAN_KM2)) %>%
#  pull(AREA_SWEPT_WINGS_MEAN_KM2) %>%
#  mean(na.rm = TRUE)

#table(dat_set$year, dat_set$SVGEAR)
dat_set <- readr::read_csv("validation_nefsc/NEFSC_spring_03282025.csv") %>%
  mutate(DT = as.POSIXlt(EST_DATE_TIME, format = "%m/%d/%Y %H:%M:%S"),
         D = as.POSIXlt(EST_DATE_TIME, format = "%m/%d/%Y")) %>%
  mutate(year = lubridate::year(DT)) %>%
  mutate(year = ifelse(is.na(year), lubridate::year(D), year)) %>%
  filter(SV_VESSEL != "DE", year >= 1983) %>%
  mutate(EFFORT = ifelse(is.na(AREA_SWEPT_WINGS_MEAN_KM2) & SV_VESSEL == "AL", 0.03995142, AREA_SWEPT_WINGS_MEAN_KM2)) %>%
  mutate(EFFORT = ifelse(is.na(AREA_SWEPT_WINGS_MEAN_KM2) & SV_VESSEL == "HB", 0.02393521, EFFORT)) %>%
  mutate(SV_VESSEL = factor(SV_VESSEL, levels = c("AL", "HB")))


# Subset for GOM and GB
dat_sf <- sf::st_as_sf(dat_set, coords = c("LON", "LAT"), crs = "NAD83")
set_intersect <- sf::st_intersects(
  area[1:2, ],
  dat_sf
) %>% unlist()


# See Miller 2010 for calibration factors from comparative fishing
#
# Table 3 for model definitions (p. 35)
# Table 23 for AIC comparisons (p. 125)
# Table 24 for rho from M_1,1 (p. 132)
# Table 31 for rho from M_2,3 (p. 152)
#
# Recommended values are in Tables 53-55 (p. 222)
comparative_rho_summary <- data.frame(
  sp = c("redfish", "yt", "whitehake"),
  SVSPP = c(155, 105, 076),
  model = c("M_1,1", "M_2,3", "M_1,1"),
  rho = c(1.4564, 2.3619, 2.2354),
  se = c(0.1301, 0.1791, 0.1726),
  lower = c(1.2224, 2.0357, 1.9214),
  upper = c(1.7351, 2.7403, 2.6008)
)


#sp <- c("redfish", "yt", "whitehake")
#sp <- c("yt", "whitehake")
sp <- "whitehake"

for (s in sp) {

  if (s == "redfish") {
    dat <- dat_set[set_intersect, ] %>%
      mutate(CATCH = ifelse(is.na(SEBASTES_F_NUM), 0, SEBASTES_F_NUM)) %>%
      mutate(NPUE = CATCH/EFFORT) %>%
      mutate(calibration = ifelse(SV_VESSEL == "AL", 1.456, 1)) # Species 155, see survdat:::apply_conversion_factors

  } else if (s == "yt") {
    dat <- dat_set[set_intersect, ] %>%
      mutate(CATCH = ifelse(is.na(LIMANDA_F_NUM), 0, LIMANDA_F_NUM)) %>%
      mutate(NPUE = CATCH/EFFORT) %>%
      mutate(calibration = ifelse(SV_VESSEL == "AL", 2.347, 1)) # Species 105, see survdat:::apply_conversion_factors

  } else {

    catch_all <- readr::read_csv("validation_nefsc/Spring/22561_UNION_FSCS_SVCAT.csv") %>%
      select(CRUISE6, STRATUM, TOW, STATION, SVSPP, EXPCATCHNUM, SCIENTIFIC_NAME) %>%
      filter(SVSPP == '076')

    dat <- dat_set[set_intersect, ] %>%
      left_join(catch_all) %>%
      mutate(CATCH = ifelse(is.na(EXPCATCHNUM), 0, EXPCATCHNUM)) %>%
      mutate(NPUE = EXPCATCHNUM/EFFORT) %>%
      mutate(calibration = ifelse(SV_VESSEL == "AL", 2.235, 1)) # Species 076, see survdat:::apply_conversion_factors
  }

  .ggsave <- function(filename, ...) {
    filename <- paste0("validation_nefsc/figures/", s, "_", filename)
    ggplot2::ggsave(filename, ...)
  }

  #### Spatial figure of set data ----
  #g <- plot_shapefile(area[1:2, ]) +
  #  geom_point(data = dat, aes(LON, LAT, colour = NPUE), size = 0.5) +
  #  scale_colour_viridis_c(trans = 'log') +
  #  facet_wrap(vars(year)) +
  #  labs(x = "Longitude", y = "Latitude", colour = "NPUE") +
  #  #theme(legend.position = "bottom", panel.spacing = unit(0, "in")) +
  #  theme(panel.spacing = unit(0, "in"))
  #.ggsave("npue_map.png", g, height = 12, width = 9)

  #### Create spatial grid ----
  dat_utm <- sdmTMB::add_utm_columns(dat, ll_names = c("LON", "LAT"), ll_crs = 4326, utm_crs = 32619, units = "km")

  grid <- expand.grid(
    X = seq(floor(min(dat_utm$X)), ceiling(max(dat_utm$X)), by = 10),
    Y = seq(floor(min(dat_utm$Y)), ceiling(max(dat_utm$Y)), by = 10)
  ) %>%
    as.data.frame() %>%
    mutate(SV_VESSEL = "HB", X_m = 1e3 * X, Y_m = 1e3 * Y, area_km2 = 100) %>%
    add_utm_columns(ll_names = c("X_m", "Y_m"), ll_crs = 32619, utm_names = c("LON", "LAT"), utm_crs = 4326, units = "m")
  #plot(grid$LON, grid$LAT)

  grid_sf <- sf::st_as_sf(grid, coords = c("LON", "LAT"), crs = "NAD83")
  grid_intersect <- sf::st_intersects(
    area[1:2, ],
    grid_sf
  ) %>% unlist()

  grid2 <- grid[grid_intersect, ]
  #plot_shapefile(area[1:2, ]) +
  #  geom_point(data = grid2, aes(LON, LAT))

  #mesh <- make_mesh(dat_utm, c("X", "Y"), type = "cutoff", cutoff = 20)

  ### Estimate calibration factor
  fit <- fit_sdmTMB(
    dat_utm,
    formula = CATCH ~ 1 + SV_VESSEL,
    offset = log(dat_utm$EFFORT)
  )
  saveRDS(fit, paste0("validation_nefsc/models/", s, "_fit_full.rds"))
  #fit <- readRDS(paste0("validation_nefsc/models/", s, "_fit_full.rds"))

  index <- make_index(fit, grid2)
  readr::write_csv(index, paste0("validation_nefsc/index/", s, "_index.csv"))
  #index <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index.csv"))

  ### Provide calibration factor from CF
  fit_rho <- fit_sdmTMB(
    dat_utm,
    formula = CATCH ~ 1,
    offset = log(dat_utm$EFFORT/dat_utm$calibration)
  )
  saveRDS(fit_rho, paste0("validation_nefsc/models/", s, "_fit_calibration.rds"))
  #fit_rho <- readRDS(paste0("validation_nefsc/models/", s, "_fit_calibration.rds"))

  index_rho <- make_index(fit_rho, grid2)
  readr::write_csv(index_rho, paste0("validation_nefsc/index/", s, "_index_calibration.csv"))
  #index_rho <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index_calibration.csv"))

  ### No calibration factor
  fit_norho <- fit_sdmTMB(
    dat_utm,
    formula = CATCH ~ 1,
    offset = log(dat_utm$EFFORT)
  )
  saveRDS(fit_norho, paste0("validation_nefsc/models/", s, "_fit_norho.rds"))
  #fit_norho <- readRDS(paste0("validation_nefsc/models/", s, "_fit_norho.rds"))

  index_norho <- make_index(fit_norho, grid2)
  readr::write_csv(index_norho, paste0("validation_nefsc/index/", s, "_index_norho.csv"))
  #index_norho <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index_norho.csv"))

  ### Plot and compare indices ----
  index_compare <- rbind(
    index_rho %>% mutate(Model = "Comparative rho"),
    index %>% mutate(Model = "Estimate rho"),
    index_norho %>% mutate(Model = "No rho")
  )
  g <- ggplot(index_compare, aes(year, est, linetype = Model, shape = Model)) +
    geom_point() +
    geom_line() +
    geom_vline(xintercept = 2008.5, linetype = 2) +
    labs(x = "Year", y = "Harry Bigelow Index") +
    scale_shape_manual(values = c(16, 1, 4)) +
    scale_linetype_manual(values = c(1, 2, 3)) +
    theme(legend.position = "bottom")
  .ggsave("index_no_comparative.png", g, height = 4, width = 6)

  dat_utm <- dat_utm %>% rename(swept.area = EFFORT)

  #### Do retrospective ----
  yret <- 2009:2014
  fit_ret <- do_retrospective(yret, dat = dat_utm, sp = s,
                              formula = CATCH ~ 1 + SV_VESSEL,
                              dir_name = "validation_nefsc")
  fit_ret <- lapply(yret, function(y) {
    readRDS(paste0("validation_nefsc/models/", s, "_ret_fit_", y, ".rds"))
  })

  index_ret <- make_index_ret(fit_ret, grid = grid2)
  readr::write_csv(index_ret, paste0("validation_nefsc/index/", s, "_index_ret.csv"))
  #index_ret <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index_ret.csv"))

  vessel_effect <- get_vessel_effect(fit_ret, par_no = 1) %>%
    mutate(year_ret = .env$yret[.data$sequence])

  log_rho_comparative <- filter(comparative_rho_summary, sp == s) %>%
    mutate(log_rho = log(rho), log_lower = log(lower), log_upper = log(upper))

  log_rho_df <- lapply(c(-Inf, Inf), function(x) mutate(log_rho_comparative, x = x)) %>%
    bind_rows()

  g <- ggplot(vessel_effect, aes(year_ret, Estimate)) +
    geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = log_lower, ymax = log_upper), inherit.aes = FALSE) +
    geom_hline(yintercept = log_rho_comparative$log_lower, linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$log_upper, linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$log_rho, linetype = 2) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    expand_limits(y = 0) +
    labs(x = "Terminal year", y = "Harry Bigelow gear effect")
  .ggsave("retrospective_vessel_effect.png", g, height = 3, width = 4)

  g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
    geom_point() +
    geom_line() +
    expand_limits(y = 0) +
    geom_vline(xintercept = 2008.5, linetype = 2) +
    #coord_cartesian(xlim = c(1985, 2000)) +
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Terminal year")
  .ggsave("index_retrospective.png", g, height = 3, width = 6)

  g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
    geom_point() +
    geom_line() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    expand_limits(y = 0) +
    facet_wrap(vars(year_ret)) +
    geom_vline(xintercept = 2008.5, linetype = 2) +
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Terminal year") +
    guides(colour = "none") +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_retrospective_ci.png", g, height = 3.5, width = 6)


  ### Fit with gap ----
  year_new_vessel <- 2009
  fit_gap <- do_retrospective_gap(gap = 1:5, dat = dat_utm, sp = s, year_new_vessel = year_new_vessel,
                                  formula = CATCH ~ 1 + SV_VESSEL, #do_fit = FALSE,
                                  dir_name = "validation_nefsc")
  #fit_gap <- lapply(1:5, function(i) {
  #  readRDS(paste0("validation/models/", s, "_gap_fit_", i, ".rds"))
  #})

  index_gap <- make_index_ret(fit_gap, grid2) %>%
    mutate(gap = gap - min(gap) + 1) %>%
    mutate(impute = year >= year_new_vessel & year < year_new_vessel + gap)
  readr::write_csv(index_gap, paste0("validation_nefsc/index/", s, "_index_gap.csv"))
  #index_gap <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index_gap.csv"))

  g <- ggplot(index_gap, aes(year, est)) +
    geom_point(aes(shape = impute, colour = factor(gap))) +
    geom_line(linetype = 2, aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = gap), aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = gap), aes(colour = factor(gap))) +
    expand_limits(y = 0) +
    geom_vline(xintercept = 2008.5, linetype = 2) +
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
    guides(shape = "none") +
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Survey gap")
  .ggsave("index_gap.png", g, height = 3, width = 6)

  g <- ggplot(index_gap, aes(year, log(est))) +
    geom_point(aes(shape = impute, colour = factor(gap))) +
    geom_line(linetype = 2, aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = gap), aes(colour = factor(gap))) +
    geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = gap), aes(colour = factor(gap))) +
    #expand_limits(y = 0) +
    geom_vline(xintercept = 2008.5, linetype = 2) +
    scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
    guides(shape = "none") +
    labs(x = "Year", y = "Harry Bigelow log(Index)", colour = "Survey gap")
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
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Terminal year") +
    guides(colour = "none", shape = "none") +
    scale_shape_manual(values = c(16, 1)) +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_gap_ci.png", g, height = 3.5, width = 6)

  vessel_effect <- get_vessel_effect(fit_gap, par_no = 1)

  g <- ggplot(vessel_effect, aes(sequence, Estimate)) +
    geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = log_lower, ymax = log_upper), inherit.aes = FALSE) +
    geom_hline(yintercept = log_rho_comparative$log_lower, linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$log_upper, linetype = 3) +
    geom_hline(yintercept = log_rho_comparative$log_rho, linetype = 2) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    labs(x = "Survey gap", y = "Harry Bigelow gear effect")
  .ggsave("gap_vessel_effect.png", g, height = 3, width = 4)

  #### NULL test ----
  grid_null <- mutate(grid2, dummy = TRUE)

  ydummy <- c(2015:2020, NA)
  dat_2009 <- dat_utm %>% filter(year >= 2009)
  fit_null <- do_null_test(ydummy, dat = dat_2009, sp = s,
                           formula = CATCH ~ 1,
                           formula_null = CATCH ~ 1 + dummy, #do_fit = FALSE,
                           dir_name = "validation_nefsc")
  #fit_null <- lapply(ydummy, function(y) {
  #  readRDS(paste0("validation/models/", s, "_null_fit_", y, ".rds"))
  #})

  index_null <- make_index_ret(fit_null, grid_null)
  readr::write_csv(index_null, paste0("validation_nefsc/index/", s, "_index_null.csv"))
  index_null <- readr::read_csv(paste0("validation_nefsc/index/", s, "_index_null.csv")) %>%
    mutate(ydummy = rep(.env$ydummy, each = length(2009:2024)))

  g <- index_null %>%
    filter(!is.na(ydummy)) %>%
    ggplot(aes(year, est)) +
    #geom_line(data = index_null %>% filter(is.na(ydummy)), linewidth = 1) +
    geom_point(aes(colour = factor(ydummy))) +
    geom_line(aes(colour = factor(ydummy))) +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Dummy year change") +
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
    labs(x = "Year", y = "Harry Bigelow Index", colour = "Dummy year change") +
    theme(panel.spacing = unit(0, "in"), strip.background = element_rect(fill = NA, colour = NA))
  .ggsave("index_null_ci.png", g, height = 3.5, width = 6)

  vessel_effect <- get_vessel_effect(fit_null[1:6], par_no = 1) %>%
    mutate(ydummy = .env$ydummy[1:6])

  g <- ggplot(vessel_effect, aes(ydummy, Estimate)) +
    geom_point() +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    geom_hline(yintercept = 0, linetype = 2) +
    labs(x = "Dummy year change", y = "Harry Bigelow gear effect")
  .ggsave("null_vessel_effect.png", g, height = 3, width = 4)
}





