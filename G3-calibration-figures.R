
library(tidyverse)

#### Load set data ----
ghl <- readr::read_csv("calibration_ghl/ghl_sets.csv")
survey_vessel <- summarise(ghl, n = n(), .by = c(year, vessel))

#### Plot NPUE ----
coast <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf")
g <- ggplot(ghl, aes(wpos, npos, colour = NPUE)) +
  geom_sf(data = coast, inherit.aes = FALSE) +
  coord_sf(xlim = c(-75, -50), ylim = c(62, 73), expand = FALSE) +
  geom_point() +
  scale_colour_viridis_c(trans = "sqrt") +
  facet_wrap(vars(year)) +
  theme(panel.spacing = unit(0, "in")) +
  labs(x = "Longitude", y = "Latitude")
ggsave("calibration_ghl/figures/NPUE.png", g, height = 8, width = 8)

#### Plot index ----
index <- rbind(
  read.csv("calibration_ghl/index/index_full.csv") %>% mutate(type = "All vessels"),
  read.csv("calibration_ghl/index/index_exHM.csv") %>% mutate(type = "ex-HM"),
  read.csv("calibration_ghl/index/index_norho.csv") %>% mutate(type = "No calibration")
) %>%
  left_join(survey_vessel)

g <- ggplot(index, aes(year, est, group = type)) +
  geom_line() +
  geom_point(data = filter(index, !is.na(vessel)), aes(shape = vessel)) +
  geom_line(aes(y = lwr), linetype = 2) +
  geom_line(aes(y = upr), linetype = 2) +
  facet_wrap(vars(type)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'))
ggsave("calibration_ghl/figures/index_compare.png", g, height = 3, width = 7)

# Design index
index_design <- readr::read_delim(
  "Data pre-2025/D4_BioAbuTOTAL_CombineGHL_99-22_6Str0A.csv",
  delim = ";"
) %>%
  left_join(survey_vessel)

g <- ggplot(index_design, aes(year, MeanAbuT, shape = vessel)) +
  geom_point() +
  geom_line() +
  geom_line(aes(y = MeanAbuT - 1.96 * sqrt(varabu)), linetype = 2) +
  geom_line(aes(y = MeanAbuT + 1.96 * sqrt(varabu)), linetype = 2) +
  geom_linerange(
    data = index_design %>% filter(year >= 2019),
    aes(ymin = MeanAbuT - 1.96 * sqrt(varabu), ymax = MeanAbuT + 1.96 * sqrt(varabu))
  ) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Design index", shape = "Vessel") +
  theme(legend.position = "bottom")
ggsave("calibration_ghl/figures/index_design.png", g, height = 4, width = 5)

# Compare with and without 0B
index <- rbind(
  read.csv(file = "calibration_ghl/index/index_full.csv") %>% mutate(type = "0AB1CD"),
  read.csv(file = "calibration_ghl/index/index_full_ex0B.csv") %>% mutate(type = "0A1CD")
) %>%
  left_join(survey_vessel)

g <- ggplot(index, aes(year, est, group = type)) +
  geom_line() +
  geom_point(data = filter(index, !is.na(vessel)), aes(shape = vessel)) +
  geom_line(aes(y = lwr), linetype = 2) +
  geom_line(aes(y = upr), linetype = 2) +
  facet_wrap(vars(type)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'))

#### Plot spatial fields ----
grid <- read.csv("calibration_ghl/grid_0AB1CD.csv") %>%
  sdmTMB::add_utm_columns(ll_names = c("X", "Y"), ll_crs = 32621,
                          utm_names = c("lon", "lat"), utm_crs = 4326, units = "m") %>%
  filter(lat < 73, lat > 61.5) %>%
  mutate(X = X/1000, Y = Y/1000, area = area/1000/1000, vessel = "TA")

source("99-functions-calibration.R")
fit <- readRDS("calibration_ghl/models/fit.rds")
pred <- make_index(fit, grid, FALSE)

coast <- rnaturalearth::ne_coastline(scale = 50, returnclass = "sf") %>%
  sf::st_crop(xmin = -75, xmax = -50, ymin = 62, ymax = 73)

g <- pred$data %>%
  filter(year == 1999) %>%
  ggplot(aes(X * 1000, Y * 1000, fill = omega_s)) +
  geom_sf(data = coast, inherit.aes = FALSE) +
  coord_sf(expand = FALSE, default_crs = 32621) +
  geom_tile(width = 10000, height = 10000) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"))
ggsave("calibration_ghl/figures/omega_full.png", g, height = 5, width = 4)

g <- pred$data %>%
  ggplot(aes(X * 1000, Y * 1000, fill = epsilon_st)) +
  geom_sf(data = coast, inherit.aes = FALSE) +
  facet_wrap(vars(year)) +
  coord_sf(expand = FALSE, default_crs = 32621) +
  geom_tile(width = 10000, height = 10000) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"))
ggsave("calibration_ghl/figures/epsilon_full.png", g, height = 8, width = 6)

g <- pred$data %>%
  ggplot(aes(X * 1000, Y * 1000, fill = est)) +
  geom_sf(data = coast, inherit.aes = FALSE) +
  facet_wrap(vars(year)) +
  coord_sf(expand = FALSE, default_crs = 32621) +
  geom_tile(width = 10000, height = 10000) +
  labs(x = "Longitude", y = "Latitude") +
  scale_fill_viridis_c()
ggsave("calibration_ghl/figures/density_full.png", g, height = 8, width = 6)


#### Retrospective indices ----
index_ret <- rbind(
  read.csv("calibration_ghl/index/index_full.csv") %>% mutate(yret = 2024),
  read.csv("calibration_ghl/index/index_ret_2023.csv") %>% mutate(yret = 2023),
  read.csv("calibration_ghl/index/index_ret_2022.csv") %>% mutate(yret = 2022)
) %>%
  left_join(survey_vessel)

g <- ggplot(index_ret, aes(year, est, group = yret, colour = factor(yret))) +
  geom_line() +
  geom_point(data = filter(index_ret, !is.na(vessel)), aes(shape = vessel)) +
  #geom_line(aes(y = lwr), linetype = 2) +
  #geom_line(aes(y = upr), linetype = 2) +
  #facet_wrap(vars(yret)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel", colour = "Terminal year") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'),
        legend.box = "vertical")
ggsave("calibration_ghl/figures/index_compare_ret.png", g, height = 3.5, width = 4)

g <- ggplot(index_ret, aes(year, est, group = yret)) +
  geom_line() +
  geom_point(data = filter(index_ret, !is.na(vessel)), aes(shape = vessel)) +
  geom_line(aes(y = lwr), linetype = 2) +
  geom_line(aes(y = upr), linetype = 2) +
  facet_wrap(vars(yret)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel", colour = "Terminal year") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'))
ggsave("calibration_ghl/figures/index_compare_ret_ci.png", g, height = 3, width = 7)

#### Gap indices ----
ygap <- 2022:2023

index_gap <- rbind(
  read.csv("calibration_ghl/index/index_full.csv") %>% mutate(type = "All vessels"),
  read.csv("calibration_ghl/index/index_gap_2023.csv") %>% mutate(type = "ex-22/23"),
  read.csv("calibration_ghl/index/index_gap_2022.csv") %>% mutate(type = "ex-2022")
) %>%
  left_join(survey_vessel)

g <- ggplot(index_gap, aes(year, est, group = type, colour = factor(type))) +
  geom_line() +
  geom_point(data = filter(index_gap, !is.na(vessel)), aes(shape = vessel)) +
  #geom_line(aes(y = lwr), linetype = 2) +
  #geom_line(aes(y = upr), linetype = 2) +
  #facet_wrap(vars(yret)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel", colour = "Data") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'),
        legend.box = "vertical")
ggsave("calibration_ghl/figures/index_compare_gap.png", g, height = 3.5, width = 4)

g <- ggplot(index_gap, aes(year, est, group = type, colour = factor(type))) +
  geom_line() +
  geom_point(data = filter(index_gap, !is.na(vessel)), aes(shape = vessel)) +
  geom_line(aes(y = lwr), linetype = 2) +
  geom_line(aes(y = upr), linetype = 2) +
  facet_wrap(vars(type)) +
  expand_limits(y = 0) +
  scale_shape_manual(values = c("HM" = 4, "TA" = 16, "PA" = 1, "NA" = NA)) +
  labs(x = "Year", y = "Tarajoq index", shape = "Vessel", colour = "Data") +
  theme(legend.position = "bottom", panel.spacing = unit(0, 'in'))
ggsave("calibration_ghl/figures/index_compare_gap_ci.png", g, height = 3, width = 7)



#### Plot rho ----
get_rho <- function(filename, fit, name) {
  if (missing(fit)) fit <- readRDS(file = filename)
  x <- sdmTMB::tidy(fit)
  if (!missing(name)) {
    x$type <- name
  }
  return(x)
}

rho_all <- get_rho("calibration_ghl/models/fit.rds", name = "All vessels")
rho_exHM <- get_rho("calibration_ghl/models/fit_exHM.rds", name = "ex-HM")

# Ret
yret <- 2022:2023
filename <- paste0("calibration_ghl/models/fit_ret_", yret, ".rds")
name_ret <- paste("Ret.", yret)
rho_ret <- Map(get_rho, filename = filename, name = name_ret) %>%
  bind_rows()

# Gap
ygap <- 2022:2023
filename2 <- paste0("calibration_ghl/models/fit_gap_", ygap, ".rds")
name_gap <- c("ex-2022", "ex-22/23")
rho_gap <- Map(get_rho, filename = filename2, name = name_gap) %>%
  bind_rows()

model_name <- c("All vessels", "ex-HM", name_ret, name_gap)
log_rho <- rbind(
  rho_all, rho_exHM, rho_ret, rho_gap
) %>%
  filter(grepl("vessel", term)) %>%
  mutate(vessel = strsplit(term, "vessel") %>% sapply(getElement, 2)) %>%
  mutate(y = paste0(vessel, " (", type, ")")) %>%
  mutate(type = factor(type, rev(model_name)))

g1 <- ggplot(log_rho, aes(type, estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = conf.low, ymax = conf.high)) +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_wrap(vars(vessel)) +
  labs(x = "Model", y = expression(widehat(log(rho)))) +
  coord_flip()

g2 <- ggplot(log_rho, aes(type, exp(estimate))) +
  geom_point() +
  geom_linerange(aes(ymin = exp(conf.low), ymax = exp(conf.high))) +
  geom_hline(yintercept = 1, linetype = 2) +
  facet_wrap(vars(vessel)) +
  labs(x = "Model", y = expression(hat(rho))) +
  expand_limits(y = 0) +
  coord_flip()

g <- ggpubr::ggarrange(g1, g2, ncol = 1)
ggsave("calibration_ghl/figures/rho_compare.png", g, height = 4, width = 4)


#### Null test - plot calibration factor ----
ytest <- 2013:2017
fit_null <- readRDS(file = "calibration_ghl/models/fit_null.rds")

rho_null <- Map(get_rho, fit = fit_null, name = ytest) %>%
  bind_rows() %>%
  filter(grepl("new_gear", term))

g <- ggplot(rho_null, aes(type, estimate, ymin = conf.low, ymax = conf.high)) +
  geom_point() +
  geom_linerange() +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(x = "Dummy year change", y = expression(widehat(log(rho))))
ggsave("calibration_ghl/figures/dummy_test.png", g, height = 3, width = 3)

#### Residuals ----
if (FALSE) {

  library(sdmTMB)

  fit <- readRDS(file = "calibration_ghl/models/fit.rds")
  fit_sims <- sdmTMB:::simulate.sdmTMB(fit, nsim = 1000)
  pred <- exp(predict(fit)$est)

  # Plot residuals

  library(DHARMa)
  res <- DHARMa::createDHARMa(fit_sims, fit$data$num, pred, integerResponse = TRUE)
  plot(res)
  plotResiduals(res)
  testOverdispersion(res)
  testOutliers(res)
  testSpatialAutocorrelation(res, x = fit$data$X, y = fit$data$Y)
  testTemporalAutocorrelation(res, time = fit$data$year)

  recalculateResiduals(res, fit$data$year) %>%
    testTemporalAutocorrelation(time = unique(fit$data$year))

  hist(res$scaledResiduals)
}
