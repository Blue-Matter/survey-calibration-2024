
library(tidyverse)
#library(NAFOdown)

#ggplot2::theme_set(theme_nafo())
theme_set(theme_bw())

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
  mutate(data.series = factor(data.series, levels = c("Engel", "Campelen")))

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

# Spatial figure of set data
nfld <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
  sf::st_crop(xmin = -59, xmax = -46, ymin = 45, ymax = 56)
g <- ggplot(dat) +
  geom_sf(data = nfld) +
  coord_sf() +
  geom_point(aes(long.start, lat.start, colour = npue), size = 0.5) +
  scale_colour_viridis_c(trans = 'log') +
  facet_wrap(vars(year)) +
  labs(x = "Longitude", y = "Latitude", colour = "NPUE") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in")) +
  ggtitle("Greenland halibut 2J+3KL")
.ggsave("npue_map.png", g, height = 12, width = 9)


g <- ggplot(dat, aes(npue, fill = data.series)) +
  facet_wrap(vars(year)) +
  geom_histogram()

g <- ggplot(dat, aes(npue, fill = data.series)) +
  #facet_grid(vars(year), vars(NAFOdiv)) +
  geom_histogram()


summarise(dat, n = n(), .by = c(year, vessel, NAFOdiv)) %>%
  arrange(year, vessel, NAFOdiv)

g <- ggplot(dat, aes(set.depth.mean, npue, colour = data.series)) +
  geom_point() +
  facet_wrap(vars(year)) +
  coord_cartesian(xlim = c(0, 800))


# Data set to fit full model
dat_fit <- filter(dat, year >= 1985, year <= 2020)

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(.data$depth > 0) %>%
  filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
  mutate(log_depth = log(depth),
         data.series = "Campelen")

source("99-functions.R")

fit <- fit_sdmTMB(dat_fit)
saveRDS(fit, "validation/models/ghl_fit_full.rds")
fit <- readRDS("validation/models/ghl_fit_full.rds")

index <- make_index(fit, grid)
readr::write_csv(index, "validation/index/ghl_index.csv")
index <- readr::read_csv("validation/index/ghl_index.csv")

fit_norho <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE))
saveRDS(fit_norho, "validation/models/ghl_fit_norho.rds")

index_norho <- make_index(fit_norho, grid)
readr::write_csv(index_norho, "validation/index/ghl_index_norho.csv")
index_norho <- readr::read_csv("validation/index/ghl_index_norho.csv")

index_best <- readr::read_csv("validation/index/ghl_index_calibration.csv")

index_compare <- rbind(
  index_best %>% mutate(Model = "Comparative rho"),
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

g <- ggplot(index_compare, aes(year, est)) +
  geom_point(shape = 1) +
  geom_line(linewidth = 0.1) +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  facet_wrap(vars(Model), ncol = 1) +
  labs(x = "Year", y = "Campelen Index") +
  scale_shape_manual(values = c(16, 1, 4)) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  theme(panel.spacing = unit(0, "in"), strip.background = element_rect(colour = NA, fill = NA))
.ggsave("index_no_comparative_ci.png", g, height = 7, width = 4)

g <- rbind(
  index %>% mutate(Model = "Estimate rho"),
  index_best %>% mutate(Model = "With comparative rho"),
  index_norho %>% mutate(Model = "No rho")
) %>%
  ggplot(aes(year, est, linetype = Model, shape = Model)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
  guides(colour = "none") +
  scale_shape_manual(values = c(1, 4, 16)) +
  scale_linetype_manual(values = c(2, 3, 1)) +
  theme(legend.position = "bottom")
.ggsave("index_no_comparative.png", g, height = 4, width = 6)

g <- ggplot(index, aes(year, est)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
  guides(colour = "none") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_no_comparative_ci.png", g, height = 3.5, width = 6)

# Do retrospective
yret <- 1995:2000
fit_ret <- do_retrospective(yret, dat = dat_fit, sp = "ghl")
fit_ret <- lapply(yret, function(y) {
  readRDS(paste0("validation/models/ghl_ret_fit_", y, ".rds"))
})
index_ret <- make_index_ret(fit_ret, grid)
readr::write_csv(index_ret, "validation/index/ghl_index_ret.csv")
index_ret <- readr::read_csv("validation/index/ghl_index_ret.csv")

vessel_effect <- get_vessel_effect(fit_ret) %>%
  mutate(year_ret = .env$yret[.data$sequence])

log_rho_comparative <- data.frame(
  x = c(-Inf, Inf),
  y1 = 0.8913 + 1.96 * 0.0616,
  y2 = 0.8913 - 1.96 * 0.0616
)

g <- ggplot(vessel_effect, aes(year_ret, Estimate)) +
  geom_ribbon(data = log_rho_comparative, fill = "grey80", aes(x, ymin = y2, ymax = y1), inherit.aes = FALSE) +
  geom_hline(yintercept = log_rho_comparative$y1[1], linetype = 3) +
  geom_hline(yintercept = log_rho_comparative$y2[1], linetype = 3) +
  geom_hline(yintercept = 0.8913, linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Terminal year", y = "Campelen gear effect")
.ggsave("retrospective_vessel_effect.png", g, height = 3, width = 4)

g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
  geom_line(data = index_best, aes(year, est), linewidth = 1, inherit.aes = FALSE) +
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
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_retrospective_ci.png", g, height = 3.5, width = 6)

# Fit with gap
fit_gap <- do_retrospective_gap(gap = 1:5, dat = dat_fit, sp = "ghl")
fit_gap <- lapply(1:5, function(i) {
  readRDS(paste0("validation/models/ghl_gap_fit_", i, ".rds"))
})
index_gap <- make_index_ret(fit_gap, grid)
readr::write_csv(index_gap, "validation/index/ghl_index_gap.csv")
index_gap <- readr::read_csv("validation/index/ghl_index_gap.csv")

g <- ggplot(index_gap, aes(year, est)) +
  geom_line(data = index_best, aes(year, est), linewidth = 1, inherit.aes = FALSE) +
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

g <- index_gap %>%
  ggplot(aes(year, log(est))) +
  geom_line(data = index_best, linewidth = 1) +
  geom_point(aes(shape = impute, colour = factor(gap))) +
  geom_line(linetype = 2, aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = gap), aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = gap), aes(colour = factor(gap))) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  coord_cartesian(xlim = c(1985, 2005)) +
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
  guides(shape = "none") +
  labs(x = "Year", y = "Campelen log(Index)", colour = "Survey gap")
.ggsave("index_gap2.png", g, height = 3, width = 6)

vessel_effect <- get_vessel_effect(fit_gap)

g <- ggplot(vessel_effect, aes(sequence, Estimate)) +
  geom_ribbon(data = log_rho_comparative, fill = "grey80", aes(x, ymin = y2, ymax = y1), inherit.aes = FALSE) +
  geom_hline(yintercept = log_rho_comparative$y1[1], linetype = 3) +
  geom_hline(yintercept = log_rho_comparative$y2[1], linetype = 3) +
  geom_hline(yintercept = 0.8913, linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Survey gap", y = "Campelen gear effect")
.ggsave("gap_vessel_effect.png", g, height = 3, width = 4)



## Design-based index ----
g <- rbind(
  dat_fit %>%
    make_design_index() %>%
    mutate(Calibrated = FALSE),
  dat_fit %>%
    mutate(number = number * ifelse(data.series == "Engel", exp(0.8913), 1)) %>%
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
