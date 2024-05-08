
library(tidyverse)
#library(NAFOdown)

#ggplot2::theme_set(theme_nafo())
theme_set(theme_bw())

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
  mutate(calibration = ifelse(data.series == "Engel", exp(0.8913), 1))

# Data set to fit full model
dat_fit <- filter(dat, year >= 1985, year <= 2020)

source("99-functions.R")
fit <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE),
                  offset = log(dat_fit$swept.area/dat_fit$calibration))
saveRDS(fit, "validation/models/cod_fit_calibration.rds")
fit <- readRDS("validation/models/cod_fit_calibration.rds")

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(.data$depth > 0) %>%
  filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
  mutate(log_depth = log(depth),
         data.series = "Campelen")

index <- make_index(fit, grid)
readr::write_csv(index, "validation/index/cod_index_calibration.csv")
index <- readr::read_csv("validation/index/cod_index_calibration.csv")

g <- ggplot(index, aes(year, est)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
  guides(colour = "none") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_comparative.png", g, height = 3.5, width = 6)

g <- ggplot(index, aes(year, est)) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
  guides(colour = "none") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_comparative_ci.png", g, height = 3.5, width = 6)
