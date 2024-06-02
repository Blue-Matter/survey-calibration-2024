
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
  mutate(data.series = factor(data.series, levels = c("Engel", "Campelen"))) %>%
  mutate(CF = FALSE)

# Data set to fit full model
dat_cf <- readr::read_csv("validation/ghl_set_cf.csv") %>%
  filter(NAFOdiv == "3L") %>%
  mutate(npue = number/swept.area) %>%
  mutate(data.series = factor(gearname, levels = c("Engel", "Campelen"))) %>%
  #mutate(long.start = -1 * long.start) %>%
  mutate(log_depth = log(set.depth.mean)) %>%
  sdmTMB::add_utm_columns(ll_names = c("long.start", "lat.start")) %>%
  mutate(CF = TRUE)

g <- ggplot(dat_cf, aes(long.start, lat.start, colour = npue)) +
  geom_point() +
  facet_wrap(vars(gearname)) +
  scale_colour_viridis_c(trans = "log")


# Data set to fit full model
dat_fit <- rbind(
  filter(dat, year >= 1985, year <= 2020) %>%
    select(year, X, Y, log_depth, data.series, number, swept.area, CF, NAFOdiv, long.start, lat.start),
  dat_cf %>%
    mutate(year = 1995) %>%
    select(year, X, Y, log_depth, data.series, number, swept.area, CF, NAFOdiv, long.start, lat.start)
)

# Compare CF data
nfld <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
  sf::st_crop(xmin = -58, xmax = -46, ymin = 45, ymax = 56)

g <- dat_fit %>%
  filter(year == 1995) %>%
  mutate(panel = ifelse(CF, paste("CF:", data.series), "Survey: Campelen")) %>%
  mutate(npue = number/swept.area) %>%
  ggplot(aes(long.start, lat.start, colour = npue)) +
  geom_sf(data = nfld, inherit.aes = FALSE) +
  coord_sf() +
  geom_point() +
  facet_wrap(vars(panel)) +
  labs(x = "Longitude", y = "Latitude", colour = "1995 NPUE") +
  theme(legend.position = "bottom",
        panel.spacing = unit(0, "in")) +
  scale_colour_viridis_c(trans = "log", breaks = c(148, 1096, 8103))
.ggsave("npue_map_cf.png", g, height = 6, width = 8)

g <- dat_fit %>%
  filter(year == 1995) %>%
  mutate(panel = ifelse(CF, paste("CF:", data.series), "Survey: Campelen")) %>%
  mutate(npue = number/swept.area) %>%
  ggplot(aes(exp(log_depth), npue)) +
  geom_point(alpha = 0.5) +
  facet_wrap(vars(panel)) +
  guides(colour = "none") +
  labs(x = "Depth (m)", y = "NPUE") +
  coord_cartesian(xlim = c(0, 1000))
.ggsave("npue_depth_cf.png", g, height = 2, width = 6)


source("99-functions.R")
fit <- fit_sdmTMB(dat_fit, formula = number ~ poly(log_depth, 2, raw = TRUE) + CF + data.series,
                  offset = log(dat_fit$swept.area))
saveRDS(fit, "validation/models/ghl_fit_cf.rds")
fit <- readRDS("validation/models/ghl_fit_cf.rds")

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(.data$depth > 0) %>%
  filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
  mutate(log_depth = log(depth),
         data.series = "Campelen",
         CF = FALSE)

index <- make_index(fit, grid)
readr::write_csv(index, "validation/index/ghl_index_cf.csv")
index <- readr::read_csv("validation/index/ghl_index_cf.csv")
