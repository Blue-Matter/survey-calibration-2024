
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

# Data set to fit full model
dat_fit <- filter(dat, year >= 1995) %>%
  mutate(after2005 = year > 2005)

source("99-functions.R")
fit_null <- fit_sdmTMB(
  dat_fit,
  number ~ poly(log_depth, 2, raw = TRUE) + after2005
)
saveRDS(fit_null, "validation/models/ghl_fit_null.rds")
fit_null <- readRDS("validation/models/ghl_fit_null.rds")

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(.data$depth > 0) %>%
  filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
  mutate(log_depth = log(depth),
         data.series = "Campelen",
         after2005 = TRUE)

index_null <- make_index(fit_null, grid)
readr::write_csv(index_null, "validation/index/ghl_index_null.csv")
index_null <- readr::read_csv("validation/index/ghl_index_null.csv")


g <- ggplot(index_null, aes(year, est)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
  guides(colour = "none") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_null.png", g, height = 3.5, width = 6)

