
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
dat_fit <- filter(dat, year >= 1995, year <= 2020)

source("99-functions.R")

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(.data$depth > 0) %>%
  filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
  mutate(log_depth = log(depth),
         data.series = "Campelen",
         dummy = TRUE)

ydummy <- c(2005:2010, NA)
fit_null <- do_null_test(ydummy, dat = dat_fit, sp = "ghl")
fit_null <- lapply(ydummy, function(y) {
  readRDS(paste0("validation/models/ghl_null_fit_", y, ".rds"))
})

index_null <- make_index_ret(fit_null, grid)

readr::write_csv(index_null, "validation/index/ghl_index_null.csv")
index_null <- readr::read_csv("validation/index/ghl_index_null.csv") %>%
  mutate(ydummy = rep(.env$ydummy, each = length(1995:2020)))


g <- index_null %>%
  filter(!is.na(ydummy)) %>%
  ggplot(aes(year, est)) +
  geom_line(data = index_null %>% filter(is.na(ydummy)), linewidth = 1) +
  geom_point(aes(colour = factor(ydummy))) +
  geom_line(aes(colour = factor(ydummy))) +
  expand_limits(y = 0) +
  labs(x = "Year", y = "Campelen Index", colour = "Dummy year change") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_null.png", g, height = 3.5, width = 6)

vessel_effect <- get_vessel_effect(fit_null[1:6]) %>%
  mutate(ydummy = .env$ydummy[1:6])

g <- ggplot(vessel_effect, aes(ydummy, Estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = "Dummy year change", y = "Campelen gear effect")
.ggsave("null_vessel_effect.png", g, height = 3, width = 4)




