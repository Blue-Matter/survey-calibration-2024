
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
  sdmTMB::add_utm_columns(ll_names = c("long.start", "lat.start"))


nfld <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
  sf::st_crop(xmin = -59, xmax = -46, ymin = 45, ymax = 56)

g <- ggplot(dat) +
  geom_sf(data = nfld) +
  coord_sf() +
  geom_point(aes(long.start, lat.start, colour = npue), size = 0.5) +
  scale_colour_viridis_c(trans = 'log', breaks = c(55, 403, 2980)) +
  facet_wrap(vars(year)) +
  labs(x = "Longitude", y = "Latitude", colour = "NPUE") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in")) +
  ggtitle("Cod 2J+3KL")
.ggsave("npue_map.png", g, height = 12, width = 9)


g <- ggplot(dat, aes(npue, fill = data.series)) +
  facet_wrap(vars(year)) +
  geom_histogram()

g <- ggplot(dat, aes(npue, fill = data.series)) +
  facet_grid(vars(year), vars(NAFOdiv)) +
  geom_histogram()


summarise(dat, n = n(), .by = c(year, vessel, NAFOdiv)) %>%
  arrange(year, vessel, NAFOdiv) %>%
  View()

g <- ggplot(dat, aes(set.depth.mean, npue, colour = data.series)) +
  geom_point() +
  facet_wrap(vars(year)) +
  coord_cartesian(xlim = c(0, 800))


# Data set to fit full model
dat_fit <- filter(dat, year >= 1985, year <= 2020)

source("99-functions.R")
fit <- fit_sdmTMB(dat_fit)
saveRDS(fit, "validation/models/cod_fit_full.rds")

grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
         depth = -1 * depth) %>%
  filter(depth > 0) %>%
  mutate(log_depth = log(depth),
         data.series = "Engel")

index <- make_index(fit, grid)
readr::write_csv(index, "validation/index/cod_index.csv")

# Do retrospective
yret <- 1995:2000
fit_ret <- do_retrospective(yret, dat = dat_fit)
index_ret <- make_index_ret(fit_ret, grid)
readr::write_csv(index_ret, "validation/index/cod_index_ret.csv")


vessel_effect <- lapply(1:length(fit_ret), function(i) {
  summary(fit_ret[[i]]$sd_report, select = "fixed") %>%
    as.data.frame() %>%
    mutate(var = rownames(.)) %>%
    mutate(year = yret[i])
}) %>%
  bind_rows() %>%
  filter(var == "b_j.3") %>%
  mutate(lwr = Estimate - 1.96 * `Std. Error`, upr = Estimate + 1.96 * `Std. Error`)

g <- ggplot(vessel_effect, aes(year, Estimate)) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Terminal year", y = "Engel vessel effect") +
  geom_hline(yintercept = 0, linetype = 2)
.ggsave("retrospective_vessel_effect.png", g, height = 3, width = 4)

g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Engel Index", colour = "Terminal year")
.ggsave("index_retrospective.png", g, height = 3, width = 6)

g <- ggplot(index_ret, aes(year, est, colour = factor(year_ret))) +
  geom_point() +
  geom_line() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  expand_limits(y = 0) +
  facet_wrap(vars(year_ret)) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Engel Index", colour = "Terminal year") +
  guides(colour = "none") +
  theme(panel.spacing = unit(0, "in"))
.ggsave("index_retrospective_ci.png", g, height = 3.5, width = 6)


## Design-based index ----
index_design <- dat_fit %>%
  filter(!is.na(strat.area)) %>%
  mutate(npue_expand = npue * strat.area) %>%
  summarise(
    index_strat = mean(npue_expand),
    var_strat = var(npue_expand),
    n = n(),
    nsamp = unique(strat.area)/sum(swept.area) * n,
    .by = c(year, strat, strat.area)
  ) %>%
  mutate(var_strat = ifelse(is.na(var_strat), 0, var_strat)) %>%
  mutate(area_total = sum(strat.area), .by = year) %>%
  summarise(
    est = sum(index_strat),
    var = sum(nsamp * (nsamp - n)/n * var_strat)/sum(nsamp)/sum(nsamp), # See SimSurvey appendix
    .by = year
  ) %>%
  mutate(
    se = sqrt(var),
    cv = se/est,
    lwr = est - 1.96 * se,
    upr = est + 1.96 * se
  )

g <- index_design %>%
  ggplot(aes(year, est)) +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  geom_point() +
  geom_line() +
  expand_limits(y = 0) +
  geom_vline(xintercept = 1994.5, linetype = 2) +
  labs(x = "Year", y = "Design-based index")
.ggsave("index_design.png", g, height = 4, width = 6)
