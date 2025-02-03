
library(tidyverse)
dat <- readRDS("validation_dogfish/wrangled-hbll-dog-sets.rds")

dat$survey_desc %>% table()

dat$survey_abbrev %>% unique()



dat %>% filter(survey_desc == "Dogfish Gear/Timing Comparison Surveys") %>% View()

dat %>%
  filter(grepl('dog', survey_abbrev)) %>%
  summarise(nsets = n(), .by = c(year, survey_abbrev)) %>%
  arrange(year)

# Calibration from 2004 survey only
comp <- filter(dat, grepl("Comparison", survey_desc), year == 2004) %>%
  mutate(catch_all = sum(catch_count), .by = fishing_event_id)

# Comparative fishing sets
comp2 <- lapply(unique(comp$fishing_event_id), function(i) {
  dfcircle <- filter(comp, fishing_event_id == i, survey_abbrev == "dog")
  dfj <- filter(comp, fishing_event_id == i, survey_abbrev == "dog-jhook")

  if (nrow(dfcircle) && nrow(dfj)) {
    data.frame(
      id = i,
      number_circle = dfcircle$catch_count,
      number_jhook = dfj$catch_count,
      offset = dfcircle$offset - dfj$offset
    )
  } else {
    data.frame()
  }
}) %>%
  bind_rows() %>%
  mutate(
    X = rnorm(nrow(.)), Y = rnorm(nrow(.)), id = factor(id)
  )

library(sdmTMB)
m <- sdmTMB(
  cbind(number_circle, number_jhook) ~ 1,
  mesh = comp2 %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
  spatial = "off",
  spatiotemporal = "off",
  data = comp2,
  offset = comp2$offset,
  family = binomial(),
  control = sdmTMBcontrol(multiphase = FALSE)
)

m2 <- sdmTMB(
  cbind(number_circle, number_jhook) ~ 1 + (1 | id),
  mesh = comp2 %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
  spatial = "off",
  spatiotemporal = "off",
  data = comp2,
  offset = comp2$offset,
  family = binomial(),
  control = sdmTMBcontrol(multiphase = FALSE)
)


# Index
# Does not have "Dogfish Gear/Timing Comparison Surveys"
# Either add the calibration factor, estimate it, or ignore the issue
log_rho <- coef(m2)["(Intercept)"]

dog <- filter(dat, grepl("Strait of Georgia", survey_desc)) %>%
  mutate(offset_rho = offset - ifelse(survey_abbrev == "dog", log_rho, 0),
         cpue = catch_count/exp(offset),
         cpue_rho = catch_count/exp(offset_rho),
         log_depth = log(depth_m)) %>%
  rename(X = UTM.lon, Y = UTM.lat) %>%
  arrange(year)

g <- ggplot(dog, aes(longitude, latitude, fill = cpue)) +
  geom_point(shape = 21) +
  facet_wrap(vars(year)) +
  scale_fill_viridis_c(trans = "log")

g <- ggplot(dog, aes(depth_m, cpue)) +
  geom_point()

# Fit index
source("99-functions.R")

grid <- as.data.frame(gfplot::dogfish_grid$grid/1000) %>%
  mutate(area_km2 = gfplot::dogfish_grid$cell_area,
         log_depth = 0,
         survey_abbrev = "dog")
#plot(grid$X, grid$Y, col = 2)
#points(dog$X, dog$Y)

# Estimated rho for circle hooks - nbinom2() family and rw spatiotemporal field
fit <- fit_sdmTMB(
  dog,
  formula = catch_count ~ 0 + poly(log_depth, 2, raw = TRUE),
  offset = dog$offset_rho,
  do_fit = TRUE,
  spatiotemporal = "off",
  spatial = "off",
  time_varying = ~ 1,
  anisotropy = FALSE,
  time_varying_type = "rw"
)
saveRDS(fit, "validation_dogfish/models/dogfish_fit_calibration.rds")

index <- make_index(fit, grid[1, ], year_new_vessel = 2004)
readr::write_csv(index, "validation_dogfish/index/dogfish_index_calibration.csv")

# Calibrate
fit <- fit_sdmTMB(
  dog,
  formula = catch_count ~ 0 + poly(log_depth, 2, raw = TRUE) + survey_abbrev,
  offset = dog$offset,
  do_fit = TRUE,
  spatiotemporal = "off",
  spatial = "off",
  time_varying = ~ 1,
  anisotropy = FALSE,
  time_varying_type = "rw"
)
saveRDS(fit, "validation_dogfish/models/dogfish_fit_full.rds")
index <- make_index(fit, grid, year_new_vessel = 2004)
readr::write_csv(index, "validation_dogfish/index/dogfish_index_full.csv")

# No rho
fit <- fit_sdmTMB(
  dog,
  formula = catch_count ~ 1 + poly(log_depth, 2, raw = TRUE),
  offset = dog$offset,
  do_fit = TRUE
)
saveRDS(fit, "validation_dogfish/models/dogfish_fit_norho.rds")

index <- make_index(fit, grid, year_new_vessel = 2004)
readr::write_csv(index, "validation_dogfish/index/dogfish_index_norho.csv")


g <- ggplot(index, aes(year, est)) +
  geom_linerange(aes(ymin = lwr, ymax = upr)) +
  geom_point(aes(shape = impute)) +
  scale_shape_manual(values = c('TRUE' = 1, 'FALSE' = 16))

