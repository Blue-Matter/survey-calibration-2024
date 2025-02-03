

library(tidyverse)
#library(glmmTMB)
library(sdmTMB)
library(DHARMa)

.ggsave <- function(filename, sp = "cod", ...) {
  filename <- file.path("validation_size", "figures", paste0(sp, "_", filename))
  ggplot2::ggsave(filename, ...)
}

setdat <- readr::read_csv("validation/CF_setdetails.csv") %>%
  select(-spec, -number, -weight)

setbio <- readxl::read_excel("validation/CF_catchdetails_438_892.xlsx") %>%
  #select(-time.mid, -bot.temp, -gear.mondevice) %>%
  #mutate(strat = as.numeric(strat)) %>%
  select(spec, number, weight, trip, set, vessel, gear)

vessel_df <- data.frame(
  vessel = c(20, 39, 30, 34),
  vesselname = c("Gadus Atlantica", "Teleost", "W. Templeman", "Alfred Needler")
)

gear_df <- data.frame(
  gear = c(61, 101, 103),
  gearname = c("Campelen", "Engel", "Engel")
)

# Cod = 438, Halibut = 892
spec_df <- data.frame(
  spec = c(438, 892),
  species = c("cod", "ghl")
)

sp <- c("cod", "ghl")

for (s in sp) {

  code <- ifelse(s == "cod", 438, 892)

  set_len <- readr::read_csv("validation_size/CF_lengthfreqs_438_892.csv") %>%
    filter(spec == code, year == 95) %>%
    filter(length.group == 1) %>%
    select(vessel, trip, set, NAFOdiv, sex, starts_with("bin")) %>%
    reshape2::melt(id.vars = c("vessel", "trip", "set", "NAFOdiv", "sex")) %>%
    rename(bin = variable) %>%
    #mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric()) %>%
    #mutate(binmid = bin) %>%
    summarise(value = sum(value),
              .by = c(vessel, trip, set, NAFOdiv, bin)) %>%
    mutate(nmeasured = sum(value), .by = c(vessel, trip, set, NAFOdiv)) %>%
    reshape2::dcast(vessel + trip + set + NAFOdiv + nmeasured ~ bin)

  set_spec <- left_join(
    setdat,
    setbio %>% filter(spec == code) # 650 rows
  ) %>%
    left_join(vessel_df) %>%
    left_join(gear_df) %>%
    filter(year == 95) %>%
    arrange(set, trip) %>%
    mutate(
      gearname = factor(gearname, levels = c("Engel", "Campelen")),
      number = ifelse(is.na(number), 0, number),
      doorspread = 1e-3 * ifelse(gearname == "Engel", 25.0403, 31.18805),
      swept.area = dist.towed * doorspread
    ) %>%
    left_join(set_len) %>%
    select(vessel, trip, set, NAFOdiv, set.depth.mean, gearname, swept.area, number, nmeasured, starts_with("bin"))

  trip_pairs <- list(c(5, 253), c(6, 254), c(7, 255))

  set_pair <- lapply(1:length(trip_pairs), function(i) {
    filter(set_spec, trip %in% trip_pairs[[i]]) %>%
      mutate(id = paste0(i, "-", set))
  }) %>%
    bind_rows() %>%
    filter(id != "2-66") %>%
    select(id, gearname, swept.area, starts_with("bin")) %>%
    reshape2::melt(id.vars = c("id", "gearname", "swept.area")) %>%
    mutate(value = ifelse(is.na(value), 0, value)) %>%
    rename(number = value)

  set_pair_pos <- set_pair %>%
    summarise(number = sum(number), .by = c(id, variable)) %>%
    filter(number > 0) %>%
    select(!number)

  set_pair <- left_join(set_pair_pos, set_pair)

  set_pair2 <- set_pair %>%
    summarise(
      number_campelen = number[1],
      number_engel = number[2],
      offset_ratio = swept.area[1]/swept.area[2],
      .by = c(id, variable)
    ) %>%
    filter(number_engel > 0 | number_campelen > 0) %>%
    mutate(p = number_campelen/(number_engel + number_campelen)) %>%
    mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric())

  # Binomial estimator, see Miller 2013, 10.1139/cjfas-2013-0136, p. 1307; Cadigan and Simpson 2023, JNWAFS 54, 1-16
  #dat_full <- lapply(1:nrow(set_pair2), function(i) {
  #  data.frame(
  #    Campelen = c(rep(1, set_pair2$number_campelen[i]), rep(0, set_pair2$number_engel[i]))
  #  ) %>%
  #    mutate(
  #      id = set_pair2$id[i],
  #      offset_ratio = set_pair2$offset_ratio[i],
  #      log_depth = log(set_pair2$depth[i]),
  #      NAFOdiv = set_pair2$NAFOdiv[i]
  #    )
  #}) %>%
  #  bind_rows()

  #m <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ 1,
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = binomial()
  #)
  #
  #m2 <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ 1,
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = betabinomial()
  #)
  #
  #m3 <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ 1 + (1 | id),
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = binomial()
  #)
  #
  #m4 <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0,
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = binomial(),
  #  REML = TRUE
  #)
  #
  #m5 <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0,
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = betabinomial(),
  #  REML = TRUE
  #)
  #
  #m6 <- glmmTMB(
  #  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0 + (1 | id),
  #  data = set_pair2,
  #  offset = log(offset_ratio),
  #  family = binomial(),
  #  REML = TRUE
  #)

  set_pair_sdm <- set_pair2 %>%
    mutate(X = rnorm(nrow(.)), Y = rnorm(nrow(.)), id = factor(id))

  m3 <- sdmTMB(
    cbind(number_campelen, number_engel) ~ 1 + (1 | id),
    mesh = set_pair_sdm %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
    spatial = "off",
    spatiotemporal = "off",
    data = set_pair_sdm,
    offset = log(set_pair_sdm$offset_ratio),
    family = binomial(),
    control = sdmTMBcontrol(multiphase = FALSE)
  )

  m6 <- sdmTMB(
    cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0 + (1 | id),
    mesh = set_pair_sdm %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
    spatial = "off",
    spatiotemporal = "off",
    data = set_pair_sdm,
    offset = log(set_pair_sdm$offset_ratio),
    family = binomial(),
    control = sdmTMBcontrol(multiphase = FALSE)
  )

  extra_bin <- setdiff(1:200, unique(set_pair_sdm$bin))
  m7 <- sdmTMB(
    cbind(number_campelen, number_engel) ~ 0,
    mesh = set_pair_sdm %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
    time = "bin",
    time_varying = ~ 1,
    time_varying_type = "rw0",
    extra_time = extra_bin,
    spatial = "off",
    spatiotemporal = "off",
    data = set_pair_sdm,
    offset = log(set_pair_sdm$offset_ratio),
    family = binomial(),
    control = sdmTMBcontrol(multiphase = FALSE)
  )

  m8 <- sdmTMB(
    cbind(number_campelen, number_engel) ~ 0 + (1 | id),
    mesh = set_pair_sdm %>% make_mesh(xy_cols = c("X", "Y"), n_knots = 10),
    time = "bin",
    time_varying = ~ 1,
    time_varying_type = "rw",
    extra_time = extra_bin,
    spatial = "off",
    spatiotemporal = "off",
    data = set_pair_sdm,
    offset = log(set_pair_sdm$offset_ratio),
    family = binomial(),
    control = sdmTMBcontrol(multiphase = FALSE)
  )

  # For cod:
  # From m3: log_rho = 2.0488, SE = 0.0811

  # For ghl:
  # From m3: log_rho = 0.9107, SE = 0.0638

  AIC(m3, m6, m7, m8)
  sapply(list(m, m2, m3, m4, m5, m6), deviance)

  #resid6 <- simulateResiduals(m6, n = 500)
  #resid6 <- DHARMa::createDHARMa(
  #  simulate(m6, n = 500),
  #  observedResponse = m6$response[, 1],
  #  fittedPredictedResponse = fitted(m6),
  #  integerResponse = TRUE
  #)
  #plot(resid6)
  #testDispersion(resid6)
  #testOutliers(resid6)
  #testResiduals(resid6)
  #testZeroInflation(resid6)

  #resid3 <- simulateResiduals(m3, n = 500)
  #plot(resid3)
  #testDispersion(resid3)
  #testResiduals(resid3)
  #testOverdispersion(resid3)

  #resid8 <- DHARMa::createDHARMa(
  #  simulate(m8, n = 500),
  #  observedResponse = m8$response[, 1],
  #  fittedPredictedResponse = fitted(m8),
  #  integerResponse = TRUE
  #)
  #plot(resid8)
  #testDispersion(resid8)
  #testResiduals(resid8)
  #testOverdispersion(resid8)

  # Make comparison figure and save rho estimates by size
  source("99-functions.R")
  g <- rbind(
    plot_rho_size(m3) %>% getElement("data") %>% mutate(Model = "Constant"),
    plot_rho_size(m6) %>% getElement("data") %>% mutate(Model = "GAM"),
    plot_rho_size(m8) %>% getElement("data") %>% mutate(Model = "RW")
  ) %>%
    mutate(rho = exp(log_rho)) %>%
    ggplot(aes(bin, colour = Model)) +
    geom_line(aes(y = log_rho)) +
    geom_line(aes(y = log_rho + 1.96 * se), linetype = 2) +
    geom_line(aes(y = log_rho - 1.96 * se), linetype = 2) +
    labs(x = "Length", y = expression(log(rho))) +
    coord_cartesian(xlim = if (s == "cod") c(0, 100) else c(0, 75),
                    ylim = if (s == "cod") c(-5, 7) else c(-4, 6))
  .ggsave("rho_comparative.png", sp = s, g, width = 6, height = 3.5)
  readr::write_csv(g$data, file = file.path("validation_size", paste0(s, "_rho_comparative.csv")))

  # Plot length comp by gear - different code compared to above
  set_len <- readr::read_csv("validation_size/CF_lengthfreqs_438_892.csv") %>%
    filter(spec == 438, year == 95) %>%
    filter(length.group == 1) %>%
    select(vessel, trip, set, NAFOdiv, sex, starts_with("bin")) %>%
    reshape2::melt(id.vars = c("vessel", "trip", "set", "NAFOdiv", "sex")) %>%
    mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric()) %>%
    #mutate(bin10 = findInterval(bin, seq(0, 200, 10))) %>%
    #mutate(binmid = paste0("bin", 10 * (bin10 - 1) + 5)) %>%
    summarise(value = sum(value),
              .by = c(vessel, bin)) %>%
    filter(bin <= ifelse(s == "cod", 100, 75)) %>%
    mutate(gear = ifelse(vessel == 20, "Engel", "Campelen"))

  g <- ggplot(set_len, aes(bin, value, colour = gear)) +
    geom_line() +
    coord_cartesian(xlim = c(0, 100)) +
    labs(x = "Length", y = "Frequency", colour = "Gear")
  .ggsave("CF_lengthcomp.png", sp = s, g, width = 6, height = 3.5)

  # plot cumulative length distribution combined from all data
  set_cumulative <- set_len %>%
    summarise(value = sum(value), .by = bin) %>%
    arrange(bin) %>%
    mutate(
      Rel = value/sum(value),
      Cumulative = cumsum(value)/sum(value)
    ) %>%
    select(!value)

  # Interpolate 33rd, 50th, and 66th percentile size
  set_cumulative %>%
    summarise(l33 = approx(Cumulative, bin, xout = 1/3)$y,
              l50 = approx(Cumulative, bin, xout = 0.5)$y,
              l66 = approx(Cumulative, bin, xout = 2/3)$y)

  g <- set_cumulative %>%
    reshape2::melt(id.vars = c("bin")) %>%
    #filter(variable == "Rel") %>%
    mutate(variable = ifelse(variable == "Rel", "Proportion", "Cumulative frequency")) %>%
    ggplot(aes(bin, value)) +
    geom_line() +
    facet_grid(vars(variable), scales = "free_y", switch = "both") +
    coord_cartesian(xlim = c(0, 100)) +
    labs(x = "Length", y = NULL, colour = "Gear") +
    guides(linetype = "none") +
    theme(
      panel.spacing = unit(0, "in"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.placement = "outside"
    ) +
    geom_vline(xintercept = ifelse(s == "cod", 42, 35), linetype = 2)


  # Plot cumulative length distribution by gear
  set_cumulative_gear <- set_len %>%
    arrange(bin) %>%
    mutate(
      Rel = value/sum(value),
      Cumulative = cumsum(value)/sum(value), .by = c(vessel, gear)
    ) %>%
    select(bin, gear, Rel, Cumulative) %>%
    reshape2::melt(id.vars = c("gear", "bin"))

  g <- set_cumulative_gear %>%
    #filter(variable == "Rel") %>%
    mutate(variable = ifelse(variable == "Rel", "Proportion", "Cumulative frequency")) %>%
    ggplot(aes(bin, value, colour = gear)) +
    geom_line() +
    facet_grid(vars(variable), scales = "free_y", switch = "both") +
    coord_cartesian(xlim = c(0, ifelse(s == "cod", 100, 75))) +
    labs(x = "Length", y = NULL, colour = "Gear") +
    guides(linetype = "none") +
    theme(
      panel.spacing = unit(0, "in"),
      strip.background = element_rect(fill = NA, colour = NA),
      strip.placement = "outside"
    )
  .ggsave("CF_cumlengthcomp.png", sp = s, g, width = 5.5, height = 5)

  # Plot set pair frequency by size
  set_pair_fn <- function(i) {
    if (i[1] > 0 && i[2] > 0) {
      return("Both")
    } else if (sum(i) == 0) {
      return("neither")
    } else if(i[1] > 0) {
      return("Campelen")
    } else {
      return("Engel")
    }
  }

  g <- set_pair %>%
    mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric()) %>%
    summarise(type = set_pair_fn(number), .by = c(id, bin)) %>%
    summarise(n = n(), .by = c(bin, type)) %>%
    ggplot(aes(bin, n, fill = type)) +
    coord_cartesian(xlim = c(0, ifelse(s == "cod", 100, 75))) +
    geom_col(width = 1, colour = "grey60", linewidth = 0.25) +
    labs(x = "Length", y = "Paired tow frequency", fill = "Presence")
  ggsave("CF_pairfrequency.png", sp = s, g, width = 6, height = 3.5)
}
