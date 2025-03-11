

library(tidyverse)
library(glmmTMB)

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

spec <- 438
species <- "cod"

for (sp in 1:2) {

  spec <- spec_df$spec[sp]
  species <- spec_df$species[sp]

  set_spec <- left_join(
    setdat,
    setbio %>% filter(.data$spec == .env$spec) # 650 rows for cod, 567 for ghl
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
    )
  readr::write_csv(set_spec, file = paste0("validation/", species, "_set_cf.csv"))

  trip_pairs <- list(c(5, 253), c(6, 254), c(7, 255))

  set_pair_key <- lapply(1:length(trip_pairs), function(i) {
    filter(set_spec, trip %in% trip_pairs[[i]]) %>%
      mutate(id = paste0(i, "-", set)) %>%
      select(id, trip, set, set.depth.mean)
  }) %>%
    bind_rows()

  set_pair <- left_join(set_spec, set_pair_key) %>%
    mutate(npue = number/swept.area) %>%
    filter(id != "2-66")
  #set_pair %>% summarise(n = n(), .by = id) %>% filter(n != 2)
  set_pair_zero <- set_pair %>%
    reshape2::dcast(id ~ gearname, value.var = "npue") %>%
    mutate(ratio = Campelen/Engel) %>%
    filter(Campelen == 0, Engel == 0) %>%
    pull(id)
  set_pair <- set_pair %>%
    filter(!id %in% set_pair_zero)

  set_pair2 <- set_pair %>%
    summarise(
      number_campelen = number[1],
      number_engel = number[2],
      offset_ratio = swept.area[1]/swept.area[2],
      depth = mean(set.depth.mean),
      .by = c(id, NAFOdiv, strat)
    ) %>%
    filter(number_engel > 0 | number_campelen > 0) %>%
    mutate(p = number_campelen/(number_engel + number_campelen))

  # Binomial estimator, see Miller 2013, 10.1139/cjfas-2013-0136, p. 1307; Cadigan and Simpson 2023, JNWAFS 54, 1-16
  dat_full <- lapply(1:nrow(set_pair2), function(i) {
    data.frame(
      Campelen = c(rep(1, set_pair2$number_campelen[i]), rep(0, set_pair2$number_engel[i]))
    ) %>%
      mutate(
        id = set_pair2$id[i],
        offset_ratio = set_pair2$offset_ratio[i],
        log_depth = log(set_pair2$depth[i]),
        NAFOdiv = set_pair2$NAFOdiv[i]
      )
  }) %>%
    bind_rows()

  # Depth effects on paired tows
  g <- ggplot(set_pair, aes(set.depth.mean, npue)) +
    geom_point(aes(colour = gearname))
  g <- ggplot(set_pair2, aes(depth, number_campelen/number_engel)) +
    geom_point()

  m <- glmmTMB(
    Campelen ~ 1,
    data = dat_full,
    offset = log(offset_ratio),
    family = binomial(link = "logit")
  )
  fixef(m)

  # There is no overdispersion!
  m2 <- glmmTMB(
    Campelen ~ 1,
    data = dat_full,
    offset = log(offset_ratio),
    family = betabinomial(link = "logit")
  )

  m3 <- glmmTMB(
    Campelen ~ 1 + poly(log_depth, 2),
    data = dat_full,
    offset = log(offset_ratio),
    family = binomial(link = "logit")
  )

  #### Cod
  # Random effect on individual pairs,
  # Effect size higher, intercept = 2.072, SE = 0.08256
  # 3L only: intercept = 2.5413, SE = 0.1367
  #
  #### GHL
  # Random effect on individual pairs,
  # Effect size higher, intercept = 0.8913, SE = 0.0616
  # 3L only, intercept = 0.9840, SE = 0.1157
  m4 <- glmmTMB(
    Campelen ~ 1 + (1|id),
    data = dat_full,
    offset = log(offset_ratio),
    family = binomial(link = "logit")
  )
  fixef(m4)

  m5 <- glmmTMB(
    Campelen ~ 1 + (1|id),
    data = dat_full,
    offset = log(offset_ratio),
    family = betabinomial(link = "logit")
  )

  AIC(m, m2, m3, m4, m5)

  library(DHARMa)
  mresid <- simulateResiduals(m4)
  residuals(mresid, quantileFunction = qnorm)
  plot(mresid)

  log_rho <- sapply(list(m, m4), function(x) fixef(x)[1] %>% as.numeric()) %>%
    as.data.frame() %>%
    rbind(0) %>%
    mutate(int = 0, Model = c("Intercept-only", "Random effects", "One-to-one")) %>%
    mutate(Model = factor(Model, levels = Model[c(2, 1, 3)]))

  g <- set_pair %>%
    reshape2::dcast(list("id", "gearname"), value.var = "npue") %>%
    ggplot(aes(Engel, Campelen)) +
    geom_point(shape = 1) +
    geom_abline(data = log_rho, aes(slope = exp(.), linetype = Model, intercept = int)) +
    labs(x = "NPUE Engel", y = "NPUE Campelen") +
    theme(legend.position = "bottom") +
    scale_linetype_manual(values = c(1, 2, 4)) +
    coord_cartesian(xlim = c(0, 5000), ylim = if (species == "cod") c(0, 15000) else c(0, 7500),
                    expand = FALSE)
  ggsave(paste0("validation/figures/", species, "_calibration.png"), g, height = 4, width = 4)
}


