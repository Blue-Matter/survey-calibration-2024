

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

# Cod = 902, Halibut = 438
spec_df <- data.frame(
  spec = c(892, 438),
  species = c("cod", "ghl")
)

set_len <- readr::read_csv("validation_size/CF_lengthfreqs_438_892.csv") %>%
  filter(spec == 892, year == 95) %>%
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
  setbio %>% filter(spec == 892) # 650 rows
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

m <- glmmTMB(
  cbind(number_campelen, number_engel) ~ 1,
  data = set_pair2,
  offset = log(offset_ratio),
  family = binomial(link = "logit")
)

m2 <- glmmTMB(
  cbind(number_campelen, number_engel) ~ 1,
  data = set_pair2,
  offset = log(set_pair2$offset_ratio),
  family = betabinomial(link = "logit")
)

m3 <- glmmTMB(
  cbind(number_campelen, number_engel) ~ 1 + (1 | id),
  data = set_pair2,
  offset = log(set_pair2$offset_ratio),
  family = binomial(link = "logit")
)

m4 <- glmmTMB(
  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0,
  data = set_pair2,
  offset = log(set_pair2$offset_ratio),
  family = binomial(link = "logit")
)

m5 <- glmmTMB(
  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0,
  data = set_pair2,
  offset = log(offset_ratio),
  family = betabinomial(link = "logit")
)

m6 <- glmmTMB(
  cbind(number_campelen, number_engel) ~ s(bin, bs = "cr") + 0 + (1 | id),
  data = set_pair2,
  offset = log(offset_ratio),
  family = binomial(link = "logit")
)

AIC(m, m2, m3, m4, m5, m6)

plot_rho <- function(m, bin = seq(0, 200, 1)) {
  est <- data.frame(bin = bin, offset_ratio = 1)
  pred <- predict(m, est, re.form = NA, se.fit = TRUE, cov.fit = TRUE)

  est$log_rho <- pred$fit
  est$se <- pred$se.fit
  g <- est %>%
    mutate(rho = exp(log_rho)) %>%
    ggplot(aes(bin)) +
    geom_line(aes(y = log_rho)) +
    geom_line(aes(y = log_rho + 1.96 * se), linetype = 2) +
    geom_line(aes(y = log_rho - 1.96 * se), linetype = 2) +
    labs(x = "Length", y = expression(log(rho)))
  g
}
g <- plot_rho(m6) +
  coord_cartesian(xlim = c(0, 75), ylim = c(-4, 3))
ggsave("validation_size/figures/cod_rho_comparative.png", g, width = 5, height = 3.5)

# Plot length comp
set_len <- readr::read_csv("validation_size/CF_lengthfreqs_438_892.csv") %>%
  filter(spec == 892, year == 95) %>%
  filter(length.group == 1) %>%
  select(vessel, trip, set, NAFOdiv, sex, starts_with("bin")) %>%
  reshape2::melt(id.vars = c("vessel", "trip", "set", "NAFOdiv", "sex")) %>%
  mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric()) %>%
  #mutate(bin10 = findInterval(bin, seq(0, 200, 10))) %>%
  #mutate(binmid = paste0("bin", 10 * (bin10 - 1) + 5)) %>%
  summarise(value = sum(value),
            .by = c(vessel, bin)) %>%
  filter(bin < 75) %>%
  mutate(gear = ifelse(vessel == 20, "Engel", "Campelen"))

g <- ggplot(set_len, aes(bin, value, colour = gear)) +
  geom_line() +
  labs(x = "Length", y = "Frequency", colour = "Gear")
ggsave("validation_size/figures/cod_CF_lengthcomp.png", g, width = 6, height = 3.5)


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
  summarise(x = set_pair_fn(number), .by = c(id, bin)) %>%
  summarise(n = n(), .by = c(bin, x)) %>%
  ggplot(aes(bin, n, colour = x)) +
  geom_line()

g <- set_pair %>%
  mutate(bin = strsplit(as.character(.data$variable), "bin") %>% sapply(getElement, 2) %>% as.numeric()) %>%
  summarise(type = set_pair_fn(number), .by = c(id, bin)) %>%
  summarise(n = n(), .by = c(bin, type)) %>%
  ggplot(aes(bin, n, fill = type)) +
  geom_col(width = 1, colour = "grey60", linewidth = 0.25) +
  labs(x = "Length", y = "Paired tow frequency", fill = "Presence")
ggsave("validation_size/figures/cod_CF_pairfrequency.png", g, width = 6, height = 3.5)




