
library(tidyverse)

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

set_spec <- left_join(
  setdat,
  setbio %>% filter(spec == 438) # 567 rows
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
  #filter((Campelen == 0 & Engel == 0) | Engel > 10000) %>%
  pull(id)
set_pair <- set_pair %>%
  filter(!id %in% set_pair_zero)

set_pair2 <- set_pair %>%
  summarise(
    number_campelen = number[1],
    number_engel = number[2],
    offset_ratio = swept.area[1]/swept.area[2],
    depth = mean(set.depth.mean),
    .by = c(id, NAFOdiv)
  ) %>%
  filter(number_engel > 0 | number_campelen > 0) %>%
  mutate(p = number_campelen/(number_engel + number_campelen))

# Binomial estimator, see Miller 2013, 10.1139/cjfas-2013-0136, p. 1307; Cadigan and Simpson 2023, JNWAFS 54, 1-16
dat_full <- lapply(1:nrow(set_pair2), function(i) {
  data.frame(
    Campelen = c(rep(1, set_pair2$number_campelen[i]), rep(0, set_pair2$number_engel[i]))
  ) %>%
    mutate(offset_ratio = set_pair2$offset_ratio[i], log_depth = log(set_pair2$depth[i]))
}) %>%
  bind_rows()

# Depth effects on paired tows
g <- ggplot(set_pair, aes(set.depth.mean, npue)) +
  geom_point(aes(colour = gearname))
g <- ggplot(set_pair2, aes(depth, number_campelen/number_engel)) +
  geom_point()

library(glmmTMB)
m <- glmmTMB(
  Campelen ~ 1,
  data = dat_full,
  offset = log(dat_full$offset_ratio),
  family = binomial(link = "logit")
)
exp(1.44256)

m2 <- glmmTMB(
  Campelen ~ 1 + poly(log_depth, 2),
  data = dat_full,
  offset = log(dat_full$offset_ratio),
  family = binomial(link = "logit")
)
AIC(m, m2)

library(DHARMa)
mresid <- simulateResiduals(m)
residuals(mresid, quantileFunction = qnorm)
plot(mresid)

#m2 <- glmmTMB(
#  Campelen ~ 1,
#  data = dat_full,
#  offset = log(dat_full$offset_ratio),
#  family = betabinomial(link = "logit")
#)

g <- set_pair %>%
  reshape2::dcast(list("id", "gearname"), value.var = "npue") %>%
  ggplot(aes(Engel, Campelen)) +
  geom_point(shape = 1) +
  geom_abline(linetype = 2, slope = 1, intercept = 0) +
  geom_abline(slope = exp(1.44256), intercept = 0) +
  geom_label(x = Inf, y = Inf, hjust = "inward", vjust = "inward", label = expression(rho==exp(1.44256))) +
  labs(x = "NPUE Engel", y = "NPUE Campelen") +
  coord_cartesian(xlim = c(0, 5000), ylim = c(0, 30000))
ggsave("validation/figures/ghl_calibration.png", g, height = 3, width = 4)
