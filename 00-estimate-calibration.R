

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


set_spec <- left_join(
  setdat,
  setbio %>% filter(spec == 892) # 650 rows
) %>%
  left_join(vessel_df) %>%
  left_join(gear_df) %>%
  filter(year == 95) %>%
  arrange(set, trip)

g <- set_spec %>%
  summarise(n = n(), .by = c(year, vesselname, gearname, NAFOdiv))
