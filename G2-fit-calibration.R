
library(tidyverse)

#### All data from 0A-0B ----
CAN <- readxl::read_excel(
  file.path("Data 2025", "AllYearsAreasGreenland HalibutData.xlsx")
)

#table(CAN$GearUsed, CAN$Year)
#table(CAN$Year, CAN$StudyArea) %>% t()

#g <- ggplot(CAN, aes(StartLongitude, StartLatitude)) +
#  geom_point(aes(colour = Count_km2)) +
#  facet_wrap(vars(Year)) +
#  scale_colour_viridis_c(trans = 'sqrt')

#### Greenland data since 2023 ----
ghl_TA <- local({
  station <- read.csv(file.path("Data 2025", "D_hauls_GHLW_2023-2024.csv.csv"))
  ghl_catch <- read.csv(file.path("Data 2025", "D_FISH_GHLW_2024-2024.csv"))

  left_join(station, select(ghl_catch, kg, num, id.2), by = "id.2") %>%
    mutate(num = ifelse(is.na(num), 0, num),
           kg = ifelse(is.na(kg), 0, kg))
})
#table(ghl_TA$year, ghl_TA$Stratum1) %>% t()

#### GHL sets prior to 2023 (use 1CD subset) ----
# Missing 2000, 2002, 2003, 2005-7, 2009, 2011, 2013
ghl_PA <- read.csv(
  file.path("Data pre-2025", "D1_CPUE_NPUE_WGHL_99-23_GHLCombine0A1CD6stra0A 35cmGearVinfo.csv.csv")
) %>%
  filter(div %in% c("1C", "1D"), year < 2023)
#table(ghl_PA$year, ghl_PA$div) %>% t()
#table(ghl_PA$year, ghl_PA$div) %>% t() %>% colSums()

# Add sets from missing years
ghl_PA2 <- local({
  # GHL sets in 1CD until 2022 (Paamiut, Helga Maria)
  # Spreadsheet is missing swept area for 2022 Tarajoq sets and catch rates
  ghl_PA2_pos <- read.csv(
    file.path("Data pre-2025", "D_hauls_GHLW_1997-2022.csv")
  ) %>%
    filter(!year %in% unique(ghl_PA$year))
  #table(ghl_PA2$year, ghl_PA$div) %>% t() %>% colSums()

  # Spreadsheet of catch rates provided by Margaret Treble, 2022
  ghl_sets <- read.csv(
    file.path("Data pre-2025", "DS3_Survey 1CD_ kg.num Densitykg Densitynum Per Station 1997-2019 StepTGadget.csv")
  ) %>%
    filter(!year %in% unique(ghl_PA$year)) %>%
    select(year, ship, trip, haul, sweptarea, kg, num) %>%
    left_join(ghl_PA2_pos %>% select(!sweptarea))
  ghl_sets
})

#### Make master spread sheet ----
ghl <- rbind(
  select(
    ghl_PA, year, id.2, sweptarea, depth, div, Stratum2, npos, wpos, num, kg, vessel
  ),
  select(
    ghl_PA2, year, id.2, sweptarea, depth, Stratum1, Stratum2, npos, wpos, num, kg, ship
  ) %>% rename(div = Stratum1, vessel = ship),
  select(
    ghl_TA, year, id.2, sweptarea, depth, Stratum1, Stratum2, npos, wpos, num, kg, ship
  ) %>% rename(div = Stratum1, vessel = ship),
  mutate(
    CAN,
    vessel = ifelse(GearUsed == "Bacalao", "TA", ifelse(Year == 2019, "HM", "PA"))
  ) %>%
    select(
      Year, TrawlID, GCRSSweptArea_km2, Depth_m, StudyArea, OrgStratum, StartLatitude,
      StartLongitude, TotalSetWt_kg, TotalSetCount, vessel
    ) %>%
    rename(
      year = Year,
      id.2 = TrawlID,
      sweptarea = GCRSSweptArea_km2,
      depth = Depth_m,
      div = StudyArea,
      Stratum2 = OrgStratum,
      npos = StartLatitude,
      wpos = StartLongitude,
      kg = TotalSetWt_kg,
      num = TotalSetCount
    )
) %>%
  mutate(NPUE = num/sweptarea, wpos = ifelse(wpos > 0, -1 * wpos, wpos),
         vessel = factor(vessel, levels = c("PA", "HM", "TA")))
#readr::write_csv(ghl, "calibration_ghl/ghl_sets.csv")

#table(ghl$year, ghl$div)

dat <- sdmTMB::add_utm_columns(ghl, ll_names = c("wpos", "npos"), utm_crs = 32621)

source('99-functions-calibration.R')

#### Prepare prediction grid ----
grid <- read.csv("calibration_ghl/grid_0AB1CD.csv") %>%
  sdmTMB::add_utm_columns(ll_names = c("X", "Y"), ll_crs = 32621,
                          utm_names = c("lon", "lat"), utm_crs = 4326, units = "m") %>%
  filter(lat < 73, lat > 61.5) %>%
  mutate(X = X/1000, Y = Y/1000, area = area/1000/1000, vessel = "TA")

grid_ex0B <- filter(grid, div != "0B")

#### Fit calibration Tarajoq + Helga Maria ----
tictoc::tic()
fit <- fit_ghl(dat, formula = num ~ 1 + vessel)
tictoc::toc()
saveRDS(fit, file = "calibration_ghl/models/fit.rds")

fit <- readRDS(file = "calibration_ghl/models/fit.rds")
index <- make_index(fit, grid)
readr::write_csv(index, file = "calibration_ghl/index/index_full.csv")

index_ex0B <- make_index(fit, grid_ex0B)
readr::write_csv(index_ex0B, file = "calibration_ghl/index/index_full_ex0B.csv")


#### Fit calibration ex. HM ----
tictoc::tic()
fit_exHM <- dat %>%
  filter(vessel != "HM") %>%
  mutate(vessel = factor(vessel, c("PA", "TA"))) %>%
  fit_ghl(formula = num ~ 1 + vessel)
tictoc::toc()
saveRDS(fit_exHM, file = "calibration_ghl/models/fit_exHM.rds")

fit_exHM <- readRDS(file = "calibration_ghl/models/fit_exHM.rds")
index <- make_index(fit_exHM, grid)
readr::write_csv(index, file = "calibration_ghl/index/index_exHM.csv")

#### Fit w/o calibration ----
tictoc::tic()
fit_norho <- fit_ghl(dat, formula = num ~ 1)
tictoc::toc()
saveRDS(fit_norho, file = "calibration_ghl/models/fit_norho.rds")

fit_norho <- readRDS(file = "calibration_ghl/models/fit_norho.rds")
index <- make_index(fit_norho, grid)
readr::write_csv(index, file = "calibration_ghl/index/index_norho.csv")


#### Retrospective ex 2024, ex 2023-24 ----
yret <- c(2022, 2023)
library(snowfall)
sfInit(TRUE, length(yret))
sfExport(list = c("fit_ghl", "make_index"))
fit_ret <- sfLapply(yret, function(y, dat, grid) {
  require(dplyr)
  fit_ret <- dat %>%
    filter(year <= y) %>%
    fit_ghl()
  saveRDS(fit_ret, file = paste0("calibration_ghl/models/fit_ret_", y, ".rds"))

  index <- make_index(fit_ret, grid)
  readr::write_csv(index, file = paste0("calibration_ghl/index/index_ret_",  y, ".csv"))
  invisible()
}, dat = dat, grid = grid)

#### Gap ex 2022 or 2023 ----
ygap <- c(2022, 2023)

fit_gap <- sfLapply(ygap, function(y, dat, grid) {
  require(dplyr)
  fit_gap <- dat %>%
    filter(!year %in% seq(2022, y)) %>%
    fit_ghl()
  saveRDS(fit_gap, file = paste0("calibration_ghl/models/fit_gap_", y, ".rds"))

  index <- make_index(fit_gap, grid)
  readr::write_csv(index, file = paste0("calibration_ghl/index/index_gap_",  y, ".csv"))
  invisible()
}, dat = dat, grid = grid)

sfStop()


#### False positive test ----
ytest <- 2013:2017
library(snowfall)
sfInit(TRUE, length(ytest))
sfExport(list = c("fit_ghl"))

fit_null <- sfLapply(ytest, function(y, dat) {
  require(dplyr)
  dat %>%
    filter(vessel == "PA") %>%
    mutate(new_gear = factor(year >= .env$y)) %>%
    fit_ghl(formula = num ~ 1 + new_gear)
}, dat = dat)
sfStop()

saveRDS(fit_null, file = "calibration_ghl/models/fit_null.rds")
