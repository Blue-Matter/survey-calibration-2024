
library(tidyverse)
theme_set(theme_bw())

Names <- c("NL Cod", "NL G. halibut", "NE Butterfish", "NE Redfish", "NE W. hake", "NE Y. flounder")

# Rho from CF
case_studies <- data.frame(
  code = c("cod", "ghl", "butterfish", "redfish", "whitehake", "yt"),
  Name = factor(Names, levels = Names),
  log_rho = c(2.072, 0.8913, log(c(1.487, 1.4564, 2.2354, 2.3619))),
  se = c(0.08256, 0.0616, 0.1475, 0.08936, 0.07723, 0.07583),
  lower = c(2.072 - 1.96 * 0.08256, 0.8913 - 1.96 * 0.0616, log(c(1.1137, 1.2224, 1.9214, 2.0357))),
  upper = c(2.072 + 1.96 * 0.08256, 0.8913 + 1.96 * 0.0616, log(c(1.9862, 1.7351, 2.6008, 2.7403)))
) %>%
  mutate(dir = ifelse(grepl("NL", Name), "validation", "validation_nefsc"),
         year_new_vessel = ifelse(grepl("NL", Name), 1995, 2009))
readr::write_csv(case_studies, "tables/rho_CF.csv")

# SE
#(log(c(1.487, 1.4564, 2.2354, 2.3619)) - log(c(1.1137, 1.2224, 1.9214, 2.0357)))/1.96
#(log(c(1.487, 1.4564, 2.2354, 2.3619)) - log(c(1.9862, 1.7351, 2.6008, 2.7403)))/1.96

log_rho_df <- lapply(1:nrow(case_studies), function(i) {
  rbind(
    case_studies[i, c("code", "Name", "log_rho", "lower", "upper")],
    case_studies[i, c("code", "Name", "log_rho", "lower", "upper")]
  ) %>%
    mutate(x = c(-Inf, Inf))
}) %>%
  bind_rows()


# Summary model objects, including rho_ST
filename <- c("validation/models/cod_fit_full.rds", "validation/models/ghl_fit_full.rds",
              "validation_nefsc/models/butterfish_fit_full.rds", "validation_nefsc/models/redfish_fit_full.rds",
              "validation_nefsc/models/whitehake_fit_full.rds", "validation_nefsc/models/yt_fit_full.rds")

mod <- lapply(filename, readRDS)

# Plot set data
NPUE <- lapply(1:length(mod), function(i) {
  if (grepl("NL", Names[i])) {
    mutate(mod[[i]][["data"]], NPUE = number/swept.area) %>%
      select(NPUE, lat.start, long.start, year) %>%
      rename(LAT = lat.start, LON = long.start) %>%
      mutate(Name = Names[i])
  } else {
    mutate(mod[[i]][["data"]], NPUE = CATCH/EFFORT) %>%
      select(NPUE, LAT, LON, year) %>%
      mutate(Name = Names[i])
  }
}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, Names))

#reshape2::dcast(NPUE, year ~ Name, fun.aggregate = length)

area <- data.frame(
  Name = Names,
  Area = c(391422, 391422, 282336, 68194, 125426, 125426)
)
readr::write_csv(area, file = "tables/survey_area.csv")

nfld <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
  sf::st_crop(xmin = -59, xmax = -46, ymin = 45, ymax = 56)
coast <- lapply(1:length(mod), function(i) {
  require(sf)
  if (grepl("NL", Names[i])) {
    nfld <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf") %>%
      st_make_valid() %>%
      sf::st_crop(xmin = -59, xmax = -46, ymin = 45, ymax = 56)
    return(nfld)
  } else {
    xrange <- range(mod[[i]]$data$LON, na.rm = TRUE) + c(-0.5, 0.5)
    yrange <- range(mod[[i]]$data$LAT, na.rm = TRUE) + c(-0.5, 0.5)

    rnaturalearth::ne_countries(scale = 10, returnclass = "sf") %>%
      st_make_valid() %>%
      sf::st_crop(xmin = xrange[1], xmax = xrange[2], ymin = yrange[1], ymax = yrange[2])
  }
})

glist <- lapply(1:length(mod), function(i) {
  g <- NPUE %>%
    filter(NPUE > 0, Name == Names[i]) %>%
    ggplot(aes(LON, LAT)) +
    geom_point(size = 0.25, data = filter(NPUE, NPUE == 0, Name == Names[i]), shape = 4, colour = "grey80") +
    geom_point(size = 0.25, aes(colour = log(NPUE))) +
    geom_sf(data = coast[[i]], inherit.aes = FALSE) +
    coord_sf() +
    scale_colour_viridis_c() +
    labs(x = NULL, y = NULL, title = Names[i],
         subtitle = parse(text = paste0("\"", format(area$Area[i], big.mark = ","), "\"~km^2")) |> eval()) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  g
})


library(patchwork)
g <- glist[[1]] + glist[[2]] + glist[[3]] + glist[[4]] + glist[[5]] + glist[[6]] +
  plot_layout(ncol = 2, nrow = 3, axis_titles = "collect")
ggsave("manuscript/figures/NPUE_casestudy.png", g, height = 8, width = 7)


#lapply(mod, sdmTMB::tidy)
rho_ind <- c(4, 4, 2, 2, 2, 2)
pars <- lapply(1:length(mod), function(i) {
  j <- summary(mod[[i]]$sd_report, select = c("fixed", "report"))
  p <- rownames(j)

  df <- j %>%
    as.data.frame() %>%
    mutate(Name = Names[i], parameter = p)

  log_rho <- filter(df, parameter == "b_j")[c(1, rho_ind[i]), ] %>%
    mutate(parameter = c("Intercept", "log(rho)"))

  rbind(
    df %>% filter(parameter != "b_j"),
    log_rho
  )

}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, Names)) %>%
  filter(parameter %in% c("sigma_O", "sigma_E", "range", "phi", "Intercept", "log(rho)")) %>%
  summarise(Estimate = mean(Estimate), `Std. Error` = mean(`Std. Error`), .by = c(Name, parameter)) %>%
  mutate(lower = Estimate - 1.96 * `Std. Error`, upper = Estimate + 1.96 * `Std. Error`)

# Write table of rho_ST
readr::write_csv(
  pars %>% filter(parameter == "log(rho)"),
  "tables/rho_ST.csv"
)

# Write table of parameters (range, phi, sigma_E, sigma_O, etc.)
readr::write_csv(pars, "tables/ST_par.csv")

# Plot comparison of log rho
rho_CF <- readr::read_csv(here::here("tables/rho_CF.csv"))
rho_ST <- readr::read_csv(here::here("tables/rho_ST.csv"))

rho <- rbind(
  rho_CF %>% select(Name, log_rho, lower, upper) %>% mutate(type = "CF"),
  rho_ST %>% select(Name, Estimate, lower, upper) %>% rename(log_rho = Estimate) %>% mutate(type = "ST")
) %>%
  mutate(Name = factor(Name, Names))

g <- ggplot(rho, aes(type, log_rho, ymin = lower, ymax = upper)) +
  facet_wrap(vars(Name), scales = "free") +
  geom_point() +
  geom_linerange() +
  coord_flip() +
  theme(strip.background = element_blank()) +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(x = NULL, y = expression(Estimate%+-%1.96~"SE")) +
  scale_x_discrete(labels = expression(log(rho[CF]), log(rho[ST])))
ggsave("manuscript/figures/rho_casestudy.png", g, height = 3, width = 6)

# Compare indices (with and without calibration)
index <- lapply(1:nrow(case_studies), function(i) {
  s <- case_studies$code[i]

  suppressMessages({
    index <- readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index.csv")))
    index_rho <- readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_calibration.csv")))
    index_norho <- readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_norho.csv")))
  })

  index_compare <- rbind(
    index_rho %>% mutate(Model = "Comparative rho"),
    index %>% mutate(Model = "Estimate rho"),
    index_norho %>% mutate(Model = "No rho")
  ) %>%
    mutate(Name = case_studies$Name[i])

  return(index_compare)
}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, levels = Names))

lab <- expression("(1) Direct"~rho[CF], "(2) Indirect"~rho[ST], "(3) None")

g <- index %>%
  #mutate(est = est, .by = Name) %>%
  ggplot(aes(year, log(est), linetype = Model)) +
  #geom_point() +
  geom_line() +
  geom_vline(data = case_studies, aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  labs(x = "Year", y = "log(Index)", linetype = "Calibration") +
  facet_wrap(vars(Name), scales = "free", ncol = 2) +
  #expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1, 4)) +
  scale_linetype_manual(values = c(1, 2, 3), labels = lab) +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_compare.png", g, height = 6, width = 6)

g <- index %>%
  #mutate(est = est, .by = Name) %>%
  ggplot(aes(year, log(est), ymin = log(lwr), ymax = log(upr), shape = Model, colour = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_line(linewidth = 0.25, linetype = 3) +
  geom_linerange(position = position_dodge(width = 0.5), linewidth = 0.25) +
  geom_vline(data = case_studies, aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  labs(x = "Year", y = "log(Index)", colour = "Calibration", shape = "Calibration") +
  facet_wrap(vars(Name), scales = "free", ncol = 2) +
  #expand_limits(y = 0) +
  scale_colour_discrete(labels = lab) +
  scale_shape_manual(values = c(16, 1, 4), labels = lab) +
  theme(legend.position = "bottom", strip.background = element_blank(), panel.grid.minor = element_blank())
ggsave("manuscript/figures/index_compare_ci.png", g, height = 6, width = 6)


# Retrospective
index_ret <- lapply(1:nrow(case_studies), function(i) {

  s <- case_studies$code[i]
  readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_ret.csv"))) %>%
    mutate(Name = case_studies$Name[i])

}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, levels = Names))

Mohn_denom <- index_ret %>%
  filter(year_ret == max(year_ret), .by = Name) %>%
  select(year, est, Name) %>%
  rename(denom = est)

Mohn_diagnostic <- index_ret %>%
  left_join(Mohn_denom) %>%
  filter(year == year_ret) %>%
  mutate(stat = est/denom - 1) %>%
  summarise(stat = mean(stat), .by = Name) %>%
  mutate(x = ifelse(Name %in% c("NL Cod", "NE W. hake"), Inf, -Inf))

g <- index_ret %>%
  select(!se) %>%
  left_join(case_studies) %>%
  mutate(npeel = year_ret - year_new_vessel + 1) %>%
  ggplot(aes(year, log_est)) +
  geom_line(aes(colour = factor(npeel))) +
  geom_vline(data = case_studies, aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  geom_label(data = Mohn_diagnostic, aes(x = x, label = round(stat, 2)),
             y = Inf, hjust = "inward", vjust = "inward") +
  facet_wrap(vars(Name), scales = "free", ncol = 2) +
  labs(x = "Year", y = "log(Index)", colour = "Elapsed years") +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_retrospective.png", g, height = 6, width = 5)


# Retrospective (gear effect)
vessel_effect_ret <- lapply(1:nrow(case_studies), function(i) {

  source("99-functions.R")

  s <- case_studies$code[i]

  if (grepl("NL", case_studies$Name[i])) {
    yret <- 1995:2000
    par_no <- 3
  } else {
    yret <- 2009:2014
    par_no <- 1
  }

  fit_ret <- lapply(yret, function(y) {
    readRDS(file.path(case_studies$dir[i], "models", paste0(s, "_ret_fit_", y, ".rds")))
  })

  vessel_effect <- get_vessel_effect(fit_ret, par_no = par_no) %>%
    mutate(year_ret = .env$yret[.data$sequence], Name = case_studies$Name[i])

  return(vessel_effect)
}) %>%
  bind_rows()


g <- vessel_effect_ret %>%
  left_join(case_studies) %>%
  mutate(npeel = year_ret - year_new_vessel + 1) %>%
  ggplot(aes(factor(npeel), Estimate)) +
  geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE) +
  #geom_hline(aes(yintercept = lower), linetype = 3) +
  #geom_hline(aes(yintercept = upper), linetype = 3) +
  geom_hline(aes(yintercept = log_rho), linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  #expand_limits(y = 0) +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  labs(x = "Elapsed years", y = expression(widehat(log(rho)))) +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_retrospective.png", g, height = 4, width = 4)

# Check if rho trend is significantly different from zero?
lapply(Names, function(i) {
  dat <- dplyr::filter(g@data, Name == i)
  fit <- lm(Estimate ~ year_ret, data = dat)
  summary(fit)
}) %>%
  structure(names = Names)


# Gap (index)
index_gap <- lapply(1:nrow(case_studies), function(i) {
  s <- case_studies$code[i]

  readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_gap.csv"))) %>%
    mutate(Name = case_studies$Name[i])
}) %>%
  bind_rows()

g <- index_gap %>%
  select(!se) %>%
  left_join(case_studies) %>%
  ggplot(aes(year, log(est))) +
  geom_point(data = dplyr::filter(index_gap, impute), shape = 1, aes(colour = factor(gap))) +
  geom_line(linetype = 2, aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = c(gap, Name)), aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = c(gap, Name)), aes(colour = factor(gap))) +
  #expand_limits(y = 0) +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  geom_vline(aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  #scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
  guides(shape = "none") +
  labs(x = "Year", y = "log(Index)", colour = "Time gap") +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_gap.png", g, height = 6, width = 5)

# Gap (gear effect)
vessel_effect_gap <- lapply(1:nrow(case_studies), function(i) {

  source("99-functions.R")

  s <- case_studies$code[i]

  if (grepl("NL", case_studies$Name[i])) {
    par_no <- 3
  } else {
    par_no <- 1
  }
  fit_gap <- lapply(1:5, function(y) {
    readRDS(file.path(case_studies$dir[i], "models", paste0(s, "_gap_fit_", y, ".rds")))
  })

  vessel_effect <- get_vessel_effect(fit_gap, par_no = par_no) %>%
    mutate(Name = case_studies$Name[i])

  return(vessel_effect)
}) %>%
  bind_rows()

g <- vessel_effect_gap %>%
  left_join(case_studies) %>%
  ggplot(aes(sequence, Estimate)) +
  geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE) +
  #geom_hline(yintercept = log_rho_comparative$log_lower, linetype = 3) +
  #geom_hline(yintercept = log_rho_comparative$log_upper, linetype = 3) +
  geom_hline(aes(yintercept = log_rho), linetype = 2) +
  geom_point() +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Time gap between surveys (years)", y = expression(widehat(log(rho)))) +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_gap.png", g, height = 4, width = 4)

# Check if rho trend is significantly different from zero?
lapply(Names, function(i) {
  dat <- dplyr::filter(g@data, Name == i)
  fit <- lm(Estimate ~ sequence, data = dat)
  summary(fit)
}) %>%
  structure(names = Names)

# Null test (index)
index_null <- lapply(1:nrow(case_studies), function(i) {
  s <- case_studies$code[i]

  if (grepl("NL", case_studies$Name[i])) {
    ydummy <- c(2005:2010, NA)
  } else {
    ydummy <- c(2015:2020, NA)
  }

  readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_null.csv"))) %>%
    mutate(Name = case_studies$Name[i],
           ydummy = rep(.env$ydummy, each = length(seq(case_studies$year_new_vessel[i], max(.data$year), 1))))
}) %>%
  bind_rows() %>%
  filter(!is.na(ydummy)) %>%
  mutate(area = grepl("NL", Name)) %>%
  mutate(ydummy2 = ydummy - min(ydummy) + 1, .by = area)

g <- index_null %>%
  ggplot(aes(year, log(est), colour = factor(ydummy2))) +
  geom_point(aes(shape = factor(year == ydummy))) +
  geom_line(aes(colour = factor(ydummy2))) +
  #expand_limits(y = 0) +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = NULL)) +
  labs(x = "Year", y = "log(Index)", colour = "Dummy year change") +
  guides(shape = "none") +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_dummy.png", g, height = 7, width = 6)

# Null test (vessel effect)
vessel_effect_null <- lapply(1:nrow(case_studies), function(i) {

  source("99-functions.R")

  s <- case_studies$code[i]

  if (grepl("NL", case_studies$Name[i])) {
    ydummy <- c(2005:2010, NA)
    par_no <- 3
  } else {
    ydummy <- c(2015:2020, NA)
    par_no <- 1
  }
  fit_null <- lapply(ydummy, function(y) {
    readRDS(file.path(case_studies$dir[i], "models", paste0(s, "_null_fit_", y, ".rds")))
  })

  vessel_effect <- get_vessel_effect(fit_null, par_no = par_no) %>%
    mutate(Name = case_studies$Name[i])

  return(vessel_effect)
}) %>%
  bind_rows()

g <- vessel_effect_null %>%
  left_join(case_studies) %>%
  mutate(ydummy = ifelse(grepl("NL", Name), 2005, 2015) + sequence - 1) %>%
  ggplot(aes(factor(ydummy), Estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Dummy year change", y = expression(widehat(log(rho)))) +
  theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5))
ggsave("manuscript/figures/rho_null.png", g, height = 4.5, width = 3.5)

# Plot diagnostics in one very large figure
g1 <- vessel_effect_null %>%
  left_join(case_studies) %>%
  mutate(ydummy = ifelse(grepl("NL", Name), 2005, 2015) + sequence - 1) %>%
  ggplot(aes(factor(ydummy), Estimate)) +
  geom_hline(yintercept = 0, linetype = 2) +
  geom_point() +
  facet_wrap(vars(Name), ncol = 1, scales = "free", strip.position = "right") +
  coord_flip() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Dummy year change", y = expression(widehat(log(rho))), subtitle = "True negative test") +
  theme(strip.background = element_blank())

g2 <- vessel_effect_ret %>%
  left_join(case_studies) %>%
  mutate(npeel = year_ret - year_new_vessel + 1) %>%
  ggplot(aes(factor(npeel), Estimate)) +
  geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE) +
  geom_hline(aes(yintercept = log_rho), linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  #expand_limits(y = 0) +
  coord_flip() +
  facet_wrap(vars(Name), ncol = 1, scales = "free", strip.position = "right") +
  labs(x = "Elapsed years", y = expression(widehat(log(rho))), subtitle = "Retrospective analysis") +
  theme(strip.background = element_blank())

g <- ggpubr::ggarrange(g1, g2, ncol = 2)
ggsave("manuscript/figures/diagnostics.png", g, height = 7, width = 4)

# Separate figure for survey gap
g <- vessel_effect_gap %>%
  left_join(case_studies) %>%
  ggplot(aes(sequence, Estimate)) +
  geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE) +
  geom_hline(aes(yintercept = log_rho), linetype = 2) +
  geom_point() +
  facet_wrap(vars(Name), ncol = 2, scales = "free", strip.position = "right") +
  coord_flip() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  labs(x = "Years between surveys", y = expression(widehat(log(rho))), subtitle = "Time gap evaluation") +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_gap.png", g, height = 4, width = 4)

# Try to use facets but doesn't seem to work well
if (FALSE) {
  df <- rbind(
    vessel_effect_ret %>%
      left_join(case_studies) %>%
      mutate(npeel = year_ret - year_new_vessel + 1) %>%
      mutate(x = npeel, column = "Elapsed years") %>%
      select(x, Estimate, lwr, upr, column, Name),
    vessel_effect_gap %>%
      left_join(case_studies) %>%
      mutate(x = sequence, column = "Years between surveys") %>%
      select(x, Estimate, lwr, upr, column, Name),
    vessel_effect_null %>%
      left_join(case_studies) %>%
      mutate(ydummy = ifelse(grepl("NL", Name), 2005, 2015) + sequence - 1) %>%
      mutate(x = factor(ydummy), column = "Dummy year change") %>%
      select(x, Estimate, lwr, upr, column, Name)
  )

  g <- df %>%
    ggplot(aes(x, Estimate)) +
    geom_point() +
    facet_grid(vars(Name), vars(column), scales = "free") +
    geom_linerange(aes(ymin = lwr, max = upr)) +
    coord_flip() +
    labs(x = NULL, y = expression(widehat(log(rho)))) +
    theme(strip.background = element_blank())
}


# Summary model objects
filename <- c("validation/models/cod_fit_full.rds", "validation/models/ghl_fit_full.rds",
              "validation_nefsc/models/butterfish_fit_full.rds", "validation_nefsc/models/redfish_fit_full.rds",
              "validation_nefsc/models/whitehake_fit_full.rds", "validation_nefsc/models/yt_fit_full.rds")

mod <- lapply(filename, readRDS)

lapply(mod, sdmTMB::tidy)
rho_ind <- c(4, 4, 2, 2, 2, 2)
pars <- lapply(1:length(mod), function(i) {
  j <- summary(mod[[i]]$sd_report, select = c("fixed", "report"))
  p <- rownames(j)

  df <- j %>%
    as.data.frame() %>%
    mutate(Name = Names[i], parameter = p)

  log_rho <- filter(df, parameter == "b_j")[c(1, rho_ind[i]), ] %>%
    mutate(parameter = c("Intercept", "log(rho)"))

  rbind(
    df %>% filter(parameter != "b_j"),
    log_rho
  )

}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, Names)) %>%
  filter(parameter %in% c("sigma_O", "sigma_E", "range", "phi", "Intercept", "log(rho)")) %>%
  summarise(Estimate = mean(Estimate), `Std. Error` = mean(`Std. Error`), .by = c(Name, parameter)) %>%
  #mutate(label = paste0(signif(Estimate, 3), " (", signif(`Std. Error`, 3), ")"))
  mutate(label = ifelse(parameter %in% c("Intercept", "log(rho)"),
                        paste0(signif(Estimate, 3), " (", signif(`Std. Error`, 3), ")"),
                        signif(Estimate, 3)))

pars %>%
  reshape2::dcast(list("parameter", "Name"), value.var = "label")

pars %>%
  reshape2::dcast(list("parameter", "Name"), value.var = "Estimate")

# Compare spatial grid
pred_grid <- lapply(1:length(mod), function(i) {

  require(sf)

  if (i < 3) {
    dat_fit <- mod[[i]]$data
    grid2 <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
      rename(X_m = X, Y_m = Y) %>%
      mutate(X = 1e-3 * X_m, Y = 1e-3 * Y_m,
             depth = -1 * depth) %>%
      filter(.data$depth > 0) %>%
      filter(.data$div %in% unique(dat_fit$NAFOdiv)) %>%
      mutate(log_depth = log(depth),
             data.series = "Campelen")
  } else {
    grid <- readr::read_csv("shapefiles/prediction_grid_nefsc.csv")

    area <- sf::st_read(dsn = system.file("extdata", "EPU.shp", package = "survdat"), quiet = TRUE)

    grid_sf <- sf::st_as_sf(grid, coords = c("LON", "LAT"), crs = "NAD83")

    if (i == 3) { # Butterfish

      grid_intersect <- sf::st_intersects(
        st_make_valid(area),
        grid_sf
      ) %>% unlist()

    } else if (i == 4) { # Redfish

      grid_intersect <- sf::st_intersects(
        area[area$EPU == "GOM", ],
        grid_sf
      ) %>% unlist()

    } else {

      grid_intersect <- sf::st_intersects(
        area[area$EPU %in% c("GOM", "GB"), ],
        grid_sf
      ) %>% unlist()

    }

    grid2 <- grid[grid_intersect, ]

  }

  pred <- predict(mod[[i]], newdata = grid2 %>% mutate(year = max(mod[[i]]$data$year))) %>%
    mutate(Name = Names[i])

}) %>%
  bind_rows()

g <- pred_grid %>%
  mutate(Name = factor(Name, Names)) %>%
  #mutate(value = omega_s) %>%
  mutate(value = pmin(omega_s, 5)) %>%
  ggplot(aes(X, Y, fill = value)) +
  geom_tile(height = 10, width = 10) +
  facet_wrap(vars(Name), scales = "free") +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = 0)

max_value <- 3
#map_range <- lapply(1:length(Names), function(i) {
#
#  val <- pred_grid %>%
#    summarise(xmin = min(X_m), xmax = max(X_m), ymin = min(Y_m), ymax = max(Y_m), .by = Name) %>%
#    filter(Name == Names[i])
#
#  crs <- ifelse(grepl("NE", Names[i]), 32619, 32622)
#
#  coast <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf") %>%
#    st_transform(crs = crs) #%>%
#    #sf::st_crop(xmin = val$xmin, xmax = val$xmax, ymin = val$ymin, ymax = val$ymax)
#
#  return(coast)
#
#})
#
#gg <- lapply(1:length(Names), function(i) {
#  crs <- ifelse(grepl("NE", Names[i]), 32619, 32622)
#
#  df <- pred_grid %>%
#    mutate(Name = factor(Name, Names)) %>%
#    mutate(value = (est - mean(est))/sd(est), .by = c(Name)) %>%
#    mutate(value = pmin(value, max_value) %>% pmax(-max_value, value)) %>%
#    filter(Name == Names[i])
#
#  coast <- rnaturalearth::ne_coastline(scale = 10, returnclass = "sf")
#
#  g <- df %>%
#    ggplot() +
#    geom_sf(data = coast) +
#    coord_sf(expand = FALSE, crs = crs, xlim = range(df$X_m), ylim = range(df$Y_m)) +
#    geom_tile(aes(X_m, Y_m, fill = value), height = 10000, width = 10000) +
#    facet_wrap(vars(Name)) +
#    scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = 0, limits = c(-1, 1) * max_value)
#  g
#})
#g <- ggpubr::ggarrange(plotlist = g, ncol = 2, nrow = 3)
#ggsave("manuscript/figures/density_casestudies.png", g, height = 10, width = 10)



max_value <- 2
g <- pred_grid %>%
  mutate(Name = factor(Name, Names)) %>%
  mutate(value = (est - mean(est))/sd(est), .by = c(Name)) %>%
  mutate(value = pmin(value, max_value) %>% pmax(-max_value)) %>%
  ggplot(aes(X, Y, fill = value)) +
  geom_tile(height = 10, width = 10) +
  facet_wrap(vars(Name), scales = "free") +
  scale_fill_gradient2(high = scales::muted("red"), low = scales::muted("blue"), mid = 0, limits = c(-1, 1) * max_value) +
  theme(panel.grid = element_blank(), legend.position = "bottom") +
  labs(fill = "Relative density")
ggsave("manuscript/figures/density_casestudies.png", g, height = 4, width = 6)

#tmap package tm_shape(FRA_sf) +
#tm_borders() +
#  tm_facets(by = "NAME_1")

