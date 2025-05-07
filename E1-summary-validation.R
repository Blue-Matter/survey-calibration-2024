
library(tidyverse)
theme_set(theme_bw())

Names <- c("NL Cod", "NL G. halibut", "NE butterfish", "NE Redfish", "NE W. hake", "NE Y. flounder")
case_studies <- data.frame(
  code = c("cod", "ghl", "butterfish", "redfish", "whitehake", "yt"),
  Name = factor(Names, levels = Names),
  log_rho = c(2.072, 0.8913, log(c(1.487, 1.4564, 2.2354, 2.3619))),
  lower = c(2.072 - 1.96 * 0.08256, 0.8913 - 1.96 * 0.0616, log(c(1.1137, 1.2224, 1.9214, 2.0357))),
  upper = c(2.072 + 1.96 * 0.08256, 0.8913 + 1.96 * 0.0616, log(c(1.9862, 1.7351, 2.6008, 2.7403)))
) %>%
  mutate(dir = ifelse(grepl("NL", Name), "validation", "validation_nefsc"),
         year_new_vessel = ifelse(grepl("NL", Name), 1995, 2009))

log_rho_df <- lapply(1:nrow(case_studies), function(i) {
  rbind(
    case_studies[i, c("code", "Name", "log_rho", "lower", "upper")],
    case_studies[i, c("code", "Name", "log_rho", "lower", "upper")]
  ) %>%
    mutate(x = c(-Inf, Inf))
}) %>%
  bind_rows()

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

g <- index %>%
  #mutate(est = est, .by = Name) %>%
  ggplot(aes(year, log(est), linetype = Model, shape = Model)) +
  #geom_point() +
  geom_line() +
  geom_vline(data = case_studies, aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  labs(x = "Year", y = "log(Index)") +
  facet_wrap(vars(Name), scales = "free", ncol = 2) +
  #expand_limits(y = 0) +
  scale_shape_manual(values = c(16, 1, 4)) +
  scale_linetype_manual(values = c(1, 2, 3)) +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_compare.png", g, height = 7, width = 6)


# Retrospective
index_ret <- lapply(1:nrow(case_studies), function(i) {

  s <- case_studies$code[i]
  readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_ret.csv"))) %>%
    mutate(Name = case_studies$Name[i])

}) %>%
  bind_rows() %>%
  mutate(Name = factor(Name, levels = Names))

g <- index_ret %>%
  left_join(case_studies) %>%
  mutate(npeel = year_ret - year_new_vessel) %>%
  #mutate(est = est/mean(est), .by = Name) %>%
  ggplot(aes(year, log(est), colour = factor(npeel))) +
  #geom_point() +
  geom_line() +
  #expand_limits(y = 0) +
  geom_vline(data = case_studies, aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  facet_wrap(vars(Name), scales = "free", ncol = 2) +
  labs(x = "Year", y = "log(Index)", colour = "Elapsed years") +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_retrospective.png", g, height = 7, width = 6)


# Retrospective (gear effect)
library(snowfall)
sfInit(TRUE, nrow(case_studies))
sfExport(list = c("case_studies"))
sfLibrary(tidyverse)
vessel_effect_ret <- sfLapply(1:nrow(case_studies), function(i) {

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
  mutate(npeel = year_ret - year_new_vessel) %>%
  ggplot(aes(npeel, Estimate)) +
  geom_ribbon(data = log_rho_df, fill = "grey80", aes(x, ymin = lower, ymax = upper), inherit.aes = FALSE) +
  #geom_hline(aes(yintercept = lower), linetype = 3) +
  #geom_hline(aes(yintercept = upper), linetype = 3) +
  geom_hline(aes(yintercept = log_rho), linetype = 2) +
  geom_point() +
  geom_linerange(aes(ymin = lwr, max = upr)) +
  #expand_limits(y = 0) +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  labs(x = "Elapsed years", y = expression(log(rho))) +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_retrospective.png", g, height = 5, width = 4)


# Gap (index)
index_gap <- lapply(1:nrow(case_studies), function(i) {
  s <- case_studies$code[i]

  readr::read_csv(file.path(case_studies$dir[i], "index", paste0(s, "_index_gap.csv"))) %>%
    mutate(Name = case_studies$Name[i])
}) %>%
  bind_rows()

g <- index_gap %>%
  left_join(case_studies) %>%
  ggplot(aes(year, log(est))) +
  geom_point(data = dplyr::filter(index_gap, impute), shape = 1, aes(colour = factor(gap))) +
  geom_line(linetype = 2, aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year < year[which(.data$impute)[1]], .by = c(gap, Name)), aes(colour = factor(gap))) +
  geom_line(data = dplyr::filter(index_gap, year > year[max(which(.data$impute))], .by = c(gap, Name)), aes(colour = factor(gap))) +
  #expand_limits(y = 0) +
  facet_wrap(vars(Name), ncol = 2, scales = "free") +
  #geom_vline(aes(xintercept = year_new_vessel - 0.5), linetype = 2) +
  #scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) +
  guides(shape = "none") +
  labs(x = "Year", y = "log(Index)", colour = "Time gap") +
  theme(legend.position = "bottom", strip.background = element_blank())
ggsave("manuscript/figures/index_gap.png", g, height = 7, width = 6)

# Gap (gear effect)
library(snowfall)
sfInit(TRUE, nrow(case_studies))
sfExport(list = c("case_studies"))
sfLibrary(tidyverse)
vessel_effect_gap <- sfLapply(1:nrow(case_studies), function(i) {

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
  labs(x = "Time gap", y = expression(log(rho))) +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_gap.png", g, height = 5, width = 4)


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
library(snowfall)
sfInit(TRUE, nrow(case_studies))
sfExport(list = c("case_studies"))
sfLibrary(tidyverse)
vessel_effect_null <- sfLapply(1:nrow(case_studies), function(i) {

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
  labs(x = "Dummy year change", y = expression(log(rho))) +
  theme(strip.background = element_blank())
ggsave("manuscript/figures/rho_null.png", g, height = 5, width = 4)

sfStop()
