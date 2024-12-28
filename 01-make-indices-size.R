
library(tidyverse)
#library(NAFOdown)

#ggplot2::theme_set(theme_nafo())
theme_set(theme_bw())

source("99-functions.R")

.ggsave <- function(filename, sp = "cod", ...) {
  filename <- file.path("validation_size", "figures", paste0(sp, "_", filename))
  ggplot2::ggsave(filename, ...)
}

sp <- c("cod", "ghl")

# "est" = use estimated calibration factors from comparative fishing
# "calibrated" = estimate calibration factor in spatiotemporal model
# "none" = no calibration factor in spatiotemporal model
type <- c("est", "calibrated", "norho")

for (s in sp) {

  formals(.ggsave)$sp <- s

  if (s == "cod") {
    # Engel = seq(7, 142, 3)
    # Camp = seq(4, 137, 1)
    #length_bin <- c(0, 30, 45)
    length_bin <- c(0, 42)

    #length_bin = 0
    dat_rho <- readr::read_csv("validation_size/cod_rho_comparative.csv") %>%
      filter(Model == "RW")

    dat_engel <- readr::read_csv("validation_size/LF_engel_438_2J3KL.csv") %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L"))

    dat_camp <- readr::read_csv("validation_size/LF_camp_438_2J3KL.csv") %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L")) %>%
      filter(year < 2020)
  } else {

    # Engel = seq(7, 142, 3)
    # Camp = seq(4, 137, 1)
    #length_bin <- c(0, 25, 45)
    length_bin <- c(0, 34)

    dat_rho <- readr::read_csv("validation_size/ghl_rho_comparative.csv") %>%
      filter(Model == "RW")

    dat_engel <- readr::read_csv("validation_size/LF_engel_892_2J3KL.csv") %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L"))

    dat_camp <- readr::read_csv("validation_size/LF_camp_892_2J3KL.csv") %>%
      filter(NAFOdiv %in% c("2J", "3K", "3L")) %>%
      filter(year < 2020)
  }

  for (i in type) {
    if (i == "est") {
      lf_engel <- expand_number_len(
        dat_lf = dat_engel,
        dat_rho = dat_rho,
        length_bin = length_bin
      )
    } else {
      lf_engel <- expand_number_len(
        dat_lf = dat_engel,
        length_bin = length_bin
      )
    }

    lf_camp <- expand_number_len(
      dat_lf = dat_camp,
      length_bin = length_bin
    )

    dat_lf <- rbind(
      lf_engel,
      lf_camp
    ) %>%
      mutate(
        long.start = -1 * long.start,
        log_depth = log(set.depth.mean),
        swept.area = dist.towed * ifelse(data.series == "Engel", 25.0403, 31.18805)
      ) %>%
      sdmTMB::add_utm_columns(ll_names = c("long.start", "lat.start")) %>%
      mutate(data.series = factor(data.series, levels = c("Engel", "Campelen"))) %>%
      mutate(npue = value.calibrate/swept.area)

    g <- ggplot(dat_lf, aes(year, npue, colour = data.series)) +
      facet_wrap(vars(length)) +
      geom_jitter()

    g <- ggplot(dat_lf, aes(long.start, lat.start, colour = npue, shape = data.series)) +
      facet_grid(vars(year), vars(length)) +
      geom_point() +
      scale_colour_viridis_c(trans = "log")

    ind <- lapply(unique(dat_lf$length), function(i) {
      dat_lf %>%
        filter(length == i) %>%
        mutate(number = npue) %>%
        #mutate(number = number.calibrate) %>%
        make_design_index() %>%
        mutate(length = i)
    }) %>%
      bind_rows()

    g <- ggplot(ind, aes(year, est)) +
      geom_line() +
      geom_point() +
      facet_wrap(vars(length))

    g <- ind %>%
      summarise(est = sum(est), .by = year) %>%
      ggplot(aes(year, est)) +
      geom_line() +
      geom_point()

    dat_lf %>%
      summarise(p = mean(npue == 0), .by = c(year, length))

    # Delta lognormal model
    formula_pos <- switch(
      i,
      "calibrated" = npue ~ 1 + poly(log_depth, 2, raw = TRUE) + data.series,
      npue ~ 1 + poly(log_depth, 2, raw = TRUE)
    )

    formula_zero <- npue ~ 1 + poly(log_depth, 2, raw = TRUE)

    tictoc::tic()
    fit_list <- do_fit_size(
      dat_lf,
      formula = list(formula_zero, formula_pos),
      length_bin = length_bin, sp = s,
      type = i,
      spatiotemporal = list("off", "rw"),
      spatial = "on",
      do_fit = TRUE
    )
    tictoc::toc()
  }

}




#tmb_obj <- fit_list[[1]]$tmb_obj
#.control <- fit_list[[1]]$nlminb_control
#
#tmb_opt <- stats::nlminb(start = tmb_obj$par, objective = tmb_obj$fn,
#                         gradient = tmb_obj$gr, #lower = lim$lower, upper = lim$upper,
#                         control = .control)
#
#fit_list <- lapply(length_bin, function(i) readRDS(paste0("validation_size/models/codfit_calibrated_", i, ".rds")))

#grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
#  mutate(X = 1e-3 * X, Y = 1e-3 * Y,
#         depth = -1 * depth) %>%
#  filter(.data$depth > 0) %>%
#  filter(.data$div %in% unique(dat_lf$NAFOdiv)) %>%
#  mutate(log_depth = log(depth),
#         data.series = "Campelen")
#
#index <- make_index_size(fit_list, length_bin, grid)
#readr::write_csv(index, "validation_size/index/cod_index_calibrated.csv")
#index <- readr::read_csv("validation_size/index/cod_index_calibrated.csv")

#g <- ggplot(index, aes(year, est)) +
#  geom_point() +
#  geom_line() +
#  expand_limits(y = 0) +
#  geom_vline(xintercept = 1994.5, linetype = 2) +
#  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
#  guides(colour = "none") +
#  theme(panel.spacing = unit(0, "in"))
#.ggsave("index_comparative.png", g, height = 3.5, width = 6)
#
#g <- ggplot(index, aes(year, est)) +
#  geom_point() +
#  geom_line() +
#  geom_linerange(aes(ymin = lwr, max = upr)) +
#  expand_limits(y = 0) +
#  geom_vline(xintercept = 1994.5, linetype = 2) +
#  labs(x = "Year", y = "Campelen Index", colour = "Terminal year") +
#  guides(colour = "none") +
#  theme(panel.spacing = unit(0, "in"))
#.ggsave("index_comparative_ci.png", g, height = 3.5, width = 6)
#
### Design-based index ----
#
#
## Compare strat area with values in prediction grid
#grid_area <- grid %>%
#  summarise(area_km2 = sum(area_km2), .by = c(strat, div))
#
#survey_area <- dat_fit %>%
#  summarise(strat.area = mean(strat.area, na.rm = TRUE), .by = c(strat, NAFOdiv)) %>%
#  rename(div = NAFOdiv) %>%
#  left_join(grid_area) %>%
#  filter(!is.na(strat), !is.na(area_km2)) %>%
#  mutate(ratio = area_km2/strat.area)
#
#g <- local({
#
#  # Expanded NPUE
#  lf_engel <- expand_number_len(
#    dat_lf = dat_engel,
#    dat_rho = dat_rho
#  )
#
#  # Nominal NPUE
#  lf_engel2 <- expand_number_len(
#    dat_lf = dat_engel
#  )
#
#  lf_camp <- expand_number_len(
#    dat_lf = dat_camp
#  )
#
#  dat_lf <- rbind(
#    lf_engel,
#    lf_camp
#  ) %>%
#    mutate(swept.area = dist.towed * ifelse(data.series == "Engel", 25.0403, 31.18805)) %>%
#    mutate(npue = number.calibrate/swept.area)
#
#  g1 <- dat_lf %>%
#    rename(number.total = number) %>%
#    rename(number = number.calibrate) %>%
#    make_design_index(simsurvey = FALSE) %>%
#    mutate(Calibrated = TRUE)
#
#  dat_lf2 <- rbind(
#    lf_engel2,
#    lf_camp
#  ) %>%
#    mutate(swept.area = dist.towed * ifelse(data.series == "Engel", 25.0403, 31.18805)) %>%
#    mutate(npue = number.calibrate/swept.area)
#
#  g2 <- dat_lf2 %>%
#    rename(number.total = number) %>%
#    rename(number = number.calibrate) %>%
#    make_design_index(simsurvey = FALSE) %>%
#    mutate(Calibrated = FALSE)
#
#  rbind(g1, g2)
#})
#
#
#gout <- g %>%
#  ggplot(aes(year, est, linetype = Calibrated, shape = Calibrated)) +
#  geom_line() +
#  geom_point() +
#  scale_shape_manual(values = c(1, 16)) +
#  scale_linetype_manual(values = c(2, 1)) +
#  expand_limits(y = 0) +
#  geom_vline(xintercept = 1994.5, linetype = 2) +
#  labs(x = "Year", y = "Design index") +
#  theme(legend.position = "bottom")
#.ggsave("index_design.png", gout, height = 4.5, width = 6)
