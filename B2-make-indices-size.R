
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
#sp <- "cod"

# "est" = use estimated calibration factors from comparative fishing
# "calibrated" = estimate calibration factor in spatiotemporal model
# "norho" = no calibration factor in spatiotemporal model
type <- c("est", "calibrated", "norho")

fit <- TRUE

# Fit models
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
      mutate(number.int = round(value.expand), offset = log(swept.area/rho.eff))

    #g <- ggplot(dat_lf, aes(year, npue, colour = data.series)) +
    #  facet_wrap(vars(length)) +
    #  geom_jitter()
    #g <- ggplot(dat_lf, aes(long.start, lat.start, colour = npue, shape = data.series)) +
    #  facet_grid(vars(year), vars(length)) +
    #  geom_point() +
    #  scale_colour_viridis_c(trans = "log")

    #ind <- lapply(unique(dat_lf$length), function(i) {
    #  dat_lf %>%
    #    filter(length == i) %>%
    #    mutate(number = npue) %>%
    #    #mutate(number = number.calibrate) %>%
    #    make_design_index() %>%
    #    mutate(length = i)
    #}) %>%
    #  bind_rows()
    #g <- ggplot(ind, aes(year, est)) +
    #  geom_line() +
    #  geom_point() +
    #  facet_wrap(vars(length))
    #g <- ind %>%
    #  summarise(est = sum(est), .by = year) %>%
    #  ggplot(aes(year, est)) +
    #  geom_line() +
    #  geom_point()
    #dat_lf %>%
    #  summarise(p = mean(npue == 0), .by = c(year, length))

    # Fit model
    if (fit) {
      formula_pos <- switch(
        i,
        "calibrated" = number.int ~ 1 + poly(log_depth, 2, raw = TRUE) + data.series,
        number.int ~ 1 + poly(log_depth, 2, raw = TRUE)
      )

      formula_zero <- number.int ~ 1 + poly(log_depth, 2, raw = TRUE)

      tictoc::tic()
      fit_list <- do_fit_size(
        dat_lf,
        formula = formula_pos,
        #formula = list(formula_zero, formula_pos),
        length_bin = length_bin,
        sp = s,
        type = i,
        #spatiotemporal = list("off", "rw"),
        spatiotemporal = "rw",
        spatial = "on",
        family = nbinom2(),
        do_fit = TRUE
      )
      tictoc::toc()

      fit_list <- lapply(length_bin, function(ell) readRDS(file.path("validation_size", "models", paste0(s, "_fit_", i, "_", ell, ".rds"))))

      grid <- readr::read_csv("shapefiles/prediction_grid.csv") %>%
        mutate(X = 1e-3 * X, Y = 1e-3 * Y,
               depth = -1 * depth) %>%
        filter(.data$depth > 0) %>%
        filter(.data$div %in% unique(dat_lf$NAFOdiv)) %>%
        mutate(log_depth = log(depth),
               data.series = "Campelen")

      index <- make_index_size(fit_list, length_bin, grid)
      readr::write_csv(index, file.path("validation_size", "index", paste0(s, "_index_", i, ".csv")))
    }
  }

  # Plot calibration factors
  # Currently there is an issue, need to properly use offset in the sdm model
  rho <- readr::read_csv(file.path("validation_size", paste0(s, "_rho_comparative.csv"))) %>%
    #filter(Model != "GAM") %>%
    filter(Model == "RW") %>%
    mutate(Model = paste("CF -", Model))

  fit_list <- lapply(length_bin, function(ell) {
    readRDS(file.path("validation_size", "models", paste0(s, "_fit_calibrated_", ell, ".rds")))
  })
  log_rho <- sapply(fit_list, function(i) coef(i)["data.seriesCampelen"])
  se <- sapply(fit_list, function(i) as.numeric(sdmTMB::tidy(i)[4, "std.error"]))

  rho_sdm <- data.frame(
    bin = unique(rho$bin),
    offset_ratio = 1,
    id = NA
  ) %>%
    mutate(log_rho = ifelse(unique(rho$bin) < length_bin[2], log_rho[1], log_rho[2]),
           se = ifelse(unique(rho$bin) < length_bin[2], se[1], se[2]),
           rho = exp(log_rho),
           Model = "Spatiotemporal"
    )

  g <- rbind(rho, rho_sdm) %>%
    ggplot(aes(bin, colour = Model)) +
    geom_line(aes(y = log_rho)) +
    geom_line(aes(y = log_rho + 1.96 * se), linetype = 2) +
    geom_line(aes(y = log_rho - 1.96 * se), linetype = 2) +
    labs(x = "Length", y = expression(log(rho))) +
    coord_cartesian(xlim = if (s == "cod") c(0, 100) else c(0, 75),
                    ylim = if (s == "cod") c(-5, 7) else c(-4, 6))
  .ggsave("rho_sdm.png", sp = s, g, width = 6, height = 3.5)


}

# Plot indices
for (s in sp) {
  index <- lapply(type, function(i) {
    readr::read_csv(file.path("validation_size", "index", paste0(s, "_index_", i, ".csv"))) %>%
      mutate(type = switch(i,
                           "est" = "Comparative rho",
                           "calibrated" = "Estimate rho",
                           "norho" = "No rho"),
             length = ifelse(length == 0,
                             ifelse(s == "cod", "Small (0-40 cm)", "Small (0-30 cm)"),
                             ifelse(s == "cod", "Large (42+ cm)", "Large (32+ cm)")
             )
      )
  }) %>%
    bind_rows()

  g <- ggplot(index, aes(year, est, linetype = type, shape = type)) +
    facet_wrap(vars(length), scales = "free_y") +
    geom_line() +
    #geom_point() +
    expand_limits(y = 0) +
    labs(x = "Year", y = "Campelen Index", linetype = "Model", shape = "Model") +
    geom_vline(xintercept = 1994, linetype = 2) +
    scale_linetype_manual(values = c(1, 2, 3)) +
    scale_shape_manual(values = c(16, 1, 4)) +
    theme(legend.position = "bottom")

  .ggsave("index_size.png", sp = s, g, height = 3, width = 6)
}
