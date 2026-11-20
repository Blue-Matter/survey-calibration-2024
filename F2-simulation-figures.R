
# Used to verify operating models and explore simulation results.

library(tidyverse)

# Generate indices of the OM ----
sim_design <- readr::read_csv("simulation_simfield/sim_design.csv")

if (FALSE) {
  nd_plot <- lapply(1:nrow(sim_design), function(i) {
    require(dplyr)
    nd_y <- readRDS(file = paste0("simulation_simfield/operating_model_", i, ".rds"))
    nsim <- length(nd_y)
    lapply(1:nsim, function(x) {
      nd_y[[x]][[1]] %>% mutate(Simulation = x)
    }) %>%
      bind_rows() %>%
      mutate(OM = paste("OM", i))
  }) %>%
    bind_rows() %>%
    mutate(
      cell_area = 0.02 * 0.02,   # Temporary! but no harm to repeat the calculation
      intercept = log(50),
      eta_pop = intercept + omega_s + epsilon_st,
      mu_pop = exp(eta_pop)
    )
  saveRDS(nd_plot, file = "simulation_simfield/simulated_fields.rds")

  index_om <- nd_plot %>%
    summarise(
      Population = sum(mu_pop * cell_area),
      Observed = sum(mu * cell_area),
      .by = c(year, OM, Simulation)
    )
  saveRDS(index_om, file = "simulation_simfield/index_om.rds")

}
# Simulated fields are large (2+ GB) and only needed for the sampling and omega_s figures
has_fields <- file.exists("simulation_simfield/simulated_fields.rds")
if (has_fields) nd_plot <- readRDS(file = "simulation_simfield/simulated_fields.rds")
index_om <- readRDS(file = "simulation_simfield/index_om.rds")


if (has_fields) {
# Manuscript figure of two fields that vary in kappa with sampling locations----
samp <- lapply(1:2, function(i) {
  readRDS(paste0("simulation_simfield/dat_list_", i, ".rds"))[[1]][["nd_y"]] %>%
    filter(year == 10) %>%
    select(!log_rho) %>%
    cbind(sim_design[i, ])
}) %>%
  bind_rows() %>%
  mutate(kappa = paste("kappa==", kappa))

pop <- nd_plot %>%
  filter(year == 10, Simulation == 1, OM %in% paste("OM", 1:2)) %>%
  left_join(sim_design, by = c("log_rho", "OM")) %>%
  mutate(kappa = paste("kappa==", kappa))

g <- pop %>%
  ggplot(aes(X, Y)) +
  geom_tile(aes(fill = eta_pop)) +
  geom_point(data = samp, shape = 21, fill = "grey80") +
  facet_wrap(vars(kappa), labeller = label_parsed) +
  coord_cartesian(expand = FALSE) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
  labs(fill = "log(Density)") +
  theme(strip.background = element_blank())
ggsave("simulation_simfield/figures/sampling.png", g, height = 2.5, width = 6)
ggsave("manuscript/figures/sampling.png", g, height = 2.5, width = 6)


# OM figures ----

# Spatial field
g <- nd_plot %>%
  filter(year == 1, Simulation == 1) %>%
  #filter(OM %in% c("OM 1", "OM 2")) %>%
  left_join(sim_design, by = "OM") %>%
  mutate(OM = factor(OM, sim_design$OM)) %>%
  ggplot(aes(X, Y, fill = omega_s)) +
  geom_tile() +
  facet_wrap(vars(OM)) +
  #facet_wrap(vars(kappa)) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
  coord_cartesian(expand = FALSE) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/omega_s.png", g, height = 8, width = 8)
}

# I don't find these very informative
if (FALSE) {
  for (i in 1:nrow(sim_design)) {
    # Spatiotemporal field
    g <- nd_plot %>%
      filter(Simulation == 1) %>%
      mutate(year = factor(paste("Year", year), paste("Year", 1:10))) %>%
      ggplot(aes(X, Y, fill = epsilon_st)) +
      geom_tile() +
      facet_wrap(vars(year)) +
      scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      coord_cartesian(expand = FALSE) +
      theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
      ggtitle(paste("OM", i))
    ggsave(paste0("simulation_simfield/figures/epsilon_", i, ".png"), g, height = 5, width = 6)

    # Abundance
    g <- nd_plot %>%
      filter(Simulation == 1) %>%
      mutate(year = factor(paste("Year", year), paste("Year", 1:10))) %>%
      ggplot(aes(X, Y, fill = eta_pop)) +
      geom_tile() +
      facet_wrap(vars(year)) +
      scale_fill_viridis_c() +
      coord_cartesian(expand = FALSE) +
      labs(fill = "log density") +
      theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
      ggtitle(paste("OM", i))
    ggsave(paste0("simulation_simfield/figures/density_", i, ".png"), g, height = 5, width = 6)
  }
}

# Plot simulated index
g <- index_om %>%
  left_join(sim_design) %>%
  mutate(Label = paste0("kappa=", kappa, ", ARrho=", ar_rho, ", rho=", log_rho, ", phi=", phi)) %>%
  mutate(OM = factor(OM, sim_design$OM)) %>%
  ggplot(aes(year, log(Population), group = Simulation)) +
  #geom_point() +
  geom_line() +
  facet_wrap(vars(OM)) +
  #geom_label(aes(label = Label), x = -Inf, y = Inf, hjust = "inward", vjust = "inward") +
  #expand_limits(y = 0) +
  labs(x = "Year", y = "log Abundance", linetype = NULL, shape = NULL) +
  #scale_shape_manual(values = c("Population" = 16, "Observed" = 1)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/index.png", g, height = 8, width = 8)

# Results ----
OM_plot <- 1:24

# Plot error in simulated vs. estimated index
# Use geometric mean for trend
get_estimated_index <- function(OM_plot, est_rho = TRUE, type = c("retro", "gap")) {
  type <- match.arg(type)

  lapply(OM_plot, function(i) {
    if (type == "retro") {
      filename <- ifelse(
        est_rho,
        paste0("simulation_simfield/results_ret_", i, ".rds"),
        paste0("simulation_simfield/results_ret_", i, "_norho.rds")
      )
    } else {
      filename <- ifelse(
        est_rho,
        paste0("simulation_simfield/results_gap_", i, ".rds"),
        paste0("simulation_simfield/results_gap_", i, "_norho.rds")
      )
    }

    if (!file.exists(filename)) {
      message("Missing ", filename)
      return(NULL)
    }
    res <- readRDS(filename)
    nsim <- length(res)
    lapply(1:nsim, function(x) {
      st_check <- res[[x]]$sd_report %>%
        filter(par_name == "sigma_E") %>%
        filter(round(Estimate, 2) > 0)

      if (nrow(st_check)) { # Check for convergence
        return(res[[x]]$index %>% mutate(Simulation = x))
      } else {
        return(data.frame())
      }
    }) %>%
      bind_rows() %>%
      mutate(OM = paste("OM", i))
  }) %>%
    bind_rows() %>%
    mutate(geom_mean = exp(mean(log(est))),
           .by = c(OM, Simulation, any_of(c("yret", "ngap")))) %>%
    mutate(est_rel = est/geom_mean)
}

index_est <- get_estimated_index(OM_plot)

index_om <- index_om %>%
  filter(OM %in% paste("OM", OM_plot))

# Plot annual error
# Normalise the OM over the same years as each retrospective peel
g <- left_join(index_est, index_om, by = c("year", "OM", "Simulation")) %>%
  mutate(Pop_rel = Population / exp(mean(log(Population))), .by = c(OM, Simulation, yret)) %>%
  filter(yret == max(yret)) %>% # full time series only
  left_join(sim_design, by = "OM") %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(Label = paste0(OM, "\nkappa=", kappa, ", AR rho=", ar_rho, ", rho=", log_rho, ", phi=", phi)) %>%
  ggplot(aes(year, log(est_rel/Pop_rel))) +
  geom_line(aes(group = Simulation), alpha = 0.1) +
  facet_wrap(vars(OM)) +
  geom_hline(yintercept = 0, linetype = 3) +
  labs(x = "Year", y = "Index error") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/index_compare.png", g, height = 7, width = 10)


# Retrospective ----

# Plot distribution of rho
get_estimated_rho <- function(OM_plot, type = c("retro", "gap")) {
  type <- match.arg(type)

  lapply(OM_plot, function(i) {
    filename <- ifelse(
      type == "retro",
      paste0("simulation_simfield/results_ret_", i, ".rds"),
      paste0("simulation_simfield/results_gap_", i, ".rds")
    )
    if (!file.exists(filename)) {
      message("Missing ", filename)
      return(NULL)
    }
    res <- readRDS(filename)
    nsim <- length(res)
    lapply(1:nsim, function(x) {
      res[[x]]$sd_report %>% filter(par_name == "b_j.1") %>%
        mutate(Simulation = x)
    }) %>%
      bind_rows() %>%
      mutate(OM = paste("OM", i))
  }) %>%
    bind_rows() %>%
    left_join(sim_design) %>%
    mutate(OM = factor(OM, levels = sim_design$OM)) %>%
    mutate(Error = Estimate - log_rho)
}

rho_ret <- get_estimated_rho(OM_plot) %>%
  mutate(y = factor(yret - 5))

g <- rho_ret %>%
  ggplot(aes(y, Estimate)) +
  geom_violin(linewidth = 0.5, draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  geom_hline(data = filter(sim_design, OM %in% paste("OM", OM_plot)), aes(yintercept = log_rho), linetype = 2) +
  labs(x = "Years after vessel change", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/rho_ret_OM.png", g, height = 6, width = 5)

# Rho coverage ----
g <- rho_ret %>%
  filter(!is.na(`Std. Error`)) %>%
  mutate(lower = Estimate - 1.96 * `Std. Error`,
         upper = Estimate + 1.96 * `Std. Error`,
         good_coverage = log_rho >= lower & log_rho <= upper) %>%
  summarise(p_coverage = mean(good_coverage), .by = c(y, log_rho, ar_rho, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(y, p_coverage, group = OM, colour = factor(kappa))) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  geom_line(position = position_dodge(0.0)) +
  geom_point(position = position_dodge(0.0), aes(group = OM, shape = factor(phi))) +
  #geom_pointrange(position = position_dodge(0.6),
  #                ymin = 0, aes(ymax = p_coverage, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0.8, 1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Coverage"~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/sim_rho_coverage_retrospective.png", g, height = 4, width = 4)
ggsave("manuscript/figures/sim_rho_coverage_retrospective.png", g, height = 4, width = 4)

# Make table for paper
p_cov_table <- g@data %>%
  filter(y == 5) %>%
  reshape2::dcast(log_rho + kappa ~ ar_rho + phi, value.var = "p_coverage")
readr::write_csv(p_cov_table, "tables/coverage_rho.csv")

# Rho error during retrospective ----
g <- rho_ret %>%
  filter(!cluster) %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(y, log_rho, ar_rho, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(y, Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Error"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/sim_rho_error_retrospective.png", g, height = 4, width = 5)
ggsave("manuscript/figures/sim_rho_error_retrospective.png", g, height = 4, width = 5)

# Rho RMSE ----
rho_ret_stat <- rho_ret %>%
  summarise(bias = mean(Error),
            RMSE = sqrt(mean(Error^2)),
            .by = c(y, OM)) %>%
  left_join(sim_design) %>%
  #mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))

g <- rho_ret %>%
  filter(!cluster) %>%
  summarise(RMSE = sqrt(mean(Error^2)),
            .by = c(y, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(y, RMSE, colour = factor(kappa))) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = RMSE, ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  #coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("RMSE"~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/rho_ret_RMSE.png", g, height = 6, width = 5)


# Trend RMSE ----
index_ret <- index_est
index_ret_norho <- get_estimated_index(OM_plot, est_rho = FALSE)

compare_index_trend_retrospective <- function(index_est, index_om) {
  om_trend <- lapply(unique(index_est$yret), function(yret) {
    index_om %>%
      filter(year <= yret) %>%
      mutate(yret = .env$yret) %>%
      filter(year > yret - 5) %>% # 5 year moving average
      summarise(slope_om = lm(log(Population) ~ year) |> coef() |> getElement(2), .by = c(yret, OM, Simulation))
  }) %>%
    bind_rows()
  index_est %>%
    filter(year > yret - 5) %>%
    summarise(slope = lm(log_est ~ year) |> coef() |> getElement(2), .by = c(yret, Simulation, OM)) %>%
    left_join(om_trend) %>%
    left_join(sim_design) %>%
    mutate(OM = factor(OM, levels = sim_design$OM)) %>%
    mutate(y = factor(yret - 5), Error = slope - slope_om)
}

trend_ret <- compare_index_trend_retrospective(index_ret, index_om) %>%
  mutate(est_rho = TRUE)
trend_ret_norho <- compare_index_trend_retrospective(index_ret_norho, index_om) %>%
  mutate(est_rho = FALSE)

# Plot the error in the log(Index) 5 year trend
g <- trend_ret %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(y, log_rho, ar_rho, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(y, Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "Error in log(Index) slope",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/slope_ret_error.png", g, height = 6, width = 5)

trend_ret_stat <- rbind(
  trend_ret, trend_ret_norho
) %>%
  mutate(Error_rel = Error/slope_om) %>%
  summarise(bias = mean(Error),
            MAE = mean(abs(Error)),
            MARE = mean(abs(Error_rel)),
            RMSE = sqrt(mean(Error^2)),
            RRMSE = sqrt(mean(Error_rel^2)),
            .by = c(y, OM, est_rho)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))

# Index trend RMSE (root mean square error)
g <- trend_ret_stat %>%
  filter(est_rho, !cluster) %>%
  ggplot(aes(y, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.8)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "RMSE in 5-year log(Index) trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_ret_stat %>% filter(!est_rho, !cluster),
                  position = position_dodge(0.5), shape = 4, linetype = 3,
                  show.legend = FALSE)
ggsave("simulation_simfield/figures/sim_trend_RMSE_retrospective.png", g, height = 4, width = 5)
ggsave("manuscript/figures/sim_trend_RMSE_retrospective.png", g, height = 4, width = 5)










# Gap ----
rho_gap <- get_estimated_rho(OM_plot, type = "gap")

# Rho error ----
g <- rho_gap %>%
  filter(!cluster) %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(ngap, Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-1.5, 1.5)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = expression("Error"~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/sim_rho_error_gap.png", g, height = 4, width = 5)
ggsave("manuscript/figures/sim_rho_error_gap.png", g, height = 4, width = 5)

# Rho coverage ----
g <- rho_gap %>%
  filter(!cluster, !is.na(`Std. Error`)) %>%
  mutate(lower = Estimate - 1.96 * `Std. Error`,
         upper = Estimate + 1.96 * `Std. Error`,
         good_coverage = log_rho >= lower & log_rho <= upper) %>%
  summarise(p_coverage = mean(good_coverage), .by = c(ngap, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(ngap, p_coverage, group = OM, colour = factor(kappa))) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  geom_line(position = position_dodge(0.0)) +
  geom_point(position = position_dodge(0.0), aes(group = OM, shape = factor(phi))) +
  #geom_pointrange(position = position_dodge(0.6),
  #                linewidth = 0,
  #                ymin = 0, aes(ymax = p_coverage, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = expression("Coverage"~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation_simfield/figures/sim_rho_coverage_gap.png", g, height = 4, width = 4)
ggsave("manuscript/figures/sim_rho_coverage_gap.png", g, height = 4, width = 4)


# Rho RMSE ----
g <- rho_gap %>%
  filter(!cluster) %>%
  summarise(RMSE = sqrt(mean(Error^2)),
            .by = c(ngap, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(ngap, Error, colour = factor(kappa))) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = RMSE, ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  #coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = expression("RMSE"~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

# Trend RMSE ----
index_gap <- get_estimated_index(OM_plot, type = "gap")
index_gap_norho <- get_estimated_index(OM_plot, est_rho = FALSE, type = "gap")

compare_index_trend_gap <- function(index_est, index_om) {
  om_trend <- index_om %>%
    filter(year > 5) %>%
    summarise(slope_om = lm(log(Population) ~ year) |> coef() |> getElement(2), .by = c(OM, Simulation))

  index_est %>%
    filter(year > 5) %>%
    summarise(slope = lm(log_est ~ year) |> coef() |> getElement(2), .by = c(ngap, Simulation, OM)) %>%
    left_join(om_trend) %>%
    left_join(sim_design) %>%
    mutate(OM = factor(OM, levels = sim_design$OM),
           Error = slope - slope_om)
}

trend_gap <- compare_index_trend_gap(index_gap, index_om) %>%
  mutate(est_rho = TRUE)
trend_gap_norho <- compare_index_trend_gap(index_gap_norho, index_om) %>%
  mutate(est_rho = FALSE)

g <- trend_gap %>%
  filter(!cluster) %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(ngap, Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "Error in 5-year log(Index) trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())


trend_gap_stat <- rbind(
  trend_gap, trend_gap_norho
) %>%
  summarise(bias = mean(Error),
            MAE = mean(abs(Error)),
            RMSE = sqrt(mean(Error^2)),
            .by = c(ngap, OM, est_rho, ar_rho)) %>%
  left_join(sim_design) %>%
  #mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))


g <- trend_gap_stat %>%
  filter(est_rho, !cluster) %>%
  ggplot(aes(ngap, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 1.5)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "RMSE in 5-year log(Index) trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_gap_stat %>% filter(!est_rho, !cluster),
                  position = position_dodge(0.5),
                  shape = 4, linetype = 3,
                  show.legend = FALSE)
ggsave("simulation_simfield/figures/sim_trend_RMSE_gap.png", g, height = 4, width = 5)
ggsave("manuscript/figures/sim_trend_RMSE_gap.png", g, height = 4, width = 5)




# Index coverage using SE's
# Absolute scale: predictions use new_vessel = FALSE, so est is on the OM scale
# Only years after the vessel change (year > 5)
calculate_index_coverage <- function(index_est, index_om) {
  index_est %>%
    filter(year > 5) %>%
    left_join(index_om, by = c("OM", "year", "Simulation")) %>%
    mutate(good_coverage = lwr <= Population & upr >= Population) %>%
    summarise(p_coverage = mean(good_coverage), .by = c(OM, yret)) %>%
    left_join(sim_design, by = "OM") %>%
    mutate(y = factor(yret - 5)) %>%
    mutate(row = paste0("log(rho)==", log_rho)) %>%
    mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))
}

coverage_ret <- calculate_index_coverage(index_ret, index_om)
coverage_ret_norho <- calculate_index_coverage(index_ret_norho, index_om)

g <- coverage_ret %>%
  ggplot(aes(y, p_coverage, group = OM, colour = factor(kappa))) +
  geom_hline(yintercept = 0.95, linetype = 3) +
  geom_line(position = position_dodge(0.0)) +
  geom_point(position = position_dodge(0.0), aes(group = OM, shape = factor(phi))) +

  geom_line(data = coverage_ret_norho, linetype = 2, position = position_dodge(0.0)) +
  geom_point(data = coverage_ret_norho, position = position_dodge(0.0), aes(group = OM, shape = factor(phi))) +

  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Coverage"~log(Index)),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("manuscript/figures/sim_index_coverage.png", g, height = 4, width = 5)
