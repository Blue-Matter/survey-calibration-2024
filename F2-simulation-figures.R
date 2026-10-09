
# Used to verify operating models and explore simulation results.

library(tidyverse)

# OM ----
sim_design <- readr::read_csv("simulation/sim_design.csv")
nd_y <- readRDS(file = "simulation/operating_model.rds")

nd_plot <- lapply(1:nrow(sim_design), function(i) {
  nd_y[[i]] %>%
    mutate(OM = paste("OM", i), cell_area = 0.02^2)
}) %>%
  bind_rows()

index_om <- nd_plot %>%
  summarise(Population = sum(mu_pop * cell_area),
            Observed = sum(mu * cell_area),
            .by = c(year, OM)) %>%
  reshape2::melt(id.vars = c("year", "OM"))

# OM figures ----

# Spatial field
g <- nd_plot %>%
  filter(year == 1) %>%
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
ggsave("simulation/figures/omega_s.png", g, height = 8, width = 8)

for (i in 1:nrow(sim_design)) {

  # Spatiotemporal field
  g <- nd_y[[i]] %>%
    mutate(year = factor(paste("Year", year), paste("Year", 1:10))) %>%
    ggplot(aes(X, Y, fill = epsilon_st)) +
    geom_tile() +
    facet_wrap(vars(year)) +
    scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
    coord_cartesian(expand = FALSE) +
    theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
    ggtitle(paste("OM", i))
  ggsave(paste0("simulation/figures/epsilon_", i, ".png"), g, height = 5, width = 6)

  # Abundance
  g <- nd_y[[i]] %>%
    mutate(year = factor(paste("Year", year), paste("Year", 1:10))) %>%
    ggplot(aes(X, Y, fill = eta_pop)) +
    geom_tile() +
    facet_wrap(vars(year)) +
    scale_fill_viridis_c() +
    coord_cartesian(expand = FALSE) +
    labs(fill = "log density") +
    theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
    ggtitle(paste("OM", i))
  ggsave(paste0("simulation/figures/density_", i, ".png"), g, height = 5, width = 6)

}

# Index
g <- index_om %>%
  left_join(sim_design) %>%
  mutate(Label = paste0("kappa=", kappa, ", ARrho=", ar_rho, ", rho=", log_rho, ", phi=", phi)) %>%
  mutate(OM = factor(OM, sim_design$OM)) %>%
  ggplot(aes(year, log(value), linetype = variable, shape = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(OM)) +
  geom_label(aes(label = Label), x = -Inf, y = Inf, hjust = "inward", vjust = "inward") +
  #expand_limits(y = 0) +
  labs(x = "Year", y = "log Abundance", linetype = NULL, shape = NULL) +
  scale_shape_manual(values = c("Population" = 16, "Observed" = 1)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/index.png", g, height = 8, width = 8)

# Results ----

# Plot index
OM_plot <- 1:24
index <- lapply(OM_plot, function(i) {
  res <- readRDS(paste0("simulation/results_ret_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {

    st_check <- res[[x]]$sd_report %>%
      filter(par_name == "sigma_E", yret == 10) %>%
      filter(round(Estimate, 2) > 0)

    if (nrow(st_check)) {
      return(res[[x]]$index %>% mutate(Simulation = x))
    } else {
      return(data.frame())
    }

  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))
}) %>%
  bind_rows() %>%
  filter(yret == 10) %>%
  mutate(geom_mean = exp(mean(log(est))), .by = c(OM, Simulation)) %>%
  mutate(value = est/geom_mean, OM = factor(OM, levels = paste("OM", OM_plot))) %>%
  left_join(sim_design, by = "OM") %>%
  mutate(Label = paste0(OM, "\nkappa=", kappa, ", AR rho=", ar_rho, ", rho=", log_rho, ", phi=", phi))

index_om2 <- index_om %>%
  filter(OM %in% paste("OM", OM_plot), variable == "Population") %>%
  rename(est = value) %>%
  mutate(geom_mean = exp(mean(log(est))), .by = OM) %>%
  mutate(value = est/geom_mean) %>%
  left_join(sim_design) %>%
  mutate(Label = paste0(OM, "\nkappa=", kappa, ", AR rho=", ar_rho, ", rho=", log_rho, ", phi=", phi))

g <- ggplot(index, aes(year, log(est))) +
  geom_line(aes(group = Simulation), alpha = 0.1) +
  geom_line(data = index_om2, colour = "red", linewidth = 1) +
  geom_point(data = index_om2, shape = 21, fill = "white", colour = "red") +
  facet_wrap(vars(Label), ncol = 4, scales = "free_y") +
  expand_limits(y = 0) +
  #coord_cartesian(ylim = c(0, 4)) +
  labs(x = "Year", y = "log(Index)") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/index_compare.png", g, height = 7, width = 10)


# Retrospective ----
OM_plot <- 1:24

# Plot distribution of rho
rho_ret <- lapply(OM_plot, function(i) {
  res <- readRDS(paste0("simulation/results_ret_", i, ".rds"))
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
  mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(Error = Estimate - log_rho)

g <- rho_ret %>%
  ggplot(aes(y, Estimate)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  geom_hline(data = sim_design, aes(yintercept = log_rho), linetype = 2) +
  labs(x = "Years after vessel change", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_ret_OM.png", g, height = 6, width = 5)

# Does CI get smaller?
g <- rho_ret %>%
  ggplot(aes(y, `Std. Error`)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  expand_limits(y = 0) +
  labs(x = "Years after vessel change", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

g <- rho_ret %>%
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
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Error"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_ret_error.png", g, height = 6, width = 5)

rho_ret_stat <- rho_ret %>%
  summarise(bias = mean(Error),
            RMSE = sqrt(mean(Error^2)),
            .by = c(y, OM)) %>%
  left_join(sim_design) %>%
  #mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))

g <- rho_ret_stat %>%
  ggplot(aes(y, RMSE, colour = factor(kappa))) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("RMSE"~~log(rho)),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_ret_RMSE.png", g, height = 6, width = 5)


# Compare trend
index_om2 <- index_om %>%
  filter(OM %in% paste("OM", OM_plot), variable == "Population") %>%
  mutate(geom_mean = exp(mean(log(value))), .by = OM) %>%
  mutate(value = value/geom_mean)

om_trend <- lapply(6:10, function(yret) {
  index_om2 %>%
    filter(year <= yret) %>%
    mutate(yret = .env$yret) %>%
    filter(year > yret - 5) %>% # 5 year moving average
    summarise(slope_om = lm(log(value) ~ year) |> coef() |> getElement(2), .by = c(yret, OM))
}) %>%
  bind_rows()
#om_trend <- index_om %>%
#  filter(year %in% c(5, 6)) %>%
#  summarise(slope_om = log(value[2]/value[1]), .by = OM)

trend_ret <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_ret_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {

    st_check <- res[[x]]$sd_report
    #%>%
    #  filter(par_name == "sigma_E") %>%
    #  filter(round(Estimate, 2) > 0.01)
    #st_check <- data.frame(yret = 6:10)

    if (nrow(st_check)) {
      #ans <- res[[x]]$index %>%
      #  filter(yret %in% st_check$yret) %>%
      #  filter(year > yret - 5) %>%
      #  summarise(slope = lm(log_est ~ year) |> coef() |> getElement(2), .by = yret) %>%
      #  mutate(Simulation = x)

      ans <- res[[x]]$index %>%
        filter(yret %in% st_check$yret) %>%
        mutate(geom_mean = exp(mean(log_est)), .by = yret) %>%
        mutate(value = est/geom_mean) %>%
        filter(year %in% c(5, 6)) %>%
        summarise(slope = log(value[2]/value[1]), .by = yret) %>%
        mutate(Simulation = x)

    } else {
      ans <- data.frame()
    }

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(y = factor(yret - 5), Error = slope - slope_om) %>%
  mutate(est_rho = TRUE)


trend_ret_norho <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_ret_", i, "_norho.rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {

    st_check <- res[[x]]$sd_report #%>%
      #filter(par_name == "sigma_E") %>%
      #filter(round(Estimate, 2) > 0.01)
    #st_check <- data.frame(yret = 6:10)

    if (nrow(st_check)) {
      #ans <- res[[x]]$index %>%
      #  filter(yret %in% st_check$yret) %>%
      #  filter(year > yret - 5) %>%
      #  summarise(slope = lm(log_est ~ year) |> coef() |> getElement(2), .by = yret) %>%
      #  mutate(Simulation = x)

      ans <- res[[x]]$index %>%
        filter(yret %in% st_check$yret) %>%
        mutate(geom_mean = exp(mean(log_est)), .by = yret) %>%
        mutate(value = est/geom_mean) %>%
        filter(year %in% c(5, 6)) %>%
        summarise(slope = log(value[2]/value[1]), .by = yret) %>%
        mutate(Simulation = x)

    } else {
      ans <- data.frame()
    }

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(y = factor(yret - 5), Error = slope - slope_om) %>%
  mutate(est_rho = FALSE)

#g <- trend_ret %>%
#  ggplot(aes(y, slope)) +
#  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
#  facet_wrap(vars(OM), ncol = 4) +
#  geom_point(data = om_trend %>% rename(slope = slope_om) %>% mutate(y = factor(yret - 5))) +
#  labs(x = "Years after vessel change", y = "Index slope") +
#  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
#ggsave("simulation/figures/slope_ret.png", g, height = 6, width = 5)

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
ggsave("simulation/figures/slope_ret_error.png", g, height = 6, width = 5)



trend_ret_stat <- rbind(
  trend_ret, trend_ret_norho
) %>%
  summarise(bias = mean(Error),
            MARE = mean(abs(Error)),
            RMSE = sqrt(mean(Error^2)),
            .by = c(y, OM, est_rho)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))

g <- trend_ret_stat %>%
  filter(est_rho) %>%
  ggplot(aes(y, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.5)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "RMSE in log(Index) slope",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_ret_stat %>% filter(!est_rho), position = position_dodge(0.5), shape = 4, linetype = 3)
ggsave("simulation/figures/slope_ret_RMSE.png", g, height = 6, width = 5)










# Gap ----
rho_gap <- lapply(OM_plot, function(i) {
  res <- readRDS(paste0("simulation/results_gap_", i, ".rds"))
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

g <- rho_gap %>%
  ggplot(aes(factor(ngap), Estimate)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  geom_hline(data = sim_design, aes(yintercept = log_rho), linetype = 2) +
  labs(x = "Time gap between surveys", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_gap.png", g, height = 6, width = 5)

# Does CI get larger?
g <- rho_gap %>%
  ggplot(aes(factor(ngap), `Std. Error`)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  expand_limits(y = 0) +
  labs(x = "Time gap between surveys", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

rho_gap_stat <- rho_gap %>%
  summarise(bias = mean(Error) %>% round(2),
            RMSE = sqrt(mean(Error^2)) %>% round(2),
            .by = c(ngap, OM)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(!cluster, "Random", "Cluster"))

g <- rho_gap %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(!cluster, "Random", "Cluster")) %>%
  ggplot(aes(factor(ngap), Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = expression("Error"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_gap_error.png", g, height = 6, width = 5)


rho_gap_stat <- rho_gap %>%
  summarise(bias = mean(Error),
            RMSE = sqrt(mean(Error^2)),
            .by = c(ngap, OM)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(!cluster, "Random", "Cluster"))

g <- rho_gap_stat %>%
  ggplot(aes(ngap, RMSE, colour = factor(kappa))) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = expression("RMSE"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_gap_RMSE.png", g, height = 6, width = 5)

trend_gap <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_gap_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {

    st_check <- res[[x]]$sd_report %>%
      filter(par_name == "sigma_E") %>%
      filter(round(Estimate, 2) > 0)

    if (nrow(st_check)) {
      ans <- res[[x]]$index %>%
        filter(year > 5) %>%
        summarise(slope = lm(log(est) ~ year) |> coef() |> getElement(2), .by = ngap) %>%
        mutate(Simulation = x)
    } else {
      ans <- data.frame()
    }

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend %>% filter(yret == 10)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(Error = slope - slope_om) %>%
  mutate(est_rho = TRUE)

trend_gap_norho <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_gap_", i, "_norho.rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {

    st_check <- res[[x]]$sd_report %>%
      filter(par_name == "sigma_E") %>%
      filter(round(Estimate, 2) > 0)

    if (nrow(st_check)) {
      ans <- res[[x]]$index %>%
        filter(year > 5) %>%
        summarise(slope = lm(log(est) ~ year) |> coef() |> getElement(2), .by = ngap) %>%
        mutate(Simulation = x)
    } else {
      ans <- data.frame()
    }

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend %>% filter(yret == 10)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(Error = slope - slope_om) %>%
  mutate(est_rho = FALSE)

#g <- trend_ret %>%
#  ggplot(aes(y, slope)) +
#  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
#  facet_wrap(vars(OM), ncol = 4) +
#  geom_point(data = om_trend %>% rename(slope = slope_om) %>% mutate(y = factor(yret - 5))) +
#  labs(x = "Years after vessel change", y = "Index slope") +
#  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
#ggsave("simulation/figures/slope_ret.png", g, height = 6, width = 5)

g <- trend_gap_norho %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(!cluster, "Random", "Cluster")) %>%
  ggplot(aes(factor(ngap), Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-0.3, 0.3)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "Error in log(Index) slope",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/slope_gap_error.png", g, height = 6, width = 5)


trend_gap_stat <- rbind(
  trend_gap, trend_gap_norho
) %>%
  summarise(bias = mean(Error),
            MARE = mean(abs(Error)),
            RMSE = sqrt(mean(Error^2)),
            .by = c(ngap, OM, est_rho, ar_rho)) %>%
  left_join(sim_design) %>%
  #mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM"))

g <- trend_gap_stat %>%
  filter(est_rho) %>%
  ggplot(aes(ngap, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "RMSE in log(Index) slope",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_gap_stat %>% filter(!est_rho),
                  position = position_dodge(0.5),
                  shape = 4, linetype = 3,
                  show.legend = FALSE)
ggsave("simulation/figures/slope_gap_RMSE.png", g, height = 6, width = 5)




# Compare sampling ----
samp <- lapply(1:2, function(i) {
  readRDS(paste0("simulation/dat_list_", i, ".rds"))[[1]] %>%
    filter(year == 1) %>%
    cbind(sim_design[i, ])
}) %>%
  bind_rows() %>%
  mutate(kappa = paste("kappa==", kappa))

pop <- lapply(1:2, function(i) {
  nd_y[[i]] %>%
    filter(year == 10) %>%
    cbind(sim_design[i, ])
}) %>%
  bind_rows() %>%
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
ggsave("manuscript/figures/sampling.png", g, height = 2.5, width = 6)


# Summary error figures
# Rho error & bias
g1 <- rho_ret %>%
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
  coord_cartesian(ylim = c(-1.1, 1.1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Error"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("manuscript/figures/sim_rho_error_retrospective.png", g1, height = 4, width = 5)

g2 <- rho_gap %>%
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
       y = expression("Error"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("manuscript/figures/sim_rho_error_gap.png", g2, height = 4, width = 5)

#g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "bottom")
#ggsave("manuscript/figures/sim_rho_error_sdmsimulate.png", g, height = 5, width = 6)

# Rho RMSE
g1 <- rho_ret %>%
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
       y = expression("RMSE"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

g2 <- rho_gap %>%
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
       y = expression("RMSE"~~widehat(log(rho))),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

#g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "bottom")
#ggsave("manuscript/figures/sim_rho_RMSE_sdmsimulate.png", g, height = 5, width = 6)

# Index trend error
g1 <- trend_ret %>%
  filter(!cluster) %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(y, log_rho, ar_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(ar_rho == 1, "RW~OM", "AR1~OM")) %>%
  ggplot(aes(y, Error, colour = factor(kappa))) +
  geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  aes(y = median, ymin = lower, ymax = upper, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "Error in trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

g2 <- trend_gap %>%
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
  coord_cartesian(ylim = c(-0.4, 0.4)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "Error in trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())

#g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "bottom")
#ggsave("manuscript/figures/sim_trend_error_sdmsimulate.png", g, height = 5, width = 6)

# Index trend RMSE
g1 <- trend_ret_stat %>%
  filter(est_rho, !cluster) %>%
  ggplot(aes(y, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "RMSE in trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_ret_stat %>% filter(!est_rho, !cluster),
                  position = position_dodge(0.5), shape = 4, linetype = 3,
                  show.legend = FALSE)
ggsave("manuscript/figures/sim_trend_RMSE_retrospective.png", g1, height = 4, width = 5)

g2 <- trend_gap_stat %>%
  filter(est_rho, !cluster) %>%
  ggplot(aes(ngap, RMSE, colour = factor(kappa), ymin = 0, ymax = RMSE, group = OM)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  #facet_wrap(vars(row), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.3)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Time gap between surveys",
       y = "RMSE in trend",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank()) +
  geom_pointrange(data = trend_gap_stat %>% filter(!est_rho, !cluster),
                  position = position_dodge(0.5),
                  shape = 4, linetype = 3,
                  show.legend = FALSE)
ggsave("manuscript/figures/sim_trend_RMSE_gap.png", g2, height = 4, width = 5)

#g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "bottom")
#ggsave("manuscript/figures/sim_trend_RMSE_sdmsimulate.png", g, height = 5, width = 6)

# log(rho) coverage
g1 <- rho_ret %>%
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
  coord_cartesian(ylim = c(0.5, 1)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = expression("Coverage"~~log(rho)),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("manuscript/figures/sim_rho_coverage_retrospective.png", g1, height = 4, width = 4)

g2 <- rho_gap %>%
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
       y = expression("Coverage"~~log(rho)),
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("manuscript/figures/sim_rho_coverage_gap.png", g2, height = 4, width = 4)

#g <- ggpubr::ggarrange(g1, g2, ncol = 1, common.legend = TRUE, legend = "bottom")
#ggsave("manuscript/figures/sim_rho_coverage_sdmsimulate.png", g, height = 5, width = 6)
