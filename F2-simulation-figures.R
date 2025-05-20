
library(tidyverse)

# OM ----
sim_design <- readr::read_csv("simulation/sim_design.csv")
nd_y <- readRDS(file = "simulation/operating_model.rds")

nd_plot <- lapply(1:nrow(sim_design), function(i) {
  nd_y[[i]] %>%
    mutate(OM = paste("OM", i))
}) %>%
  bind_rows() %>%
  mutate(OM = factor(OM, sim_design$OM))

index_om <- nd_plot %>%
  summarise(Population = sum(pred_pop),
            Observed = sum(pred),
            .by = c(year, OM)) %>%
  reshape2::melt(id.var = c("year", "OM")) %>%
  mutate(OM = factor(OM, levels = paste("OM", 1:nrow(sim_design))))

# OM figures ----

# Spatial field
g <- nd_plot %>%
  filter(year == 1) %>%
  ggplot(aes(X, Y, fill = omega_s)) +
  geom_tile() +
  facet_wrap(vars(OM)) +
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
    mutate(year = paste("Year", year)) %>%
    ggplot(aes(X, Y, fill = log_pred_pop)) +
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
  ggplot(aes(year, log(value), linetype = variable, shape = variable)) +
  geom_point() +
  geom_line() +
  facet_wrap(vars(OM)) +
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
    res[[x]]$index %>% mutate(Simulation = x)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))
}) %>%
  bind_rows() %>%
  filter(yret == 10) %>%
  mutate(geom_mean = exp(mean(log(est))), .by = c(OM, Simulation)) %>%
  mutate(value = est/geom_mean, OM = factor(OM, levels = paste("OM", OM_plot)))

index_om2 <- index_om %>%
  filter(OM %in% paste("OM", OM_plot), variable == "Population") %>%
  mutate(geom_mean = exp(mean(log(value))), .by = OM) %>%
  mutate(value = value/geom_mean)

g <- ggplot(index, aes(year, value)) +
  geom_line(aes(group = Simulation), alpha = 0.1) +
  geom_line(data = index_om2, colour = "red", linewidth = 1) +
  geom_point(data = index_om2, shape = 21, fill = "white", colour = "red") +
  facet_wrap(vars(OM), ncol = 4) +
  expand_limits(y = 0) +
  coord_cartesian(ylim = c(0, 4)) +
  labs(x = "Year", y = "Index") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/index_compare.png", g, height = 7, width = 6)


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
  mutate(Error = Estimate - log_rho) %>%
  mutate(coverage = log_rho < Estimate + 1.96 * `Std. Error` & log_rho > Estimate - 1.96 * `Std. Error`) %>%
  mutate(coverage_zero = Estimate + 1.96 * `Std. Error` > 0 & Estimate - 1.96 * `Std. Error` < 0)

g <- rho_ret %>%
  ggplot(aes(y, Estimate)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  geom_hline(data = sim_design, aes(yintercept = log_rho), linetype = 2) +
  labs(x = "Years after vessel change", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_ret_OM.png", g, height = 6, width = 5)

g <- rho_ret %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(y, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
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
  mutate(column = ifelse(cluster, "Random", "Cluster"))

g <- rho_ret_stat %>%
  ggplot(aes(y, RMSE, colour = factor(kappa))) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "RMSE",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_ret_RMSE.png", g, height = 6, width = 5)


#g <- rho_ret %>%
#  summarise(p_cover = mean(coverage), .by = c(y, OM, kappa, log_rho, cluster, phi)) %>%
#  mutate(row = paste0("log(rho)==", log_rho)) %>%
#  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
#  ggplot(aes(y, p_cover, colour = factor(kappa), shape = factor(phi))) +
#  #geom_point(position = position_dodge(0.6)) +
#  #geom_hline(yintercept = 0, linetype = 3) +
#  geom_pointrange(fatten = 3.5, #shape = 1,
#                  position = position_dodge(0.6),
#                  ymin = 0,
#                  aes(ymax = p_cover, group = OM, shape = factor(phi))) +
#  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
#  #coord_cartesian(ylim = c(-1.5, 1.5)) +
#  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
#  labs(x = "Years after vessel change",
#       y = expression("Coverage"~~widehat(log(rho))),
#       colour = expression(kappa),
#       shape = expression(phi)) +
#  coord_cartesian(ylim = c(0, 1)) +
#  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
#ggsave("simulation/figures/log_rho_ret_coverage.png", g, height = 6, width = 5)
#
#g <- rho_ret %>%
#  summarise(p_cover = mean(coverage_zero), .by = c(y, OM, kappa, log_rho, cluster, phi)) %>%
#  mutate(row = paste0("log(rho)==", log_rho)) %>%
#  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
#  ggplot(aes(y, p_cover, colour = factor(kappa), shape = factor(phi))) +
#  #geom_point(position = position_dodge(0.6)) +
#  #geom_hline(yintercept = 0, linetype = 3) +
#  geom_pointrange(fatten = 3.5, #shape = 1,
#                  position = position_dodge(0.6),
#                  ymin = 0,
#                  aes(ymax = p_cover, group = OM, shape = factor(phi))) +
#  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
#  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
#  labs(x = "Years after vessel change",
#       y = expression("Coverage-zero"~~widehat(log(rho))),
#       colour = expression(kappa),
#       shape = expression(phi)) +
#  coord_cartesian(ylim = c(0, 1)) +
#  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
#ggsave("simulation/figures/log_rho_ret_coverage_zero.png", g, height = 6, width = 5)


# Compare trend
index_om2 <- index_om %>%
  filter(OM %in% paste("OM", OM_plot), variable == "Population") %>%
  mutate(geom_mean = exp(mean(log(value))), .by = OM) %>%
  mutate(value = value/geom_mean)

om_trend <- lapply(6:10, function(yret) {
  index_om2 %>%
    filter(year <= yret) %>%
    summarise(slope_om = lm(log(value) ~ year) |> coef() |> getElement(2), .by = OM) %>%
    mutate(yret = yret)
}) %>%
  bind_rows()



trend_ret <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_ret_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {
    ans <- res[[x]]$index %>%
      summarise(slope = lm(log(est) ~ year) |> coef() |> getElement(2), .by = yret) %>%
      mutate(Simulation = x)

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(y = factor(yret - 5), Error = slope - slope_om)

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
            .by = c(y, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
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



trend_ret_stat <- trend_ret %>%
  summarise(bias = mean(Error),
            RMSE = sqrt(mean(Error^2)),
            .by = c(y, OM)) %>%
  left_join(sim_design) %>%
  #mutate(y = factor(yret - 5)) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster"))

g <- trend_ret_stat %>%
  ggplot(aes(y, RMSE, colour = factor(kappa))) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(position = position_dodge(0.5), aes(ymin = 0, ymax = RMSE, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  coord_cartesian(ylim = c(0, 0.2)) +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "RMSE in log(Index) slope",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/slope_ret_RMSE.png", g, height = 6, width = 5)





coverage_ret <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_ret_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {
    res[[x]]$index %>% mutate(Simulation = x)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(index_om2) %>%
  mutate(geom_mean_est = exp(mean(log(est))), .by = c(yret, Simulation, OM)) %>%
  mutate(est = est/geom_mean_est, lwr = lwr/geom_mean_est, upr = upr/geom_mean_est) %>%
  mutate(coverage = value < upr & value > lwr)

g <- coverage_ret %>%
  summarise(coverage_all_years = all(coverage), .by = c(yret, Simulation, OM)) %>%
  summarise(p_cover = mean(coverage_all_years), .by = c(yret, OM)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(y = factor(yret - 5)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
  ggplot(aes(y, p_cover, colour = factor(kappa), shape = factor(phi))) +
  #geom_point(position = position_dodge(0.6)) +
  #geom_hline(yintercept = 0, linetype = 3) +
  geom_pointrange(fatten = 3.5, #shape = 1,
                  position = position_dodge(0.6),
                  ymin = 0,
                  aes(ymax = p_cover, group = OM, shape = factor(phi))) +
  facet_grid(vars(row), vars(column), labeller = "label_parsed") +
  scale_shape_manual(values = c("0.2" = 16, "4" = 1)) +
  labs(x = "Years after vessel change",
       y = "Index coverage",
       colour = expression(kappa),
       shape = expression(phi)) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/coverage_ret.png", g, height = 6, width = 5)







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
  mutate(Error = Estimate - log_rho) %>%
  mutate(coverage = log_rho < Estimate + 1.96 * `Std. Error` & log_rho > Estimate - 1.96 * `Std. Error`) %>%
  mutate(coverage_zero = Estimate + 1.96 * `Std. Error` > 0 & Estimate - 1.96 * `Std. Error` < 0)

g <- rho_gap %>%
  ggplot(aes(factor(ngap), Estimate)) +
  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
  facet_wrap(vars(OM), ncol = 4) +
  geom_hline(data = sim_design, aes(yintercept = log_rho), linetype = 2) +
  labs(x = "Time gap between surveys", y = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/rho_gap.png", g, height = 6, width = 5)

#rho_gap_stat <- rho_gap %>%
#  summarise(bias = mean(Error) %>% round(2),
#            RMSE = sqrt(mean(Error^2)) %>% round(2),
#            .by = c(ngap, OM)) %>%
#  left_join(sim_design) %>%
#  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
#  mutate(row = paste0("log(rho)==", log_rho)) %>%
#  mutate(column = ifelse(cluster, "Random", "Cluster"))

g <- rho_gap %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
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


trend_gap <- lapply(OM_plot, function(i) {

  res <- readRDS(paste0("simulation/results_gap_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {
    ans <- res[[x]]$index %>%
      summarise(slope = lm(log(est) ~ year) |> coef() |> getElement(2), .by = ngap) %>%
      mutate(Simulation = x)

    return(ans)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))

}) %>%
  bind_rows() %>%
  left_join(om_trend %>% filter(yret == 10)) %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM)) %>%
  mutate(Error = slope - slope_om)

#g <- trend_ret %>%
#  ggplot(aes(y, slope)) +
#  geom_violin(draw_quantiles = c(0.025, 0.5, 0.975)) +
#  facet_wrap(vars(OM), ncol = 4) +
#  geom_point(data = om_trend %>% rename(slope = slope_om) %>% mutate(y = factor(yret - 5))) +
#  labs(x = "Years after vessel change", y = "Index slope") +
#  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
#ggsave("simulation/figures/slope_ret.png", g, height = 6, width = 5)

g <- trend_gap %>%
  summarise(median = median(Error),
            lower = quantile(Error, 0.025),
            upper = quantile(Error, 0.975),
            .by = c(ngap, log_rho, cluster, OM, kappa, phi)) %>%
  mutate(row = paste0("log(rho)==", log_rho)) %>%
  mutate(column = ifelse(cluster, "Random", "Cluster")) %>%
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




# Compare sampling ----
samp <- lapply(c(1:2, 7:8), function(i) {
  readRDS(paste0("simulation/dat_list_", i, ".rds"))[[1]] %>%
    filter(year == 1) %>%
    cbind(sim_design[i, ])
}) %>%
  bind_rows() %>%
  mutate(design = ifelse(cluster, "Cluster", "Random"), kappa = paste("kappa = ", kappa))

pop <- lapply(c(1:2, 7:8), function(i) {
  nd_y[[i]] %>%
    filter(year == 1) %>%
    cbind(sim_design[i, ])
}) %>%
  bind_rows() %>%
  mutate(design = ifelse(cluster, "Cluster", "Random"), kappa = paste("kappa = ", kappa))


g <- ggplot(pop, aes(X, Y)) +
  geom_tile(aes(fill = log_pred_pop)) +
  geom_point(data = samp, shape = 21, fill = "grey80") +
  facet_grid(vars(kappa), vars(design)) +
  coord_cartesian(expand = FALSE) +
  scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
  labs(fill = "log Density") +
  theme(strip.background = element_blank())
ggsave("simulation/figures/sampling.png", g, height = 4, width = 6)
