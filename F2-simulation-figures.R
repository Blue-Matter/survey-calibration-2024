
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
  res <- readRDS(paste0("simulation/results_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {
    res[[x]]$index %>% mutate(Simulation = x)
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))
}) %>%
  bind_rows() %>%
  mutate(geom_mean = exp(mean(log(est))), .by = OM) %>%
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

# Plot log rho
rho_est <- lapply(OM_plot, function(i) {
  res <- readRDS(paste0("simulation/results_", i, ".rds"))
  nsim <- length(res)
  lapply(1:nsim, function(x) {
    data.frame(
      est = res[[x]]$sd_report[2, 1],
      Simulation = x
    )
  }) %>%
    bind_rows() %>%
    mutate(OM = paste("OM", i))
}) %>%
  bind_rows() %>%
  left_join(sim_design) %>%
  mutate(OM = factor(OM, levels = sim_design$OM))

rho_stat <- rho_est %>%
  summarise(
    bias = mean(est - log_rho) %>% round(2),
    RMSE = sqrt(mean((est - log_rho)^2)) %>% round(2),
    .by = OM
  ) %>%
  left_join(sim_design)

g <- ggplot(rho_est) +
  geom_histogram(aes(x = est), linewidth = 0.1, colour = "black", fill = "grey80") +
  facet_wrap(vars(OM), ncol = 4) +
  geom_vline(data = rho_stat, aes(xintercept = log_rho), linetype = 2) +
  geom_label(data = rho_stat, aes(label = bias), x = Inf, y = Inf, hjust = "inward", vjust = "inward") +
  geom_label(data = rho_stat, aes(label = RMSE), x = -Inf, y = Inf, hjust = "inward", vjust = "inward") +
  labs(x = expression(log(rho))) +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/log_rho.png", g, height = 6, width = 5)


# Compare trend
index_trend <- index %>%
  summarise(slope = lm(log(value) ~ year) |> coef() |> getElement(2),
            .by = c(OM, Simulation))

om_trend <- index_om2 %>%
  summarise(slope = lm(log(value) ~ year) |> coef() |> getElement(2),
            .by = OM)

g <- ggplot(index_trend) +
  geom_histogram(aes(x = slope), linewidth = 0.1, colour = "black", fill = "grey80") +
  facet_wrap(vars(OM), ncol = 4) +
  geom_vline(data = om_trend, aes(xintercept = slope), linetype = 2) +
  #geom_label(data = rho_stat, aes(label = bias), x = Inf, y = Inf, hjust = "inward", vjust = "inward") +
  #geom_label(data = rho_stat, aes(label = RMSE), x = -Inf, y = Inf, hjust = "inward", vjust = "inward") +
  labs(x = "Index slope") +
  theme(legend.position = "bottom", panel.spacing = unit(0, "in"), strip.background = element_blank())
ggsave("simulation/figures/index_slope.png", g, height = 6, width = 5)
