rm(list = ls())
library(tidyverse)
source("simulationParameters.R")
source("figure_save_parameters.R")
source("useful_plot_functions.R")
# One obs, several starts -------------------------------------------------
title_size <- 25
axis_size <- 20
foo_plot <- function(name, sel = T,  legendPosition = "none"){
  thetaEst <- read.table(paste0("simulation_results/est_", name, ".txt"),
                         h = F, sep = ";")
  thetaEst <- thetaEst[, order(thetaEst[1, ])]
  theta_gg <- do.call(rbind.data.frame,
                      lapply(1:ncol(thetaEst), function(i){
                        data.frame(Estimate = thetaEst[, i], Time = times,
                                   Start = factor(rep(i, nrow(thetaEst)),
                                                  levels = 1:ncol(thetaEst)))[sel, ]
                      }))
  thetaMoy <- apply(thetaEst, 2, smoothEst)
  theta_gg_smooth <- do.call(rbind.data.frame,
                             lapply(1:ncol(thetaMoy), function(i){
                               data.frame(Estimate = thetaMoy[, i], Time = times,
                                          Start = factor(rep(i, nrow(thetaMoy)),
                                                         levels = 1:ncol(thetaEst)))[sel, ]
                             }))
  savePdf(name)
  print(foo(theta_gg, legendPosition))
  dev.off()
  embedPdf(name)
  name_smooth <- paste(name, "smooth", sep = "_")
  savePdf(name_smooth)
  print(foo(theta_gg_smooth, legendPosition))
  dev.off()
  embedPdf(name_smooth)
}

names <- c("oneObs_severalStarts", "severalObs_oneStart")
foo_plot("oneObs_severalStarts", sel = seq(1, 5000, by = 5))
foo_plot("severalObs_oneStart", sel = seq(1, 5000, by = 20))
foo_legend(name = "oneObs_oneStart_severalGrads")

# Gradient Steps ----------------------------------------------------------

gradient_df <- data.frame(GradientStep = get_grad_steps(gradient_power = 0.6, cst = 8)[, 1],
                          Time = times)

savePdf("gradient_steps")
ggplot(data = gradient_df, aes(x = Time, y = GradientStep)) + 
  geom_line(col = "gray") + geom_point(size = 0.8) +
  scale_y_continuous(limits = c(0, 0.5)) + labs(y = "Gradient step") + 
  theme_bw() + 
  theme(legend.position = "none",
        axis.title = element_text(size = title_size),
        axis.text = element_text(size = axis_size))
dev.off()
embedPdf("gradient_steps")

# Data --------------------------------------------------------------------

seed <- 122 
simulated_POD <- read.table(paste0("simulated_data/simul_data_seed", 
                                   seed, ".txt"),
                            sep = ";", header = T)
observations <- simulated_POD[, "observations"]

data_df <- data.frame(Time = rep(times, 2), 
                      Value = c(simulated_POD[, "hiddenStates"], 
                                simulated_POD[, "observations"]),
                      Nature = factor(rep(c("True state", "Observation"), 
                                          rep(nrow(simulated_POD), 2))))

savePdf("SINE_data")
ggplot(data = data_df, aes(x = Time, y = Value, colour = Nature)) +
  geom_point(size = 0.2) +
  labs(y = "Value") + 
  guides(colour = guide_legend(override.aes = list(size = 6))) +
  theme_bw() + 
  theme(legend.position = c(0.5, 0.2), legend.background = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = axis_size), 
        axis.title = element_text(size = title_size),
        axis.text = element_text(size = axis_size))
dev.off()
embedPdf("SINE_data")


# Density -----------------------------------------------------------------

foo_dens <- function(name){
  thetaEst <- read.table(paste0("simulation_results/est_", name, ".txt"),
                         h = F, sep = ";")
  thetaEst <- thetaEst[, order(thetaEst[1, ])]
  thetaMoy <- apply(thetaEst, 2, smoothEst)
  theta_gg_smooth <- do.call(rbind.data.frame,
                             lapply(1:ncol(thetaMoy), function(i){
                               data.frame(Estimate = thetaMoy[, i], Time = times,
                                          Start = factor(rep(i, nrow(thetaMoy)),
                                                         levels = 1:ncol(thetaMoy)))
                             }))
  last_pts <- do.call(rbind.data.frame, lapply(split(theta_gg_smooth, 
                                              theta_gg_smooth$Start, drop = T),
                                              function(x) x[nrow(x),]))
  p <- ggplot(last_pts, aes(Estimate)) + labs(y = "Density") +
    geom_density(col = "white", fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = pi / 4, col = "red", linetype = 2) +
    theme_bw() + 
    theme(axis.title = element_text(size = title_size),
          axis.text = element_text(size = axis_size))
  savePdf(paste0(name, "_smooth_dens"))
  print(p)
  dev.off()
  embedPdf(paste0(name, "_smooth_dens"))
}

foo_dens("severalObs_oneStart")



# Version JCGS ------------------------------------------------------------

my_theme <- function(base_size = 8){
  theme_bw()  %+replace%
    theme(
      panel.border = element_rect(colour = "black", 
                                  fill = rgb(0, 0, 0, 0)),
      # plot.background = element_rect(fill = "white"),# bg around panel
      legend.background = element_blank(), 
      text = element_text(family = "LM Roman 10", size = base_size),
      axis.title = element_text(size = rel(1)),
      legend.text = element_text(size = rel(1)),
      legend.title = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      strip.background = element_rect(fill = "lightgoldenrod1",
                                      color = "black"),
      plot.subtitle = element_text(hjust = 0.5, size = rel(1)),
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5))
}
theme_set(my_theme(base_size = 14))


# Figure 4 ----------------------------------------------------------------

seed <- 122
data_plot <- read.table(paste0("simulated_data/simul_data_seed", 
                               seed, ".txt"),
                        sep = ";", header = TRUE) %>% 
  rename(Time = simulationTimes) %>% 
  pivot_longer(-Time, names_to = "Nature", values_to = "Value") %>% 
  mutate(Nature = factor(Nature, levels = c("hiddenStates", "observations"),
                         labels = c("True state", "Observation"))) %>% 
  ggplot(aes(x = Time, y = Value, colour = Nature)) +
  geom_point(size = 0.2) +
  labs(y = "Process value") + 
  scale_color_viridis_d(option = "D") +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  theme(legend.position = c(0.75, 0.9), legend.background = element_blank(),
        legend.title = element_blank())
grad_steps_plot <- data.frame(GradientStep = get_grad_steps(gradient_power = 0.6, cst = 8)[, 1],
                              Time = times) %>% 
  ggplot(aes(x = Time, y = GradientStep)) + 
  geom_line(col = "gray") + geom_point(size = 0.8) +
  scale_y_continuous(limits = c(0, 0.5)) + 
  labs(y = "Gradient step")
estimate_plot <- read.table(paste0("simulation_results/est_oneObs_severalStarts_IS.txt"),
                            header = FALSE, sep = ";") %>% 
  mutate(Time = times) %>% # Creation of times
  filter((Time %% 20) == 0) %>% # Select only part of estimates
  pivot_longer(-Time, names_to = "Start", values_to = "Estimate", names_prefix = "V") %>% 
  mutate(Start = factor(Start) %>% 
           fct_reorder(Estimate, .fun = function(x) x[1])) %>% 
  ggplot() +
  aes(x = Time, y = Estimate, color = Start) +
  scale_color_viridis_d() + 
  geom_path() +
  theme(legend.position = "none") +
  geom_hline(yintercept = trueTheta, color = "red")
gridExtra::grid.arrange(data_plot, grad_steps_plot, estimate_plot,
                        nrow = 1) %>% 
  ggsave(filename = "../../21-backwardIS/JCGS/Revision/Figure4.pdf", plot = ., device = cairo_pdf, 
         width = 12, height = 4)


# Figure 5 ----------------------------------------------------------------

estimation_df <- read.table("simulation_results/est_severalObs_oneStart.txt",
                            sep = ";") %>% 
  mutate(Time = times) %>% # Creation of times
  filter((Time %% 20) == 0 | Time == 4999) %>% # Select only part of estimates
  pivot_longer(-Time, names_to = "DataSet", values_to = "Estimate", names_prefix = "V") %>% 
  mutate(DataSet = factor(DataSet) %>% 
           fct_reorder(Estimate, .fun = function(x) x[1])) %>% 
  group_by(DataSet) %>% 
  mutate(Smooth = smoothEst(Estimate)) %>% 
  ungroup()

estimate_plot <- estimation_df %>% 
  ggplot() +
  aes(x = Time, y = Estimate, color = DataSet) +
  scale_color_viridis_d() + 
  geom_path() +
  theme(legend.position = "none") +
  geom_hline(yintercept = trueTheta, color = "red")
smoothed_plot <- estimation_df %>% 
  ggplot() +
  aes(x = Time, y = Smooth, color = DataSet) +
  scale_color_viridis_d() + 
  geom_path() +
  theme(legend.position = "none") +
  labs(y = "Estimate") +
  geom_hline(yintercept = trueTheta, color = "red")
density_plot <- estimation_df %>% 
  filter(Time == 4999) %>% 
  ggplot() +
  aes(x = Smooth) +
  geom_density(color = "white", fill = "steelblue", alpha = .5) +
  labs(x = "Estimate", y = "Density") +
  geom_vline(xintercept = trueTheta, col = "red", linetype = 2)
gridExtra::grid.arrange(estimate_plot, smoothed_plot, density_plot,
                        nrow = 1) %>% 
  ggsave(filename = "../../21-backwardIS/JCGS/Figure5.pdf", 
         plot = ., device = cairo_pdf, 
         width = 12, height = 4)


# Figure 6 ----------------------------------------------------------------


estimation_df <- read.table("simulation_results/est_oneObs_oneStart_severalGrads.txt",
                            sep = ";") %>% 
  mutate(Time = times) %>% # Creation of times
  filter((Time %% 20) == 0 | Time == 4999) %>% # Select only part of estimates
  pivot_longer(-Time, names_to = "kappa", 
               values_to = "Estimate", names_prefix = "V") %>% 
  mutate(kappa = factor(kappa, levels = paste(1:6), 
                        labels = paste(seq(.5, 1, by = .1)))) %>% 
  group_by(kappa) %>% 
  mutate(Smooth = smoothEst(Estimate)) %>% 
  ungroup()
estimate_plot <- estimation_df %>% 
  ggplot() +
  aes(x = Time, y = Estimate, color = kappa) +
  scale_color_viridis_d() + 
  geom_path() +
  theme(legend.position = "none") +
  geom_hline(yintercept = trueTheta, color = "red")
smoothed_plot <- estimation_df %>% 
  ggplot() +
  aes(x = Time, y = Smooth, color = kappa) +
  scale_color_viridis_d() + 
  geom_path() +
  theme(legend.position = c(.8, .7)) +
  labs(color = expression(kappa), y = "Estimate") +
  geom_hline(yintercept = trueTheta, color = "red")
gridExtra::grid.arrange(estimate_plot, smoothed_plot,
                        nrow = 1) %>% 
  ggsave(filename = "../../21-backwardIS/JCGS/Figure6.pdf", 
         plot = ., device = cairo_pdf, 
         width = 12, height = 4)
