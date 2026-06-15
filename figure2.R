library(ggplot2)
library(patchwork)

#function for Equation 1
calc_r_tau <- function(tau, N_ratio) {
  (1 + N_ratio * (exp(tau) - 1)) * exp(-tau)
}

#define the 9 decline scenarios (specify ancestral population sizes + decline magnitudes)
scenarios <- expand.grid(
  N0 = c(100000, 1000000, 10000000),
  FoldDecline = c(10, 100, 1000)
)

scenarios$N <- scenarios$N0 / scenarios$FoldDecline
scenarios$Ratio <- scenarios$N / scenarios$N0

#generate data for time in absolute generations; long-term
t_seq <- seq(0, 10000, length.out = 300)
plot_data_absolute <- data.frame()

for (i in 1:nrow(scenarios)) {
  tau_scenario <- t_seq / (2 * scenarios$N[i])
  temp_df <- data.frame(
    t = t_seq,
    r_tau = calc_r_tau(tau_scenario, scenarios$Ratio[i]),
    N0 = as.factor(scenarios$N0[i]),
    FoldDecline = as.factor(scenarios$FoldDecline[i])
  )
  plot_data_absolute <- rbind(plot_data_absolute, temp_df)
}

#zoom in on early decline phase
t_seq <- seq(0, 200, length.out = 300)
plot_data_zoomed <- data.frame()

for (i in 1:nrow(scenarios)) {
  tau_scenario <- t_seq / (2 * scenarios$N[i])
  temp_df <- data.frame(
    t = t_seq,
    r_tau = calc_r_tau(tau_scenario, scenarios$Ratio[i]),
    N0 = as.factor(scenarios$N0[i]),
    FoldDecline = as.factor(scenarios$FoldDecline[i])
  )
  plot_data_zoomed <- rbind(plot_data_zoomed, temp_df)
}


#plot A: absolute time; long-term
p1 <- ggplot(plot_data_absolute, aes(x = t, y = r_tau, color = FoldDecline, linetype = N0)) +
  annotate("rect", xmin = 0, xmax = 200, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.4) +
  geom_line(linewidth = 1, alpha = 0.9) +
  scale_linetype_manual(values = c("1e+05" = "solid", "1e+06" = "dashed", "1e+07" = "dotted")) +
  scale_color_manual(values = c("10" = "#56B4E9", "100" = "#E69F00", "1000" = "#D55E00")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_x_continuous(n.breaks = 5, expand = c(0, 5), limits = c(0, 10000)) +
  labs(x = "Generations (t)", y = expression(r[tau] ~ "Relative Diversity")) +
  theme(
  	panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  ) 

#plot B: zoomed 
p2 <- ggplot(plot_data_zoomed, aes(x = t, y = r_tau, color = FoldDecline, linetype = N0)) +
  annotate("rect", xmin = 0, xmax = 200, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.4) +
  geom_line(linewidth = 1, alpha = 0.9) +
  scale_linetype_manual(values = c("1e+05" = "solid", "1e+06" = "dashed", "1e+07" = "dotted")) +
  scale_color_manual(values = c("10" = "#56B4E9", "100" = "#E69F00", "1000" = "#D55E00")) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.25), limits = c(0, 1), expand = c(0.01, 0.01)) +
  scale_x_continuous(n.breaks = 5, expand = c(0, 0), limits = c(0, 200)) +
  labs(x = "Generations (t)", y = NULL, color = "Fold Decrease", linetype = "Ancestral N") +
  theme(
  	panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key.width = unit(1, "cm")
  )

combined_plot <- p1 + p2 +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(theme = theme(plot.title = element_text(face = "bold", size = 16, hjust = 0.5)))

ggsave("figure2.pdf", plot = combined_plot, width = 13, height = 6)