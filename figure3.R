library(ggplot2)
library(dplyr)
library(ggrepel)
library(ape)
library(caper)
library(geiger)
library(ggpubr)
library(tidyr)
library(scales)

IUCN_colours <- c("least_concern" = "#60C659", "near_threatened" = "#CCE226", 
                 "vulnerable" = "#F9E814", "endangered" = "#FC7F3F", 
                 "critically_endangered" = "#D81E05", "not_assessed" = "#D1D1C6")

theme_minimal_axes <- theme(
  panel.background = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line(color = "black", linewidth = 0.4),
  axis.ticks = element_line(color = "black", linewidth = 0.4),
  legend.position = "top",
  text = element_text(size = 12)
)

#read in and prepare phylogeny
tree <- read.tree("clean_data/phylo_final_Nc.nwk")
tree <- multi2di(tree, random = TRUE)
is.binary(tree) 
tree <- chronos(tree)
is.ultrametric(tree) 

#read in data
data_Nc <- read.csv("clean_data/Nc_final.csv", header = TRUE)
data_combined <- data_Nc %>%
  mutate(
    log10_Ne = log10(Ne),
    log10_Nc_estimate = log10(Nc_estimate),
    log_Ne_Nc_ratio = log10(Ne / Nc_estimate),
    log10_gen_time_yr = log10(gen_time_yr)
  )

#match data to tree and create comparative object for PGLS
tree_species <- tree$tip.label
data_combined <- data_combined %>% filter(species %in% tree_species)
data_combined <- as.data.frame(data_combined)
rownames(data_combined) <- data_combined$species

comp_data <- comparative.data(
  phy = tree, 
  data = data_combined, 
  names.col = "species", 
  vcv = TRUE, 
  warn.dropped = TRUE,
  na.omit = FALSE)

all_species_pred <- comp_data$data
all_species_pred$species <- rownames(all_species_pred)
all_species_pred <- all_species_pred %>% drop_na(log10_gen_time_yr)
n_sims <- 10000


#METHOD 1: Ne/Nc ~ Generation Time 

ne_nc_data <- comp_data$data
ne_nc_data$species <- rownames(ne_nc_data)

ne_nc_data <- ne_nc_data %>%
  filter(IUCN_status %in% c("least_concern", "not_assessed"),
         species != "Homo_sapiens") %>%
  drop_na(log10_gen_time_yr, log_Ne_Nc_ratio)

comp_ne_nc <- comparative.data(phy = comp_data$phy, data = ne_nc_data, 
                               names.col = "species", vcv = TRUE, 
                               warn.dropped = TRUE, na.omit = FALSE)

#fit PGLS model
ne_nc_model <- pgls(log_Ne_Nc_ratio ~ log10_gen_time_yr, data = comp_ne_nc, lambda = "ML")

#extract summary statistics 
summary_ne_nc   <- summary(ne_nc_model)
r2_ne_nc        <- round(summary_ne_nc$adj.r.squared, 3)
pval_ne_nc      <- signif(summary_ne_nc$coefficients[2, 4], 3)
lambda_ne_nc    <- round(ne_nc_model$param["lambda"], 3)
coef_mean_ne_nc <- coef(ne_nc_model)

coef_vcov_ne_nc <- as.numeric(ne_nc_model$RMS) * solve(t(ne_nc_model$x) %*% solve(ne_nc_model$Vt, ne_nc_model$x))
slope_se_ne_nc  <- sqrt(diag(coef_vcov_ne_nc))["log10_gen_time_yr"]

set.seed(42)
sim_slopes_ne_nc <- rnorm(n_sims, mean = coef_mean_ne_nc["log10_gen_time_yr"], sd = slope_se_ne_nc)

comp_ne_nc$data$fitted <- fitted(ne_nc_model)

#ribbon data for Supplementary Figure 2
x_seq_ne_nc <- seq(min(comp_ne_nc$data$log10_gen_time_yr), max(comp_ne_nc$data$log10_gen_time_yr), length.out = 100)
exact_preds_ne_nc <- coef_mean_ne_nc[1] + coef_mean_ne_nc[2] * x_seq_ne_nc
sim_preds_ne_nc   <- outer(x_seq_ne_nc, sim_slopes_ne_nc, FUN = function(x, s) coef_mean_ne_nc[1] + s * x)
pred_sd_ne_nc     <- apply(sim_preds_ne_nc, 1, sd)

ribbon_data_ne_nc <- data.frame(
  log10_gen_time_yr = x_seq_ne_nc,
  fitted            = exact_preds_ne_nc,
  ci_lower          = exact_preds_ne_nc - 1.96 * pred_sd_ne_nc,
  ci_upper          = exact_preds_ne_nc + 1.96 * pred_sd_ne_nc
)

#Supplementary Figure 2
pdf("predicting_Nc_figs/Supp_Fig_2_Ne_Nc_vs_gen_time.pdf", width = 10, height = 8)
ggplot(comp_ne_nc$data, aes(x = log10_gen_time_yr)) +
  geom_ribbon(data = ribbon_data_ne_nc,
              aes(x = log10_gen_time_yr, ymin = ci_lower, ymax = ci_upper),
              fill = "blue", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = ribbon_data_ne_nc,
            aes(x = log10_gen_time_yr, y = fitted),
            color = "blue", linewidth = 1, inherit.aes = FALSE) +
  geom_point(aes(y = log_Ne_Nc_ratio), size = 3, alpha = 0.7, color = "#333333") +
  labs(x = "Log10 Generation Time (years)",
       y = "Log10(Ne / Nc)",
       subtitle = paste("PGLS: log10(Ne/Nc) ~ log10(gen_time), R² =",
                        r2_ne_nc, ", p =", pval_ne_nc)) +
  theme_minimal_axes
dev.off()


#predict ancestral Ne/Nc & ancestral Nc using METHOD 1
all_species_pred$predicted_log10_Ne_Nc_ancestral <- 
  coef_mean_ne_nc[1] + coef_mean_ne_nc[2] * all_species_pred$log10_gen_time_yr

all_species_pred$predicted_log10_Nc_method1 <- log10(all_species_pred$Ne) - all_species_pred$predicted_log10_Ne_Nc_ancestral
all_species_pred$predicted_Nc_method1 <- 10^all_species_pred$predicted_log10_Nc_method1

#check if current Ne/Nc is above or below the regression line for plotting
all_species_pred <- all_species_pred %>%
  mutate(
    Ne_Nc_above_line = log_Ne_Nc_ratio > predicted_log10_Ne_Nc_ancestral,
    species_label = gsub("^([A-Z])[a-z]*_(.*)", "\\1. \\2", species)
  )

#Figure 3: current Ne/Nc vs predicted ancestral Ne/Nc for all species
pdf("predicting_Nc_figs/Figure_3_predicted_vs_observed_NeNc.pdf", width = 11, height = 8)
plot_fig3 <- ggplot(all_species_pred, aes(x = log10_gen_time_yr)) +
  geom_point(aes(y = predicted_log10_Ne_Nc_ancestral, color = IUCN_status), 
             size = 3.5, shape = 1, stroke = 1, alpha = 1) +
  geom_point(aes(y = log_Ne_Nc_ratio, color = IUCN_status, fill = IUCN_status, shape = Ne_Nc_above_line), 
             size = 3.5, stroke = 0.5) +
  geom_segment(aes(x = log10_gen_time_yr, y = predicted_log10_Ne_Nc_ancestral,
                   xend = log10_gen_time_yr, yend = log_Ne_Nc_ratio, color = IUCN_status),
               linetype = "dashed", alpha = 1) +
  geom_line(data = ribbon_data_ne_nc,
            aes(x = log10_gen_time_yr, y = fitted),
            color = "blue", linewidth = 1, inherit.aes = FALSE) +
  geom_ribbon(data = ribbon_data_ne_nc,
              aes(x = log10_gen_time_yr, ymin = ci_lower, ymax = ci_upper),
              fill = "blue", alpha = 0.15, inherit.aes = FALSE) +
  geom_text_repel(data = subset(all_species_pred, IUCN_status %in% c("endangered", "critically_endangered", "vulnerable")),
                  aes(y = log_Ne_Nc_ratio, label = species_label), 
                  size = 3, fontface = "italic", box.padding = 0.8, point.padding = 0.5, min.segment.length = 0, segment.color = "grey50", segment.size = 0.3, max.overlaps = Inf) +
  scale_color_manual(values = IUCN_colours) +
  scale_fill_manual(values = IUCN_colours, guide = "none") +
  scale_shape_manual(values = c("FALSE" = 25, "TRUE" = 24), guide = "none") +
  labs(x = "Log10 Generation Time (years)",
       y = "Log10(Ne/Nc)",
       color = "IUCN Status") +
  theme_minimal_axes
print(plot_fig3)
dev.off()


#METHOD 2: Nc ~ Generation Time

lc_na <- comp_data$data
lc_na$species <- rownames(lc_na)

lc_na <- lc_na %>%
  filter(IUCN_status %in% c("least_concern", "not_assessed"),
         species != "Homo_sapiens") %>%
  drop_na(log10_gen_time_yr, Nc_estimate)

comp_lc_na <- comparative.data(phy = comp_data$phy, data = lc_na,
                               names.col = "species", vcv = TRUE,
                               warn.dropped = TRUE, na.omit = FALSE)

nc_model <- pgls(log10_Nc_estimate ~ log10_gen_time_yr, data = comp_lc_na, lambda = "ML")

nc_model_summary <- summary(nc_model)
r2_nc        <- round(nc_model_summary$adj.r.squared, 3)
pval_nc      <- signif(nc_model_summary$coefficients[2, 4], 3)
lambda_nc    <- round(nc_model$param["lambda"], 3)
coef_mean_nc <- coef(nc_model)

coef_vcov_nc <- as.numeric(nc_model$RMS) * solve(t(nc_model$x) %*% solve(nc_model$Vt, nc_model$x))
slope_se_nc  <- sqrt(diag(coef_vcov_nc))["log10_gen_time_yr"]

set.seed(21)
sim_slopes_nc <- rnorm(n_sims, mean = coef_mean_nc["log10_gen_time_yr"], sd = slope_se_nc)

comp_lc_na$data$fitted <- fitted(nc_model)

x_seq_nc <- seq(min(comp_lc_na$data$log10_gen_time_yr), max(comp_lc_na$data$log10_gen_time_yr), length.out = 100)
exact_preds_nc <- coef_mean_nc[1] + coef_mean_nc[2] * x_seq_nc
sim_preds_nc   <- outer(x_seq_nc, sim_slopes_nc, FUN = function(x, s) coef_mean_nc[1] + s * x)
pred_sd_nc     <- apply(sim_preds_nc, 1, sd)

ribbon_data_nc <- data.frame(
  log10_gen_time_yr = x_seq_nc,
  fitted            = exact_preds_nc,
  ci_lower          = exact_preds_nc - 1.96 * pred_sd_nc,
  ci_upper          = exact_preds_nc + 1.96 * pred_sd_nc
)


#Supplementary Figure 3
pdf("predicting_Nc_figs/Supp_Fig_3_Nc_vs_gen_time.pdf", width = 10, height = 8)
ggplot(comp_lc_na$data, aes(x = log10_gen_time_yr)) +
  geom_ribbon(data = ribbon_data_nc,
              aes(x = log10_gen_time_yr, ymin = ci_lower, ymax = ci_upper),
              fill = "blue", alpha = 0.15, inherit.aes = FALSE) +
  geom_point(aes(y = log10_Nc_estimate), size = 3, alpha = 0.7, color = "#333333") +
  geom_line(data = ribbon_data_nc,
            aes(x = log10_gen_time_yr, y = fitted),
            color = "blue", linewidth = 1, inherit.aes = FALSE) +
  labs(x = "Log10 Generation Time (years)",
       y = "Log10 Census Population Size (Nc)",
       subtitle = paste("PGLS: log10(Nc) ~ log10(gen_time_yr), R² =", r2_nc,
                        ", p =", pval_nc)) +
  theme_minimal_axes
dev.off()


#predict ancestral Nc using METHOD 2
all_species_pred$predicted_log10_Nc_method2 <- coef_mean_nc[1] + coef_mean_nc[2] * all_species_pred$log10_gen_time_yr
all_species_pred$predicted_Nc_method2       <- 10^all_species_pred$predicted_log10_Nc_method2


#Supplementary Figure 4
pdf("predicting_Nc_figs/Supp_Fig_4_predicted_vs_observed_Nc.pdf", width = 12, height = 8)
plot_supp4 <- ggplot(all_species_pred, aes(x = log10_gen_time_yr)) +  
  geom_ribbon(data = ribbon_data_nc,
              aes(x = log10_gen_time_yr, ymin = ci_lower, ymax = ci_upper),
              fill = "blue", alpha = 0.15, inherit.aes = FALSE) +
  geom_line(data = ribbon_data_nc,
            aes(x = log10_gen_time_yr, y = fitted),
            color = "blue", linewidth = 1, inherit.aes = FALSE) +
  geom_point(aes(y = predicted_log10_Nc_method2, color = IUCN_status), size = 3, shape = 1) +
  geom_point(aes(y = log10_Nc_estimate,
                 shape = log10_Nc_estimate > predicted_log10_Nc_method2,
                 color = IUCN_status, fill = IUCN_status),
             size = 3.5, stroke = 0.05) +
  geom_segment(aes(x = log10_gen_time_yr, y = log10_Nc_estimate,
                   xend = log10_gen_time_yr, yend = predicted_log10_Nc_method2,
                   color = IUCN_status),
               linetype = "dashed", alpha = 1) +
  scale_color_manual(values = IUCN_colours) +
  scale_fill_manual(values = IUCN_colours, guide = "none") +
  scale_shape_manual(values = c(25, 24), guide = "none") +
  labs(x = "Log10 Generation Time (years)",
       y = "Log10 Census Population Size (Nc)",
       color = "IUCN Status") +
  guides(colour = guide_legend(override.aes = list(shape = 16))) +
  theme_minimal_axes +
  theme(legend.position = "top", plot.subtitle = element_text(size = 10, color = "gray30"))
print(plot_supp4)
dev.off()


#calculate fold Nc declines for histograms
comparison_df <- all_species_pred %>%
  filter(IUCN_status %in% c("endangered", "critically_endangered", "vulnerable")) %>%
  dplyr::select(species, IUCN_status, log10_gen_time_yr, log10_Nc_estimate,
                predicted_log10_Nc_method1, 
                predicted_log10_Nc_method2) %>%
  mutate(
    method1_fold = 10^(predicted_log10_Nc_method1 - log10_Nc_estimate),
    method2_fold = 10^(predicted_log10_Nc_method2 - log10_Nc_estimate)
  )

#Figure 4: histogram of Nc declines (METHOD 1) 
pdf("predicting_Nc_figs/Figure_4_histogram_Nc_decline_method1.pdf", width = 7, height = 6)
hist_method1 <- ggplot(comparison_df, aes(x = log10(method1_fold), fill = IUCN_status)) +
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.9) +
  scale_fill_manual(values = IUCN_colours) +
  scale_x_continuous(labels = function(x) comma(10^x)) +
  labs(x = "Ancestral Nc / Current Nc (Fold Decline - Method 1)", 
       y = "Number of species", fill = "IUCN Status") +
  theme_minimal_axes +
  theme(legend.position = "top")
print(hist_method1)
dev.off()

#Supplementary Figure 5
pdf("predicting_Nc_figs/Supp_Fig_5_histogram_Nc_decline_method2.pdf", width = 7, height = 7)
hist_method2 <- ggplot(comparison_df, aes(x = log10(method2_fold), fill = IUCN_status)) +
  geom_histogram(binwidth = 0.5, color = "black", alpha = 0.9) +
  scale_fill_manual(values = IUCN_colours) +
  scale_x_continuous(labels = function(x) comma(10^x)) +
  labs(x = "Ancestral Nc / Current Nc (Fold Decline - Method 2)", 
       y = "Number of species", fill = "IUCN Status") +
  theme_minimal_axes +
  theme(legend.position = "top")
print(hist_method2)
dev.off()

#correlation between METHOD 1 and METHOD 2 decline magnitude estimates
comparison_clean <- comparison_df %>%
 filter(!is.na(predicted_log10_Nc_method1) & !is.na(predicted_log10_Nc_method2))

pearson_res <- cor.test(comparison_clean$predicted_log10_Nc_method1, 
                        comparison_clean$predicted_log10_Nc_method2, 
                        method = "pearson")

#Supplementary Figure 6
pdf("predicting_Nc_figs/Supp_Fig_6_method1_vs_method2_correlation.pdf", width = 7, height = 7)
plot_supp6 <- ggplot(comparison_clean, aes(x = predicted_log10_Nc_method2, y = predicted_log10_Nc_method1)) +
  geom_point(aes(color = IUCN_status), size = 3) +
  geom_smooth(method = "lm", color = "blue", se = TRUE, alpha = 0.1) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top") +
  scale_color_manual(values = IUCN_colours) +
  labs(x = "Method 2: Log10 Predicted Nc", y = "Method 1: Log10 Predicted Nc") +
  theme_minimal_axes
print(plot_supp6)
dev.off()
