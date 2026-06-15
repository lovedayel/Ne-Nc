library(ggplot2)
library(ggpubr)
library(ape)
library(caper)
library(dplyr)
library(geiger)
library(phylolm) 

#read in Newick file, check for and resolve polytomies, make sure tree is ultrametric
tree = read.tree("clean_data/phylo_final_Nc.nwk")
resolved_tree = multi2di(tree, random = TRUE)
resolved_tree = chronos(resolved_tree)

#read in data
data = read.csv("clean_data/Nc_final.csv", h = TRUE, row.names = 1)
data$species <- rownames(data)

#filter and log transform data
filtered_data <- data %>%
  filter(!is.na(Ne), !is.na(Ne_lower), !is.na(Ne_upper), !is.na(Nc_estimate)) %>%
  mutate(
    log_Ne = log10(Ne),
    log_Nc = log10(Nc_estimate),
    log_Ne_Nc_ratio = log10(Ne / Nc_estimate)
  ) %>%
  select(-Nc_lower, -Nc_upper, -Nc_source)

#filter data to only include species in the tree
filtered_data <- filtered_data[filtered_data$species %in% resolved_tree$tip.label, ]

#set colour scheme
IUCN_colors <- c(
  "least_concern" = "#60C659",
  "near_threatened" = "#CCE226",
  "vulnerable" = "#F9E814",
  "endangered" = "#FC7F3F",
  "critically_endangered" = "#D81E05",
  "not_assessed" = "#D1D1C6"
)

#Fig 1A: Ne vs. Nc with PGLS correction

fig1A_data <- filtered_data
fig1A_data$IUCN_status[fig1A_data$IUCN_status == ""] <- "not_assessed"

#create comparative data object for caper
comp_data = comparative.data(phy = resolved_tree, data = fig1A_data, 
                             names.col = "species", vcv = TRUE, warn.dropped = TRUE)

#fit PGLS model
pgls_model = pgls(log_Ne ~ log_Nc, data = comp_data, lambda = "ML")
model_summary = summary(pgls_model)

#get PGLS coefficient and p-value
pgls_coef_slope = round(model_summary$coefficients[2, "Estimate"], 2)
pgls_p = signif(model_summary$coefficients[2, "Pr(>|t|)"], 2)
r2 = round(model_summary$adj.r.squared, 2)

plot_data = data.frame(
  log_Nc = comp_data$data$log_Nc, 
  observed_values = comp_data$data$log_Ne, 
  fitted_values = fitted(pgls_model),
  IUCN_status = factor(comp_data$data$IUCN_status, 
                       levels = c("critically_endangered", "endangered",
                                  "vulnerable", "near_threatened", 
                                  "least_concern", "not_assessed"))
)

#plot Figure 1A
pdf("fig1A.pdf", width = 8, height = 4.5)
fig1A = ggplot(plot_data, aes(x = log_Nc, y = observed_values, color = IUCN_status)) +
  geom_point(size = 3) +  
  geom_line(aes(y = fitted_values), colour = "blue", linewidth = 1) +  
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey", linewidth = 1) +  
  scale_x_continuous(labels = function(x) parse(text = paste0("10^", x))) + 
  scale_y_continuous(labels = function(y) parse(text = paste0("10^", y))) + 
  scale_color_manual(
    name = "IUCN Red List Status",
    values = IUCN_colors,
    labels = c("critically_endangered" = "Critically Endangered",
               "endangered" = "Endangered", "vulnerable" = "Vulnerable",
               "near_threatened" = "Near Threatened", "least_concern" = "Least Concern",
               "not_assessed" = "Not Assessed")
  ) + 
  labs(
    x = expression("Log10 Census population size (" * italic(N)[c] * ")"), 
    y = expression("Log10 Effective population size (" * italic(N)[e] * ")")
  ) +
  coord_cartesian(ylim = c(NA, 10)) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key = element_blank(),
    legend.position = "right"
  ) +
  annotate("text", x = -Inf, y = Inf, 
           label = paste0("      PGLS Coef = ", pgls_coef_slope, ", p = ", pgls_p, "\n      R² = ", r2), 
           hjust = 0, vjust = 1.7, size = 4.5)

print(fig1A)
dev.off()


#Fig 1B: relationship between Ne/Nc ratios and IUCN status

IUCN_order = c("least_concern", "near_threatened", "vulnerable", "endangered", "critically_endangered")

data_IUCN <- filtered_data %>%
  filter(!is.na(IUCN_status) & IUCN_status != "" & IUCN_status != "not_assessed") %>%
  mutate(IUCN_status = factor(IUCN_status, levels = IUCN_order))

rownames(data_IUCN) <- data_IUCN$species

obj_anova <- name.check(resolved_tree, data_IUCN)
tree_anova <- drop.tip(resolved_tree, obj_anova$tree_not_data)
anova_data <- data_IUCN[rownames(data_IUCN) %in% tree_anova$tip.label, ]

#perform one-way ANOVA for all categories
anova_result <- aov(log_Ne_Nc_ratio ~ IUCN_status, data = anova_data)
cat("\n--- ANOVA Results ---\n")
print(summary(anova_result))

#proceed with pairwise comparisons if ANOVA is significant
run_pairwise_phylolm <- function(status1, status2, full_data, tree) {
  sub_data <- full_data %>%
    filter(IUCN_status %in% c(status1, status2)) %>%
    mutate(comparison_group = factor(IUCN_status))
  
  chk <- name.check(tree, sub_data)
  sub_tree <- drop.tip(tree, chk$tree_not_data)
  sub_data <- sub_data[rownames(sub_data) %in% sub_tree$tip.label, ]
  
  model <- phylolm(log_Ne_Nc_ratio ~ comparison_group, data = sub_data, phy = sub_tree, model = "BM")
  return(model)
}

#run pairwise comparisons
pgls_EN_vs_LC <- run_pairwise_phylolm("endangered", "least_concern", data_IUCN, resolved_tree)
pgls_CR_vs_LC <- run_pairwise_phylolm("critically_endangered", "least_concern", data_IUCN, resolved_tree)

cat("\n--- Phylolm: Endangered vs Least Concern ---\n")
print(summary(pgls_EN_vs_LC))
cat("\n--- Phylolm: Critically Endangered vs Least Concern ---\n")
print(summary(pgls_CR_vs_LC))

#plot fig 1B
IUCN_abbreviations <- c(
  "least_concern" = "LC",
  "near_threatened" = "NT",
  "vulnerable" = "VU",
  "endangered" = "EN",
  "critically_endangered" = "CR"
)

pdf("fig1B.pdf", width = 8, height = 4.5)
fig1B <- ggplot(data_IUCN, aes(x = IUCN_status, y = log_Ne_Nc_ratio, fill = IUCN_status)) +
  geom_boxplot() +
  labs(
    x = "IUCN Red List Status", 
    y = expression("Log10 " * italic(N)[e] * "/" * italic(N)[c] * " Ratio")
  ) +
  scale_x_discrete(labels = IUCN_abbreviations) +  
  scale_fill_manual(
    name = "IUCN Red List Status",
    values = IUCN_colors,
    labels = c("least_concern" = "Least Concern", "near_threatened" = "Near Threatened",
               "vulnerable" = "Vulnerable", "endangered" = "Endangered",
               "critically_endangered" = "Critically Endangered")
  ) +
  theme(
    panel.background = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.key = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(fig1B)
dev.off()