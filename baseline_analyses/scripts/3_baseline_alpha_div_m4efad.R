## Code summary:
# This script calculates alpha diversity metrics (Shannon, Richness, Evenness) from species-level data,
# compares them between MAM and Healthy using Wilcoxon tests,
# summarizes the results (mean ± SD),
# and visualizes all three metrics in a single faceted plot.





# load packages  ----------------------------------------------------------

library(tidyverse) # required for data wrangling and ggplot
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(Maaslin2) # required for differential pathway analysis
library(smplot2) # fits line to scatter plots and adds pearson correlation
library(scales) # required for label_comma()
library(data.table) # require for gene table cleaning
library(dplyr)
library(stringr)
library(readxl)
library(reshape2)
library(igraph)
library(pheatmap)
library(stats)
library(broom)




# load data  --------------------------------------------------------------

# baseline metadata 
load("data/processed/baseline_data.Rdata")
# baseline species - raw 
load("data/processed/baseline_species.Rdata")


# Check if baseline_species has row names, and if so, remove them
baseline_species <- baseline_species %>% 
  column_to_rownames("Seq_ID")

#calculate the Richness, Shannon and evenness of species 
alpha_div_mam_healthy <- tibble(
  Seq_ID = rownames(baseline_species),
  Richness = vegan::specnumber(baseline_species),
  Shannon = vegan::diversity(baseline_species, index = "shannon"),
  Evenness = Shannon / log(Richness)
) %>%
  inner_join(baseline_data)



# summary (mean + sd) -----------------------------------------------------

alpha_div_summary <- alpha_div_mam_healthy %>%
  group_by(Condition) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    mean_richness = mean(Richness, na.rm = TRUE),
    sd_richness = sd(Richness, na.rm = TRUE),
    mean_evenness = mean(Evenness, na.rm = TRUE),
    sd_evenness = sd(Evenness, na.rm = TRUE),
    .groups = "drop"
  )




# stats -------------------------------------------------------------------

# Perform the Wilcoxon rank-sum test
wilcox_mam_healthy <- wilcox.test(Shannon ~ Condition, data = alpha_div_mam_healthy) # p val = 0.2909

wilcox_mam_healthy_richness <- wilcox.test(Richness ~ Condition, data = alpha_div_mam_healthy) # p = 0.3181

wilcox_mam_healthy_eveness <- wilcox.test(Evenness ~ Condition, data = alpha_div_mam_healthy) # p = 0.2279


# plots  ------------------------------------------------------------------

# all the three metrics of interest are in a single facetted plot 

# Reshape to long format
alpha_div_long <- alpha_div_mam_healthy %>%
  select(Seq_ID, Condition, Shannon, Richness, Evenness) %>%
  pivot_longer(cols = c(Shannon, Richness, Evenness),
               names_to = "Alpha_metric",
               values_to = "Value")

# Set order of metrics
alpha_div_long$Alpha_metric <- factor(
  alpha_div_long$Alpha_metric,
  levels = c("Shannon", "Richness", "Evenness")
)

# Plot
alpha_div_combined_plot <- ggplot(alpha_div_long, aes(x = Condition, y = Value, fill = Condition)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Condition)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3, position = position_dodge(width = 0.75)) +
  facet_wrap(~ Alpha_metric, scales = "free_y") +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Alpha Diversity Metric")


# pdf
ggsave(file = "plots/chp1/clean_plots/alpha_div_combined_plot.pdf", 
       plot = alpha_div_combined_plot, units = "cm", width=10, height=9)






