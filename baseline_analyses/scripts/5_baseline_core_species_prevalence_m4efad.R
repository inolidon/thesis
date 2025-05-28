## code summary 
## identify core species present in MAM v well-nourished 
## data can be directly loaded by using the files in the "save datasets" section 



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
library(AMR)



# load data  --------------------------------------------------------------

# baseline metadata 
load("data/processed/baseline_data.Rdata")
# baseline species - raw 
load("data/processed/baseline_species.Rdata")



# Core species - Prevalence of >85% ------------------------------------------------------

# Join with metadata to get Condition
species_long <- baseline_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  left_join(baseline_data %>% select(Seq_ID, Condition), by = "Seq_ID")

# Calculate prevalence by Condition group
core_species_by_group <- species_long %>%
  group_by(Condition, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.85)

# Clean species names
core_species_by_group <- core_species_by_group %>%
  mutate(Taxa = gsub("^s__", "", Taxa),        # remove s__ prefix
         Taxa = gsub("_", " ", Taxa))          # replace _ with space

# View core species per group
core_species_by_group %>% group_by(Condition) %>% summarise(n_core = n_distinct(Taxa))


# Plot --------------------------------------------------------------------


core_prevalent_species <- ggplot(core_species_by_group, aes(x = prevalence, y = reorder(Taxa, prevalence), fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  labs(
    x = "Prevalence",
    y = "",
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    plot.title = element_text(size = 9),
    legend.position = "none"
  )





# pdf
ggsave(file = "plots/chp1/clean_plots/core_prevalent_species.pdf", 
       plot = core_prevalent_species, units = "cm", width=12, height=9)



# Rothi + Strep - Prevalence  ---------------------------------------------

# Define species of interest
target_species <- c("s__Rothia_mucilaginosa", "s__Streptococcus_salivarius")

# Filter and calculate prevalence by Condition
species_prevalence_mam_healthy <- baseline_species %>%
  select(Seq_ID, all_of(target_species)) %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  inner_join(baseline_data %>% select(Seq_ID, Condition), by = "Seq_ID") %>%
  group_by(Condition, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# View result
species_prevalence_mam_healthy


# plot --------------------------------------------------------------------

# Prepare labels nicely
species_prevalence_mam_healthy <- species_prevalence_mam_healthy %>%
  mutate(
    Taxa = gsub("^s__", "", Taxa),
    Taxa = gsub("_", " ", Taxa)
  )

# Plot
prevalence_plot_rothia_strep <- ggplot(species_prevalence_mam_healthy, aes(x = Taxa, y = prevalence, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.75), color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Prevalence")

# Save PDF
ggsave(file = "plots/chp1/clean_plots/prevalence_plot_rothia_strep.pdf", 
       plot = prevalence_plot_rothia_strep, units = "cm", width = 9, height = 9)






