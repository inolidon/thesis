## code summary 
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
library(patchwork)


# load data ---------------------------------------------------------------
# metadata
# refeed_data: contains paired data for 12mo and 15mo
load("data/processed/refeed_data.Rdata")

# contains recovery data information 
load("data/processed/recovery_data.Rdata")

# contains group related data 
load("data/processed/group_data.Rdata")

# species - raw for each variable 
load("data/processed/refeed_species.Rdata")
load("data/processed/recovery_species.Rdata")
load("data/processed/group_species.Rdata")


# Core species - Prevalence of >85% ------------------------------------------------------


# 1: 12 v 15 --------------------------------------------------------------

# Join with metadata to get Condition
species_long <- refeed_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  left_join(refeed_data %>% select(Seq_ID, Age_months), by = "Seq_ID")

# Calculate prevalence by Condition group
core_species_by_age <- species_long %>%
  group_by(Age_months, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.85)

# Define the core species across age groups
core_species_age <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Escherichia_coli",
  "s__Streptococcus_salivarius",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Prevotella_copri"
)

# Recalculate prevalence across age groups
prevalence_by_age <- species_long %>%
  group_by(Age_months, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter only those core species
prevalence_filtered_age <- prevalence_by_age %>%
  filter(Taxa %in% core_species_age) %>%
  arrange(Taxa, Age_months)


# Clean species names
core_species_by_age <- prevalence_filtered_age %>%
  mutate(Taxa = gsub("^s__", "", Taxa),        # remove s__ prefix
         Taxa = gsub("_", " ", Taxa))          # replace _ with space

# View core species per group
core_species_by_age %>% group_by(Age_months) %>% summarise(n_core = n_distinct(Taxa))



# 2: grA v grB ---------------------------------------------------------------

# Pivot longer and join with metadata to get Group info
species_long_group <- group_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  left_join(group_data %>% select(Seq_ID, Group), by = "Seq_ID")

# Calculate prevalence by Group
core_species_by_group <- species_long_group %>%
  group_by(Group, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.85)

# Define the core species across groups
core_species_group <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Escherichia_coli",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Collinsella_aerofaciens",
  "s__Faecalibacterium_prausnitzii",
  "s__Lactobacillus_ruminis",
  "s__Prevotella_copri"
)

# Recalculate prevalence across groups
prevalence_by_group <- species_long_group %>%
  group_by(Group, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter only those core species
prevalence_filtered_group <- prevalence_by_group %>%
  filter(Taxa %in% core_species_group) %>%
  arrange(Taxa, Group)


# Clean species names
core_species_by_group <- prevalence_filtered_group %>%
  mutate(Taxa = gsub("^s__", "", Taxa),
         Taxa = gsub("_", " ", Taxa))

# View core species per group
core_species_by_group %>%
  group_by(Group) %>%
  summarise(n_core = n_distinct(Taxa))



# 3: R v nR ---------------------------------------------------------------

# Pivot longer and join with metadata to get Recovery info
species_long_recovery <- recovery_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  left_join(recovery_data %>% select(Seq_ID, Recovery), by = "Seq_ID")

# Calculate prevalence by Recovery status
core_species_by_recovery <- species_long_recovery %>%
  group_by(Recovery, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.85)

# Define core species across recovery status
core_species_recovery <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Escherichia_coli",
  "s__Faecalibacterium_prausnitzii",
  "s__Prevotella_copri"
)

# Recalculate prevalence across recovery status
prevalence_by_recovery <- species_long_recovery %>%
  group_by(Recovery, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter for defined core species
prevalence_filtered_recovery <- prevalence_by_recovery %>%
  filter(Taxa %in% core_species_recovery) %>%
  arrange(Taxa, Recovery)


# Clean species names
core_species_by_recovery <- prevalence_filtered_recovery %>%
  mutate(Taxa = gsub("^s__", "", Taxa),
         Taxa = gsub("_", " ", Taxa))

# View core species per recovery status
core_species_by_recovery %>%
  group_by(Recovery) %>%
  summarise(n_core = n_distinct(Taxa))



# Plot --------------------------------------------------------------------


# 1: 12 v 15 --------------------------------------------------------------

core_species_by_age_plot <- ggplot(core_species_by_age, aes(x = prevalence, y = reorder(Taxa, prevalence), fill = Age_months)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  labs(x = "Prevalence", y = "") +
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
ggsave(file = "plots/chp2/clean_plots/core_species_by_age_plot.pdf", 
       plot = core_species_by_age_plot, units = "cm", width=10, height=9)


# 2: grA v grB ------------------------------------------------------------
core_species_by_group_plot <- ggplot(core_species_by_group, aes(x = prevalence, y = reorder(Taxa, prevalence), fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
  labs(x = "Prevalence", y = "") +
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
ggsave(file = "plots/chp2/clean_plots/core_species_by_group_plot.pdf", 
       plot = core_species_by_group_plot, units = "cm", width=10, height=9)

# 3: R v nR ---------------------------------------------------------------

core_species_by_recovery_plot <- ggplot(core_species_by_recovery, aes(x = prevalence, y = reorder(Taxa, prevalence), fill = Recovery)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  labs(x = "Prevalence", y = "") +
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
ggsave(file = "plots/chp2/clean_plots/core_species_by_recovery_plot.pdf", 
       plot = core_species_by_recovery_plot, units = "cm", width=10, height=9)






# Test bay  ---------------------------------------------------------------

# 12 v 15
# Define the core species across age groups
core_species_age <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Escherichia_coli",
  "s__Streptococcus_salivarius",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Prevotella_copri"
)

# Recalculate prevalence across age groups
prevalence_by_age <- species_long %>%
  group_by(Age_months, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter only those core species
prevalence_filtered_age <- prevalence_by_age %>%
  filter(Taxa %in% core_species_age) %>%
  arrange(Taxa, Age_months)

# View result
prevalence_filtered_age


# grA v grB
# Define the core species across groups
core_species_group <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Escherichia_coli",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Collinsella_aerofaciens",
  "s__Faecalibacterium_prausnitzii",
  "s__Lactobacillus_ruminis",
  "s__Prevotella_copri"
)

# Recalculate prevalence across groups
prevalence_by_group <- species_long_group %>%
  group_by(Group, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter only those core species
prevalence_filtered_group <- prevalence_by_group %>%
  filter(Taxa %in% core_species_group) %>%
  arrange(Taxa, Group)

# View result
prevalence_filtered_group

# R v nR
# Define core species across recovery status
core_species_recovery <- c(
  "s__Bifidobacterium_bifidum",
  "s__Bifidobacterium_breve",
  "s__Bifidobacterium_longum",
  "s__Bifidobacterium_pseudocatenulatum",
  "s__Escherichia_coli",
  "s__Faecalibacterium_prausnitzii",
  "s__Prevotella_copri"
)

# Recalculate prevalence across recovery status
prevalence_by_recovery <- species_long_recovery %>%
  group_by(Recovery, Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Filter for defined core species
prevalence_filtered_recovery <- prevalence_by_recovery %>%
  filter(Taxa %in% core_species_recovery) %>%
  arrange(Taxa, Recovery)

# View result
prevalence_filtered_recovery




