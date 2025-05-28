## code summary  
## Data can be directly loaded using the files in the "save datasets" section  
## This script performs baseline gut microbiome analyses (12-month timepoint)  
## comparing MAM vs Healthy infants using MetaPhlAn3 taxonomic outputs.

## Step 1: Load required packages and all taxonomic/metadata datasets  
## Step 2: Calculate and compare the proportion of UNKNOWN taxa between groups  
## Step 3: Prepare relative abundance (RA) tables at all taxonomic levels for baseline samples  
## Step 4: Plot overall RA distribution by taxonomic level and condition  
## Step 5: Run Wilcoxon tests for RA differences across taxonomic levels  
## Step 6: Count and plot number of unique taxa per level and condition  
## Step 7: Filter phylum, genus, and species tables for prevalence >10%  
## Step 8: Plot phylum-level RA distributions and run Wilcoxon tests on filtered phyla  
## Step 9: Calculate, plot, and test Firmicutes/Bacteroidetes (F/B) ratio  
## Step 10: Calculate, plot, and test Prevotella/Bacteroides (P/B) ratio  
## Step 11: Investigate relationship between Firmicutes abundance and WLZ/WHZ scores using two modeling frameworks:  
##          a) Firmicutes as predictor of WLZ/WHZ  
##          b) WLZ/WHZ as predictor of Firmicutes abundance  
## Step 12: Plot both modeling approaches, stratified by Condition (MAM vs Healthy)  
## Step 13: Save all processed and filtered datasets for downstream analysis  




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




# load data ---------------------------------------------------------------

# all samples - all taxa
load("data/processed/metaphlan3_taxonomy.Rdata") # loads cleaned file called "taxonomy"

# taxa with RA for each level 
load("data/processed/species_all_samples.Rdata")
load("data/processed/genus_all_samples.Rdata")
load( "data/processed/family_all_samples.Rdata")
load("data/processed/phylum_all_samples.Rdata")
load("data/processed/kingdom_all_samples.Rdata")
load("data/processed/class_all_samples.Rdata")
load("data/processed/order_all_samples.Rdata")

# baseline metadata 
load("data/processed/baseline_data.Rdata")
# baseline meta with anthropometrics 
load("data/processed/all_meta_baseline.Rdata")




# UNKNOWN taxa proportions ------------------------------------------------

# There is a large prportion of UNKNOWN taxa in each sample, 
# we want to. know if this differs between MAM and Healthy

# Filter the taxonomy table for baseline data 
taxonomy_baseline <- baseline_data %>%
  left_join(taxonomy, by = "Seq_ID")


# check Seq_IDs without Uknown
seq_ids_without_unknown <- taxonomy_baseline %>%
  group_by(Seq_ID) %>%
  summarise(has_unknown = any(Taxa == "UNKNOWN")) %>%
  filter(!has_unknown) %>%
  pull(Seq_ID)


check <- taxonomy_baseline %>%
  filter(Taxa == "UNKNOWN") 


# Categorize Taxa as Labeled vs Unknown
taxonomy_baseline <- taxonomy_baseline %>%
  mutate(Taxa_Label = case_when(
    Taxa == "UNKNOWN" ~ "Unknown",
    TRUE ~ "Labeled"
  ))



# Calculate the correct proportion of UNKNOWN based on RA
unknown_RA_per_sample <- taxonomy_baseline %>%
  group_by(Seq_ID) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA
  )


# Now calculate the mean and standard deviation of UNKNOWN RA across samples
unknown_RA_summary <- unknown_RA_per_sample %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )


# by Condition
# First, calculate UNKNOWN RA per sample and keep Condition
unknown_RA_per_sample_condition <- taxonomy_baseline %>%
  group_by(Seq_ID, Condition) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA,
    .groups = "drop"
  )

# Now summarise by Condition
unknown_RA_summary_by_condition <- unknown_RA_per_sample_condition %>%
  group_by(Condition) %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# Wilcoxon rank-sum test comparing prop_unknown between MAM and Healthy
wilcox_unknown_RA <- wilcox.test(
  prop_unknown ~ Condition,
  data = unknown_RA_per_sample_condition
)









# Taxa distribution - MAM v Healthy ---------------------------------------
# Q: What is the distribution of taxa between MAM and Healthy at 12 months?
# plot - bar plot 

# modify taxa datasets and filter for baseline data 
# we turn these into objects to be used downstream

species_all_samples_edit <- species_all_samples %>% 
  rownames_to_column("Seq_ID")

genus_all_samples_edit <- genus_all_samples %>% 
  rownames_to_column("Seq_ID")

family_all_samples_edit <- family_all_samples %>% 
  rownames_to_column("Seq_ID")

phylum_all_samples_edit <- phylum_all_samples %>% 
  rownames_to_column("Seq_ID")

kingdom_all_samples_edit <- kingdom_all_samples %>% 
  rownames_to_column("Seq_ID")

order_all_samples_edit <- order_all_samples %>% 
  rownames_to_column("Seq_ID")

class_all_samples_edit <- class_all_samples %>% 
  rownames_to_column("Seq_ID")


# Filter to only baseline samples 

baseline_species <- species_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_genus <- genus_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_family <- family_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_phylum <- phylum_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_kingdom <- kingdom_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_order <- order_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

baseline_class <- class_all_samples_edit %>%
  semi_join(baseline_data, by = "Seq_ID") 

# Convert to table format 
baseline_species_table <- baseline_species %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_genus_table <- baseline_genus %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_family_table <- baseline_family %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_phylum_table <- baseline_phylum %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_kingdom_table <- baseline_kingdom %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_order_table <- baseline_order %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")

baseline_class_table <- baseline_class %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data, by = "Seq_ID")


# Combine baseline *table datasets with a Level column
combined_taxa_baseline <- bind_rows(
  baseline_kingdom_table %>% mutate(Level = "Kingdom"),
  baseline_phylum_table %>% mutate(Level = "Phylum"),
  baseline_class_table %>% mutate(Level = "Class"),
  baseline_order_table %>% mutate(Level = "Order"),
  baseline_family_table %>% mutate(Level = "Family"),
  baseline_genus_table %>% mutate(Level = "Genus"),
  baseline_species_table %>% mutate(Level = "Species")
)


# remove all RA = 0 rows
combined_taxa_baseline <- combined_taxa_baseline %>%
  filter(RA != 0)


# Plot RA by Taxonomic Level and Condition
all_taxa_plot <- combined_taxa_baseline %>%
  ggplot(aes(x = Level, y = RA, fill = Condition)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Condition)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  scale_y_log10() +
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
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Relative Abundance (log10-transformed)") 


# Run Wilcoxon test for each taxonomic level
wilcoxon_results <- combined_taxa_baseline %>%
  group_by(Level) %>%
  summarise(
    p_value = wilcox.test(RA ~ Condition)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adjusted = p.adjust(p_value, method = "fdr")  # Adjust for multiple comparisons
  )

# View results
print(wilcoxon_results)



# plot - unique taxa ------------------------------------------------------

# Count unique taxa per Level and Condition
unique_taxa_counts <- combined_taxa_baseline %>%
  group_by(Level, Condition) %>%
  summarise(unique_taxa_count = n_distinct(Taxa), .groups = "drop")

# Set factor levels to enforce taxonomic order
unique_taxa_counts$Level <- factor(
  unique_taxa_counts$Level,
  levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
)


# Plot the bar plot with ordered x-axis
unique_taxa_plot <- ggplot(unique_taxa_counts, aes(x = Level, y = unique_taxa_count, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  labs(
    x = "",
    y = "Number of Unique Taxa"
  ) +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


# pdf
ggsave(file = "plots/chp1/clean_plots/unique_taxa_plot.pdf", 
       plot = unique_taxa_plot, units = "cm", width=12, height=12)



# Prevalence filtering  -------------------------------------------------------------
# this section does the following:

# filters the phylum dataset for prev >10% only 
# filters the species dataset for prev > 10% only 
# filteres the genus dataset for prev >10%
# !!! these are the datasets used in downstream analyses 


# phylum ------------------------------------------------------------------

# Calculate prevalence of each phylum (proportion of samples with RA > 0)
baseline_phylum_prevalence <- baseline_phylum %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# prevalence of each phylum present and RA
phylum_prevalence_detailed <- baseline_phylum %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  mutate(
    Count_present = sum(RA > 0),
    Prevalence = Count_present / n()
  ) %>%
  ungroup()


# Filter for 10%
filtered_baseline_phylum <- baseline_phylum %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.1) %>%
  pull(Taxa) %>%
  { baseline_phylum[, c("Seq_ID", .)] }




# genus -------------------------------------------------------------------

# Calculate prevalence of each genus
baseline_genus_prevalence <- baseline_genus %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Prevalence of each genus present and RA
genus_prevalence_detailed <- baseline_genus %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  mutate(
    Count_present = sum(RA > 0),
    Prevalence = Count_present / n()
  ) %>%
  ungroup()

# Filter for 10%
filtered_baseline_genus <- baseline_genus %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.1) %>%
  pull(Taxa) %>%
  { baseline_genus[, c("Seq_ID", .)] }



# species -----------------------------------------------------------------
# Calculate prevalence of each species
baseline_species_prevalence <- baseline_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop")

# Prevalence of each species present and RA
species_prevalence_detailed <- baseline_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  mutate(
    Count_present = sum(RA > 0),
    Prevalence = Count_present / n()
  ) %>%
  ungroup()

# Filter for 10%
filtered_baseline_species <- baseline_species %>%
  pivot_longer(-Seq_ID, names_to = "Taxa", values_to = "RA") %>%
  group_by(Taxa) %>%
  summarise(prevalence = mean(RA > 0), .groups = "drop") %>%
  filter(prevalence >= 0.1) %>%
  pull(Taxa) %>%
  { baseline_species[, c("Seq_ID", .)] }
















# Phylum profiles  --------------------------------------------------------

# convert the data table to a different format - use the prevalence filter added phylum data
p_baseline <- filtered_baseline_phylum %>% 
  gather(Taxa, RA, -Seq_ID) %>% 
  mutate(Taxa = gsub("^p__", "", Taxa)) %>%  # remove 'p__' prefix
  inner_join(baseline_data, by = "Seq_ID")

# plot 
phylum_plot <- ggplot(p_baseline, aes(x = Taxa, y = RA, fill = Condition)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Condition)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  scale_y_log10() +
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
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  xlab("") +
  ylab("Relative Abundance (log10-transformed)") 

# pdf
ggsave(file = "plots/chp1/clean_plots/phylum_plot.pdf", 
       plot = phylum_plot, units = "cm", width=9, height=9)


# convert the data table to a different format - use the prevalence filter added phylum data
p_baseline <- filtered_baseline_phylum %>% 
  gather(Taxa, RA, -Seq_ID) %>% 
  mutate(Taxa = gsub("^p__", "", Taxa)) %>%  # remove 'p__' prefix
  inner_join(baseline_data, by = "Seq_ID")


# Wilcoxon 

# Run Wilcoxon rank-sum test for each phylum and adjust p-values
wilcox_phylum_results <- p_baseline %>%
  group_by(Taxa) %>%
  summarise(
    wilcox_test = list(wilcox.test(RA ~ Condition)),
    .groups = "drop"
  ) %>%
  mutate(
    p_value = map_dbl(wilcox_test, ~ .x$p.value),
    stat = map_dbl(wilcox_test, ~ .x$statistic),
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  select(Taxa, stat, p_value, p_adj)

# no sig difference in phyla 



# Phylum bar plot ---------------------------------------------------------

# Bar plot of Phyla at a glance for MAM v Healthy after a 10% prevalence filter 
bar_phylum <- filtered_baseline_phylum %>% 
  # rownames_to_column("Seq_ID") %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(baseline_data) 

# Calculate relative abundance (RA) as percentages within each Condition
bar_phylum <- within(bar_phylum, {
  RA_percent <- RA / sum(RA) * 100
})

phylum_bar_plot <- ggplot(bar_phylum, aes(x = Condition, y = RA_percent, fill = Taxa)) +
  geom_bar(stat = "identity", position = "fill", color = NA, linewidth = 0) +
  scale_fill_brewer(palette = "Dark2") +
  labs(
    x = "",
    y = "Relative Abundance (Proportion)"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 9)
  )

# pdf
ggsave(file = "plots/chp1/clean_plots/phylum_bar_plot.pdf", 
       plot = phylum_bar_plot, units = "cm", width=12, height=9)


# F/B ratio  --------------------------------------------------------------

phylum <- phylum_all_samples

FBratio_mam_healthy <- phylum %>%
  select(p__Firmicutes, p__Bacteroidetes) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    p__Firmicutes = p__Firmicutes + 1e-6,
    p__Bacteroidetes = p__Bacteroidetes + 1e-6,
    FBratio = p__Firmicutes / p__Bacteroidetes
  ) %>%
  inner_join(baseline_data)


# plot FB -----------------------------------------------------------------

FB_plot <- ggplot(FBratio_mam_healthy, aes(x = Condition, y = FBratio, fill = Condition)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Condition)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) + # add threshold line (>1 = Firmicutes dominance)
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
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
    legend.position = "none",
    axis.text.x = element_text()
  ) +
  xlab("") +
  ylab("F/B Ratio (log10-transformed)")

# pdf
ggsave(file = "plots/chp1/clean_plots/FB_plot.pdf", 
       plot = FB_plot, units = "cm", width=6, height=9)



# stat - F/B --------------------------------------------------------------------

# Wilcoxon test
wilcox_result_mam_healthy_f_b <- wilcox.test(FBratio ~ Condition, data = FBratio_mam_healthy)
# p-value = 0.02798

# Group-level mean-based F/B ratio (phylum-level)
FB_summary <- FBratio_mam_healthy %>%
  group_by(Condition) %>%
  summarise(
    mean_firmicutes = mean(p__Firmicutes),
    mean_bacteroidetes = mean(p__Bacteroidetes),
    F_B_ratio = mean_firmicutes / mean_bacteroidetes,
    .groups = "drop"
  )



# P/B ratio ---------------------------------------------------------------

genus <- genus_all_samples

PBratio_mam_healthy <- genus %>% 
  select(g__Prevotella, g__Bacteroides) %>% #change to genus level 
  rownames_to_column("Seq_ID") %>% 
  mutate(g__Prevotella = g__Prevotella+1e-6, # add small pseudocount to account for zero abundance
         g__Bacteroides = g__Bacteroides+1e-6, # add small pseudocount to account for zero abundance
         PBratio = g__Prevotella/g__Bacteroides) %>% # we are diving P by B so above 1 would indicate higher F than B. Zero would indicate equl amounts of F and B. 
  inner_join(baseline_data)


# plot - PB ---------------------------------------------------------------

PB_plot <- ggplot(PBratio_mam_healthy, aes(x = Condition, y = PBratio, fill = Condition)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Condition)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) + # add threshold line (1 = Prevotella dominance)
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
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
    legend.position = "none",
    axis.text.x = element_text()
  ) +
  xlab("") +
  ylab("P/B Ratio (log10-transformed)")

# pdf
ggsave(file = "plots/chp1/clean_plots/PB_plot.pdf", 
       plot = PB_plot, units = "cm", width=6, height=9)


# stat - PB --------------------------------------------------------------------

# Statistic check
wilcox_result_mam_healthy_p_b <- wilcox.test(PBratio_mam_healthy$PBratio ~ PBratio_mam_healthy$Condition)
# p val = 0.0147

# Calculate the total (or mean) abundance of each genus per condition
PB_summary <- PBratio_mam_healthy %>%
  group_by(Condition) %>%
  summarise(
    mean_prevotella = mean(g__Prevotella - 1e-6),  # remove pseudocount
    mean_bacteroides = mean(g__Bacteroides - 1e-6),
    P_B_ratio = mean_prevotella / mean_bacteroides
  )




# Firmicutes ~ WLZ/WHZ ----------------------------------------------------

# Firmicutes phylum has been associated with weight gain and we see higher Firmicutes in 
# healthy samples, therefore, we want to identify a correlation between the WLZ_WHZ ratio and Firmicutes abundance
# we choose WLZ_WHZ ratio as this is the Z score for weight length and weight height which standardizes the 
# spread of weight over lenght and height for each infant thus making it comparable

# Prepare dataset
firmicutes_data <- baseline_phylum %>%
  select(Seq_ID, p__Firmicutes) %>%
  inner_join(baseline_data, by = "Seq_ID") %>%
  mutate(p__Firmicutes = log10(p__Firmicutes + 1e-6)) # log-transform Firmicutes


# Q1 ----------------------------------------------------------------------

# Q: Does Firmicutes abundance explain variation in WLZ/WHZ scores, accounting for Condition (MAM vs Healthy) and Sex?”
# This treats Firmicutes as a potential predictor, and WLZ/WHZ as the response,

# Fit linear model
model_firmicutes <- lm(base_WLZ_WHZ ~ p__Firmicutes + Condition + Sex, data = firmicutes_data)

# Summary
summary(model_firmicutes)
# p = 0.0791 for firmicutes 

## Run on each Condition 

# Split data by Condition
mam_data <- firmicutes_data %>% filter(Condition == "MAM")
healthy_data <- firmicutes_data %>% filter(Condition == "Healthy")

# Linear model for MAM only
model_mam <- lm(base_WLZ_WHZ ~ p__Firmicutes + Sex, data = mam_data)
summary(model_mam)

# Linear model for Healthy only
model_healthy <- lm(base_WLZ_WHZ ~ p__Firmicutes + Sex, data = healthy_data)
summary(model_healthy)





# Q2 ----------------------------------------------------------------------

# Q2: As infants become more malnourished (lower WLZ/WHZ), does the abundance of Firmicutes decline?

# Whole dataset
model_firmicutes_response <- lm(p__Firmicutes ~ base_WLZ_WHZ + Condition + Sex, data = firmicutes_data)
summary(model_firmicutes_response)

# MAM only
model_firmicutes_mam <- lm(p__Firmicutes ~ base_WLZ_WHZ + Sex, data = mam_data)
summary(model_firmicutes_mam)

# Healthy only
model_firmicutes_healthy <- lm(p__Firmicutes ~ base_WLZ_WHZ + Sex, data = healthy_data)
summary(model_firmicutes_healthy)



# plot - Q1 & Q2 ----------------------------------------------------------

# Plot 1: WLZ/WHZ as response ~ log-transformed Firmicutes
fimicutes_Q1_plot <- ggplot(firmicutes_data, aes(x = p__Firmicutes, y = base_WLZ_WHZ, color = Condition)) +
  geom_point(alpha = 0.7) +
  #geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Condition) +
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(),
    axis.text.x = element_text(size = 7),  # x-axis text size
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 9),
    legend.position = "none"  # remove legend
  ) +
  xlab("Relative Abundance (log10-transformed)") +
  ylab("WLZ/WHZ Score") 

# pdf
ggsave(file = "plots/chp1/clean_plots/fimicutes_Q1_plot.pdf", 
       plot = fimicutes_Q1_plot, units = "cm", width=8, height=6)


# Plot 2: Firmicutes as response ~ WLZ/WHZ
fimicutes_Q2_plot <- ggplot(firmicutes_data, aes(x = base_WLZ_WHZ, y = p__Firmicutes, color = Condition)) +
  geom_point(alpha = 0.7) +
  #geom_smooth(method = "lm", se = TRUE, color = "black") +
  facet_wrap(~ Condition) +
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(),
    axis.text.x = element_text(size = 7),  # x-axis text size
    axis.text.y = element_text(size = 9),  # y-axis text size
    axis.title = element_text(size = 9),
    strip.text = element_text(size = 9),
    legend.position = "none"  # remove legend
  ) +
  xlab("WLZ/WHZ Score") +
  ylab("Relative Abundance (log10-transformed)") 


# pdf
ggsave(file = "plots/chp1/clean_plots/fimicutes_Q2_plot.pdf", 
       plot = fimicutes_Q2_plot, units = "cm", width=8, height=6)







# save datasets  ----------------------------------------------------------

# Save each baseline dataset
save(baseline_species, file = "data/processed/baseline_species.Rdata")
save(baseline_genus, file = "data/processed/baseline_genus.Rdata")
save(baseline_family, file = "data/processed/baseline_family.Rdata")
save(baseline_phylum, file = "data/processed/baseline_phylum.Rdata")
save(baseline_kingdom, file = "data/processed/baseline_kingdom.Rdata")
save(baseline_order, file = "data/processed/baseline_order.Rdata")
save(baseline_class, file = "data/processed/baseline_class.Rdata")

# Prevalence filtered
save(filtered_baseline_phylum, file = "data/processed/filtered_baseline_phylum.Rdata") 
save(filtered_baseline_genus, file = "data/processed/filtered_baseline_genus.Rdata")
save(filtered_baseline_species, file = "data/processed/filtered_baseline_species.Rdata")


