# This script compares the phyla of each test group 
# Test groups included in the script:

# - 12mo v 15mo 
# - GrA v GrB
# - R v nR 


# Load Packages  ----------------------------------------------------------
library(tidyverse) # required for data wrangling and ggplot
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(dplyr)
library(patchwork)

# There are sections for each analysis and each dataset required for it is added 
# under the Load Data section for each analysis

# load data ---------------------------------------------------------------
# refeed_data: contains paired data for 12mo and 15mo
load("data/processed/refeed_data.Rdata")
load("data/processed/group_data.Rdata")
load("data/processed/recovery_data.Rdata")

# Taxa data:
# phylum 
# the following datasets only contain phyla with a prevalence > 10% 
load("data/processed/filtered_refeed_phylum.Rdata") 
load("data/processed/filtered_group_phylum.Rdata")
load("data/processed/filtered_recovery_phylum.Rdata") 

# taxa with RA for each level 
load("data/processed/species_all_samples.Rdata")
load("data/processed/genus_all_samples.Rdata")
load( "data/processed/family_all_samples.Rdata")
load("data/processed/phylum_all_samples.Rdata")
load("data/processed/kingdom_all_samples.Rdata")
load("data/processed/class_all_samples.Rdata")
load("data/processed/order_all_samples.Rdata")


# Phylum Profiles ---------------------------------------------------------

# PLOTS -------------------------------------------------------------------

# 1: 12 v 15  -------------------------------------------------------

# stucture the data for the plot 
p_refeed <- filtered_refeed_phylum %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(refeed_data) 

# Violin plot with mean and log-scaled y-axis - not used in thesis 
ggplot(p_refeed, aes(x = Taxa, y = RA, fill = Age_months)) +
  geom_quasirandom(dodge.width = 0.75, shape = 21) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  labs(
    x = "Phyla",
    y = "Relative Abundance (log10 scale)"
  ) +
  scale_y_log10() +  # Add log scale to the y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("plots/chp2/Phyla-refeed.jpg", h = 6, w = 8)
ggsave("plots/chp2/Phyla-refeed.pdf", h = 6, w = 8)



# 2: grA v grB -----------------------------------------------------------
# there are no significant differences between Gr A v GrB

# structure data for the plot
p_group <- filtered_refeed_phylum %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(group_data) 

# plot not used in thesis
ggplot(p_group, aes(x = Taxa, y = RA, fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(aes(group = Group), position = position_dodge(width = 0.75), size = 0.001) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 0.001,
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
  labs(
    title = "Comparison of Phylum Relative Abundance between Refeed 1 and Refeed 2",
    x = "Phyla",
    y = "Relative Abundance (log10 scale)"
  ) +
  scale_y_log10() +  # Add log scale to the y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 



# 3: R v nR ------------------------------------------------------------------

# structure the data for the plot
p_recovery <- filtered_recovery_phylum %>%
  gather(Taxa, RA, -Seq_ID) %>%
  mutate(Taxa = str_remove(Taxa, "^p__")) %>%
  group_by(Taxa) %>%
  inner_join(recovery_data, by = "Seq_ID")


phylum_plot_recovery <- ggplot(p_recovery, aes(x = Taxa, y = RA, fill = Recovery)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Recovery)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  scale_y_log10() +
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_color_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
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
ggsave(file = "plots/chp2/clean_plots/phylum_plot_recovery.pdf", 
       plot = phylum_plot_recovery, units = "cm", width=9, height=9)


# STATS -------------------------------------------------------------------

# 1: 12 v 15  ----------------------------------------------------------------

## make the dataset match up in timepoints
p_refeed <- filtered_refeed_phylum %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(refeed_data) 

p_refeed <- p_refeed %>%
  group_by(Subject_ID) %>%
  filter(all(c("12", "15") %in% Age_months) && n_distinct(Age_months) == 2) %>% 
  ungroup()

# p value adjusted code 
pvals_refeed <- p_refeed %>%
  group_by(Taxa) %>%
  summarise(
    wilcox_test = list(wilcox.test(RA ~ Age_months, paired = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(
    stat = map_dbl(wilcox_test, ~ .x$statistic),
    p_value = map_dbl(wilcox_test, ~ .x$p.value),
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  select(Taxa, stat, p_value, p_adj)




# 2: grA v grB ---------------------------------------------------------------

pvals_group<- p_group %>%
  #filter(Taxa %in% c("p__Bacteroidetes", "p__Firmicutes", "p__Actinobacteria", "p__Proteobacteria", "p__Spirochaetes", "p__Verrucomicrobia", "p__Viruses_unclassified")) %>%
  group_by(Taxa) %>%
  summarise(p_value = wilcox.test(RA ~ Group)$p.value)

# p val adjusted code 
pvals_group <- p_group %>%
  group_by(Taxa) %>%
  summarise(
    wilcox_test = list(wilcox.test(RA ~ Group)),
    .groups = "drop"
  ) %>%
  mutate(
    stat = map_dbl(wilcox_test, ~ .x$statistic),
    p_value = map_dbl(wilcox_test, ~ .x$p.value),
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  select(Taxa, stat, p_value, p_adj)


# 3: R v nR ------------------------------------------------------------------

pavls_recovery <- p_recovery %>%
  #filter(Taxa %in% c("p__Bacteroidetes", "p__Firmicutes", "p__Actinobacteria", "p__Proteobacteria", "p__Spirochaetes", "p__Verrucomicrobia", "p__Viruses_unclassified")) %>%
  group_by(Taxa) %>%
  summarise(p_value = wilcox.test(RA ~ Recovery)$p.value)


# p val adjusted code 
pvals_recovery <- p_recovery %>%
  group_by(Taxa) %>%
  summarise(
    wilcox_test = list(wilcox.test(RA ~ Recovery)),
    .groups = "drop"
  ) %>%
  mutate(
    stat = map_dbl(wilcox_test, ~ .x$statistic),
    p_value = map_dbl(wilcox_test, ~ .x$p.value),
    p_adj = p.adjust(p_value, method = "fdr")
  ) %>%
  select(Taxa, stat, p_value, p_adj)


### Calculate the mean RA value for the phylum Proteobacteria as this is the only 
# phylum that has a sig diff 

mean_ra <- p_recovery %>%
  filter(Taxa %in% c("p__Bacteroidetes",
                     "p__Proteobacteria",
                     "p__Firmicutes",
                     "p__Viruses_unclassified")) %>%
  group_by(Taxa, Recovery) %>%
  summarise(
    Mean_RA = mean(RA, na.rm = TRUE),
    SD_RA = sd(RA, na.rm = TRUE),
    .groups = "drop"
  )

# View the results
print(mean_ra)



# F/B ratio  --------------------------------------------------------------


# 1: 12 v 15 --------------------------------------------------------------
phylum <- phylum_all_samples

FBratio_refeed <- phylum %>%
  select(p__Firmicutes, p__Bacteroidetes) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    p__Firmicutes = p__Firmicutes + 1e-6,
    p__Bacteroidetes = p__Bacteroidetes + 1e-6,
    FBratio = p__Firmicutes / p__Bacteroidetes
  ) %>%
  inner_join(refeed_data, by = "Seq_ID")


# plot 
FB_plot_refeed <- ggplot(FBratio_refeed, aes(x = Age_months, y = FBratio, fill = Age_months)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Age_months)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  scale_color_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
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

# Save to PDF
ggsave(file = "plots/chp2/clean_plots/FB_plot_refeed.pdf", 
       plot = FB_plot_refeed, units = "cm", width = 6, height = 9)

# stat
# samples need to be paired to run paired wilcoxon
# Filter to keep only subjects with both timepoints
FBratio_paired <- FBratio_refeed %>%
  filter(Age_months %in% c("12", "15")) %>%
  group_by(Subject_ID) %>%
  filter(n() == 2) %>%
  ungroup()

# Reshape to wide format for paired test
FBratio_wide <- FBratio_paired %>%
  select(Subject_ID, Age_months, FBratio) %>%
  pivot_wider(names_from = Age_months, values_from = FBratio)

# Run paired Wilcoxon test
wilcox_result_refeed <- wilcox.test(FBratio_wide$`12`, FBratio_wide$`15`, paired = TRUE)
wilcox_result_refeed
# p-value = 0.01827

# Group-level mean-based F/B ratio (phylum-level)
FB_summary_refeed <- FBratio_refeed %>%
  group_by(Age_months) %>%
  summarise(
    mean_firmicutes = mean(p__Firmicutes),
    mean_bacteroidetes = mean(p__Bacteroidetes),
    F_B_ratio = mean_firmicutes / mean_bacteroidetes,
    .groups = "drop"
  )



# 2: grA v grB ------------------------------------------------------------

# F/B Ratio calculation
FBratio_group <- phylum_all_samples %>%
  select(p__Firmicutes, p__Bacteroidetes) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    p__Firmicutes = p__Firmicutes + 1e-6,
    p__Bacteroidetes = p__Bacteroidetes + 1e-6,
    FBratio = p__Firmicutes / p__Bacteroidetes
  ) %>%
  inner_join(group_data, by = "Seq_ID")

# Plot
FB_plot_group <- ggplot(FBratio_group, aes(x = Group, y = FBratio, fill = Group)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Group)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
  scale_color_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
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

ggsave(file = "plots/chp2/clean_plots/FB_plot_group.pdf",
       plot = FB_plot_group, units = "cm", width = 6, height = 9)

# Wilcoxon test
wilcox_result_group <- wilcox.test(FBratio ~ Group, data = FBratio_group)
# p-value = 0.8026

# Summary
FB_summary_group <- FBratio_group %>%
  group_by(Group) %>%
  summarise(
    mean_firmicutes = mean(p__Firmicutes),
    mean_bacteroidetes = mean(p__Bacteroidetes),
    F_B_ratio = mean_firmicutes / mean_bacteroidetes,
    .groups = "drop"
  )





# 3: R v nR ---------------------------------------------------------------

# F/B Ratio calculation
FBratio_recovery <- phylum_all_samples %>%
  select(p__Firmicutes, p__Bacteroidetes) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    p__Firmicutes = p__Firmicutes + 1e-6,
    p__Bacteroidetes = p__Bacteroidetes + 1e-6,
    FBratio = p__Firmicutes / p__Bacteroidetes
  ) %>%
  inner_join(recovery_data, by = "Seq_ID")

# Plot
FB_plot_recovery <- ggplot(FBratio_recovery, aes(x = Recovery, y = FBratio, fill = Recovery)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 16,
    size = 1,
    alpha = 0.7,
    aes(color = Recovery)
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(width = 0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_color_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
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

ggsave(file = "plots/chp2/clean_plots/FB_plot_recovery.pdf",
       plot = FB_plot_recovery, units = "cm", width = 6, height = 9)

# Wilcoxon test
wilcox_result_recovery <- wilcox.test(FBratio ~ Recovery, data = FBratio_recovery)
# p-value = 0.02378

# Summary
FB_summary_recovery <- FBratio_recovery %>%
  group_by(Recovery) %>%
  summarise(
    mean_firmicutes = mean(p__Firmicutes),
    mean_bacteroidetes = mean(p__Bacteroidetes),
    F_B_ratio = mean_firmicutes / mean_bacteroidetes,
    .groups = "drop"
  )




# Combined FB plot  -------------------------------------------------------

# Combine plots horizontally using patchwork
combined_FB_plot_horizontal <- FB_plot_refeed | FB_plot_group | FB_plot_recovery +
  plot_layout(nrow = 1) +
  plot_annotation(title = "Firmicutes-to-Bacteroidetes (F/B) Ratio Across Conditions")

# Save to PDF
ggsave("plots/chp2/clean_plots/combined_FB_plot_horizontal.pdf",
       plot = combined_FB_plot_horizontal, width = 12, height = 9, units = "cm")


# P/B ratio  --------------------------------------------------------------

genus <- genus_all_samples


# 1: 12 v 15 --------------------------------------------------------------

PBratio_refeed <- genus_all_samples %>%
  select(g__Prevotella, g__Bacteroides) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    g__Prevotella = g__Prevotella + 1e-6,
    g__Bacteroides = g__Bacteroides + 1e-6,
    PBratio = g__Prevotella / g__Bacteroides
  ) %>%
  inner_join(refeed_data, by = "Seq_ID")

PB_plot_refeed <- ggplot(PBratio_refeed, aes(x = Age_months, y = PBratio, fill = Age_months)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Age_months)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  scale_color_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.text = element_text(), axis.line = element_line(),
        axis.title = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.position = "none", axis.text.x = element_text()) +
  xlab("") + ylab("P/B Ratio (log10-transformed)")

ggsave("plots/chp2/clean_plots/PB_plot_refeed.pdf", plot = PB_plot_refeed, units = "cm", width = 6, height = 9)

# stats
# Paired Wilcoxon for Age_months
PB_paired <- PBratio_refeed %>%
  filter(Age_months %in% c("12", "15")) %>%
  group_by(Subject_ID) %>%
  filter(n() == 2) %>%
  ungroup()

PB_wide <- PB_paired %>%
  select(Subject_ID, Age_months, PBratio) %>%
  pivot_wider(names_from = Age_months, values_from = PBratio)

wilcox_result_refeed_pb <- wilcox.test(PB_wide$`12`, PB_wide$`15`, paired = TRUE)
# p-value = 0.1083

# Summary
PB_summary_refeed <- PBratio_refeed %>%
  group_by(Age_months) %>%
  summarise(
    mean_prevotella = mean(g__Prevotella - 1e-6),
    mean_bacteroides = mean(g__Bacteroides - 1e-6),
    P_B_ratio = mean_prevotella / mean_bacteroides,
    .groups = "drop"
  )

# 2: grA v grB ------------------------------------------------------------

PBratio_group <- genus_all_samples %>%
  select(g__Prevotella, g__Bacteroides) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    g__Prevotella = g__Prevotella + 1e-6,
    g__Bacteroides = g__Bacteroides + 1e-6,
    PBratio = g__Prevotella / g__Bacteroides
  ) %>%
  inner_join(group_data, by = "Seq_ID")

PB_plot_group <- ggplot(PBratio_group, aes(x = Group, y = PBratio, fill = Group)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Group)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
  scale_color_manual(values = c("Local RUSF (A)" = "azure3", "ERUSF (B)" = "black")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.text = element_text(), axis.line = element_line(),
        axis.title = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.position = "none", axis.text.x = element_text()) +
  xlab("") + ylab("P/B Ratio (log10-transformed)")

ggsave("plots/chp2/clean_plots/PB_plot_group.pdf", plot = PB_plot_group, units = "cm", width = 6, height = 9)

# Stat
# Wilcoxon for Group
wilcox_result_group_pb <- wilcox.test(PBratio_group$PBratio ~ PBratio_group$Group)
# p-value = 0.7877

# Summary
PB_summary_group <- PBratio_group %>%
  group_by(Group) %>%
  summarise(
    mean_prevotella = mean(g__Prevotella - 1e-6),
    mean_bacteroides = mean(g__Bacteroides - 1e-6),
    P_B_ratio = mean_prevotella / mean_bacteroides,
    .groups = "drop"
  )




# R v nR ------------------------------------------------------------------

PBratio_recovery <- genus_all_samples %>%
  select(g__Prevotella, g__Bacteroides) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    g__Prevotella = g__Prevotella + 1e-6,
    g__Bacteroides = g__Bacteroides + 1e-6,
    PBratio = g__Prevotella / g__Bacteroides
  ) %>%
  inner_join(recovery_data, by = "Seq_ID")

PB_plot_recovery <- ggplot(PBratio_recovery, aes(x = Recovery, y = PBratio, fill = Recovery)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Recovery)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(limits = c(1e-6, 200), breaks = c(1e-6, 1e-4, 1e-2, 1, 100)) +
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_color_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid = element_blank(),
        axis.text = element_text(), axis.line = element_line(),
        axis.title = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.position = "none", axis.text.x = element_text()) +
  xlab("") + ylab("P/B Ratio (log10-transformed)")

ggsave("plots/chp2/clean_plots/PB_plot_recovery.pdf", plot = PB_plot_recovery, units = "cm", width = 6, height = 9)

# stat
# Wilcoxon for Recovery
wilcox_result_recovery_pb <- wilcox.test(PBratio_recovery$PBratio ~ PBratio_recovery$Recovery)
# p-value = 0.4724

# Summary
PB_summary_recovery <- PBratio_recovery %>%
  group_by(Recovery) %>%
  summarise(
    mean_prevotella = mean(g__Prevotella - 1e-6),
    mean_bacteroides = mean(g__Bacteroides - 1e-6),
    P_B_ratio = mean_prevotella / mean_bacteroides,
    .groups = "drop"
  )


