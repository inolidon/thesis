




# load data ---------------------------------------------------------------
# refeed_data: contains paired data for 12mo and 15mo
load("data/processed/refeed_data.Rdata")
load("data/processed/group_data.Rdata")
load("data/processed/recovery_data.Rdata")

# Taxa data:
# species 
load("data/processed/refeed_species.Rdata")
load("data/processed/group_species.Rdata")
load("data/processed/recovery_species.Rdata")



# 1: 12 v 15 -----------------------------------------------------------------

# refeed_species must not have a column for the identifier, convert it to a row
refeed_species <- refeed_species %>% 
  column_to_rownames("Seq_ID")

# Refeed Alpha Diversity
alpha_div_refeed <- tibble(
  Seq_ID = rownames(refeed_species),
  Richness = vegan::specnumber(refeed_species),
  Shannon = vegan::diversity(refeed_species, index = "shannon"),
  Evenness = Shannon / log(Richness)
) %>%
  inner_join(refeed_data, by = "Seq_ID")

# Summary
alpha_div_refeed_summary <- alpha_div_refeed %>%
  group_by(Age_months) %>%
  summarise(
    mean_shannon = mean(Shannon), sd_shannon = sd(Shannon),
    mean_richness = mean(Richness), sd_richness = sd(Richness),
    mean_evenness = mean(Evenness), sd_evenness = sd(Evenness),
    .groups = "drop"
  )

# Stats (paired)
alpha_paired <- alpha_div_refeed %>%
  group_by(Subject_ID) %>%
  filter(n() == 2) %>%
  select(Subject_ID, Age_months, Shannon, Richness, Evenness) %>%
  pivot_wider(names_from = Age_months, values_from = c(Shannon, Richness, Evenness))

wilcox_shannon <- wilcox.test(alpha_paired$`Shannon_12`, alpha_paired$`Shannon_15`, paired = TRUE) # p-value = 0.1352
wilcox_richness <- wilcox.test(alpha_paired$`Richness_12`, alpha_paired$`Richness_15`, paired = TRUE) # p-value = 0.02718
wilcox_evenness <- wilcox.test(alpha_paired$`Evenness_12`, alpha_paired$`Evenness_15`, paired = TRUE) # p-value = 0.2633



# 2: grA v grB ------------------------------------------------------------

# group_species must not have a column for the identifier, convert it to a row
group_species <- group_species %>% 
  column_to_rownames("Seq_ID")

# Group Alpha Diversity
alpha_div_group <- tibble(
  Seq_ID = rownames(group_species),
  Richness = vegan::specnumber(group_species),
  Shannon = vegan::diversity(group_species, index = "shannon"),
  Evenness = Shannon / log(Richness)
) %>%
  inner_join(group_data, by = "Seq_ID")

# Summary
alpha_div_group_summary <- alpha_div_group %>%
  group_by(Group) %>%
  summarise(
    mean_shannon = mean(Shannon), sd_shannon = sd(Shannon),
    mean_richness = mean(Richness), sd_richness = sd(Richness),
    mean_evenness = mean(Evenness), sd_evenness = sd(Evenness),
    .groups = "drop"
  )

# Stats
wilcox_group_shannon <- wilcox.test(Shannon ~ Group, data = alpha_div_group) # p-value = 0.03692
wilcox_group_richness <- wilcox.test(Richness ~ Group, data = alpha_div_group) # p-value = 0.322
wilcox_group_evenness <- wilcox.test(Evenness ~ Group, data = alpha_div_group) # p-value = 0.05437



# 3: R v nR ---------------------------------------------------------------

# recovery_species must not have a column for the identifier, convert it to a row
recovery_species <- recovery_species %>% 
  column_to_rownames("Seq_ID")


# Recovery Alpha Diversity
alpha_div_recovery <- tibble(
  Seq_ID = rownames(recovery_species),
  Richness = vegan::specnumber(recovery_species),
  Shannon = vegan::diversity(recovery_species, index = "shannon"),
  Evenness = Shannon / log(Richness)
) %>%
  inner_join(recovery_data, by = "Seq_ID")

# Summary
alpha_div_recovery_summary <- alpha_div_recovery %>%
  group_by(Recovery) %>%
  summarise(
    mean_shannon = mean(Shannon), sd_shannon = sd(Shannon),
    mean_richness = mean(Richness), sd_richness = sd(Richness),
    mean_evenness = mean(Evenness), sd_evenness = sd(Evenness),
    .groups = "drop"
  )

# Stats
wilcox_recovery_shannon <- wilcox.test(Shannon ~ Recovery, data = alpha_div_recovery) # p-value = 0.5816
wilcox_recovery_richness <- wilcox.test(Richness ~ Recovery, data = alpha_div_recovery) # p-value = 0.7376
wilcox_recovery_evenness <- wilcox.test(Evenness ~ Recovery, data = alpha_div_recovery) # p-value = 0.4057




# plots  ------------------------------------------------------------------

# 1: 12 v 15 --------------------------------------------------------------
# Reshape to long format
alpha_div_long_refeed <- alpha_div_refeed %>%
  select(Seq_ID, Age_months, Shannon, Richness, Evenness) %>%
  pivot_longer(cols = c(Shannon, Richness, Evenness),
               names_to = "Alpha_metric",
               values_to = "Value")

# Set order
alpha_div_long_refeed$Alpha_metric <- factor(
  alpha_div_long_refeed$Alpha_metric,
  levels = c("Shannon", "Richness", "Evenness")
)

# Plot
alpha_div_plot_refeed <- ggplot(alpha_div_long_refeed, aes(x = Age_months, y = Value, fill = Age_months)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Age_months)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Alpha_metric, scales = "free_y") +
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
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Alpha Diversity Metric")

ggsave(file = "plots/chp2/clean_plots/alpha_div_plot_refeed.pdf",
       plot = alpha_div_plot_refeed, units = "cm", width = 10, height = 9)




# 2: grA v grB ------------------------------------------------------------

alpha_div_long_group <- alpha_div_group %>%
  select(Seq_ID, Group, Shannon, Richness, Evenness) %>%
  pivot_longer(cols = c(Shannon, Richness, Evenness),
               names_to = "Alpha_metric",
               values_to = "Value")

alpha_div_long_group$Alpha_metric <- factor(
  alpha_div_long_group$Alpha_metric,
  levels = c("Shannon", "Richness", "Evenness")
)

alpha_div_plot_group <- ggplot(alpha_div_long_group, aes(x = Group, y = Value, fill = Group)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Group)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Alpha_metric, scales = "free_y") +
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
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Alpha Diversity Metric")

ggsave(file = "plots/chp2/clean_plots/alpha_div_plot_group.pdf",
       plot = alpha_div_plot_group, units = "cm", width = 10, height = 9)





# 3: R v nR ---------------------------------------------------------------


alpha_div_long_recovery <- alpha_div_recovery %>%
  select(Seq_ID, Recovery, Shannon, Richness, Evenness) %>%
  pivot_longer(cols = c(Shannon, Richness, Evenness),
               names_to = "Alpha_metric",
               values_to = "Value")

alpha_div_long_recovery$Alpha_metric <- factor(
  alpha_div_long_recovery$Alpha_metric,
  levels = c("Shannon", "Richness", "Evenness")
)

alpha_div_plot_recovery <- ggplot(alpha_div_long_recovery, aes(x = Recovery, y = Value, fill = Recovery)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Recovery)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3,
               position = position_dodge(width = 0.75)) +
  facet_wrap(~ Alpha_metric, scales = "free_y") +
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
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Alpha Diversity Metric")

ggsave(file = "plots/chp2/clean_plots/alpha_div_plot_recovery.pdf",
       plot = alpha_div_plot_recovery, units = "cm", width = 10, height = 9)





















