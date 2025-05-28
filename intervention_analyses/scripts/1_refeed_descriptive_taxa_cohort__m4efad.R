# This script contains all the code for generating the descriptive stats for the
# 15 months time point, i.e., the refeed 



# Load packages  ----------------------------------------------------------
library(tidyverse)
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(dplyr)


# Load Data ---------------------------------------------------------------

# the data we will use is the 
# 1. metadata: 
# refeed_data: contains paired data for 12mo and 15mo
load("data/processed/refeed_data.Rdata") 

# recovery_data: contains information about refeed group and recovery status 
load("data/processed/recovery_data.Rdata")

# taxonomy table with all the taxa - for UKNOWN taxa 
load("data/processed/metaphlan3_taxonomy.Rdata") # loads cleaned file called "taxonomy"

# 2: Taxa data:
# species 
load("data/processed/refeed_species.Rdata")

# genus
load("data/processed/refeed_genus.Rdata")

# phylum 
load("data/processed/refeed_phylum.Rdata")

# class 
load("data/processed/refeed_class.Rdata")

# order 
load("data/processed/refeed_order.Rdata")

# kingdom 
load("data/processed/refeed_kingdom.Rdata")

# family
load("data/processed/refeed_family.Rdata")


# DESCRIPTIVE STATS -------------------------------------------------------


# 1: baseline (12mo) vs refeed (15mo) -------------------------------------

# the 12mo and 15mo samples have been paired to only include the participants who
# had a baseline sample and a post refeed sample
# this is because paired statistical tests require equal matched pairs, eg: Wilcoxon signed rank test

# Sample numbers in each time point
age_sample_no_table <- refeed_data %>%
  group_by(Age_months) %>%
  summarize(Count = n(), .groups = 'drop') # removes grouping to preserve the data structure
# n = 148 in each age group


# samples without pairs were dropped 
# here we identify how many samples did not have pairs

# load the main meta file 
load("data/processed/meta_arm.Rdata")

# select the refeed timepoint variables without removing any points 
all_refeed_data <- meta_arm %>% 
  filter(Arm %in% c("MAM_12", "MAM_15"))%>%
  arrange(Subject_ID)

# identify the unmatched samples, i.e., ones without a second timepoint
unpaired_subjects <- all_refeed_data %>%
  group_by(Subject_ID) %>%
  filter(n_distinct(Age_months) == 1) # 9 samples have only one time point, unmatched

# saving this as an excel for supplementary data use
library(writexl)

# Save dataset as Excel file
write_xlsx(unpaired_subjects, "outputs/tables/unpaired_subjects_chp2.xlsx")


# 2: refeed type 1 vs type 2 ----------------------------------------------

# only the 15mo MAM samples are being used here as they are the only time point 
# at which the intervention was tested 

# Sample numbers that received type A v type B refeed 

refeed_type_sample_no_table <- recovery_data %>%
  group_by(Group.x) %>%
  summarize(Count = n(), .groups = 'drop')
# Local RUSF (A) = 75
# ERUSF (B) = 73


# 3: recovered vs not recovered -------------------------------------------

# the 15mo samples are stratified into the outcome of the intervention, 
# either recovered or not recovered in 3 months. 

# Sample numbers that Recovered v not recovered

recovery_type_sample_no_table <- recovery_data %>%
  group_by(Recovery) %>%
  summarize(Count = n(), .groups = 'drop')
# FALSE = 82
# TRUE = 66


# 4: stratified groups and recovery ----------------------------------------------------

# the numbers within each category of refeed type stratified by recovery 
group_by_recovery_sample_no_table <- recovery_data %>%
  group_by(Recovery, Group.x) %>%
  summarize(Count = n(), .groups = 'drop')

# Recovery Group          Count
# FALSE    ERUSF (B)         41
# FALSE    Local RUSF (A)    41
# TRUE     ERUSF (B)         32
# TRUE     Local RUSF (A)    34


# UNKNOWN taxa proportions ------------------------------------------------

# There is a large prportion of UNKNOWN taxa in each sample, 
# we want to. know if this differs between baseline and post refeed
# Q: did refeed introduce more unknown taxa? - change in the microbial profile 


# 1: 12 v 15 --------------------------------------------------------------


# Filter the taxonomy table for baseline data 
taxonomy_refeed <- refeed_data %>%
  left_join(taxonomy, by = "Seq_ID")


# check Seq_IDs without Uknown
seq_ids_without_unknown <- taxonomy_refeed %>%
  group_by(Seq_ID) %>%
  summarise(has_unknown = any(Taxa == "UNKNOWN")) %>%
  filter(!has_unknown) %>%
  pull(Seq_ID)

# how many samples have unknown 
check <- taxonomy_refeed %>%
  filter(Taxa == "UNKNOWN") 


# Categorize Taxa as Labeled vs Unknown
taxonomy_refeed <- taxonomy_refeed %>%
  mutate(Taxa_Label = case_when(
    Taxa == "UNKNOWN" ~ "Unknown",
    TRUE ~ "Labeled"
  ))



# Calculate the correct proportion of UNKNOWN based on RA
unknown_RA_per_sample_refeed <- taxonomy_refeed %>%
  group_by(Seq_ID) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA
  )


# Now calculate the mean and standard deviation of UNKNOWN RA across samples
unknown_RA_summary_refeed <- unknown_RA_per_sample_refeed %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )


# by Age_months
# First, calculate UNKNOWN RA per sample and keep Age group
unknown_RA_per_sample_age_group <- taxonomy_refeed %>%
  group_by(Seq_ID, Age_months) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA,
    .groups = "drop"
  )

# Now summarise by Age_months
unknown_RA_summary_by_age_group <- unknown_RA_per_sample_age_group %>%
  group_by(Age_months) %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# Step 1: Add Subject_ID and Age_months to the unknown RA summary
id_map <- refeed_data %>% select(Seq_ID, Subject_ID, Age_months)

paired_unknown_RA <- unknown_RA_per_sample_refeed %>%
  left_join(id_map, by = "Seq_ID")

# Step 2: Keep only subjects with both 12 and 15 month values
paired_wide <- paired_unknown_RA %>%
  select(Subject_ID, Age_months, prop_unknown) %>%
  distinct() %>%
  pivot_wider(names_from = Age_months, values_from = prop_unknown) %>%
  drop_na(`12`, `15`)  # Ensure both values are present

# Step 3: Run paired Wilcoxon test on matched pairs
wilcox.test(paired_wide$`12`, paired_wide$`15`, paired = TRUE)


# Plot - as diff significant  ---------------------------------------------

# Prepare long-format data for plotting
paired_long <- paired_wide %>%
  pivot_longer(cols = c(`12`, `15`), names_to = "Age_months", values_to = "prop_unknown")


# Box plot with quasirandom points
plot_unknown_box <- ggplot(paired_long, aes(x = Age_months, y = prop_unknown, fill = Age_months, color = Age_months)) +
  geom_quasirandom(
    dodge.width = 0.9,
    shape = 21,
    size = 1,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.6,
    position = position_dodge(0.9),
    alpha = 0.5,
    outlier.colour = NA,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  scale_color_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  labs(
    y = "Proportion of Unknown Taxa",
    x = NULL
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )


# pdf
ggsave(file = "plots/chp2/clean_plots/plot_unknown_box.pdf", 
       plot = plot_unknown_box, units = "cm", width=5, height=7)


# 2: grA v grB ------------------------------------------------------------
# Filter the taxonomy table using group_data
taxonomy_group <- group_data %>%
  left_join(taxonomy, by = "Seq_ID")

# check Seq_IDs without UNKNOWN
seq_ids_without_unknown <- taxonomy_group %>%
  group_by(Seq_ID) %>%
  summarise(has_unknown = any(Taxa == "UNKNOWN")) %>%
  filter(!has_unknown) %>%
  pull(Seq_ID)

# how many samples have unknown
check <- taxonomy_group %>%
  filter(Taxa == "UNKNOWN") 

# Categorize Taxa as Labeled vs Unknown
taxonomy_group <- taxonomy_group %>%
  mutate(Taxa_Label = case_when(
    Taxa == "UNKNOWN" ~ "Unknown",
    TRUE ~ "Labeled"
  ))

# Calculate the correct proportion of UNKNOWN based on RA
unknown_RA_per_sample_group <- taxonomy_group %>%
  group_by(Seq_ID) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA
  )

# Now calculate the mean and standard deviation of UNKNOWN RA across samples
unknown_RA_summary_group <- unknown_RA_per_sample_group %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# by Group
# First, calculate UNKNOWN RA per sample and keep Group
unknown_RA_per_sample_group_grouped <- taxonomy_group %>%
  group_by(Seq_ID, Group) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA,
    .groups = "drop"
  )

# Now summarise by Group
unknown_RA_summary_by_group <- unknown_RA_per_sample_group_grouped %>%
  group_by(Group) %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# Wilcoxon rank-sum test comparing prop_unknown between groups
wilcox_unknown_RA_group <- wilcox.test(
  prop_unknown ~ Group,
  data = unknown_RA_per_sample_group_grouped
)



# 3: R v nR ---------------------------------------------------------------
# Filter the taxonomy table using recovery_data
taxonomy_recovery <- recovery_data %>%
  left_join(taxonomy, by = "Seq_ID")

# check Seq_IDs without UNKNOWN
seq_ids_without_unknown <- taxonomy_recovery %>%
  group_by(Seq_ID) %>%
  summarise(has_unknown = any(Taxa == "UNKNOWN")) %>%
  filter(!has_unknown) %>%
  pull(Seq_ID)

# how many samples have unknown
check <- taxonomy_recovery %>%
  filter(Taxa == "UNKNOWN") 

# Categorize Taxa as Labeled vs Unknown
taxonomy_recovery <- taxonomy_recovery %>%
  mutate(Taxa_Label = case_when(
    Taxa == "UNKNOWN" ~ "Unknown",
    TRUE ~ "Labeled"
  ))

# Calculate the correct proportion of UNKNOWN based on RA
unknown_RA_per_sample_recovery <- taxonomy_recovery %>%
  group_by(Seq_ID) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA
  )

# Now calculate the mean and standard deviation of UNKNOWN RA across samples
unknown_RA_summary_recovery <- unknown_RA_per_sample_recovery %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# by Recovery
# First, calculate UNKNOWN RA per sample and keep Recovery status
unknown_RA_per_sample_recovery_grouped <- taxonomy_recovery %>%
  group_by(Seq_ID, Recovery) %>%
  summarise(
    total_RA = sum(RA, na.rm = TRUE),
    unknown_RA = sum(RA[Taxa == "UNKNOWN"], na.rm = TRUE),
    prop_unknown = unknown_RA / total_RA,
    .groups = "drop"
  )

# Now summarise by Recovery
unknown_RA_summary_by_recovery <- unknown_RA_per_sample_recovery_grouped %>%
  group_by(Recovery) %>%
  summarise(
    mean_unknown_RA = mean(prop_unknown, na.rm = TRUE),
    sd_unknown_RA = sd(prop_unknown, na.rm = TRUE)
  )

# Wilcoxon rank-sum test comparing prop_unknown between recovery states
wilcox_unknown_RA_recovery <- wilcox.test(
  prop_unknown ~ Recovery,
  data = unknown_RA_per_sample_recovery_grouped
)




# Counting Taxa -----------------------------------------------------------

# Kingdom
kg_rf <- refeed_kingdom %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
ph_rf <- refeed_phylum %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
cl_rf <- refeed_class %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
or_rf <- refeed_order %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
fm_rf <- refeed_family %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
gn_rf <- refeed_genus %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 

# Species 
sp_rf <- refeed_species %>%
  gather(Taxa, RA, -Seq_ID) %>%  # Reshape dataset, excluding Seq_ID
  distinct(Taxa) 




# Test bay  ---------------------------------------------------------------

# Step 1: Add Subject_ID and Age_months to the unknown RA summary
id_map <- refeed_data %>% select(Seq_ID, Subject_ID, Age_months)

paired_unknown_RA <- unknown_RA_per_sample_refeed %>%
  left_join(id_map, by = "Seq_ID")

# Step 2: Keep only subjects with both 12 and 15 month values
paired_wide <- paired_unknown_RA %>%
  select(Subject_ID, Age_months, prop_unknown) %>%
  distinct() %>%
  pivot_wider(names_from = Age_months, values_from = prop_unknown) %>%
  drop_na(`12`, `15`)  # Ensure both values are present

# Step 3: Run paired Wilcoxon test on matched pairs
wilcox.test(paired_wide$`12`, paired_wide$`15`, paired = TRUE)


# Prepare long-format data for plotting
paired_long <- paired_wide %>%
  pivot_longer(cols = c(`12`, `15`), names_to = "Age_months", values_to = "prop_unknown")

# Box plot with quasirandom points
plot_unknown_box <- ggplot(paired_long, aes(x = Age_months, y = prop_unknown, fill = Age_months, color = Age_months)) +
  geom_quasirandom(
    dodge.width = 0.9,
    shape = 21,
    size = 1,
    alpha = 0.7
  ) +
  geom_boxplot(
    width = 0.6,
    position = position_dodge(0.9),
    alpha = 0.5,
    outlier.colour = NA,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  scale_color_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  labs(
    y = "Proportion of Unknown Taxa (log10 scale)",
    x = NULL
  ) +
  scale_y_log10() +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )








