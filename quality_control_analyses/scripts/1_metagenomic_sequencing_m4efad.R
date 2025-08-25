## code summary 
# This script focuses on quality control for the metagenomic sequences. 
# The focus is on the results from Kneaddata workflow 

# load Packages -----------------------------------------------------------
library(tidyverse) # required for data wrangling and ggplot
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(scales) # required for scale_y_continuous() which prevents y axis on plot from defaulting to scientific notation
library(dplyr)
library(vegan)
library(readxl)
library(car) # for Levene's Test - testing whether variances are equal between groups
library(rstatix)
library(ggplot2)
library(writexl)


# load Data ---------------------------------------------------------------

# raw kneaddata output table available in data folder - has been converted to Rdata as below
# kneaddata <- read_tsv("data/raw/m4efad_kneaddata_read_counts_april2024.tsv")

# the kneaddata Rdata files, this contains the following files:
# kneaddata_raw = whole dataset includes 36 months samples + controls, 
# kneaddata_3yr_removed = 36 months samples and controls removed, 
# kneaddata_controls_only = controls only,
# kneaddata_no_controls = controls removed but contains 36 months samples
# kneaddata_metadata_combined = 36 months samples and controls removed and Sample_ID (Seq_ID) combined to metadata - main file used here
load("data/processed/kneaddata.Rdata")

# contains species richness and Shannon diversity values 
load("data/processed/sp_metrics.Rdata")

# contains species taxa RA for each sample 
load("data/processed/sp_ra.Rdata")

# Number of reads at each qc stage  ---------------------------------------
# key stages: 
# trimming = removal of adaptor sequences and low quality sequences 
# human removal = removal of contaminant host DNA and other contaminating DNA or RNA (example, enviromental duirng extraction)

# Combine read pairs/orphans (unmatched single reads) to calculate total number of reads before, during, and after QC steps
read_count_table <- kneaddata_metadata_combined %>% 
  mutate(Raw = `raw pair1`+`raw pair2`,
         Trimmed = `trimmed pair1`+`trimmed pair2`+`trimmed orphan1`+`trimmed orphan2`,
         Human_removal = `decontaminated hg37dec_v0.1 pair1`+`decontaminated hg37dec_v0.1 pair2`+`decontaminated hg37dec_v0.1 orphan1`+`decontaminated hg37dec_v0.1 orphan2`,
         Final = `final pair1`+`final pair2`+`final orphan1`+`final orphan2`,
         Lost_reads = Raw-Final)


# Calculate the average + sd Raw the whole dataset
averages_whole_dataset <- read_count_table %>%
  summarize(
    avg_Raw = mean(Raw, na.rm = TRUE),
    sd_Raw = sd(Raw, na.rm = TRUE),
    avg_Final = mean(Final, na.rm = TRUE),
    sd_Final = sd(Final, na.rm = TRUE)
  )

# Calculate the average + sd Raw and Final counts for each Condition
averages <- read_count_table %>%
  group_by(Condition) %>%
  summarize(
    avg_Raw = mean(Raw, na.rm = TRUE),
    sd_Raw = sd(Raw, na.rm = TRUE),
    avg_Final = mean(Final, na.rm = TRUE),
    sd_Final = sd(Final, na.rm = TRUE)
  )

# Test distribution of the Raw and Final read counts for normality using the Shapiro-Wilk test
# Final counts 
shapiro_final <- shapiro.test(read_count_table$Final)
print(shapiro_final) # p < 2.2e-16, not a normal distribution 

# Raw counts 
shapiro_raw <- shapiro.test(read_count_table$Raw)
print(shapiro_raw) # p < 2.2e-16, not a normal distribution 

# plot of read counts at stage
# create a table with each stage of the qc process of interest 
qc_steps <- read_count_table %>% 
  select(Seq_ID, Raw, Trimmed, Final) %>% 
  gather(Steps, Reads, -Seq_ID) %>% 
  mutate(Steps = factor(Steps, levels = c("Raw", "Trimmed", "Final")))

# plot 
whole_dataset <- ggplot(qc_steps, aes(x = Steps, y = Reads)) +
  geom_quasirandom(dodge.width = 0.75, shape = 21, size = 0.1) +  # Edited size of the points
  geom_boxplot(outlier.colour = NA, alpha = 0.5, fill = "slategray3") +
  theme_bw() +
  theme(
    axis.text = element_text(),  # Edited text size for axis labels, the x and y 
    axis.title = element_text(),  # Edited text size for axis titles
    legend.text = element_text(),  # Edited text size for legend text
    legend.title = element_text()  # Edited text size for legend title
  ) +
  xlab("") +
  ylab("Read count")

# pdf
# ggsave(file = "plots/qc/clean_plots/whole_dataset.pdf", 
# plot = whole_dataset, units = "cm", width=8, height=8)

# svg
# ggsave(file = "plots/qc/clean_plots/whole_dataset.svg", 
#       plot = whole_dataset, units = "cm", width=8, height=8)


# investigate if Final read counts differ between Conditions  -------------

# Final quality control reads
# Perform a Wilcox test (non parametric test) to compare Final counts between MAM and Healthy
wilcox_test_result <- wilcox.test(Final ~ Condition, data = read_count_table)
print(wilcox_test_result)
# p = 0.7877 no significant difference in the final read counts

# Trimmed stage
# Perform a Wilcox test (non parametric test) to compare Trimmed counts between MAM and Healthy
wilcox_test_result_trim <- wilcox.test(Trimmed ~ Condition, data = read_count_table)
print(wilcox_test_result_trim)
# p = 0.051 no significant difference in the final read counts

# Raw stage
# Perform a Wilcox test (non parametric test) to compare Trimmed counts between MAM and Healthy
wilcox_test_result_raw <- wilcox.test(Raw ~ Condition, data = read_count_table)
print(wilcox_test_result_raw)
# p = 0.02567 significant difference in the final read counts

# Arrange dataset for plotting: convert from wide to long format
qc_steps_condition <- read_count_table %>%
  select(Seq_ID, Raw, Trimmed, Final, Condition) %>%
  pivot_longer(cols = c(Raw, Trimmed, Final), 
               names_to = "Steps", 
               values_to = "Reads") %>%
  mutate(Steps = factor(Steps, levels = c("Raw", "Trimmed", "Final"))) %>% 
  filter(Condition %in% c("MAM", "Healthy")) %>% 
  group_by(Steps)

# Filter and recode Condition if necessary (here we assume that Conditions are already "MAM" and "Healthy")
qc_steps_condition_filtered <- qc_steps_condition %>%
  filter(Condition %in% c("MAM", "Healthy")) %>% 
  group_by(Steps)

# plot 
qc_steps_condition_filtered_plot <- ggplot(qc_steps_condition, aes(x = Steps, y = Reads, fill = Condition)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Condition)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +  # Boxplot fill colors
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +  # Point colors
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 9),  # Edited text size for axis titles
    axis.title.y = element_text(size = 9),  # Ensure y-axis label matches axis title size
    legend.text = element_text(),  # Edited text size for legend text
    legend.title = element_text()  # Edited text size for legend title
  ) +
  xlab("") +
  ylab("Sequencing depth")


# pdf
ggsave(file = "plots/qc/clean_plots/qc_steps_condition_filtered_plot.pdf", 
 plot = qc_steps_condition_filtered_plot, units = "cm", width=12, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/qc_steps_condition_filtered_plot.svg", 
      plot = qc_steps_condition_filtered_plot, units = "cm", width=12, height=8)



# investigate the proportion of contaminant reads (host and other DNA)  ----------------------------------------
# Q: What proportion of total reads were contaminant reads in samples and the two types of controls?

# step1: create s new column in sample and control data sets

# dataset with only samples 
kneaddata_samples <- kneaddata_metadata_combined %>% 
  mutate(type = "Sample")

# separate the microbial controls from the DNA controls as they are used for different purposes 
# i.e., they are not equal in what they control (Seq bias v extraction bias)
kneaddata_controls <- kneaddata_controls_only %>%
  mutate(type = case_when(
    str_starts(Subject_ID, "P") ~ "Control_Microbial",
    str_starts(Subject_ID, "D") ~ "Control_DNA"
   ))

# step2: row bind two datasets
kneaddata_samples_controls <- bind_rows(kneaddata_samples, 
                                        kneaddata_controls)


# calculate proportion of contaminant reads 
read_count_table <- kneaddata_samples_controls %>% 
  mutate(Raw = `raw pair1`+`raw pair2`,
         Trimmed = `trimmed pair1`+`trimmed pair2`+`trimmed orphan1`+`trimmed orphan2`,
         Human_removal = `decontaminated hg37dec_v0.1 pair1`+`decontaminated hg37dec_v0.1 pair2`+`decontaminated hg37dec_v0.1 orphan1`+`decontaminated hg37dec_v0.1 orphan2`,
         Final = `final pair1`+`final pair2`+`final orphan1`+`final orphan2`,
         Contaminants = Trimmed-Final,
         Lost_reads = Raw-Final,
         Contaminant_prop = Contaminants/Raw)

# save(read_count_table, file = "data/processed/read_count_table.Rdata")

# Plot the contaminant proportion - all data 
contaminant_plot <- ggplot(read_count_table, aes(x = type, y = Contaminant_prop)) +
  geom_quasirandom(dodge.width = 0.75, shape = 21, size = 0.1) +  # Edited size of the points
  geom_boxplot(outlier.colour = NA, alpha = 0.5, fill = "slategray3") +
  theme_bw() +
  theme(
    axis.text = element_text(),  # Edit text size for axis labels
    axis.title = element_text(),  # Edit text size for axis titles
    legend.text = element_text(),  # Edit text size for legend text
    legend.title = element_text()  # Edit text size for legend title
  ) +
  xlab("") +
  ylab("Contaminant Proportion") 

# pdf
# ggsave(file = "plots/qc/clean_plots/contaminant_plot.pdf", 
# plot = contaminant_plot, units = "cm", width=9, height=8)

# svg
# ggsave(file = "plots/qc/clean_plots/contaminant_plot.svg", 
#      plot = contaminant_plot, units = "cm", width=9, height=8)


# Summarise the results into a table 
contaminant_prop_summary <- read_count_table %>%
  group_by(type) %>%
  summarize(
    `>50%` = sum(Contaminant_prop > 0.5, na.rm = TRUE),
    `>75%` = sum(Contaminant_prop > 0.75, na.rm = TRUE),
    `<10%` = sum(Contaminant_prop < 0.1, na.rm = TRUE),
    `<1%` = sum(Contaminant_prop < 0.01, na.rm = TRUE)
  )

# percentage of samples & controls that have <10% contamination
contaminants_less_than_10_percent <- read_count_table %>%
  group_by(type) %>%
  summarize(
    total_samples = n(),
    samples_less_than_0_1 = sum(Contaminant_prop < 0.1, na.rm = TRUE)
  ) %>%
  mutate(percentage = (samples_less_than_0_1 / total_samples) * 100)
#  82.3% of samples, 100% of controls


# contaminants between MAM v Healthy 

# calculate proportion of contaminant reads 
read_count_table_conditions <- kneaddata_samples_controls %>% 
  mutate(Raw = `raw pair1`+`raw pair2`,
         Trimmed = `trimmed pair1`+`trimmed pair2`+`trimmed orphan1`+`trimmed orphan2`,
         Human_removal = `decontaminated hg37dec_v0.1 pair1`+`decontaminated hg37dec_v0.1 pair2`+`decontaminated hg37dec_v0.1 orphan1`+`decontaminated hg37dec_v0.1 orphan2`,
         Final = `final pair1`+`final pair2`+`final orphan1`+`final orphan2`,
         Contaminants = Trimmed-Final,
         Lost_reads = Raw-Final,
         Contaminant_prop = Contaminants/Raw) %>% 
  filter(Condition %in% c("MAM", "Healthy")) %>% 
  group_by(Condition)

# Final quality control reads
# Perform a Wilcox test (non parametric test) to compare Final counts between MAM and Healthy
wilcox_test_result_contaminats_conditions <- wilcox.test(Contaminants ~ Condition, data = read_count_table_conditions)
print(wilcox_test_result_contaminats_conditions)
# p = 2.381e-06 significant difference in the final read counts


# Calculate the average + sd Contaminants each Condition
average_contaminants <- read_count_table_conditions %>%
  group_by(Condition) %>%
  summarize(
    avg_Contaminants = mean(Contaminants, na.rm = TRUE),
    sd_Contaminants = sd(Contaminants, na.rm = TRUE)
  )


# Plot the contaminant proportion - MAM v Healthy
contaminant_conditions_plot <- ggplot(read_count_table_conditions, aes(x = Condition, y = Contaminant_prop, fill = Condition)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Condition)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +  # Boxplot fill colors
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +  # Point colors
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 10),  # Edited text size for axis titles
    axis.title.y = element_text(size = 10),  # Ensure y-axis label matches axis title size
    legend.position = "none"  # Removes the legend
  ) +
  xlab("") +
  ylab("Contaminant Proportion")


# pdf
ggsave(file = "plots/qc/clean_plots/contaminant_conditions_plot.pdf", 
plot = contaminant_conditions_plot, units = "cm", width=8, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/contaminant_conditions_plot.svg", 
     plot = contaminant_conditions_plot, units = "cm", width=8, height=8)



# investigate the proportion of reads below target sequencing depth -------
# Q: What proportion of total reads were below the sequencing depth (40mil paired)

# Filter and count reads in the Final column below 40 million, separated by type
reads_below_40_million <- read_count_table %>%
  filter(Final < 40000000) %>%
  group_by(type) 
# 175 samples had reads below 40mil seq depth - all controls above target seq depth

# Q: what proportion of reads had less than half the target seq depth? (40mil target, therefore 20mil is half)
reads_below_20_million <- read_count_table %>%
  filter(Final < 20000000) %>%
  group_by(type) 
# 22 samples had reads below 20 million reads - - all controls above target seq depth

# proportion calculation
# below 40 million 
proportion_below_40_million <- read_count_table %>%
  group_by(type) %>%
  summarize(
    total_reads = n(),
    reads_below_40_million = sum(Final < 40000000)
  ) %>%
  mutate(proportion = (reads_below_40_million / total_reads) * 100)

proportion_below_40_million_condition <- read_count_table %>%
  group_by(Condition) %>%
  summarize(
    total_reads = n(),
    reads_below_40_million = sum(Final < 40000000)
  ) %>%
  mutate(proportion = (reads_below_40_million / total_reads) * 100)

# below 20 million 
proportion_below_20_million <- read_count_table %>%
  group_by(type) %>%
  summarize(
    total_reads = n(),
    reads_below_20_million = sum(Final < 20000000)
  ) %>%
  mutate(proportion = (reads_below_20_million / total_reads) * 100)

# below 20 million - Conditions 
proportion_below_20_million_condition <- read_count_table %>%
  group_by(Condition) %>%
  summarize(
    total_reads = n(),
    reads_below_20_million = sum(Final < 20000000)
  ) %>%
  mutate(proportion = (reads_below_20_million / total_reads) * 100)


# plot of read counts of final stage
qc_final_only <- read_count_table %>% 
  select(Seq_ID, Final) %>% 
  gather(Steps, Reads, -Seq_ID) %>% 
  mutate(Steps = factor(Steps, levels = c("Final"))) %>% 
  rename("Post quality control sequences" = Steps)

# plot 
final_reads_only <- ggplot(qc_final_only, aes(x = "Post quality control sequences", y = Reads)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, color = "slategray3") +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", fill = "slategray3", linewidth = 0.3) +  # Thinner box plot lines
  geom_hline(yintercept = 20000000, linetype = "dashed", color = "red", size = 0.7) +  # Dashed red line at 20M reads
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 9),  # Edited text size for axis titles
    axis.title.y = element_text(size = 9),  # Ensure y-axis label matches axis title size
    legend.position = "none"  # Removes the legend
  ) +
  xlab("") +
  ylab("Sequencing depth")


# pdf
ggsave(file = "plots/qc/clean_plots/final_reads_only.pdf", 
       plot = final_reads_only, units = "cm", width=8, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/final_reads_only.svg", 
       plot = final_reads_only, units = "cm", width=8, height=8)



# investigate the relationship between DNA quality and low sequenc --------

# Q: Do the low seq depth samples have lower quality DNA?
# we will use the qc_report from the manufacturer to check this

# load data
load("data/processed/qc_reports.Rdata")

# combine the qc report to the low seq depth samples 
reads_below_20_million <- reads_below_20_million %>%
  inner_join(qc_reports, by = "Seq_ID") %>%  # Merge qc_reports with reads_below_20_million based on Seq_ID
  mutate(Grade = factor(Grade, levels = unique(Grade)))  # Convert Grade to a factor with custom level order

# Define colors for each grade
grade_colors <- c("A" = "darkseagreen",  # High quality, high quantity 
                  "B" = "azure2",  # Medium quality, sufficient quantity 
                  "C" = "azure3",  # Low quality, low quantity 
                  "D" = "azure4")  # Very poor quality, very low quantity 

# Plot stacked bar plot 
low_seq_sample_grades_plot <- ggplot(reads_below_20_million, aes(x = Grade, fill = Grade)) +
  geom_bar(aes(color = Condition), position = "stack", linewidth = 2) +  # Stacked bars with condition-based border color
  scale_fill_manual(values = grade_colors,
                    labels = c(
                      "A" = "High quality, high quantity",
                      "B" = "Medium quality, sufficient quantity",
                      "C" = "Low quality, low quantity",
                      "D" = "Very poor quality, very low quantity")) +  # Fill bars by Grade
  scale_color_manual(values = c("MAM" = "darkgoldenrod2", "Healthy" = "cyan3")) +  # Border color by Condition
  theme_bw() +
  ylab("Number of Samples") +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.text = element_text(),  # Edit text size for axis labels
    axis.title = element_text(),  # Edit text size for axis titles
    legend.text = element_text(size = 10),  # Edit text size for legend text
    legend.title = element_text(size = 10)  # Edit text size for legend title
  )

# pdf
ggsave(file = "plots/qc/clean_plots/low_seq_sample_grades_plot.pdf", 
       plot = low_seq_sample_grades_plot, units = "cm", width=16, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/low_seq_sample_grades_plot.svg", 
       plot = low_seq_sample_grades_plot, units = "cm", width=16, height=8)


# investigate Final read counts between sequencing batches ------------------------------------------------------------
# Is there a difference in the final read counts between Seq_batches?  

# number of samples and controls in each batch 
seq_batch_read_count_summary <- read_count_table %>%
  group_by(Seq_batch, type, Condition) %>%
  summarise(count = n(), .groups = "drop")

# save the file for sample numbers
write_xlsx(seq_batch_read_count_summary, "data/processed/seq_batch_read_count_summary.xlsx")

# Samples
# Calculate mean final reads for each sequencing batch
mean_reads_samples <- read_count_table %>%
  filter(type == "Sample") %>%
  group_by(Seq_batch, Condition) %>%
  summarize(
    mean_Final = mean(Final, na.rm = TRUE),
    sd_Final = sd(Final, na.rm = TRUE),
    .groups = "drop"
  )

# save the file for sample numbers
write_xlsx(mean_reads_samples, "data/processed/mean_reads_samples.xlsx")

# Is there a difference in the read counts between the same type of controls?
# i.e., between the 5 microbial controls 
# between the 5 DNA controls 

# MB controls
mean_reads_control_mb <- read_count_table %>%
  filter(type == "Control_Microbial") %>%
  group_by(Seq_batch) %>%
  summarize(
    mean_Final = mean(Final, na.rm = TRUE),
    sd_Final = sd(Final, na.rm = TRUE),
    .groups = "drop"
  )

# DNA controls
mean_reads_control_dna <- read_count_table %>%
  filter(type == "Control_DNA") %>%
  group_by(Seq_batch) %>%
  summarize(
    mean_Final = mean(Final, na.rm = TRUE),
    sd_Final = sd(Final, na.rm = TRUE),
    .groups = "drop"
  )

# output summary of means
mean_read_summary_table <- mean_reads_samples %>%
  full_join(mean_reads_control_mb, by = "Seq_batch") %>%
  full_join(mean_reads_control_dna, by = "Seq_batch")

# rename appropriately
mean_read_summary_table <- mean_read_summary_table %>% 
  rename(Sample_mean = mean_Final.x,
         Sample_sd = sd_Final.x,
         Control_Microbial_mean = mean_Final.y,
         Control_Microbial_sd = sd_Final.y,
         Control_DNA_mean = mean_Final,
         Control_DNA_sd = sd_Final)

# Save the dataframe as an Excel file
write_xlsx(mean_read_summary_table, "data/processed/mean_read_summary_table.xlsx")

# run statistical tests to compare the differences between batches  --------------------------------------------------------------

# *for Samples ------------------------------------------------------------

# Prepare the data by filtering only samples
samples_data <- read_count_table %>%
  filter(type == "Sample")

# Perform Kruskal-Wallis test - as the distribution is not normal 
kruskal_result <- kruskal.test(Final ~ Seq_batch, data = samples_data)
# p-value = 0.003184 - significant differences in Final read counts between sample batches 

# for raw
kruskal_result2 <- kruskal.test(Raw ~ Seq_batch, data = samples_data)

# Pairwise comparisons: Dunn's test with BH correction
pairwise_results <- samples_data %>%
  dunn_test(Final ~ Seq_batch, p.adjust.method = "BH")
print(pairwise_results)

# b2 v b4, b3 v b5 and b4 v b5 differ significantly in their final read counts 

# Raw - Pairwise comparisons: Dunn's test with BH correction
pairwise_results2 <- samples_data %>%
  dunn_test(Raw ~ Seq_batch, p.adjust.method = "BH")
print(pairwise_results2)


# *for Controls -----------------------------------------------------------

# microbial ---------------------------------------------------------------
# Prepare the data by filtering only mb controls
microbial_data <- read_count_table %>%
  filter(type == "Control_Microbial")

# Perform Kruskal-Wallis test
kruskal_result_mb <- kruskal.test(Final ~ Seq_batch, data = microbial_data)
# p = 0.3023 


# dna ---------------------------------------------------------------------
# Prepare the data by filtering only dna controls
dna_data <- read_count_table %>%
  filter(type == "Control_DNA")

# Perform Kruskal-Wallis test
kruskal_result_dna <- kruskal.test(Final ~ Seq_batch, data = dna_data)
# p = 0.4461


# Plot the final read counts for each Seq_batch 
samples_by_seq_batch_plot <- ggplot(samples_data, aes(x = Seq_batch, y = Final, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border, colored by batch
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines with black border
  scale_fill_manual(values = setNames(scales::hue_pal()(length(unique(samples_data$Seq_batch))), unique(samples_data$Seq_batch))) +  # Fill colors for batches
  scale_color_manual(values = setNames(scales::hue_pal()(length(unique(samples_data$Seq_batch))), unique(samples_data$Seq_batch))) +  # Match point colors to batches
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 9),  # Edited text size for axis titles
    axis.title.y = element_text(size = 9),  # Ensure y-axis label matches axis title size
    legend.text = element_text(),  # Edited text size for legend text
    legend.title = element_text(),  # Edited text size for legend title
    legend.position = "none"  # Removes the legend
  ) +
  xlab("Sequencing batch") +
  ylab("Sequencing depth")



# pdf
 ggsave(file = "plots/qc/clean_plots/samples_by_seq_batch_plot.pdf", 
    plot = samples_by_seq_batch_plot, units = "cm", width=10, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/samples_by_seq_batch_plot.svg", 
     plot = samples_by_seq_batch_plot, units = "cm", width=10, height=8)




# Plot the final read counts for each Seq_batch 
samples_by_seq_batch_plot2 <- ggplot(samples_data, aes(x = Seq_batch, y = Raw, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border, colored by batch
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines with black border
  scale_fill_manual(values = setNames(scales::hue_pal()(length(unique(samples_data$Seq_batch))), unique(samples_data$Seq_batch))) +  # Fill colors for batches
  scale_color_manual(values = setNames(scales::hue_pal()(length(unique(samples_data$Seq_batch))), unique(samples_data$Seq_batch))) +  # Match point colors to batches
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 9),  # Edited text size for axis titles
    axis.title.y = element_text(size = 9),  # Ensure y-axis label matches axis title size
    legend.text = element_text(),  # Edited text size for legend text
    legend.title = element_text(),  # Edited text size for legend title
    legend.position = "none"  # Removes the legend
  ) +
  xlab("Sequencing batch") +
  ylab("Sequencing depth")



# Optional plot 
# combine the datasets by row bind 
control_data <- bind_rows(microbial_data, dna_data)

# plot
controls_by_seq_batch_plot <- # Create the plot
  ggplot(control_data, aes(x = Seq_batch, y = Final)) +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  geom_quasirandom(aes(color = Seq_batch), width = 0.2, alpha = 0.5) +  # Add quasirandom points with color
  scale_y_continuous(labels = scales::label_comma()) +
  theme_bw() +
  ylab("Read count") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +  # Rotate x-axis labels for better readability
  theme(
    axis.text = element_text(),  # Edit text size for axis labels, the x and y 
    axis.title = element_text(),  # Edit text size for axis titles
    legend.text = element_text(),  # Edit text size for legend text
    legend.title = element_text()  # Edit text size for legend title
  )



# investigate the influence of final read counts on Species richness  ---------------------------------
# Q: Did the sequencing depth influence species richness?

# check normality assumption to determine test 
shapiro_final <- shapiro.test(sp_metrics$Final) 
shapiro_richness <- shapiro.test(sp_metrics$Richness)

# Correlation test:
# Use Pearson if normal, Spearman if not
if(shapiro_final$p.value > 0.05 & shapiro_richness$p.value > 0.05) {
  cor_test <- cor.test(sp_metrics$Final, sp_metrics$Richness, method = "pearson")
} else {
  cor_test <- cor.test(sp_metrics$Final, sp_metrics$Richness, method = "spearman")
}
print(cor_test)

# p = 0.0001171
# rho = 0.1567495 
# there is a significant positive relationship between the final read counts and species richness

# plot of read counts against species richness 
richness_read_count_plot <- ggplot(sp_metrics, aes(x = Final, y = Richness)) +
  geom_point(color = "darkblue", size = 1) +  # Plot individual points
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the full border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.line = element_line()  # Keeps only x and y axis lines
  ) +
  labs(x = "Post quality control sequence depth",
       y = "Species Richness")



# pdf
ggsave(file = "plots/qc/clean_plots/richness_read_count_plot.pdf", 
      plot = richness_read_count_plot, units = "cm", width=10, height=8)

# svg
ggsave(file = "plots/qc/clean_plots/richness_read_count_plott.svg", 
     plot = richness_read_count_plot, units = "cm", width=10, height=8)




