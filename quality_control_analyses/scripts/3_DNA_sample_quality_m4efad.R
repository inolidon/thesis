## code summary 
# This script focuses on DNA yield, DNA purity and microbial load 
# The focus is on the results from lab extractions 

# ---- UNIT NOTES ----
# These values were provided by the sequencing company **prior to sequencing** 
# after their **quality control (QC) assessment**.
#
# Conc (Concentration) column: Measured in **ng/µL** (nanograms per microliter)
#   - Represents the DNA concentration in the sample.
#
# Volume column: Measured in **µL** (microliters)
#   - Represents the total volume of the DNA sample.
#
# Amount column: Now converted to **ng** (nanograms) after multiplying by 1000.
#   - Originally provided in **µg (micrograms)** and converted for consistency.
#   - Formula used: Amount (ng) = Concentration (ng/µL) × Volume (µL).
#
# These values were **QC-verified** by the sequencing provider before proceeding 
# with library preparation and sequencing.

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
library(ggpubr)
library(corrplot)
library(Hmisc)
library(FSA)  # For Dunn test

# load Data ---------------------------------------------------------------

# contains the DNA purity info
load("data/processed/dna_purity_quantity.Rdata")

# the kneaddata Rdata files, this contains the following files:
# kneaddata_raw = whole dataset includes 36 months samples + controls, 
# kneaddata_3yr_removed = 36 months samples and controls removed, 
# kneaddata_controls_only = controls only,
# kneaddata_no_controls = controls removed but contains 36 months samples
# kneaddata_metadata_combined = 36 months samples and controls removed and Sample_ID (Seq_ID) combined to metadata - main file used here
load("data/processed/kneaddata.Rdata")
load("data/processed/read_count_table.Rdata")

# select the columns of interest 
read <- read_count_table %>% 
  select(Seq_ID, type, Raw, Trimmed, Human_removal, Final, Contaminants, Lost_reads, Contaminant_prop)


# Combine the DNA_purity_quantity dataset to the read_count_table
reads_dna_purity <- read %>%
  left_join(dna_purity_quantity, 
            by = "Seq_ID")

# convert the Amount column to ng as its currently in micro grams
reads_dna_purity <- reads_dna_purity %>%
  mutate(Amount = Amount * 1000)



# investigate the difference in the quality of DNA samples between Conditions --------

# calculate average values for columns of interest

# List of columns to summarize
columns_to_summarize <- c("A260280", "A260230", "Conc", "Amount")

# Replace NA values with zero and calculate mean & SD for each Condition
reads_dna_summary <- reads_dna_purity %>%
  mutate(across(all_of(columns_to_summarize), ~replace_na(.x, 0))) %>%
  group_by(Condition) %>%
  summarise(
    across(all_of(columns_to_summarize), list(
      mean = ~mean(.x, na.rm = TRUE),
      sd = ~sd(.x, na.rm = TRUE)
    ), .names = "{col}_{fn}")
  )

# Display the summary table
print(reads_dna_summary)

# filter out control samples 
reads_dna_samples <- reads_dna_purity %>%
  filter(!is.na(Condition) & Condition != "control")


# stat test  --------------------------------------------------------------

wilcox_A260280 <- wilcox.test(A260280 ~ Condition, data = reads_dna_samples, exact = FALSE) # p-value = 0.306
wilcox_A260230 <- wilcox.test(A260230 ~ Condition, data = reads_dna_samples, exact = FALSE) # p-value = 2.507e-08
wilcox_conc <- wilcox.test(Conc ~ Condition, data = reads_dna_samples, exact = FALSE) # p-value = 0.2324
wilcox_amount <- wilcox.test(Amount ~ Condition, data = reads_dna_samples, exact = FALSE) # p-value = 0.9445

# plots -------------------------------------------------------------------

# plot - A260280 for Conditions 
qc_A260280_condition_plot <- ggplot(reads_dna_samples, aes(x = Condition, y = A260280, fill = Condition)) +
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
    legend.position = "none",  # Removes the legend
    axis.text.x = element_blank()  # Removes x-axis labels
  ) +
  xlab("") +
  ylab("A260/280 Ratio")  # Y-axis label


# plot - A260230 for Conditions 
qc_A260230_condition_plot <- ggplot(reads_dna_samples, aes(x = Condition, y = A260230, fill = Condition)) +
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
    legend.position = "none",  # Removes the legend
    axis.text.x = element_blank()  # Removes x-axis labels
  ) +
  xlab("") +
  ylab("A260/230 Ratio")  # Y-axis label




# plot - Conc for Conditions 
qc_conc_condition_plot <- ggplot(reads_dna_samples, aes(x = Condition, y = Conc, fill = Condition)) +
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
    legend.position = "none",  # Removes the legend
    axis.text.x = element_blank()  # Removes x-axis labels
  ) +
  xlab("") +
  ylab(expression("DNA concentration (ng/" * mu * "L)"))  # Y-axis label



# plot - Amount for Conditions 
qc_amount_condition_plot <- ggplot(reads_dna_samples, aes(x = Condition, y = Amount, fill = Condition)) +
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
    legend.position = "none",  # Removes the legend
    axis.text.x = element_blank()  # Removes x-axis labels
  ) +
  xlab("") +
  ylab(expression("Amount of DNA (ng)"))  # Y-axis label


# save plots  -------------------------------------------------------------
# all plots saved in PDF format 

# A260280
ggsave(file = "plots/qc/clean_plots/qc_A260280_condition_plot.pdf", 
       plot = qc_A260280_condition_plot, units = "cm", width=5, height=8)

# A260230
ggsave(file = "plots/qc/clean_plots/qc_A260230_condition_plot.pdf", 
       plot = qc_A260230_condition_plot, units = "cm", width=5, height=8)


# Conc
ggsave(file = "plots/qc/clean_plots/qc_conc_condition_plot.pdf", 
       plot = qc_conc_condition_plot, units = "cm", width=5, height=8)


# Amount
ggsave(file = "plots/qc/clean_plots/qc_amount_condition_plot.pdf", 
       plot = qc_amount_condition_plot, units = "cm", width=5, height=8)




# investigate which DNA metric influences Final seq depth  ----------------

# Amount (DNA yield) is used and not Conc as it is more deterministic of Seq depth than 
# conc as low yield results in insufficient amount of DNA to sequence 

# Q: Can the DNA metrics predict Raw and Final sequence depth?

# Fit the multiple regression model
lm_model <- lm(Final ~ A260280 + A260230 + Amount, data = reads_dna_samples)
lm_model2 <- lm(Raw ~ A260280 + A260230 + Amount, data = reads_dna_samples)

# Get summary of the model
summary(lm_model2)


# Q: which metrics are correlated with Final Seq depth?

# Select relevant columns
columns_to_check <- c("A260280", "A260230", "Conc", "Amount", "Final")

# Compute Spearman correlation matrix
cor_results <- rcorr(as.matrix(reads_dna_samples[, columns_to_check]), type = "spearman")

# Extract correlation coefficients and p-values
cor_matrix <- cor_results$r  # Correlation values
p_values <- cor_results$P    # P-values for significance

# Print results
print("Correlation Matrix:")
print(cor_matrix)

print("P-value Matrix:")
print(p_values)

# Convert p-values matrix to a data frame while keeping row and column names
p_values_df <- as.data.frame(p_values)

# Ensure row names are properly set
p_values_df <- cbind(Comparison = rownames(p_values_df), p_values_df)

# Convert cor_matrix matrix to a data frame while keeping row and column names
cor_values_df <- as.data.frame(cor_matrix)

# Ensure row names are properly set
cor_values_df <- cbind(Comparison = rownames(cor_values_df), cor_values_df)

## Which associations were significant?

# Define significance threshold
significance_threshold <- 0.05

# Convert all columns (except 'Comparison') to numeric for filtering
p_values_df[, -1] <- lapply(p_values_df[, -1], as.numeric)

# Replace non-significant values (p >= 0.05) with NA
significant_p_values_df <- p_values_df
significant_p_values_df[, -1] <- apply(significant_p_values_df[, -1], 2, function(x) ifelse(x < significance_threshold, x, NA))

# Print only significant values
print(significant_p_values_df)




# investigate the differences of DNA metrics between sequence batc --------

# Final seq depth differed between batches 
# Final seq depth is associated with all the DNA metrics 
# therefore, run Kruskal wallis for all metrics across batches 


# Define the columns to test
columns_to_test <- c("A260280", "A260230", "Amount", "Conc")

# Create an empty data frame to store results
kruskal_results_df <- data.frame(
  Variable = character(),
  Chi_Square = numeric(),
  P_Value = character(),  # Store p-values as characters for scientific notation
  Significance = character(),
  stringsAsFactors = FALSE
)

# Loop through each column and run Kruskal-Wallis test
for (col in columns_to_test) {
  test_result <- kruskal.test(reads_dna_samples[[col]] ~ as.factor(reads_dna_samples$Seq_batch))
  
  # Convert small p-values to scientific notation if below threshold
  formatted_p_value <- ifelse(test_result$p.value < 0.001, formatC(test_result$p.value, format = "e", digits = 3), round(test_result$p.value, 4))
  
  # Determine significance
  significance <- ifelse(test_result$p.value < 0.05, "Significant", "Not Significant")
  
  # Store results in the data frame
  kruskal_results_df <- rbind(kruskal_results_df, data.frame(
    Variable = col,
    Chi_Square = round(test_result$statistic, 4),
    P_Value = formatted_p_value,
    Significance = significance
  ))
}

# Print the final results dataset
print(kruskal_results_df)


# all the metrics differ significantly across batches 

# therefore, post hoc Pairwise comparisons: Dunn's test with BH correction

# A260280
A260280_pairwise <- reads_dna_samples %>%
  dunn_test(A260280 ~ Seq_batch, p.adjust.method = "BH")
print(A260280_pairwise)

# A260230
A260230_pairwise <- reads_dna_samples %>%
  dunn_test(A260230 ~ Seq_batch, p.adjust.method = "BH")
print(A260230_pairwise)

# Conc
Conc_pairwise <- reads_dna_samples %>%
  dunn_test(Conc ~ Seq_batch, p.adjust.method = "BH")
print(Conc_pairwise)

# Amount
Amount_pairwise <- reads_dna_samples %>%
  dunn_test(Amount ~ Seq_batch, p.adjust.method = "BH")
print(Amount_pairwise)








# make plots for DNA metrics across seq batches ---------------------------

# A260280 
qc_A260280_seqbatch_plot <- ggplot(reads_dna_samples, aes(x = Seq_batch, y = A260280, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
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
  ylab("A260/280 Ratio") 

# A260230 
qc_A260230_seqbatch_plot <- ggplot(reads_dna_samples, aes(x = Seq_batch, y = A260230, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
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
  ylab("A260/230 Ratio") 


# Conc 
qc_conc_seqbatch_plot <- ggplot(reads_dna_samples, aes(x = Seq_batch, y = Conc, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
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
  ylab(expression("DNA conc. (ng/" * mu * "L)"))  # Y-axis label

# Amount 
qc_amount_seqbatch_plot <- ggplot(reads_dna_samples, aes(x = Seq_batch, y = Amount, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
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
  ylab(expression("Amount of DNA (ng)"))  # Y-axis label




# save plots --------------------------------------------------------------

# A260280
ggsave(file = "plots/qc/clean_plots/qc_A260280_seqbatch_plot.pdf", 
       plot = qc_A260280_seqbatch_plot, units = "cm", width=4, height=4)

# A260230
ggsave(file = "plots/qc/clean_plots/qc_A260230_seqbatch_plot.pdf", 
       plot = qc_A260230_seqbatch_plot, units = "cm", width=4, height=4)

# Conc
ggsave(file = "plots/qc/clean_plots/qc_conc_seqbatch_plot.pdf", 
       plot = qc_conc_seqbatch_plot, units = "cm", width=4, height=4)

# Amount
ggsave(file = "plots/qc/clean_plots/qc_amount_seqbatch_plot.pdf", 
       plot = qc_amount_seqbatch_plot, units = "cm", width=4.3, height=4)
















# TESTING -----------------------------------------------------------------





library(ggplot2)
library(ggpubr)

# List of columns to analyze
columns_to_check <- c("A260280", "A260230", "Conc", "Amount")

# Function to plot distribution and normality check
plot_distribution_normality <- function(data, column) {
  p1 <- ggplot(data, aes_string(x = column)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "blue", alpha = 0.5) +
    geom_density(color = "red", size = 1) +
    ggtitle(paste("Distribution of", column)) +
    theme_minimal()
  
  p2 <- ggqqplot(data[[column]], title = paste("Q-Q Plot of", column))
  
  print(ggarrange(p1, p2, ncol = 2, nrow = 1))
  
  # Perform Shapiro-Wilk normality test
  shapiro_test <- shapiro.test(data[[column]])
  print(paste("Shapiro-Wilk Test for", column, ": W =", round(shapiro_test$statistic, 4), ", p =", round(shapiro_test$p.value, 4)))
  
  if (shapiro_test$p.value > 0.05) {
    print(paste(column, "data appears to be normally distributed (p > 0.05)."))
  } else {
    print(paste(column, "data does not appear to be normally distributed (p ≤ 0.05)."))
  }
}

# Apply function to each column
for (col in columns_to_check) {
  plot_distribution_normality(reads_dna_purity, col)
}




# Function to create box plots comparing Condition
plot_boxplot_condition <- function(data, column) {
  p <- ggplot(data, aes(x = Condition, y = !!sym(column), fill = Condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = 16, outlier.size = 2) +
    geom_jitter(width = 0.2, alpha = 0.5, size = 1.5) + # Adds individual points
    ggtitle(paste("Boxplot of", column, "by Condition")) +
    ylab(column) +
    theme_minimal() +
    theme(legend.position = "none") # Remove legend since Condition is on x-axis
  print(p)
}

# Generate box plots for each column
for (col in columns_to_check) {
  plot_boxplot_condition(reads_dna_purity, col)
}



# Function to perform Wilcoxon Rank Sum Test
run_wilcoxon_test <- function(data, column) {
  # Filter data for MAM and Healthy conditions
  filtered_data <- data %>% filter(Condition %in% c("MAM", "Healthy"))
  
  # Perform Wilcoxon test
  test_result <- wilcox.test(filtered_data[[column]] ~ filtered_data$Condition, exact = FALSE)
  
  # Print results
  print(paste("Wilcoxon Rank Sum Test for", column))
  print(paste("W =", round(test_result$statistic, 4), ", p =", round(test_result$p.value, 4)))
  
  if (test_result$p.value < 0.05) {
    print("Significant difference (p < 0.05)\n")
  } else {
    print("No significant difference (p >= 0.05)\n")
  }
}

# Run Wilcoxon test for each column
for (col in columns_to_check) {
  run_wilcoxon_test(reads_dna_purity, col)
}



# Select relevant columns
columns_to_check <- c("A260280", "A260230", "Conc", "Amount", "Raw", "Final", "Contaminants")

# Filter dataset to remove missing values
reads_dna_purity_clean <- reads_dna_purity %>% select(all_of(columns_to_check)) %>% na.omit()

# Compute Spearman correlation matrix
cor_matrix <- cor(reads_dna_purity_clean, method = "spearman")

# Plot correlation heatmap
corrplot(cor_matrix, method = "color", type = "upper", tl.col = "black", tl.srt = 45)

# Function to create scatter plots
plot_relationship <- function(data, x_col, y_col) {
  ggplot(data, aes_string(x = x_col, y = y_col)) +
    geom_point(color = "blue", alpha = 0.6) +
    geom_smooth(method = "lm", color = "red", se = FALSE) +
    ggtitle(paste("Scatter plot of", x_col, "vs", y_col)) +
    theme_minimal()
}

# Generate scatter plots for each combination
plot_list <- list()
for (x_col in c("A260280", "A260230", "Conc", "Amount")) {
  for (y_col in c("Raw", "Final", "Contaminants")) {
    plot_list[[paste(x_col, y_col, sep = "_")]] <- plot_relationship(reads_dna_purity_clean, x_col, y_col)
  }
}

# Arrange scatter plots in a grid
ggarrange(plotlist = plot_list, ncol = 3, nrow = 4)


# Filter dataset to remove missing values
reads_dna_purity_clean <- reads_dna_purity %>% select(all_of(columns_to_check)) %>% na.omit()

# Compute correlation coefficients and p-values
cor_results <- rcorr(as.matrix(reads_dna_purity_clean), type = "spearman")

# Extract correlation matrix
cor_matrix <- cor_results$r

# Extract p-value matrix
p_value_matrix <- cor_results$P

# Print correlation coefficients
print("Correlation Coefficients:")
print(cor_matrix)

# Print p-values
print("P-values:")
print(p_value_matrix)

# Convert p-value matrix to a data frame
p_values_df <- as.data.frame(as.table(p_value_matrix))

# Filter for p-values < 0.05
significant_p_values <- p_values_df %>%
  filter(Freq < 0.05)

# Print significant p-values
print("Significant p-values (p < 0.05):")
print(significant_p_values)


# test plot 
testplot <- ggplot(reads_dna_samples, aes(x = Condition, y = A260280, fill = Condition)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Condition)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4"), name = "Condition") +  # Boxplot fill colors
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4"), name = "Condition") +  # Point colors
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 10),  # Edited text size for axis titles
    axis.title.y = element_text(size = 10),  # Ensure y-axis label matches axis title size
    legend.position = "right",  # **Adds the legend back**
    legend.title = element_text(size = 10),  # Legend title size
    legend.text = element_text(size = 9),  # Legend text size
    axis.text.x = element_blank()  # Removes x-axis labels
  ) +
  xlab("") +
  ylab("A260/280 Ratio")  # Y-axis label

# A260280
ggsave(file = "plots/qc/clean_plots/TEST.pdf", 
       plot = testplot, units = "cm", width=5, height=8)


# Updated Plot - Amount of DNA by Seq_batch with Legend
testplot2seqbatches <- ggplot(reads_dna_samples, aes(x = Seq_batch, y = Amount, fill = Seq_batch)) +
  geom_quasirandom(dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7, aes(color = Seq_batch)) +  # Points without black border
  geom_boxplot(outlier.colour = NA, alpha = 0.5, color = "black", linewidth = 0.3) +  # Thinner box plot lines
  theme_bw() +
  theme(
    panel.border = element_blank(),  # Removes the border around the plot
    panel.grid = element_blank(),  # Removes all grid lines
    axis.text = element_text(),  # Edited text size for axis labels
    axis.line = element_line(),  # Keeps only x and y axis lines
    axis.title = element_text(size = 9),  # Edited text size for axis titles
    axis.title.y = element_text(size = 9),  # Ensure y-axis label matches axis title size
    legend.position = "right",  # **Adds the legend back**
    legend.title = element_text(size = 9),  # Set legend title size
    legend.text = element_text(size = 8)  # Set legend text size
  ) +
  xlab("Seq Batch") +
  ylab(expression("Amount of DNA (ng)")) +  # Y-axis label
  scale_fill_discrete(name = "Seq Batch") +  # Ensures legend title appears for fill
  scale_color_discrete(name = "Seq Batch")  # Ensures legend title appears for color

# A260280
ggsave(file = "plots/qc/clean_plots/SeqbatchesTEST.pdf", 
       plot = testplot2seqbatches, units = "cm", width=5, height=8)
