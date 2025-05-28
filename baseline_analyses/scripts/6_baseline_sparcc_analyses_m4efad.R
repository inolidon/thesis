## code summary 
# This script uses the SparCC algorithm to obtain covarience. 

# SparCC is used as it accounts for the compositionality of the gut microbiome dataset
# it also deals with the sparsity of the data

# Input data must be normalized to be relative, ie, sum to 1. 

# No prevalence filter must be added, data set must be in the raw format



# load packages  ----------------------------------------------------------

library(tidyverse) # required for data wrangling and ggplot
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(data.table) # require for gene table cleaning
library(dplyr)
library(stringr)
library(readxl)
library(reshape2) # has the melt funct  ion, which I use to wrangle data
library(igraph)


# load data ---------------------------------------------------------------

# use the baseline_species dataset to find co-occurences of species 
# this dataset does not have a prevalence filter 
# subsetted the data into health conditions 

load("data/processed/baseline_species.Rdata")
load("data/processed/baseline_data.Rdata") # load metadata 


# load tidy output files  -------------------------------------------------
# load the modified outputs directly here:
load("data/processed/mam_cor_matrix.Rdata") # raw matrix - MAM
load("data/processed/mam_pvals_matrix.Rdata") # raw pvalue matrix - MAM
load("data/processed/healthy_cor_matrix.Rdata") # raw matrix - Healthy
load("data/processed/healthy_pvals_matrix.Rdata") # raw pvalue matrix - Healthy
load("data/processed/mam_cor_table.Rdata") # correlation table format - MAM
load("data/processed/mam_pvals_table.Rdata") # p value table format - MAM
load("data/processed/healthy_cor_table.Rdata") # correlation table format - Healthy
load("data/processed/healthy_pvals_table.Rdata")  # p value table format - Healthy
load("data/processed/format_sig_mam_cor_table.Rdata") # cor table filtered for significant pairwise - MAM
load("data/processed/format_sig_healthy_cor_table.Rdata") # cor table filtered for significant pairwise - Healthy


# Prepare Data for SparCC INPUT  -----------------------------------------------------------
# We will be producing independent networks for each variable
# In this script, the two variables are the two conditions: MAM and Healthy
# We will filter the main metadata file and the species dataframe for each condition

# (make sure the baseline_data is in the right format - has a column called "Seq_ID")

# run this code if the baseline data is not in the right format 
baseline_data <- baseline_data %>% 
  rownames_to_column("Seq_ID")


# - create individual data  -----------------------------------------------

# create individual metadata
M_data <- baseline_data %>% 
  filter(Condition == "MAM")

H_data <- baseline_data %>% 
  filter(Condition == "Healthy")

# create individual species data for each condition
MAM_sp <- baseline_species %>% 
  semi_join(M_data, by = "Seq_ID")

Healthy_sp <- baseline_species %>% 
  semi_join(H_data, by = "Seq_ID")

# convert datasets to have the Seq_ID as a rowname

MAM_sp <- MAM_sp %>% 
  column_to_rownames("Seq_ID")

Healthy_sp <- Healthy_sp %>% 
  column_to_rownames("Seq_ID")

# SparCC requires the datasets to be in txt format 
# Rows = Species names, Columns = samples

# we need to tranpose the data and remove the Seq_ID column

# MAM
MAM_sp_transposed <- t(MAM_sp)

# MAM_sp_transposed<- is.data.frame(MAM_sp_transposed)

# Add a label for the row names column
MAM_sp_sparcc <- MAM_sp_transposed %>% 
  as.data.frame() %>% 
  rownames_to_column("species") 

# write table
write.table(MAM_sp_sparcc, "data/processed/sparcc_data/MAM_sp_sparcc.txt", sep = "\t", 
            row.names = F, col.names = TRUE, quote = FALSE)



# Healthy

# MAM
Healthy_sp_transposed <- t(Healthy_sp)

# Add a label for the row names column
Healthy_sp_sparcc <- Healthy_sp_transposed %>% 
  as.data.frame() %>% 
  rownames_to_column("species") 

# write table 
write.table(Healthy_sp_sparcc, "data/processed/sparcc_data/Healthy_sp_sparcc.txt", sep = "\t", 
            row.names = F, col.names = TRUE, quote = FALSE)




# Formatting OUTPUT ------------------------------------------------------

## There are two files from SparCC output:
# 1: cor_sparcc.csv - this is the correlation matrix that contains the rho for each pairwise species interaction 

# 2: pvals_one_sided.csv - this contains the pvalues for the rhos between each pairwise species interaction 

# on local PC 
# the data path is: python/SparCC/project/output/ ... mam or healthy 
# p value file is within the pvals folder inside each mam or healthy folder 


# load the raw outputs

# cor matrix --------------------------------------------------------------

# MAM
raw_m_cor_sparcc <- read_csv("python/SparCC/project/output/mam/m_cor_sparcc.csv")

# Healthy
raw_h_cor_sparcc <- read_csv("python/SparCC/project/output/healthy/h_cor_sparcc.csv")

# NOTE: the first column of the dataset is numbering and not needed - this will be removed

# p matix ----------------------------------------------------------------

# MAM
raw_m_pvals_one_sided <- read_csv("python/SparCC/project/output/mam/pvals/m_pvals_one_sided.csv")

# Healthy
raw_h_pvals_one_sided <- read_csv("python/SparCC/project/output/healthy/pvals/h_pvals_one_sided.csv")

# NOTE: the first column of the dataset is numbered and not needed - this will be removed

# FORMATTING DATASETS -----------------------------------------------------

# The output from Sparcc produces a matrix where the row names, and the column names i.e,
# the name of the species is lost. Therefore, we will rename the rows and columns using the 
# RA dataset used as the input for SparCC

# we will start by creating new datasets by renaming the original datasets 
# the purpose of creating new datasets is to preserve the original datasets to re run/check anything if needed

# dataset input: (ra = relative abundance)
mam_ra <- MAM_sp_sparcc  
healthy_ra <- Healthy_sp_sparcc

# outputs to format
# cor matrix
mam_cor <- raw_m_cor_sparcc
healthy_cor <- raw_h_cor_sparcc

# p val matrix 
mam_pvals <- raw_m_pvals_one_sided
healthy_pvals <- raw_h_pvals_one_sided


# We will add row names and column names to the outputs by matching to 
# the species names in the RA datasets 

# Cor - MAM  --------------------------------------------------------------

### MAM 

# Cor Table

# Removing the first column from the mam_cor dataset
# check that the number of rows and columns are now equal (1015 x 1015)
mam_cor <- mam_cor[ , -1]

# add the species columm from the ra data 
mam_cor <- cbind(species = mam_ra$species, mam_cor)

#to add the column names as the species names from the RA dataset, we must first 
# change the column "species" to rownames. This adjusts the columns correctly for 
# addition of column names- species 
mam_cor <- mam_cor %>% 
  column_to_rownames("species")

# change the column names using the ra dataset
colnames(mam_cor) <- mam_ra$species

save(mam_cor, file = "data/processed/mam_cor_matrix.Rdata")

# P vals - MAM  --------------------------------------------------------------

# P vals Table

# Removing the first column from the mam_pvals dataset
# check that the number of rows and columns are now equal (1015 x 1015)
mam_pvals <- mam_pvals[ , -1]

# add the species columm from the ra data 
mam_pvals <- cbind(species = mam_ra$species, mam_pvals)

# to add the column names as the species names from the RA dataset, we must first 
# change the column "species" to rownames. This adjusts the columns correctly for 
# addition of column names- species 
mam_pvals <- mam_pvals %>% 
  column_to_rownames("species")

# change the column names using the ra dataset
colnames(mam_pvals) <- mam_ra$species

save(mam_pvals, file = "data/processed/mam_pvals_matrix.Rdata")

# Cor - Healthy  --------------------------------------------------------------

### HEALTHY

# Cor Table

# Removing the first column from the mam_cor dataset
# check that the number of rows and columns are now equal (1015 x 1015)
healthy_cor <- healthy_cor[ , -1]

# add the species columm from the ra data 
healthy_cor <- cbind(species = healthy_ra$species, healthy_cor)

# to add the column names as the species names from the RA dataset, we must first 
# change the column "species" to rownames. This adjusts the columns correctly for 
# addition of column names- species 
healthy_cor <- healthy_cor %>% 
  column_to_rownames("species")

# change the column names using the ra dataset
colnames(healthy_cor) <- healthy_ra$species

save(healthy_cor, file = "data/processed/healthy_cor_matrix.Rdata")

# P vals - Healthy  --------------------------------------------------------------

# P vals Table

# Removing the first column from the mam_pvals dataset
# check that the number of rows and columns are now equal (1015 x 1015)
healthy_pvals <- healthy_pvals[ , -1]

# add the species columm from the ra data 
healthy_pvals <- cbind(species = healthy_ra$species, healthy_pvals)

# to add the column names as the species names from the RA dataset, we must first 
# change the column "species" to rownames. This adjusts the columns correctly for 
# addition of column names- species 
healthy_pvals <- healthy_pvals %>% 
  column_to_rownames("species")

# change the column names using the ra dataset
colnames(healthy_pvals) <- healthy_ra$species

save(healthy_pvals, file = "data/processed/healthy_pvals_matrix.Rdata")


# RESHAPE MATRIX - TABLE --------------------------------------------------

# Next, we will reshape the matrix to tables 

# we have 4 datasets to reshape
# 2 cor matrix 
# 2 p val matrix 

# data matrix names: 
# mam_cor, healthy_cor
# mam_pvals, healthy_pvals

# the matrix is formatted to a table with two columns Var 1 and Var 2 and value for 
# the rho or pval - the format resembles an edge list 

# MAM ---------------------------------------------------------------------
# Make the rownames into a column

# Cor table 
mam_cor_table <- mam_cor %>% 
  rownames_to_column("Var1")

mam_cor_table <- mam_cor_table %>% 
  gather(key = "Var2", value = "value", -Var1)

save(mam_cor_table, file = "data/processed/mam_cor_table.Rdata")

# P vals table 

mam_pvals_table <- mam_pvals %>% 
  rownames_to_column("Var1")

mam_pvals_table <- mam_pvals_table %>% 
  gather(key = "Var2", value = "value", -Var1)

save(mam_pvals_table, file = "data/processed/mam_pvals_table.Rdata")


# Healthy -----------------------------------------------------------------

# Cor table 

healthy_cor_table <- healthy_cor %>% 
  rownames_to_column("Var1")

healthy_cor_table <- healthy_cor_table %>% 
  gather(key = "Var2", value = "value", -Var1)

save(healthy_cor_table, file = "data/processed/healthy_cor_table.Rdata")

# P vals table 

healthy_pvals_table <- healthy_pvals %>% 
  rownames_to_column("Var1")

healthy_pvals_table <- healthy_pvals_table %>% 
  gather(key = "Var2", value = "value", -Var1)

save(healthy_pvals_table, file = "data/processed/healthy_pvals_table.Rdata")



# SIGNIFICANCE FILTER -----------------------------------------------------
# The coefficients are filtered based on significant p values 
# for this, the p vals table is first filtered to obtain p < 0.05 values 
# this table is matched to the coef table 

# Significant p values 
## MAM
sig_mam_pvals_table <- mam_pvals_table %>% 
  filter(value < 0.05) # n = 754

## Healthy
sig_healthy_pvals_table <- healthy_pvals_table %>% 
  filter(value < 0.05) # n = 864

# Now, using the significant pairwise interactions, filter the coef table 

# Cor values 
## MAM
sig_mam_cor_table <- mam_cor_table %>%
  semi_join(sig_mam_pvals_table, by = c("Var1", "Var2")) 


# manually filtering the dataset for:

# format the tables
# 1: remove self duplicates 
# 2: remove bi directional pairwise

# Filter to keep only one instance of each species pair, regardless of order
format_sig_mam_cor_table <- sig_mam_cor_table %>%
  # Ensure pairs are ordered alphabetically to avoid duplicate (A-B and B-A) scenarios
  mutate(Var1 = pmin(Var1, Var2),
         Var2 = pmax(Var1, Var2)) %>%
  # Remove duplicates after reordering
  distinct(Var1, Var2, value, .keep_all = TRUE) %>% 
  filter(Var1 != Var2)
# n = 377

save(format_sig_mam_cor_table, file = "data/processed/format_sig_mam_cor_table.Rdata")

## Healthy
sig_healthy_cor_table <- healthy_cor_table %>%
  semi_join(sig_healthy_pvals_table, by = c("Var1", "Var2")) 


# format the tables
# 1: remove self duplicates 
# 2: remove bi directional pairwise

# Filter to keep only one instance of each species pair, regardless of order
format_sig_healthy_cor_table <- sig_healthy_cor_table %>%
  # Ensure pairs are ordered alphabetically to avoid duplicate (A-B and B-A) scenarios
  mutate(Var1 = pmin(Var1, Var2),
         Var2 = pmax(Var1, Var2)) %>%
  # Remove duplicates after reordering
  distinct(Var1, Var2, value, .keep_all = TRUE) %>% 
  filter(Var1 != Var2)
# n = 432

save(format_sig_healthy_cor_table, file = "data/processed/format_sig_healthy_cor_table.Rdata")



# Output Analyses  --------------------------------------------------------

summary(format_sig_mam_cor_table$value)
summary(format_sig_healthy_cor_table$value)

range(format_sig_mam_cor_table$value)
range(format_sig_healthy_cor_table$value)

# the correlation coefficient values are too low to be meaningful 


# plots -------------------------------------------------------------------


# Healthy SparCC correlation distribution
healthy_sparcc_plot <- ggplot(format_sig_healthy_cor_table, aes(x = value)) +
  geom_histogram(binwidth = 0.01, fill = "cyan4", color = "black", linewidth = 0.3) +
  labs(x = "SparCC Correlation Coefficient", y = "Frequency") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp1/clean_plots/healthy_sparcc_plot.pdf", 
       plot = healthy_sparcc_plot, units = "cm", width=6, height=8)

# MAM SparCC correlation distribution
mam_sparcc_plot <- ggplot(format_sig_mam_cor_table, aes(x = value)) +
  geom_histogram(binwidth = 0.01, fill = "darkgoldenrod1", color = "black", linewidth = 0.3) +
  labs(x = "SparCC Correlation Coefficient", y = "Frequency") +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp1/clean_plots/mam_sparcc_plot.pdf", 
       plot = mam_sparcc_plot, units = "cm", width=6, height=8)















