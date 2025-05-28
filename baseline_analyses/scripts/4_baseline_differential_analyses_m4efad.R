## code summary 
## - Running MaAsLin2 to detect differentially abundant species and pathways
## - Annotating species with oxygen tolerance (aerobic vs anaerobic)
## - Plotting high-effect species and significant species (Rothia, Strep)
## - Visualizing significant pathway shifts at the community level
## - Running MaAsLin2 on species-stratified pathways
## - Linking differentially abundant pathways to specific species



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
# has a prevalence filter of 10%
load("data/processed/filtered_baseline_species.Rdata")
# pathway abundance 
load("data/processed/pathways.Rdata")
# species stratified pathways 
load("data/processed/pathways_strat.Rdata")


# Maaslin on species abundance  -------------------------------------------

# Rownames need to be in the same order for both metadata and phyla data
baseline_data <- baseline_data[sort(rownames(baseline_data)),]

all(rownames(baseline_data) == rownames(baseline_species)) # check they match up

# Run maaslin2
Maaslin2(input_data = baseline_species, 
         input_metadata = baseline_data,
         output = "outputs/maaslin2/species_baseline", # specify directory to save output files
         fixed_effects = "Condition", # testing variable (can add multiple variables)
         transform = "log", # log transform species counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.10, # species needs to be found in atleast 10% of samples to be included
         cores = 4) # speeds things up


# Read in output files
all_results_species <- read_tsv("outputs/maaslin2/species_baseline/all_results.tsv")
sig_results_species <- read_tsv("outputs/maaslin2/species_baseline/significant_results.tsv")

n_distinct(all_results_species$feature) # 128 species differentially abundant 
n_distinct(sig_results_species$feature) # 2 species significantly differentially abundant

# we can look at the coefficient value to work out whether species increased/decreased in relative abundance
n_distinct(filter(all_results_species, coef >0)) # 60 species increased
n_distinct(filter(all_results_species, coef <0)) # 68 species decreased

# top 10 species - only two incresaed their RA in MAM
head(arrange(sig_results_species, desc(coef)), 10) 

# how many species are reduced in MAM
n_distinct(filter(all_results_species, coef <0)) # 68 species are reduced 



# plots -------------------------------------------------------------------

# Ordered PLOT
all_results_species %>% ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("cyan4", "darkgoldenrod1")) +
  labs(x = "Effect size and direction of change in abundance",
       y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust=0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 6)) +
  guides(fill = FALSE)  # Remove the legend


# Plot - coef > 1 | coef < -1 ---------------------------------------------

# select the species that had an effect size (coef) of > 1 and less than -1
high_coef_sp <- all_results_species %>% 
  filter(coef > 1 | coef < -1) %>% 
  mutate(feature = sub("^s__", "", feature))

# plot
high_coef_sp %>% ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("cyan4", "darkgoldenrod1")) +
  labs(x = "Effect size",
       y = "Species") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, vjust=0.5, hjust=1, size = 10),
        axis.text.y = element_text(size = 10)) +
  guides(fill = FALSE)  # Remove the legend


# assign aerobic/anaerobic to species  ------------------------------------

highsp <- high_coef_sp

highsp <- highsp %>%
  mutate(
    species = gsub("_", " ", feature),   # Replace underscores with spaces in the 'feature' column
    anaerobic = mo_is_anaerobic(species) # Apply the mo_is_anaerobic function on the 'species' column
  )


highsp <- highsp %>%
  mutate(
    o_toll = mo_oxygen_tolerance(species) # Apply the mo_oxygen_tolerance function on the 'species' column
  )


# showing anaerobic - adds a dot to the bars 
diff_species_plot <- highsp %>% 
  ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) + 
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_point(
    aes(shape = factor(anaerobic)),
    position = position_nudge(x = 0),
    size = 4,
    color = "red"
  ) +  # Add point to anaerobic species
  scale_shape_manual(values = c(NA, 16)) +  # Point only for anaerobic species
  scale_fill_manual(values = c("cyan4", "darkgoldenrod1")) +
  labs(
    x = "Effect size",
    y = "Species"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )


# pdf
ggsave(file = "plots/chp1/clean_plots/diff_species_plot.pdf", 
       plot = diff_species_plot, units = "cm", width=12, height=9)



# Plot - significant species ----------------------------------------------

# Plot the RA for the two sig spp
# Rothia + Strep 

# create a list of the significant species of MAM v Healthy at 12 months
mam_12_sig_list <- c("s__Rothia_mucilaginosa",
                     "s__Streptococcus_salivarius")

# We will use the "filtered_baseline_species" dataset

# Filter out the significant species 
spRoth_spStrep <- filtered_baseline_species %>%
  gather(Taxa, RA, -Seq_ID) %>% 
  inner_join(baseline_data, by = "Seq_ID") %>% 
  filter(Taxa %in% mam_12_sig_list) %>% 
  group_by(Subject_ID, Taxa) %>%
  ungroup()%>%
  group_by(Taxa) %>% 
  mutate(Taxa = sub("^s__", "", Taxa))



# Plot
rothia_strep_plot <- ggplot(spRoth_spStrep, aes(x = Taxa, y = RA, fill = Condition)) +
  geom_quasirandom(
    dodge.width = 0.9,
    shape = 21,
    size = 1,
    alpha = 0.7,
    aes(color = Condition)
  ) +
  geom_boxplot(
    width = 0.6,
    position = position_dodge(0.9),
    alpha = 0.5,
    outlier.colour = NA,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  scale_y_log10() +
  labs(
    y = "Relative Abundance (log10 scale)",
    x = NULL
  ) +
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
ggsave(file = "plots/chp1/clean_plots/rothia_strep_plot.pdf", 
       plot = rothia_strep_plot, units = "cm", width=9, height=9)

# stat --------------------------------------------------------------------

# run stat
pval_roth_strep <- spRoth_spStrep%>%
  summarise(p_value = wilcox.test(RA ~ Condition)$p.value)

# Rothia_mucilaginosa      0.00105 
# Streptococcus_salivarius 0.000686


# Maaslin on pathway abundance  -------------------------------------------

# Prepare metadata - check that the baseline_data set has Seq_ID as rows
baseline_data <- baseline_data %>% 
  column_to_rownames("Seq_ID")

# Prepare matching pathways data
baseline_pathways <- filter(pathways, rownames(pathways) %in% rownames(baseline_data))

# Rownames need to be in the same order for both metadata and phyla data
baseline_data <- baseline_data[sort(rownames(baseline_data)),]

baseline_pathways <- baseline_pathways[sort(rownames(baseline_pathways)),]

all(rownames(baseline_data) == rownames(baseline_pathways)) # check they match up

# Run maaslin2
Maaslin2(input_data = baseline_pathways, 
         input_metadata = baseline_data,
         output = "outputs/maaslin2/pathways_baseline", # specify directory to save output files
         fixed_effects = "Condition", # testing variable (can add multiple variables)
         transform = "log", # log transform pathway counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.05, # pathway needs to be found in atleast 5% of samples to be included
         cores = 4) # speeds things up

# Read in output files
all_results_pathways <- read_tsv("outputs/maaslin2/pathways_baseline/all_results.tsv")
sig_results_pathways <- read_tsv("outputs/maaslin2/pathways_baseline/significant_results.tsv")

n_distinct(all_results_pathways$feature) # 467 pathways included in testing
n_distinct(sig_results_pathways$feature) # 5 pathways found to be significantly different between MAM and Healthy

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(sig_results_pathways, coef >0)) # 5 pathways increased
n_distinct(filter(sig_results_pathways, coef <0)) # 0 pathways decreased

# top 10 pathways that increased following FMT
head(arrange(sig_results_pathways, desc(coef)), 10) 




# Plot - significant pathways  -----------------------------------------

## Ordered plot - all the results 

# plot the pathways that include the 5 sig pathways
pathway_plot <- sig_results_pathways %>%
  ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("darkgoldenrod1")) +
  labs(
    x = "Effect Size",
    y = "Community-Level Pathways"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp1/clean_plots/pathway_plot.pdf", 
       plot = pathway_plot, units = "cm", width=15, height=9)


# Pathways associated with species  -----------------------------------------------------

# Plot the significant species and their associated significant pathways
# Compare the community pathway changes and and the cumulative changes 

# First, we need to run Maaslin on the species stratified pathways table
# this is to obtain the differential pathways


# * Maaslin - Pathways Stratified  ------------------------------------------

# we wil use the pathways_strat file as the format us suitable for Maasline

# Prepare metadata
baseline_data <- baseline_data

basep_strat <- pathways_strat %>% 
  column_to_rownames(var = "Seq_ID")

# Prepare matching pathways data
baseline_pathways_strat <- filter(basep_strat, rownames(basep_strat) %in% rownames(baseline_data))

# the two datasets did not match in the number of rows to we filter to match
# basep_strat <- semi_join(basep_strat, baseline_data, by = "Seq_ID")

# convert the column Seq_ID back to a row for data entry into Maaslin,
# an error pops up sayong that there are rownames already existing, so we remove them as follows:

row.names(baseline_data) <- NULL
row.names(basep_strat) <- NULL
baseline_data <- baseline_data %>%  column_to_rownames(var = "Seq_ID")
basep_strat <- basep_strat %>%  column_to_rownames(var = "Seq_ID")

# Rownames need to be in the same order for both metadata and phyla data
baseline_data <- baseline_data[sort(rownames(baseline_data)),]

baseline_pathways_strat <- baseline_pathways_strat[sort(rownames(baseline_pathways_strat)),]

all(rownames(baseline_data) == rownames(baseline_pathways_strat)) # check they match up

# Run maaslin2
Maaslin2(input_data = baseline_pathways_strat, 
         input_metadata = baseline_data,
         output = "outputs/maaslin2/pathways_strat_baseline", # specify directory to save output files
         fixed_effects = "Condition", # testing variable (can add multiple variables)
         #random_effects = "Subject_ID", # have to add given we have multiple samples from the same individual (i.e. non-independent)
         transform = "log", # log transform pathway counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.05, # pathway needs to be found in atleast 5% of samples to be included
         cores = 4) # speeds things up


# Read in output files
all_results_baseline_pathways_strat <- read_tsv("outputs/maaslin2/pathways_strat_baseline/all_results.tsv")
sig_results_baseline_pathways_strat <- read_tsv("outputs/maaslin2/pathways_strat_baseline/significant_results.tsv")

n_distinct(all_results_baseline_pathways_strat$feature) # 7112 pathways included in testing
n_distinct(sig_results_baseline_pathways_strat$feature) # 66 pathways found to be significantly different between MAM and Healthy

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(sig_results_baseline_pathways_strat, coef >0)) # 66 pathways increased
n_distinct(filter(sig_results_baseline_pathways_strat, coef <0)) # 0 pathways decreased

# top 10 pathways that increased following refeed
head(arrange(sig_results_baseline_pathways_strat, desc(coef)), 10) 

# top 10 pathways that decreased
head(arrange(sig_results_baseline_pathways_strat, (coef)), 20) 


# * Associate species to pathways -----------------------------------------

# Data we will use

# species stratfied differential results 
sig_results_baseline_pathways_strat <- read_tsv("outputs/maaslin2/pathways_strat_baseline/significant_results.tsv")

# community level pathways differential 
sig_results_pathways <- read_tsv("outputs/maaslin2/pathways_baseline/significant_results.tsv")

# species differential results
sig_results_species <- read_tsv("outputs/maaslin2/species_baseline/significant_results.tsv")

# Rename the datasets to allow modifications 

diff_paths_sp_strat <- sig_results_baseline_pathways_strat

diff_paths_comm <- sig_results_pathways 

diff_spp <- sig_results_species

### Step 1 - MODIFICATIONS TO DATASETS 
# 1: Differential community pathways dataset

# The diff_paths_sp_strat has the species in the same label as the feature or pathway
# we must separate this into genus and species 

diff_paths_sp_strat <- diff_paths_sp_strat %>%
  separate(feature, into = c("feature", "Genus"), sep = "\\.\\.g__|\\.g__", extra = "merge") %>%
  mutate(Species = ifelse(grepl("\\.unclassified$", feature), "unclassified", Genus))

# to separate the genus column into species
diff_paths_sp_strat <- diff_paths_sp_strat %>% 
  separate(Genus, into = c("Genus", "Species"), sep = "\\.s__", extra = "merge", fill = "right")

# rename the coef column
diff_paths_sp_strat <- diff_paths_sp_strat %>%
  rename(coef_paths_strat = coef) 

# select the columns of interest - here we want feature, Species, coef
diff_paths_sp_strat <- diff_paths_sp_strat %>% 
  select (feature, Species, coef_paths_strat)

# 2: Differential community pathways 

# rename the codef column
diff_paths_comm <- diff_paths_comm %>%
  rename(coef_paths_comm = coef) %>% 
  select (feature,coef_paths_comm)


# 3: Differential species

# rename the codef column
diff_spp <- diff_spp %>%
  rename(coef_spp = coef) %>% 
  select (feature,coef_spp) %>% 
  rename(Species = feature) 

diff_spp$Species <- gsub("^s__", "", diff_spp$Species)

### Step 2 - Filter the datasets 

# Filter the diff_paths_sp_strat by diff_paths_comm 

comm_diff_paths_sp_strat <- diff_paths_sp_strat %>%
  filter(feature %in% diff_paths_comm$feature)


# Filter the diff_strat_paths to only have the diff_spp

spp_diff_paths_sp_strat <- diff_paths_sp_strat %>%
  filter(Species %in% diff_spp$Species)

# this was a test run - came back with one pathway from the community paths that also was sig in strat paths
pp1 <- spp_diff_paths_sp_strat %>%
  filter(feature %in% diff_paths_comm$feature)


#SPECIES WITH THE GREEATEST CHANGE

# this has all the species that were tested 
all_results_species <- read_tsv("outputs/maaslin2/species_baseline/all_results.tsv")


# We have filtred the spp with coef change
# Select species with coeff > 1 and > -1 from all_sig dataset
species_coef <- all_results_species %>%
  filter(coef > 1 | coef < -1)

# make a vector of the 13 species
species_coef_test <- species_coef %>% 
  select(feature, coef) %>% 
  rename (coef_sp_all = coef) %>% 
  rename(Species = feature)

# modify the dataset
species_coef_test$Species <- gsub("^s__", "", species_coef_test$Species)


# Filter the diff_paths_sp_strat for the 13 spp that had highest effect size
high_effect_spp_diff_paths_sp_strat <- diff_paths_sp_strat %>%
  filter(Species %in% species_coef_test$Species)

# INFERENCE: the only species that is si different in MAM and is associated with a sig pathway is 
# Streptococcus salivarius and the pathway is valine biosynthesis, both in the community and when 
# associated with Strep, the pathway increases in MAM














