# code summary 
# The script focuses on bias control of extraction and sequencing 

# Positive controls were included with samples, these are: 
# 1 - Microbial Community Standard (microbial_controls) - to control for extraction bias 
# 2 - Microbial Community DNA Standard (dna_controls) - to control for sequencing bias 

# The ZYMOBIOMICS Microbial Community and Microbial Community DNA standards were used. 
# The manufacturer has provided expected yields/abundance of the defined set of 
# microorganisms. These are compared against the control samples. 


# load Packages -----------------------------------------------------------
library(tidyverse) # required for data wrangling and ggplot
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(scales) # required for scale_y_continuous() which prevents y axis on plot from defaulting to scientific notation
library(dplyr)
library(vegan)
library(readxl)
library(purrr) # used for stats
library(readr)
library(pals) # contains color palettes 
library(Polychrome) # allows creation of colur palettes with seed colors 
library("AMR") # to get microorganism properties we will use a package called AMR


# load data ---------------------------------------------------------------

# this contains 3 files, all_controls, microbial_controls and dna_controls
load("data/processed/controls.Rdata")

# read in the reference standard - both the microbial and the DNA have the same 
# proportion of the reference species - this file was created in the data_processing.R script
load("data/processed/reference_standard.Rdata")

# species dataset for comparison
load("data/processed/species_all_samples.Rdata")


# 1) Microbial Standard ---------------------------------------------------

# the microbial standard is used to test for extraction bias 

# Q: What species are present in the microbial standard? 

# match the Seq_ID of the microbial_standard dataset to the species dataset

# create the Seq_ID column in the species dataset,
# gather the species into a Taxa column,
# inner_join to the microbial_controls dataset
# convert RA into a percentage/proportion for easy comparison
# remove zero RA species

mb_sp <- species_all_samples %>% 
  rownames_to_column("Seq_ID") %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(microbial_controls) %>% 
  mutate(RA = RA * 100,
         Taxa_plot = ifelse(RA > 1, Taxa, "Other (< 1%)")) %>% 
  select(Taxa, RA, Subject_ID) %>% 
  filter(RA > 0)

# Make the table wider and make all NA values 0
wide_mb_sp <- mb_sp %>%
  drop_na() %>% 
  pivot_wider(names_from = Subject_ID, values_from = RA) %>%
  mutate_all(~ifelse(is.na(.), 0, .))

# combine with the reference 
combined_mb_ref <- left_join(wide_mb_sp, reference_standard, by = "Taxa") %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) # Replace all the N/A values with zero for RA profile building

# how many species are present in the microbial standard that is not present in the reference
sum(combined_mb_ref$theoretical_RA == 0)
# There are 23 species present in addition to the 10 in the reference 

# prepare for the plot
# convert to long format for plotting 
long_combined_mb_ref <- reshape2::melt(combined_mb_ref, id.vars = "Taxa")

# Calculate the proportion of RA within each column
long_combined_mb_ref <- long_combined_mb_ref %>%
  group_by(variable) %>%
  mutate(RA_proportion = value / sum(value) * 100)

# create a new column called batch to make the plot more elegant 
long_combined_mb_ref <- long_combined_mb_ref %>%
  mutate(batch = case_when(
    variable == "PC1" ~ "B1",
    variable == "PC1_EH" ~ "B1-2",
    variable == "PC2" ~ "B2",
    variable == "PC3" ~ "B3",
    variable == "PC4" ~ "B4",
    variable == "PC5" ~ "B5",
    variable == "theoretical_RA" ~ "Ref",
    TRUE ~ NA_character_
  )) %>% 
  arrange(batch)


# obtain the proportions of the listed species in each sample from batches 

# Remove rows where theoretical_RA is 0
combined_mb_ref_listed_only <- combined_mb_ref %>%
  filter(theoretical_RA != 0)

# save the file for sample numbers
write_xlsx(combined_mb_ref_listed_only, "data/processed/combined_mb_ref_listed_only.xlsx")

long_combined_mb_ref_listed_only <- reshape2::melt(combined_mb_ref_listed_only, id.vars = "Taxa")

long_combined_mb_ref_listed_only <- long_combined_mb_ref_listed_only %>%
  group_by(variable) %>%
  mutate(RA_proportion = value / sum(value) * 100)


# plot the actual species profile -----------------------------------------

# plot 1  ---------

# plot 1 - Any species present in addition to expected species?
# the grey regions show the proportion of species present that were added on to the 
# microbial standard during the extraction, sequencing and data processing 

# Existing code to create the named vector of colors for the selected taxa
# list of taxa from reference 
selected_taxa <- c(
  "s__Pseudomonas_aeruginosa_group",
  "s__Escherichia_coli",
  "s__Salmonella_enterica",
  "s__Lactobacillus_fermentum",
  "s__Enterococcus_faecalis",
  "s__Staphylococcus_aureus",
  "s__Listeria_monocytogenes",
  "s__Bacillus_subtilis_group",
  "s__Saccharomyces_cerevisiae",
  "s__Cryptococcus_neoformans"
)

pastel_colors <- c(
  "thistle", # light red
  "#FFCC99", # light orange
  "#FFFF99", # light yellow
  "#99FF99", # light green
  "#99CCFF", # light blue
  "#CC99FF", # light purple
  "#FF99CC", # light pink
  "#99FFFF", # light cyan
  "lightsalmon", # light lemon
  "#FFCCFF"  # light magenta
)

taxa_colors <- setNames(pastel_colors, selected_taxa)


# Create a new ggplot with the specific and default colors
mb_species_expected_plot <- ggplot(long_combined_mb_ref, aes(x = batch, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Expected Species",
       x = "Microbial Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = taxa_colors, na.value = "azure3") +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )


# pdf
ggsave(file = "plots/qc/clean_plots/mb_species_expected_plot.pdf", 
plot = mb_species_expected_plot, units = "cm", width=15, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/mb_species_expected_plot.svg", 
 plot = mb_species_expected_plot, units = "cm", width=15, height=10)


# plot of listed species only 
ggplot(long_combined_mb_ref_listed_only, aes(x = variable, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Expected Species",
       x = "Microbial Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  #scale_fill_manual(values = taxa_colors, na.value = "azure3") +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )






# plot 2  -----------------------------------------------------------------

# plot 2 - What are the additional species?
# plot the whole species profile without adding a grey filter

# using the Polychrome package we create a colur palette of 40 colurs using seed colors 
# Define the seed colors
seed_colors <- c("#CC99FF", "#FF99CC", "#99FFFF", "#FFCC99")

# Create a 40-color palette based on the seed colors
P40 <- createPalette(40, seed_colors)

# Display the palette using the swatch function
swatch(P40)

# Assign grey to selected_taxa and colors to others
taxa_colors_2 <- setNames(ifelse(unique(long_combined_mb_ref$Taxa) %in% selected_taxa, 
                                 "grey", P40[1:length(unique(long_combined_mb_ref$Taxa))]), 
                          unique(long_combined_mb_ref$Taxa))

# Create the ggplot with the specific and default colors
mb_species_additional_plot <- ggplot(long_combined_mb_ref, aes(x = batch, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Observed Species",
       x = "Microbial Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = taxa_colors_2, na.value = "azure3") +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )

# pdf
ggsave(file = "plots/qc/clean_plots/mb_species_additional_plot.pdf", 
      plot = mb_species_additional_plot, units = "cm", width=23, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/mb_species_additional_plot.svg", 
      plot = mb_species_additional_plot, units = "cm", width=23, height=10)

# types of unlisted species 
unexp <- combined_mb_ref %>%
  filter(theoretical_RA == 0) %>% 
  mutate(Taxa = gsub("^s__", "", Taxa)) %>% 
  select(Taxa) %>% 
  mutate(type = case_when(
    grepl("phage", Taxa, ignore.case = TRUE) ~ "Phage",
    grepl("virus", Taxa, ignore.case = TRUE) ~ "Virus",
    TRUE ~ "Other"
  ))

# Calculate the proportion of each type
type_proportions <- unexp %>%
  group_by(type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count) * 100)



# Mean RA of each species  ------------------------------------------------

# Calculate the mean relative abundance for each species
mean_ra_per_species <- wide_mb_sp %>%
  rowwise() %>%
  mutate(mean_RA = mean(c_across(-Taxa), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Taxa, mean_RA)

# check that the means total to 100 as RA = 1 (100%)
total_mean_RA <- sum(as.numeric(mean_ra_per_species$mean_RA))

# convert the numbers to 3 didgits for readability 
mean_ra_per_species <- mean_ra_per_species %>%
  mutate(mean_RA = formatC(mean_RA, format = "f", digits = 3))

# Combine the theoretical_RA column from combined_mb_ref with mean_ra_per_species by Taxa
mean_ra_per_species <- mean_ra_per_species %>%
  left_join(combined_mb_ref %>% select(Taxa, theoretical_RA), by = "Taxa")  %>% 
  mutate(Taxa = gsub("^s__", "", Taxa))


# Gramness bias  ----------------------------------------------------------

# The gramness of the detected species is relevant for the microbial standard as this is 
# extracted along with other samples. Therefore, we will use a package called AMR to assign
# Gram-ness and check if there is any bias towards a particular Gram type 

# create a new dataset for ease
gram <- mean_ra_per_species %>%
  mutate(gram_stain = mo_gramstain(Taxa))


# prepare plot 3 ------------------------------------------------------------------
# reassign a new name to the Taxa column as we need to ammend the string names 
gram_taxa <- gram %>% 
  mutate(Taxa = paste0("s__", Taxa))

# to plot the grammness we need to have the gram function run on the ling format, 
# therefore, we will run the fram function again. on the long format data 

gram_plot <- long_combined_mb_ref %>% 
  mutate(sp = gsub("^s__", "", Taxa))

# run the gram stain 
gram_plot <- gram_plot %>% 
  mutate(gram_stain = mo_gramstain(sp))


# we need to make some adjustments to the gram_plot dataset
# we will separate the phages so that they are not confused with the gram stain values 
gram_plot <- gram_plot %>%
  mutate(type = case_when(
    grepl("phage", Taxa, ignore.case = TRUE) | grepl("virus", Taxa, ignore.case = TRUE) ~ "phage/virus",
    TRUE ~ "bacteria"
  ))

# create a new column to be used in the plot 
gram_plot <- gram_plot %>%
  mutate(gram = ifelse(type == "bacteria", gram_stain, "phage/virus"))

# we will fill in some entries manually that didnt have gram stains identified
gram_plot <- gram_plot %>%
  mutate(gram = case_when(
    Taxa == "s__Cryptococcus_neoformans" ~ "Eukaryote-no_stain",
    Taxa == "s__Saccharomyces_cerevisiae" ~ "Eukaryote-no_stain",
    type == "bacteria" ~ gram_stain,
    type == "phage/virus" ~ "phage/virus",
    TRUE ~ gram  # Preserve existing values for other entries
  ))


# to calculte proportions, first remove all entries with zero RA values 
gram_plot2 <- gram_plot %>% 
  filter(value!= 0)


# Reclassify batch column: Combine all "B*" entries into "Samples" and keep "Ref" as "Reference"
gram_plot_cleaned <- gram_plot2 %>%
  mutate(batch_group = ifelse(grepl("^B", batch), "Samples", 
                              ifelse(batch == "Ref", "Reference", batch)))

# Count the total number of samples in each batch_group
sample_counts <- gram_plot_cleaned %>%
  group_by(batch_group) %>%
  summarise(total_samples = n_distinct(Taxa), .groups = "drop")  # Count unique samples

# Calculate proportions within each batch_group and include total counts
gram_proportions <- gram_plot_cleaned %>%
  group_by(batch_group, gram) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(batch_group) %>%
  mutate(total_count = sum(count),  # Total count within each batch_group
         proportion = (count / total_count) * 100) %>%  # Calculate proportion
  ungroup()

# save the file for sample numbers
write_xlsx(gram_proportions, "data/processed/gram_proportions.xlsx")


# Testing the proportions when phages are removed
# Remove "phage/virus" entries and recalculate proportions
gram_proportions_phage_removed <- gram_plot_cleaned %>%
  filter(gram != "phage/virus") %>%  # Exclude phage/virus entries
  group_by(batch_group, gram) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(batch_group) %>%
  mutate(total_count = sum(count),  # Total count within each batch_group after removing phage/virus
         proportion = (count / total_count) * 100) %>%  # Recalculate proportion
  ungroup()

# Calculate proportions within each batch - this is to compare proportions between batches 
# Remove "Ref" entries and calculate proportions within each batch
gram_proportions_batch <- gram_plot2 %>%
  filter(batch != "Ref") %>%  # Exclude Ref batch
  group_by(batch, gram) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(batch) %>%
  mutate(total_count = sum(count),  # Total count within each batch after removing Ref
         proportion = (count / total_count) * 100) %>%  # Recalculate proportion
  ungroup()


# Remove "Ref" batch and filter out unwanted gram categories
gram_proportions_batch_only <- gram_plot2 %>%
  filter(batch != "Ref", !gram %in% c("phage/virus", "Eukaryote-no_stain")) %>%  # Exclude specified categories
  group_by(batch, gram) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(batch) %>%
  mutate(total_count = sum(count),  # Total count within each batch after filtering
         proportion = (count / total_count) * 100) %>%  # Recalculate proportion
  ungroup()


# Create a grouped bar plot - contains phages and eukaryotes 
ggplot(gram_proportions_batch, aes(x = gram, y = proportion, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +  # Position dodge for side-by-side bars
  labs(title = "Proportion of Gram Categories Across Batches",
       x = "Gram Category",
       y = "Proportion (%)",
       fill = "Batch") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels for readability

# Define a custom color palette
custom_colors <- c("B1" = "#E69F00", "B2" = "#56B4E9", "B3" = "#009E73", 
                   "B4" = "#F0E442", "B5" = "#D55E00", "B6" = "#CC79A7", "B1-2" = "magenta4")


# Create the grouped bar plot with the new color scheme - no phages and aukaryotes 
between_batches_plot <- ggplot(gram_proportions_batch_only, aes(x = gram, y = proportion, fill = batch)) +
  geom_bar(stat = "identity", position = "dodge") +  # Side-by-side bars
  labs(x = "",
       y = "RA Proportion (%)",
       fill = "Batch") +
  theme_bw() +  # White background theme for consistency
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),  # Centered x-axis labels
        panel.border = element_blank(),  # Removes border
        panel.grid = element_blank(),  # Removes grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        axis.text = element_text(size = 9),  # Adjusted text size for labels
        axis.title = element_text(size = 9),  # Adjusted text size for axis titles
        legend.text = element_text(size = 9),  # Adjusted legend text size
        legend.title = element_text(size = 10)  # Adjusted legend title size
  ) + scale_fill_manual(values = custom_colors)  # Apply custom colors

# pdf
ggsave(file = "plots/qc/clean_plots/between_batches_plot.pdf", 
       plot = between_batches_plot, units = "cm", width=10, height=9)

# svg
ggsave(file = "plots/qc/clean_plots/between_batches_plot.svg", 
       plot = between_batches_plot, units = "cm", width=10, height=9)

# plot 3 ------------------------------------------------------------------

# gram stain plot 

gram_stain_plot <- ggplot(gram_plot, aes(x = batch, y = RA_proportion, fill = gram)) +
  geom_bar(stat = "identity") +
  labs(title = "Gram Stain for Microbial Standard",
       x = "Microbial Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Gram stain") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = c("Gram-negative" = "hotpink1", "Gram-positive" = "mediumpurple3", 
                               "Eukaryote-no_stain" = "seagreen", "phage/virus" = "navajowhite2")) +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )

# pdf
 ggsave(file = "plots/qc/clean_plots/gram_stain_plot.pdf", 
      plot = gram_stain_plot, units = "cm", width=12.8, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/gram_stain_plot.svg", 
      plot = gram_stain_plot, units = "cm", width=12.8, height=10)



# difference in RA between samples v ref ----------------------------------

# We want to know what the difference in RA is between the samples overall and the reference
# for this, we will group the species into their respective gram stain groups and look at the 
# difference in their relative abundance 

# filter for theoretical ref != 0, ie, listed samples only 

listed_sp <- gram %>% 
  filter(theoretical_RA != 0) %>% 
  rename(ref_RA = theoretical_RA) 

# Convert mean_RA and ref_RA to numeric if they are not already
listed_sp <- listed_sp %>%
  mutate(mean_RA = as.numeric(mean_RA),
         ref_RA = as.numeric(ref_RA))

# Calculate the percentage change and add a new column
listed_sp <- listed_sp %>%
  mutate(percentage_change_RA = ((mean_RA - ref_RA) / ref_RA) * 100)

# there seems to be a bias in the extraction method where gram-negative bacteria are favored



# statistical difference in species between batches  ----------------------

# Is there a statistical difference in the species abundance between batches?

## There is only one sample per species per batch and this is insufficient to perform a stat test on



# 2) DNA Standard  --------------------------------------------------------

# Q: What species are present in the DNA standard? 
# as the DNA standard is pure DNA and no extraction process occured, we would
# expect to only see the 10 species that is listed in the reference standard 

# match the Seq_ID of the dna_standard dataset to the species dataset

# create the Seq_ID column in the species dataset,
# gather the species into a Taxa column,
# inner_join to the dna_controls dataset
# convert RA into a percentage/proportion for easy comparison
# remove zero RA species

dna_sp <- species_all_samples %>% 
  rownames_to_column("Seq_ID") %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(dna_controls) %>% 
  mutate(RA = RA * 100,
         Taxa_plot = ifelse(RA > 1, Taxa, "Other (< 1%)")) %>% 
  select(Taxa, RA, Subject_ID) %>% 
  filter(RA > 0)

# Make the table wider and make all NA values 0
wide_dna_sp <- dna_sp %>%
  drop_na() %>% 
  pivot_wider(names_from = Subject_ID, values_from = RA) %>%
  mutate_all(~ifelse(is.na(.), 0, .))


# combine with the reference 
combined_dna_ref <- left_join(wide_dna_sp, reference_standard, by = "Taxa") %>% 
  mutate_all(~ifelse(is.na(.), 0, .)) # Replace all the N/A values with zero for %RA profile building

# how many species are present in the microbial standard that is not present in the reference
sum(combined_dna_ref$theoretical_RA == 0)
# There are 14 species present in addition to the 10 in the reference 

# prepare for the plot
# convert to long format for plotting 
long_combined_dna_ref <- reshape2::melt(combined_dna_ref, id.vars = "Taxa")

# Calculate the proportion of RA within each column
long_combined_dna_ref <- long_combined_dna_ref %>%
  group_by(variable) %>%
  mutate(RA_proportion = value / sum(value) * 100)

# create a new column called batch to make the plot more elegant 
long_combined_dna_ref <- long_combined_dna_ref %>%
  mutate(batch = case_when(
    variable == "D1" ~ "B1",
    variable == "D2" ~ "B2",
    variable == "D3" ~ "B3",
    variable == "D4" ~ "B4",
    variable == "D5_1" ~ "B5-1",
    variable == "D5_2" ~ "B5-2",
    variable == "theoretical_RA" ~ "Ref",
    TRUE ~ NA_character_
  )) %>% 
  arrange(batch)


# obtain the proportions of the listed species in each sample from batches 

# Remove rows where theoretical_RA is 0
combined_dna_ref_listed_only <- combined_dna_ref %>%
  filter(theoretical_RA != 0)

# save the file for sample numbers
write_xlsx(combined_dna_ref_listed_only, "data/processed/combined_dna_ref_listed_only.xlsx")

long_combined_dna_ref_listed_only <- reshape2::melt(combined_dna_ref_listed_only, id.vars = "Taxa")

long_combined_dna_ref_listed_only <- long_combined_dna_ref_listed_only %>%
  group_by(variable) %>%
  mutate(RA_proportion = value / sum(value) * 100)


# plot the actual species profile -----------------------------------------


# plot 1  ---------

# plot 1 - Any species present in addition to expected species?
# the grey regions show the proprtion of species present that were added on to the 
# microbial standard during the extraction, sequencing and data processing 

# Existing code to create the named vector of colors for the selected taxa
selected_taxa <- c(
  "s__Pseudomonas_aeruginosa_group",
  "s__Escherichia_coli",
  "s__Salmonella_enterica",
  "s__Lactobacillus_fermentum",
  "s__Enterococcus_faecalis",
  "s__Staphylococcus_aureus",
  "s__Listeria_monocytogenes",
  "s__Bacillus_subtilis_group",
  "s__Saccharomyces_cerevisiae",
  "s__Cryptococcus_neoformans"
)

pastel_colors <- c(
  "thistle", # light red
  "#FFCC99", # light orange
  "#FFFF99", # light yellow
  "#99FF99", # light green
  "#99CCFF", # light blue
  "#CC99FF", # light purple
  "#FF99CC", # light pink
  "#99FFFF", # light cyan
  "lightsalmon", # light lemon
  "#FFCCFF"  # light magenta
)

taxa_colors <- setNames(pastel_colors, selected_taxa)

# Create a new ggplot with the specific and default colors
dna_species_expected_plot <- ggplot(long_combined_dna_ref, aes(x = batch, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Expected Species",
       x = "DNA Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = taxa_colors, na.value = "azure3") +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )


# pdf
ggsave(file = "plots/qc/clean_plots/dna_species_expected_plot.pdf", 
plot = dna_species_expected_plot, units = "cm", width=15, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/dna_species_expected_plot.svg", 
 plot = dna_species_expected_plot, units = "cm", width=15, height=10)



# plot of listed species only 
ggplot(long_combined_dna_ref_listed_only, aes(x = variable, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Expected Species",
       x = "DNA Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
 # scale_fill_manual(values = taxa_colors, na.value = "azure3") +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )


# plot 2  -----------------------------------------------------------------

# plot 2 - What are the additional species?
# for, this we plot the whole species profile without adding a grey filter

# using the Polychrome package we create a colur palette of 40 colurs using seed colors 
# Define the seed colors
seed_colors <- c("#CC99FF", "#FF99CC", "#99FFFF", "#FFCC99")

# Create a 40-color palette based on the seed colors
P40 <- createPalette(40, seed_colors)

# Display the palette using the swatch function
swatch(P40)


# Create a named vector of colors for the taxa

# Assign grey to selected_taxa and colors to others
taxa_colors_2 <- setNames(ifelse(unique(long_combined_dna_ref$Taxa) %in% selected_taxa, 
                                 "grey", P40[1:length(unique(long_combined_dna_ref$Taxa))]), 
                          unique(long_combined_dna_ref$Taxa))

# Create the ggplot with the specific and default colors
dna_species_additional_plot <- ggplot(long_combined_dna_ref, aes(x = batch, y = RA_proportion, fill = Taxa)) +
  geom_bar(stat = "identity") +
  labs(title = "Observed Species",
       x = "DNA Standard in each batch",
       y = "RA Proportion (%)",
       fill = "Species") +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  scale_fill_manual(values = taxa_colors_2, na.value = "azure3") +
  theme_bw() +
  theme_bw() +
  theme(panel.border = element_blank(),  # Removes the border around the plot
        panel.grid = element_blank(),  # Removes all grid lines
        axis.line = element_line(),  # Keeps only x and y axis lines
        plot.title = element_text(),  # Added text size for plot title
        axis.text = element_text(),  # Edited text size for axis labels, the x and y 
        axis.title = element_text(),  # Edited text size for axis titles
        legend.text = element_text(size = 10),  # Edited text size for legend text
        legend.title = element_text()  # Edited text size for legend title
  )

# pdf
ggsave(file = "plots/qc/clean_plots/dna_species_additional_plot.pdf", 
       plot = dna_species_additional_plot, units = "cm", width=23, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/dna_species_additional_plot.svg", 
       plot = dna_species_additional_plot, units = "cm", width=23, height=10)


# types of unlisted species 
unexp_dna <- combined_dna_ref %>%
  filter(theoretical_RA == 0) %>% 
  mutate(Taxa = gsub("^s__", "", Taxa)) %>% 
  select(Taxa) %>% 
  mutate(type = case_when(
    grepl("phage", Taxa, ignore.case = TRUE) ~ "Phage",
    grepl("virus", Taxa, ignore.case = TRUE) ~ "Virus",
    TRUE ~ "Other"
  ))

# Calculate the proportion of each type
type_proportions <- unexp_dna %>%
  group_by(type) %>%
  summarise(count = n()) %>%
  mutate(proportion = count / sum(count) * 100)


# Mean RA of each species  ------------------------------------------------

# Calculate the mean relative abundance for each species
mean_ra_per_species_dna <- wide_dna_sp %>%
  rowwise() %>%
  mutate(mean_RA = mean(c_across(-Taxa), na.rm = TRUE)) %>%
  ungroup() %>%
  select(Taxa, mean_RA)

# check that the means total to 100 as RA = 1 (100%)
total_mean_RA <- sum(as.numeric(mean_ra_per_species_dna$mean_RA))

# convert the numbers to 3 didgits for readability 
mean_ra_per_species_dna <- mean_ra_per_species_dna %>%
  mutate(mean_RA = formatC(mean_RA, format = "f", digits = 3))

# Combine the theoretical_RA column from combined_mb_ref with mean_ra_per_species by Taxa
mean_ra_per_species_dna <- mean_ra_per_species_dna %>%
  left_join(combined_dna_ref %>% select(Taxa, theoretical_RA), by = "Taxa")  %>% 
  mutate(Taxa = gsub("^s__", "", Taxa))


# difference in RA between samples v ref ----------------------------------

# We want to know what the difference in RA is between the samples overall and the reference
# for this, we will look at the difference in their relative abundance 

# filter for theoretical ref != 0, ie, listed samples only 

listed_sp_dna <- mean_ra_per_species_dna %>% 
  filter(theoretical_RA != 0) %>% 
  rename(ref_RA = theoretical_RA) 

# Convert mean_RA and ref_RA to numeric if they are not already
listed_sp_dna <- listed_sp_dna %>%
  mutate(mean_RA = as.numeric(mean_RA),
         ref_RA = as.numeric(ref_RA))

# Calculate the percentage change and add a new column
listed_sp_dna <- listed_sp_dna %>%
  mutate(percentage_change_RA = ((mean_RA - ref_RA) / ref_RA) * 100)

# all the dna standard species are reduced in comparison to the reference values 



# 3) RA of all the batches  -----------------------------------------------

# This section generates a table to observe the RA of species across batches, i,e.,
# comparing the RA of any particular species between batches 

# mb ----------------------------------------------------------------------

batches_mb_ref <- combined_mb_ref %>%
  mutate(across(-Taxa, ~ as.numeric(formatC(as.numeric(.), format = "f", digits = 3))))

# save as csv for Supplememtary data
write.csv(batches_mb_ref, "outputs/tables/batches_mb_ref.csv", row.names = FALSE)

# dna ---------------------------------------------------------------------

# Convert the values within combined_dna_ref to numeric values with 4 digits
batches_dna_ref <- combined_dna_ref %>%
  mutate(across(-Taxa, ~ as.numeric(formatC(as.numeric(.), format = "f", digits = 3))))

# save as csv for Supplememtary data
write.csv(batches_dna_ref, "outputs/tables/batches_dna_ref.csv", row.names = FALSE)



# 4) Prevalence of all taxa  ----------------------------------------------

# Function to calculate prevalence
calculate_prevalence <- function(df) {
  sample_columns <- setdiff(names(df), c("Taxa", "theoretical_RA"))  # Identify sample ID columns
  df$prevalence <- (rowSums(!is.na(df[, sample_columns]) & df[, sample_columns] > 0) / length(sample_columns)) * 100
  return(df)
}


# Apply the function to both datasets
combined_mb_ref_prev <- calculate_prevalence(combined_mb_ref)
combined_dna_ref_prev <- calculate_prevalence(combined_dna_ref)

# modify data for plotting

# Function to classify organism type
classify_organism <- function(df) {
  df$organism <- ifelse(grepl("virus", df$Taxa, ignore.case = TRUE), "virus",
                        ifelse(grepl("phage", df$Taxa, ignore.case = TRUE), "phage", "bacteria"))
  
  df$Taxa <- gsub("^s__", "", df$Taxa)  # Remove "s__" from Taxa names
  
  return(df)
}


# Apply function to dataset
combined_mb_ref_prev <- classify_organism(combined_mb_ref_prev)
combined_dna_ref_prev <- classify_organism(combined_dna_ref_prev)

# plot

# Micrbial control
prevalence_plot_mb <- ggplot(combined_mb_ref_prev, aes(x = prevalence, y = reorder(Taxa, prevalence), color = organism)) +
  geom_point(size = 3) +  # Points instead of bars
  labs(x = "Prevalence (%)",
       y = "",
       color = "Organism") +
  scale_color_manual(values = c("bacteria" = "azure3", "virus" = "red2", "phage" = "red2")) +  # Custom colors
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),  # Removes all grid lines (both major and minor)
    axis.line = element_line(),
    plot.title = element_text(),
    axis.text = element_text(),
    axis.title = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# pdf
ggsave(file = "plots/qc/clean_plots/prevalence_plot_mb.pdf", 
       plot = prevalence_plot_mb, units = "cm", width=14, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/prevalence_plot_mb.svg", 
       plot = prevalence_plot_mb, units = "cm", width=14, height=10)

# DNA control
prevalence_plot_dna <- ggplot(combined_dna_ref_prev, aes(x = prevalence, y = reorder(Taxa, prevalence), color = organism)) +
  geom_point(size = 3) +  # Points instead of bars
  labs(x = "Prevalence (%)",
       y = "",
       color = "Organism") +
  scale_color_manual(values = c("bacteria" = "azure3", "virus" = "red2", "phage" = "red2")) +  # Custom colors
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),  # Removes all grid lines (both major and minor)
    axis.line = element_line(),
    plot.title = element_text(),
    axis.text = element_text(),
    axis.title = element_text(),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 12)
  )

# pdf
ggsave(file = "plots/qc/clean_plots/prevalence_plot_dna.pdf", 
       plot = prevalence_plot_dna, units = "cm", width=14, height=10)

# svg
ggsave(file = "plots/qc/clean_plots/prevalence_plot_dna.svg", 
       plot = prevalence_plot_dna, units = "cm", width=14, height=10)
