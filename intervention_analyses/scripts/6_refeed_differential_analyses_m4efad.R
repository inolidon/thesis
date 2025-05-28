# This script contains the results of the Maaslin run and the relevant plots 
# Load Packages  ----------------------------------------------------------
library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(DescTools)
library(lme4)
library(lmerTest)
library(circlize)
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(vegan)
library(lme4)
library(broom.mixed)
library(knitr)
library(car)
library(AMR)


# 12 v 15 -----------------------------------------------------------------


# Species  ----------------------------------------------------------------


# load data ---------------------------------------------------------------

# the script post_refeed.R script in the following path contains all the Maaslin runs 
# path to script: ~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/scripts

# path to output files: ~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2

# load data from the output path

# Read in output files
all_results_species_refeed <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/species_refeed/all_results.tsv")
sig_results_species_refeed <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/species_refeed/significant_results.tsv")

n_distinct(all_results_species_refeed$feature) # 139 species included in testing
n_distinct(sig_results_species_refeed$feature) # 88 species found to be significantly different 

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(sig_results_species_refeed, coef >0)) # 62 species increased
n_distinct(filter(sig_results_species_refeed, coef <0)) # 26 species decreased

# top 10 species that increased following re feed
head(arrange(sig_results_species_refeed, desc(coef)), 10) 

# how many species are reduced post refeed
n_distinct(filter(sig_results_species_refeed, coef <0)) # 26 species are reduced 


#  Plot Species - differential abundance in MAM ---------------------------

sig_results_species_refeed <- sig_results_species_refeed %>%
  mutate(feature = gsub("^s__", "", feature),                 # remove s__ prefix
         feature = sub("_", " ", feature))                    # replace first _ with a space


# Differentially abundant species in 15 months (increase in 15 magenta)
ordered_plot <- sig_results_species_refeed %>%
  ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("darkgoldenrod1", "darkmagenta")) +
  labs(
    x = "Effect Size",
    y = "Species"
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
ggsave(file = "plots/chp2/clean_plots/ordered_plot.pdf", 
       plot = ordered_plot, units = "cm", width=15, height=20)


# identify Anaerobic species  ---------------------------------------------

# prepare the data 
sp_refeed <- sig_results_species_refeed %>%
  mutate(
    species = gsub("_", " ", feature),   # Replace underscores with spaces in the 'feature' column
    anaerobic = mo_is_anaerobic(species) # Apply the mo_is_anaerobic function on the 'species' column
  )

sp_refeed <- sp_refeed %>%
  mutate(
    o_toll = mo_oxygen_tolerance(species) # Apply the mo_oxygen_tolerance function on the 'species' column
  )




# grA v grB ---------------------------------------------------------------

# laod the data from the previous Maaslin run 
# use correct path: ~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2

# Read in output files
all_results_species_refeed_group <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/species_refeed_group/all_results.tsv")
sig_results_species_refeed_group <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/species_refeed_group/significant_results.tsv")

n_distinct(all_results_species_refeed_group$feature) # 148 species included in testing
n_distinct(sig_results_species_refeed_group$feature) # 6 species found to be significantly different 

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(sig_results_species_refeed_group, coef >0)) # 3 species increased
n_distinct(filter(sig_results_species_refeed_group, coef <0)) # 3 species decreased

# top 10 species that increased following re feed
head(arrange(all_results_species_refeed_group, desc(coef)), 10) 

# how many species are reduced in group 2
n_distinct(filter(all_results_species_refeed_group, coef <0)) # 81 species are reduced 



# Plot  -------------------------------------------------------------------
sig_results_species_refeed_group <- sig_results_species_refeed_group %>%
  mutate(feature = gsub("^s__", "", feature),                 # remove s__ prefix
         feature = sub("_", " ", feature))                    # replace first _ with a space


sig_results_species_refeed_group_plot <- sig_results_species_refeed_group %>%
  ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("azure3", "black")) +
  labs(
    x = "Effect Size",
    y = "Species"
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
ggsave(file = "plots/chp2/clean_plots/sig_results_species_refeed_group_plot.pdf", 
       plot = sig_results_species_refeed_group_plot, units = "cm", width=15, height=5)


# R v nR ------------------------------------------------------------------

# There are no differentially abundant species in the recovered group compared to the unrecovered
# There is no plot for this section 

# Prepare metadata
recovery_data <- recovery_data %>% 
  filter(Recovery %in% c("FALSE", "TRUE")) %>% 
  column_to_rownames("Seq_ID")

# Prepare matching species data
recovery_species <- filter(species, rownames(species) %in% rownames(recovery_data))

# there is an additonal data point in the recovery_data set, therefore, we remove it to match
recovery_data <- recovery_data %>% rownames_to_column("Seq_ID")
recovery_species <- recovery_species %>% rownames_to_column("Seq_ID")
recovery_data <- semi_join(recovery_data, recovery_species, by = "Seq_ID")

# convert the column Seq_ID back to a row for data entry into Maaslin
recovery_data <- recovery_data %>%  column_to_rownames("Seq_ID")
recovery_species <- recovery_species %>%  column_to_rownames("Seq_ID")

# Rownames need to be in the same order for both metadata and phyla data
recovery_data <- recovery_data[sort(rownames(recovery_data)),]

recovery_species <- recovery_species[sort(rownames(recovery_species)),]
save(recovery_species, file = "data/recovery_species.Rdata")

all(rownames(recovery_data) == rownames(recovery_species)) # check they match up

# Run maaslin2
Maaslin2(input_data = recovery_species, 
         input_metadata = recovery_data,
         output = "data/maaslin2/species_recovery", # specify directory to save output files
         fixed_effects = c("Recovery"), # testing variable (can add multiple variables)
         #random_effects = "Feed_type", # have to add given we have multiple samples from the same individual (i.e. non-independent)
         transform = "clr", # log transform pathway counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.1, # species needs to be found in atleast 10%%% of samples to be included
         cores = 4) # speeds things up



# Read in output files
testall_results_species_recovery <- read_tsv("data/maaslin2/species_recovery/all_results.tsv")
testsig_results_species_recovery <- read_tsv("data/maaslin2/species_recovery/significant_results.tsv")

n_distinct(all_results_species_recovery$feature) # 150 species included in testing
n_distinct(sig_results_species_recovery$feature) # 0 species found to be significantly different 

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(all_results_species_recovery, coef >0)) # 79 species increased
n_distinct(filter(all_results_species_recovery, coef <0)) # 71 species decreased

# top 10 species that increased in Group 2
head(arrange(all_results_species_recovery, desc(coef)), 10) 

# how many species are reduced in Group 2
n_distinct(filter(all_results_species_recovery, coef <0)) # 71 species are reduced 



# Baseline differential species comparison --------------------------------

# Plot for chosen species 
# create a list of the species of interest
sp_interest_list <- c("s__Rothia_mucilaginosa",
                      "s__Streptococcus_salivarius")

sp_interest <- long_spp_MAM_whole %>%
  filter(Taxa %in% sp_interest_list) %>%
  ungroup() %>%
  mutate(Taxa = gsub("^s__", "", Taxa),
         Taxa = sub("_", " ", Taxa)) %>%
  group_by(Taxa)


# 1: 12 v 15 -------------------------------------------------------------

# for the difference in age between recovery & refeed group
sig_sp_in_12v15 <- ggplot(sp_interest, aes(x = Taxa, y = RA, fill = Age_months, color = Age_months)) +
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
  scale_y_log10() +
  labs(
    y = "Relative Abundance (log10 scale)",
    x = NULL
  ) +
  facet_grid(Recovery ~ Group) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp2/clean_plots/sig_sp_in_12v15t.pdf", 
       plot = sig_sp_in_12v15, units = "cm", width=9, height=9)

# 2: R v nR ---------------------------------------------------------------

# for difference in recovery between age & group
sig_sp_in_RvnR <- ggplot(sp_interest, aes(x = Taxa, y = RA, fill = Recovery, color = Recovery)) +
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
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_color_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_y_log10() +
  labs(
    y = "Relative Abundance (log10 scale)",
    x = NULL
  ) +
  facet_grid(Group ~ Age_months) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp2/clean_plots/sig_sp_in_RvnR.pdf", 
       plot = sig_sp_in_RvnR, units = "cm", width=9, height=9)

# 3: grA v grB ------------------------------------------------------------

# for difference in group between age & recovery
sig_sp_in_Gr1v2 <- ggplot(sp_interest, aes(x = Taxa, y = RA, fill = Group, color = Group)) +
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
  scale_fill_manual(values = c("1" = "azure3", "2" = "black")) +
  scale_color_manual(values = c("1" = "azure3", "2" = "black")) +
  scale_y_log10() +
  labs(
    y = "Relative Abundance (log10 scale)",
    x = NULL
  ) +
  facet_grid(Recovery ~ Age_months) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp2/clean_plots/sig_sp_in_Gr1v2.pdf", 
       plot = sig_sp_in_Gr1v2, units = "cm", width=9, height=9)

# Wilcoxon Tests  ---------------------------------------------------------

## Prepare datasets 
# Main dataset used is spp_interest 
# Divide it into 12 months and 15 months 

### 12 months 
sp_12mo <- sp_interest %>% 
  filter(Age_months == "12")

### 15 months 
unedit_sp_15mo <- sp_interest %>% 
  filter(Age_months == "15") #this has 880 datapoints, we further editing required

# 15mo datasets have duplicated rows 
# these will be removed before proceeding

duplicated_subjects_15mo <- unedit_sp_15mo %>%
  filter(duplicated(Subject_ID) | duplicated(Subject_ID, fromLast = TRUE))

# this dataset contains one value for Rothia and one for Strep for each Subject ID
unedit_15mo_unique <- unedit_sp_15mo %>%
  distinct(Subject_ID, .keep_all = TRUE)

# rerun all the subsets to have the right values
# rename the dataset for ease

sp_15mo <- unedit_15mo_unique


## Subset the data further into relevant test factor

### A) Spp in Recovery

# sp_12mo Subsetting 

sp_12_R <- sp_12mo %>% 
  filter(Recovery == "TRUE") %>% 
  group_by(Taxa)

sp_12_nR <- sp_12mo %>% 
  filter(Recovery == "FALSE") %>% 
  group_by(Taxa)


# sp_15mo Subsetting 

sp_15_R <- sp_15mo %>% 
  filter(Recovery == "TRUE") %>% 
  group_by(Taxa)

sp_15_nR <- sp_15mo %>% 
  filter(Recovery == "FALSE") %>% 
  group_by(Taxa)




### B) Spp in Group

# sp_12mo Subsetting 

sp_12_gr1 <- sp_12mo %>% 
  filter(Group == "1") %>% 
  group_by(Taxa)

sp_12_gr2 <- sp_12mo %>% 
  filter(Group == "2") %>% 
  group_by(Taxa)


# sp_15mo Subsetting 

sp_15_gr1 <- sp_15mo %>% 
  filter(Group == "1") %>% 
  group_by(Taxa)

sp_15_gr2 <- sp_15mo %>% 
  filter(Group == "2") %>% 
  group_by(Taxa)


# C) For Age 

sp_gr1 <- sp_interest %>% 
  filter(Group == "1") %>%
  distinct(Seq_ID, .keep_all = TRUE)

sp_gr2 <- sp_interest %>% 
  filter(Group == "2") %>%
  distinct(Seq_ID, .keep_all = TRUE)


# subset the data further by each recovery state

# Group 1 
sp_gr1_R <- sp_gr1 %>% 
  filter(Recovery == "TRUE")

sp_gr1_nR <- sp_gr1 %>% 
  filter(Recovery == "FALSE")


# Group 2 
sp_gr2_R <- sp_gr2 %>% 
  filter(Recovery == "TRUE")

sp_gr2_nR <- sp_gr2 %>% 
  filter(Recovery == "FALSE")




## Run Wilcoxon for each Test factor 

# For Recovery - testing the dif between recovered v not recovered
# we will use the subsets of the groups to find each p val within each Group for each age

# quadrant - 12 mo Gr1
wil_recov_12gr1 <- sp_12_gr1 %>% 
  summarise(p_value = wilcox.test(RA ~ Recovery)$p.value)
# Rothia mucilaginosa        0.453
# Streptococcus salivarius   0.261

# quadrant - 12 mo Gr2
wil_recov_12gr2 <- sp_12_gr2 %>% 
  summarise(p_value = wilcox.test(RA ~ Recovery)$p.value)
# Rothia mucilaginosa        0.184
# Streptococcus salivarius   0.951

# quadrant - 15 mo Gr1
wil_recov_15gr1 <- sp_15_gr1 %>% 
  summarise(p_value = wilcox.test(RA ~ Recovery)$p.value)
# Rothia mucilaginosa        0.386
# Streptococcus salivarius   0.432

# quadrant - 15 mo Gr2
wil_recov_15gr2 <- sp_15_gr2 %>% 
  summarise(p_value = wilcox.test(RA ~ Recovery)$p.value)
# Rothia mucilaginosa        0.271
# Streptococcus salivarius   0.838


# For Groups - testing the diff between gr1 and gr 2 
# we will use the subsets of the recovery to find each p val within each Recovery for each age

# quadrant - 12 mo R
wil_recov_12_R <- sp_12_R %>% 
  summarise(p_value = wilcox.test(RA ~ Group)$p.value)
# Rothia mucilaginosa        0.581
# Streptococcus salivarius   0.157

# quadrant - 12 mo nR
wil_recov_12_nR <- sp_12_nR %>% 
  summarise(p_value = wilcox.test(RA ~ Group)$p.value)
# Rothia mucilaginosa        0.224
# Streptococcus salivarius   0.436

# quadrant - 15 mo R
wil_recov_15_R <- sp_15_R %>% 
  summarise(p_value = wilcox.test(RA ~ Group)$p.value)
# Rothia mucilaginosa       0.0762
# Streptococcus salivarius  0.719

# quadrant - 15 mo nR
wil_recov_15_nR <- sp_15_nR %>% 
  summarise(p_value = wilcox.test(RA ~ Group)$p.value)
# Rothia mucilaginosa        0.956
# Streptococcus salivarius   0.294


# For Age - testing the diff between 12mo and 15mo 
# we will use the subsets of the groups to find each p val within each group for each recovery state

# quadrant - gr1 R
wil_recov_gr1_R <- sp_gr1_R %>% 
  summarise(p_value = wilcox.test(RA ~ Age_months, paired = TRUE)$p.value)
# Rothia mucilaginosa       0.0101
# Streptococcus salivarius  0.317 

# quadrant - gr1 nR #### there 4 extra rows that did not match and were removed,
# check the testing area at the bottom to see how the new sp_gr1_nR was formed
wil_recov_gr1_nR <- sp_gr1_nR %>% 
  summarise(p_value = wilcox.test(RA ~ Age_months, paired = TRUE)$p.value)
# Rothia mucilaginosa        0.746
# Streptococcus salivarius   0.695

# quadrant - gr2 R
wil_recov_gr2_R <- sp_gr2_R %>% 
  summarise(p_value = wilcox.test(RA ~ Age_months, paired = TRUE)$p.value)
# Rothia mucilaginosa        0.533
# Streptococcus salivarius   0.961

# quadrant - gr2 nR
wil_recov_gr2_nR <- sp_gr2_nR %>% 
  summarise(p_value = wilcox.test(RA ~ Age_months, paired = TRUE)$p.value)
# pairing mismatch - resolved below 

paired_sp_gr2_nR <- sp_gr2_nR %>%
  group_by(Subject_ID, Taxa) %>%
  filter(n_distinct(Age_months) == 2) %>%  # keep only if both 12 & 15 months exist
  ungroup()

paired_sp_gr2_nR %>%
  group_by(Taxa) %>%
  summarise(p_value = wilcox.test(RA ~ Age_months, paired = TRUE)$p.value)

# Rothia mucilaginosa        0.639
# Streptococcus salivarius   0.481


# Pathways  ---------------------------------------------------------------

# the Maaslin output from the script results_2023 is being used here 
# change the path to make it suitable 

# Read in output files
all_results_refeed_pathways <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/pathways_refeed/all_results.tsv")
sig_results_refeed_pathways <- read_tsv("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/Wellcome_Leap/Data_Analysis/m4efad/data/maaslin2/pathways_refeed/significant_results.tsv")

n_distinct(all_results_refeed_pathways$feature) # 466 pathways included in testing
n_distinct(sig_results_refeed_pathways$feature) # 383 pathways found to be significantly different between MAM and Healthy

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(sig_results_refeed_pathways, coef >0)) # 75 pathways increased
n_distinct(filter(sig_results_refeed_pathways, coef <0)) # 308 pathways decreased

# top 10 pathways that increased following refeed
head(arrange(sig_results_refeed_pathways, desc(coef)), 10) 

# top 10 pathways that decreased
head(arrange(sig_results_refeed_pathways, (coef)), 20) 


