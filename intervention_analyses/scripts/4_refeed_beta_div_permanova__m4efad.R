# This script focusses on:
#1: Beta diversity 
#2: PERMANOVA - on several variables

# Variables included in the script:

# - 12mo v 15mo 
# - GrA v GrB
# - R v nR 

# Bray-Curtis for 15mo samples only → for cross-sectional PERMANOVA at 15mo
# Use when testing if microbial composition at 15mo differs by recovery, refeed, WLZ, etc.

# Bray-Curtis within-subject (12mo vs 15mo) → for longitudinal change analysis
# Use when testing if *degree of change* in microbiome relates to outcomes at 15mo

# Full Bray-Curtis matrix (all samples) → for clustering or ordination (PCoA, NMDS)
# Use for visualizing sample similarity patterns across all subjects/timepoints


# Load Packages  ----------------------------------------------------------
library(tidyverse) # required for data wrangling and ggplot
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(dplyr)
library(stats)
library(broom)
library(readxl)



# load data ---------------------------------------------------------------
# refeed_data: contains paired data for 12mo and 15mo
load("data/processed/refeed_data.Rdata")

# contains recovery data information 
load("data/processed/recovery_data.Rdata")

# contains group related data 
load("data/processed/group_data.Rdata")

# Taxa data:
# species 
load("data/processed/refeed_species.Rdata")

# contains anthro data upto 15 months
load("data/processed/mam_anthro_measures.Rdata")

# contains all confounding factors
load("data/processed/all_meta_baseline.Rdata")

# LOAD this for complete 15mo metadata with confounding factors
load("data/processed/meta_15mo_full.Rdata")

# contains a combination of both timepoints and both group and recovery
load("data/processed/MAM_meta_whole_dataset.Rdata")

# contains the species for each of the above in a table format - used in facet plots
load('data/processed/long_spp_MAM_whole.Rdata')


# Within subject species change  --------------------------
# A distance matrix is constructed on the 12 months and 15 months timepoints for each sample
# only matched samples are used
# within sample species composition changes are generated 
# this is checked for correlation with the change in WHZ_WLZ score for the 13 weeks they were measured (some samples dropped before 13 weeks)

# Data prep 

# add timepoint and subject_ID info to the species 
refeed_species_meta <- refeed_species %>%
  rownames_to_column("Seq_ID") %>%
  inner_join(refeed_data %>% select(Seq_ID, Subject_ID, Age_months), by = "Seq_ID")


# separate samples for 12mo and 15mo 
species_12 <- refeed_species_meta %>%
  filter(Age_months == "12") %>%
  arrange(Subject_ID)

species_15 <- refeed_species_meta %>%
  filter(Age_months == "15") %>%
  arrange(Subject_ID)


# make sure matrices match
# Ensure same subjects in same order
common_subjects <- intersect(species_12$Subject_ID, species_15$Subject_ID)

# Keep only subjects that have both 12mo and 15mo samples
species_12 <- species_12 %>% filter(Subject_ID %in% common_subjects)
species_15 <- species_15 %>% filter(Subject_ID %in% common_subjects)

# Select species abundance columns (starting with 's__') and convert to matrix
mat_12 <- species_12 %>% select(starts_with("s__")) %>% as.matrix()
mat_15 <- species_15 %>% select(starts_with("s__")) %>% as.matrix()


# calculate Bray-Curtis dissimilarity between 12mo vs 15mo per Subject (ie, within individuals)
bray_within <- vegan::vegdist(rbind(mat_12, mat_15), method = "bray")

# Get diagonal from 12 vs 15
n <- nrow(mat_12)
bray_diag <- diag(as.matrix(bray_within)[1:n, (n+1):(2*n)])

bray_df <- tibble(
  Subject_ID = species_12$Subject_ID,
  Bray_Curtis = bray_diag
)


# WLZ change ~ intra subject MB changes -----------------------------------

# Get baseline (1.00)
wlz_baseline <- mam_anthro_measures %>%
  filter(An_Time == 1.00) %>%
  select(Subject_ID, WLZ_baseline = WLZ_WHZ)

# Get each subject's *last available* follow-up WLZ
wlz_final <- mam_anthro_measures %>%
  filter(An_Time > 1.00) %>%
  group_by(Subject_ID) %>%
  slice_max(An_Time) %>%
  ungroup() %>%
  select(Subject_ID, WLZ_final = WLZ_WHZ)

# Merge and calculate change
wlz_change_df <- wlz_baseline %>%
  inner_join(wlz_final, by = "Subject_ID") %>%
  mutate(WLZ_change = WLZ_final - WLZ_baseline)

# remove lines with 99.99 - ie, dropped from study entirely
wlz_change_df <- wlz_change_df %>%
  filter(WLZ_baseline != 99.99, WLZ_final != 99.99)


bray_wlz_merge <- bray_df %>%
  inner_join(wlz_change_df, by = "Subject_ID")


# test correlation
cor_test <- cor.test(bray_wlz_merge$Bray_Curtis, bray_wlz_merge$WLZ_change, method = "spearman")
print(cor_test)


# plot relationship
ggplot(bray_wlz_merge, aes(x = Bray_Curtis, y = WLZ_change)) +
  geom_point(alpha = 0.8) +
  #geom_smooth(method = "lm", color = "darkred", se = TRUE) +
  theme_minimal() +
  labs(
    x = "Bray–Curtis Dissimilarity (12 vs 15mo)",
    y = "Change in WLZ",
    title = "Microbiome Shifts vs WLZ Recovery"
  )






# Linear Model  --------------------------------------------------------

# Data structuring 
# select final anthropometric measures 

# Data structuring --------------------------------------------------------

# Get each subject's *last available* follow-up anthro
# HC removed as there is only one value for most subjects,ie., not accurate of the 15mo timepoint
anthro_final_15mo <- mam_anthro_measures %>%
  filter(An_Time > 1.00) %>%
  group_by(Subject_ID) %>%
  slice_max(An_Time) %>%
  ungroup() %>%
  select(Subject_ID, 
         WLZ_final = WLZ_WHZ,
         Weight_final = Weight,
         Length_final = Length,
        MUAC_final = MUAC) # n = 159

# remove lines with 99.99 - ie, dropped from study entirely
anthro_final_15mo <- anthro_final_15mo %>%
  filter(WLZ_final != 99.99) # n = 144

# Check for any 99.99 values in the dataset
anthro_final_15mo %>%
  summarise(across(everything(), ~ sum(. == 99.99)))

# pull other variables suh as sex, group, recovery status form relevant datasets

meta1_15mo <- recovery_data %>% 
  select(Subject_ID, Seq_ID, Recovery, base_Weight, base_Length,
         base_MUAC, base_WLZ_WHZ, Seq_batch, Sample_Name)

meta2_15mo <- group_data %>% 
  select(Subject_ID, Seq_ID, Group)

meta_15mo <- meta1_15mo %>%
  inner_join(meta2_15mo, by = c("Subject_ID", "Seq_ID"))


confound_factors <- all_meta_baseline %>% 
  select(Subject_ID, Sex, Place_of_birth, 
         Breastfeeding_duration_months, Delivery_Mode)

meta_15mo_full <- meta_15mo %>%
  inner_join(confound_factors, by = c("Subject_ID"))

meta_15mo_full <- meta_15mo %>%
  inner_join(confound_factors, by = "Subject_ID") %>%
  inner_join(anthro_final_15mo, by = "Subject_ID")


save(meta_15mo_full, file = "data/processed/meta_15mo_full.Rdata")



# running  ------------------------------------------------------------

# Combine Bray-Curtis dissimilarity values with full 15mo metadata by Subject_ID
lm_df <- bray_df %>%
  inner_join(meta_15mo_full, by = "Subject_ID")

# Run linear model
lm_result <- lm(
  Bray_Curtis ~ Recovery + Group + Sex + base_WLZ_WHZ + WLZ_final + 
    Place_of_birth + Delivery_Mode + Breastfeeding_duration_months,
  data = lm_df
)

summary(lm_result)

# univariate models 
# Run univariate linear models
lm_recovery <- lm(Bray_Curtis ~ Recovery, data = lm_df)
lm_group <- lm(Bray_Curtis ~ Group, data = lm_df)
lm_sex <- lm(Bray_Curtis ~ Sex, data = lm_df)
lm_base_wlz <- lm(Bray_Curtis ~ base_WLZ_WHZ, data = lm_df)
lm_final_wlz <- lm(Bray_Curtis ~ WLZ_final, data = lm_df)
lm_birthplace <- lm(Bray_Curtis ~ Place_of_birth, data = lm_df)
lm_delivery <- lm(Bray_Curtis ~ Delivery_Mode, data = lm_df)
lm_breastfeed <- lm(Bray_Curtis ~ Breastfeeding_duration_months, data = lm_df)

# View summaries
summary(lm_recovery)
summary(lm_group)
summary(lm_sex)
summary(lm_base_wlz)
summary(lm_final_wlz)
summary(lm_birthplace)
summary(lm_delivery)
summary(lm_breastfeed)



# Permanova  --------------------------------------------------------------

# Q: which combination of variables best explain the variation in species composition?

# BC run all the whole species dataset which contains all the combinations of samples and variables 

# Run Bray Curtis on the whole species dataset
# Calculate Bray-Curtis Dissimilarity Index (vegan package is used)
BC_dis <- vegdist(species, method = "bray") 

# Data re-structuring 
BC_dis_long <- BC_dis %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample1") %>%
  gather(Sample2, BC_dis, -Sample1) %>%
  mutate(ID = apply(cbind(as.character(Sample1), as.character(Sample2)), 1, function(x) {str_c(sort(x), collapse = ":")})) %>%
  distinct(ID, .keep_all = T) %>% # removes duplicates
  filter(!(Sample1 == Sample2)) %>% # removes rows where samples are compared to themselves
  select(-ID)

# For visualisation of the Brays-Curtis dissimilarity index - NMDS
# Perform Ordination using Nonmetric Multidimensional Scaling (nMDS)
set.seed(280) # set a seed to ensure reproducible results (applicable to methods involving randomness)

# Running mds on the Bray-Curtis distance
mds <- metaMDS(BC_dis)

# join to the MAM table 
mds_data_all_variables <- mds$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Seq_ID") %>% 
  inner_join(long_spp_MAM_whole)



# Plots  ------------------------------------------------------------------

# Differences in Group for Age and Recovery
plot_12_15_beta <- ggplot(mds_data_all_variables, aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = Age_months, color = Age_months),
    shape = 21,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  scale_color_manual(values = c("12" = "darkgoldenrod1", "15" = "darkmagenta")) +
  facet_grid(Recovery ~ Group) +
  theme_bw() +
  theme(
   panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    aspect.ratio = 1
  ) +
  xlab("NMDS1") +
  ylab("NMDS2") 

# pdf
ggsave(file = "plots/chp2/clean_plots/plot_12_15_beta.pdf", 
       plot = plot_12_15_beta, units = "cm", width=11, height=11)


# Differences in Recovery for Age and Group
plot_recovery <- ggplot(mds_data_all_variables, aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = Recovery, color = Recovery),
    shape = 21,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  scale_color_manual(values = c("TRUE" = "darkkhaki", "FALSE" = "darkblue")) +
  facet_grid(Group ~ Age_months) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    aspect.ratio = 1
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")

# pdf
ggsave(file = "plots/chp2/clean_plots/plot_recovery.pdf", 
       plot = plot_recovery, units = "cm", width=11, height=11)


# Differences in Group for Age and Recovery
plot_group <- ggplot(mds_data_all_variables, aes(x = MDS1, y = MDS2)) +
  geom_point(
    aes(fill = Group, color = Group),
    shape = 21,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("1" = "azure3", "2" = "black")) +
  scale_color_manual(values = c("1" = "azure3", "2" = "black")) +
  facet_grid(Recovery ~ Age_months) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    aspect.ratio = 1
  ) +
  xlab("NMDS1") +
  ylab("NMDS2") 


# pdf
ggsave(file = "plots/chp2/clean_plots/plot_group.pdf", 
       plot = plot_group, units = "cm", width=11, height=11)

# Stats -------------------------------------------------------------------

# We will be running 8 permanovas in total
# These will be checking the effect of species on Recovery within each Refeed Type and,
# the effect of Refeed Type within each recovery status on species composition divided up by Age 

# First, we will create subsets of data

# We will use the MAM_meta_whole_dataset as this contains the Arms for each age group


# Prepare data ------------------------------------------------------------

### STEP 1: Subset the data into 12mo and 15mo (MAM_meta_whole_dataset, n = 732)

# 12mo dataset "mo12"
# n = 292

mo12 <- MAM_meta_whole_dataset %>% 
  filter(Age_months == "12")


# 15mo dataset "mo15"
# n = 440

mo15 <- MAM_meta_whole_dataset %>% 
  filter(Age_months == "15")


### STEP 2: Remove "Pre-recovery" data points from mo12 as this is redudant

mo12 <- mo12 %>% 
  filter(Recovery != "Pre-recovery")
# n = 146

### STEP 3: There will be further subsetting based on the test factor 


# A) Spp ~ Recovery in Groups ------------------------------------------------
##  Spp ~ Recovery

# Subset mo12 and mo15 for into each Refeed Group - 4 datasets in total

# mo12_gr1, mo12_gr2 
# mo15_gr1, mo15_gr2

# mo12 Subsetting 

mo12_gr1 <- mo12 %>% 
  filter(Group == "1")

mo12_gr2 <- mo12 %>% 
  filter(Group == "2")


# mo15 Subsetting 

mo15_gr1 <- mo15 %>% 
  filter(Group == "1")

mo15_gr2 <- mo15 %>% 
  filter(Group == "2")

## MO15 datasets have duplicated rows 
# these will be removed before proceeding

duplicated_subjects_mo15 <- mo15 %>%
  filter(duplicated(Subject_ID) | duplicated(Subject_ID, fromLast = TRUE))

mo15_unique <- mo15 %>%
  distinct(Subject_ID, .keep_all = TRUE)

## rerun all the subsets to have the right values
## rename the dataset for ease

mo15 <- mo15_unique

# mo15 Subsetting 

mo15_gr1 <- mo15 %>% 
  filter(Group == "1")

mo15_gr2 <- mo15 %>% 
  filter(Group == "2")


### STEP 4: filter the spp_whole_dataset for each age_group combo

# convert the spp_whole_dataset rows into column by Seq_ID

spp_whole_dataset <- spp_whole_dataset %>% 
  rownames_to_column(var = "Seq_ID")

### MO12
## mo12_gr1_spp

mo12_gr1_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo12_gr1$Seq_ID)

## mo12_gr2_spp

mo12_gr2_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo12_gr2$Seq_ID)

### MO15
## mo15_spp 

#### Create a test mo dataset to combine to spp 

# mo15_gr1_spp
mo15_gr1_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo15_gr1$Seq_ID)

# mo15_gr2_spp
mo15_gr2_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo15_gr2$Seq_ID)


# remove the Seq_ID column in all the spp datasets

mo12_gr1_spp <- mo12_gr1_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo12_gr2_spp <- mo12_gr2_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo15_gr1_spp <- mo15_gr1_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo15_gr2_spp <- mo15_gr2_spp %>% 
  column_to_rownames(var = "Seq_ID")



# Run Permanova -----------------------------------------------------------

### STEP 5: run the permanova for each group 

## MO12
# Recovery ~ mo12gr1
#species against recovery (y ~ x), 
# Is there a significant difference in the species between 
# recovered v not recovered in the 12 months infants in the refeed group 1 and 2?

permmo12gr1 <- adonis2(mo12_gr1_spp ~ Recovery, 
                       data = mo12_gr1,
                       by = "margin", 
                       permutations = 999,
                       method = "bray")
# p = 0.847
# R2 = 0.009016266


# Recovery ~ mo12gr2

permmo12gr2 <- adonis2(mo12_gr2_spp ~ Recovery, 
                       data = mo12_gr2,
                       by = "margin", 
                       permutations = 999,
                       method = "bray")

# p = 0.175
# R2 = 0.01836217


## MO15
# Recovery ~ mo15gr1
#species against recovery (y ~ x), 
# Is there a significant difference in the species between 
# recovered v not recovered in the 15 months infants in the refeed group 1 and 2?

permmo15gr1 <- adonis2(mo15_gr1_spp ~ Recovery, 
                       data = mo15_gr1,
                       by = "margin", 
                       permutations = 999,
                       method = "bray")

# p = 0.763
# R2 = 0.009669809

# Recovery ~ mo12gr2
permmo15gr2 <- adonis2(mo15_gr2_spp ~ Recovery, 
                       data = mo15_gr2,
                       by = "margin", 
                       permutations = 999,
                       method = "bray")

# p = 0.118
# R2 = 0.02082594




# B) Spp ~ Groups in R v nR --------------------------------------------------
##  Spp ~ Groups 

# Subset mo12 and mo15 for into each recovery state - 4 datasets in total

# mo12_R, mo12_nR 
# mo15_R, mo15_nR

# mo12 Subsetting 

mo12_R <- mo12 %>% 
  filter(Recovery == "TRUE")

mo12_nR <- mo12 %>% 
  filter(Recovery == "FALSE")


# mo15 Subsetting 

mo15_R <- mo15 %>% 
  filter(Recovery == "TRUE")

mo15_nR <- mo15 %>% 
  filter(Recovery == "FALSE")


## Combine to form matching spp dataset
# convert the spp_whole_dataset rows into column by Seq_ID

spp_whole_dataset <- spp_whole_dataset %>% 
  rownames_to_column(var = "Seq_ID")

### MO12
## mo12_R_spp

mo12_R_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo12_R$Seq_ID)

## mo12_gr2_spp

mo12_nR_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo12_nR$Seq_ID)

### MO15
## mo15_spp 

#### Create a test mo dataset to combine to spp 

## mo15_R_spp

mo15_R_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo15_R$Seq_ID)

## mo15_nR_spp

mo15_nR_spp <- spp_whole_dataset %>% 
  filter(Seq_ID %in% mo15_nR$Seq_ID)


# remove the Seq_ID column in all the spp datasets

mo12_R_spp <- mo12_R_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo12_nR_spp <- mo12_nR_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo15_R_spp <- mo15_R_spp %>% 
  column_to_rownames(var = "Seq_ID")

mo15_nR_spp <- mo15_nR_spp %>% 
  column_to_rownames(var = "Seq_ID")



#  Run Permanova ----------------------------------------------------------

# Group ~ mo12_R
permmo12R <- adonis2(mo12_R_spp ~ Group, 
                     data = mo12_R,
                     by = "margin", 
                     permutations = 999,
                     method = "bray")
# p = 0.936
# R2 = 0.008457898


# Group ~ mo12_nR

permmo12nR <- adonis2(mo12_nR_spp ~ Group, 
                      data = mo12_nR,
                      by = "margin", 
                      permutations = 999,
                      method = "bray")

# p = 0.174
# R2 = 0.01656208


## MO15
# Recovery ~ mo15gr1
#species against recovery (y ~ x), 
# Is there a significant difference in the species between 
# recovered v not recovered in the 15 months infants in the refeed group 1 and 2?

permmo15R <- adonis2(mo15_R_spp ~ Group, 
                     data = mo15_R,
                     by = "margin", 
                     permutations = 999,
                     method = "bray")

# p = 0.047
# R2 = 0.02776802

# Recovery ~ mo12gr2
permmo15nR <- adonis2(mo15_nR_spp ~ Group, 
                      data = mo15_nR,
                      by = "margin", 
                      permutations = 999,
                      method = "bray")

# p = 0.455
# R2 = 0.01183985























