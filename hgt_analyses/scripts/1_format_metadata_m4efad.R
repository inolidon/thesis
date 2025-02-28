## code summary:
# format metadata of m4efad study to match the datasets of Gut-Bugs and FOCUS
# make group names consistent
# subset samples by timepoint and conditions (MAM v Healthy, 12, 15 and 24 months)
# check metadata sample IDs consistent with waafle and metaphlan sample IDs (very important)
# quantify samples with waafle data and metaphlan data


# load libraries  ---------------------------------------------------------

library(tidyverse)

## set working directory

# setwd("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad")


# load master m4efad metadata file  ---------------------------------------

# use on local PC
# load("data/processed/updated_meta.Rdata")

# modify metadata to match other studies  ---------------------------------

# format column names and remove additional columns 
# formatted_metadata_m4efad <- updated_meta %>%
#  rename(Sample_ID = Seq_ID,                 # renaming the column names to match
#         Participant_ID = Subject_ID,
#         Group = Condition,
#         Timepoint = Age_months) %>%
#  select(-Seq_batch, -Sample_Name)  # removing Seq_batch and Sample_name column 
# n = 681 


# create generalised timepoints for plotting ## 12 v15 v 24 timepoints (MAM and Healthy)
# formatted_metadata_m4efad <- formatted_metadata_m4efad %>%
#  mutate(Timepoint_general = case_when(Timepoint == "12" & Group == "MAM" ~ "MAM_12", # make generalised timepoints
#                                       Timepoint == "15" & Group == "MAM" ~ "MAM_15",
#                                       Timepoint == "24" & Group == "MAM" ~ "MAM_24",
#                                       Timepoint == "12" & Group == "Healthy" ~ "Healthy_12",
#                                       Timepoint == "24" & Group == "Healthy" ~ "Healthy_24")) %>%
#  drop_na(Timepoint_general) # remove the rows that have Timepoint_general == N/A
# these are the controls (standards and 36 months timepoint)
# n = 588 



# Start Here  -------------------------------------------------------------


# load functions ----------------------------------------------------------

# subset samples by generalised timepoint ## subset to a specific timepoint and compare conditions
subset_samples <- function(formatted_metadata, timepoint_subset){
  ## formatted_metadata: metadata containing at least columns 'Sample_ID', 'Participant_ID', 'Group', 'Timepoint', 'Timepoint_general' (df)
  ## timepoint_subset: timepoint from 'Timepoint_general' column to subset (str)
  formatted_metadata %>% filter(Timepoint_general == timepoint_subset)
}

# load formatted metaphlan4 GTDB or SGB species output
load_formatted_m4_spp <- function(formatted_metaphlan){
  ## formatted_metaphlan: path to metaphlan data, data contains the object 'species' where rownames = sample IDs, colnames = species (str)
  load(formatted_metaphlan)
  species %>%
    rownames_to_column("Sample_ID") %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) # remove underscores from metaphlan sample ID to match metadata
}



# load data  ----------------------------------------------

# 1.metadata --------------------------------------------------------------
# use on local PC
# load("data/processed/updated_meta.Rdata")

# use on github repo 
load("hgt_analyses/data/updated_meta.Rdata")

# modify metadata to match other studies  ---------------------------------

# format column names and remove additional columns 
formatted_metadata_m4efad <- updated_meta %>%
  rename(Sample_ID = Seq_ID,                 # renaming the column names to match
         Participant_ID = Subject_ID,
         Group = Condition,
         Timepoint = Age_months) %>%
  select(-Seq_batch, -Sample_Name)  # removing Seq_batch and Sample_name column 
# n = 681 


# create generalised timepoints for plotting ## 12 v15 v 24 timepoints (MAM and Healthy)
formatted_metadata_m4efad <- formatted_metadata_m4efad %>%
  mutate(Timepoint_general = case_when(Timepoint == "12" & Group == "MAM" ~ "MAM_12", # make generalised timepoints
                                       Timepoint == "15" & Group == "MAM" ~ "MAM_15",
                                       Timepoint == "24" & Group == "MAM" ~ "MAM_24",
                                       Timepoint == "12" & Group == "Healthy" ~ "Healthy_12",
                                       Timepoint == "24" & Group == "Healthy" ~ "Healthy_24")) %>%
  drop_na(Timepoint_general) # remove the rows that have Timepoint_general == N/A
# these are the controls (standards and 36 months timepoint)
# n = 588 



# 2. waafle ---------------------------------------------------------------
# this dataset contains the main hgt dataframe, the directed hgt, the undirected hgt
# and the species level filtered hgts 

# on local PC
# load("data/processed/hgt/hgt_dataframes.Rdata") # use the hgt_df (main - whole output)

# rename dataset & column in hgt_df to match other studies + formatted_metadata

# hgt_df_m4efad <- hgt_df %>% 
#  rename(Sample_ID = Seq_ID)


# running -----------------------------------------------------------------

# subset samples and quantify

# MAM_12
MAM_12_samples <- subset_samples(formatted_metadata_m4efad, "MAM_12") # n = 157 

# MAM_15
MAM_15_samples <- subset_samples(formatted_metadata_m4efad, "MAM_15") # n = 148 

# MAM_24
MAM_24_samples <- subset_samples(formatted_metadata_m4efad, "MAM_24") #n = 144

# Healthy_12
Healthy_12_samples <- subset_samples(formatted_metadata_m4efad, "Healthy_12") # n = 74

# Healthy_24
Healthy_24_samples <- subset_samples(formatted_metadata_m4efad, "Healthy_24") # n = 65



# Do not run this section ---------------------------------------------------------------

## TESTS

# check formatting of sample IDs is the same as in waafle
length(which(formatted_metadata_m4efad$Sample_ID %in% hgt_df_m4efad$Sample_ID)) # 577 = 11 samples missing 
# (samples did not have any hgts)


# load formatted metaphlan4 GTDB or SGB species output
# note: load formatted metaphlan output into function as they all have objects with the same names

# m4 sgb
m4_spp_sgb_m4efad <- load_formatted_m4_spp("data/processed/clean_outputs/metaphlan4_sgb_m4efad.Rdata") 

# m4 gtdb
m4_spp_gtdb_m4efad <- load_formatted_m4_spp("data/processed/clean_outputs/metaphlan4_gtdb_m4efad.Rdata") 


# !! 24 months samples missing  -------------------------------------------

# check formatting of sample IDs is the same as in metaphlan species
# note: underscores have been removed from metaphlan sample IDs, must do this again when using the data
length(which(formatted_metadata_m4efad$Sample_ID %in% m4_spp_sgb_m4efad$Sample_ID)) # 378 
length(which(formatted_metadata_m4efad$Sample_ID %in% m4_spp_gtdb_m4efad$Sample_ID)) # 378 

missing_samples <- formatted_metadata_m4efad$Sample_ID[!formatted_metadata_m4efad$Sample_ID %in% m4_spp_sgb_m4efad$Sample_ID]

missing_samples_data <- formatted_metadata_m4efad %>% 
  filter(Sample_ID %in% missing_samples)


matching_samples_24 <- formatted_metadata_m4efad %>% 
  filter(Sample_ID %in% m4_spp_sgb_m4efad$Sample_ID, Timepoint == 24)



# Continue here  ----------------------------------------------------------


# save data files ---------------------------------------------------------

# save(formatted_metadata_m4efad, file = "data/processed/clean_outputs/metadata_m4efad.RData")
# save(MAM_12_samples, MAM_15_samples, MAM_24_samples, Healthy_12_samples, Healthy_24_samples, 
  #   file = "data/processed/clean_outputs/sample_subsets.RData")

# waafle file
# save(hgt_df_m4efad, file = "data/processed/clean_outputs/hgt_df_m4efad.RData")




