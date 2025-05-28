## This script contains the codes to process raw data files such as the 
# Species etc the outputs of Metaphlan and any other downstream processed 
# data table including baseline specific taxonomy files, edited metadata files etc. 


# Packages  ---------------------------------------------------------------

# load or install packages as needed, below are the packages used for this script

library(tidyverse) # required for data wrangling and ggplot
library(scales) # required for label_comma()
library(data.table) # require for gene table cleaning
library(readxl) # to read in data from excel files such as metadata


# PART 1: MetaPhlAn3 - Data cleaning -----------------------------------------------------------

## When the data updates, run this section again to updata all taxonony files 
# This is data as at OCT 2023 (seqbx4)
# updated with seqbx5 - April 2024 

taxonomy_raw <- read_tsv("data/raw/m4efad_metaphlan3_profiles_april2024.tsv", skip = 1)

# Generate a taxa list and summarise the number of distinct taxa within each rank
taxa_list <- taxonomy_raw %>% 
  select(clade_name) %>% 
  filter(grepl("s__", clade_name), !grepl("t__", clade_name)) %>% 
  distinct() %>% 
  separate(clade_name, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "\\|", remove = T)

taxa_list %>% 
  gather(Rank, clade_name) %>% 
  mutate(Rank = factor(Rank, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))) %>% 
  group_by(Rank) %>% 
  summarise(n_distinct = n_distinct(clade_name))

table(taxa_list$Kingdom)
table(taxa_list$Phylum)

# Create an output table to identify the distinct phyla and the number of samples that have that particle phylum

output_table <- taxa_list %>% 
  group_by(Phylum) %>%
  summarise(Number_of_Occurrences = n()) %>%
  arrange(desc(Number_of_Occurrences))

save(taxa_list, file = "data/processed/metaphlan3_taxa_list.Rdata") 


# Clean up taxonomic table
taxonomy <- taxonomy_raw  %>% 
  select(-NCBI_tax_id) %>% 
  rename(Taxa = clade_name) %>% 
  gather(Seq_ID, RA, -Taxa) %>% # covert to long format so it's easier to edit the Sample IDs
  mutate(Seq_ID = str_remove(Seq_ID, ".metaphlan"), # remove suffix from Sample ID
         Taxa = str_remove(Taxa, ".*\\|")) # simplify taxa name (removes everything before the last | delimiter)

save(taxonomy, file = "data/processed/metaphlan3_taxonomy.Rdata")


# isolate taxonomic ranks
species <- filter(taxonomy, grepl("s__", Taxa)) %>% spread(Seq_ID, RA)
genus <- filter(taxonomy, grepl("g__", Taxa)) %>% spread(Seq_ID, RA)
family <- filter(taxonomy, grepl("f__", Taxa)) %>% spread(Seq_ID, RA)
phylum <- filter(taxonomy, grepl("p__", Taxa)) %>% spread(Seq_ID, RA)
kingdom <- filter(taxonomy, grepl("k__", Taxa)) %>% spread(Seq_ID, RA)
class <- filter(taxonomy, grepl("c__", Taxa)) %>% spread(Seq_ID, RA)
order <- filter(taxonomy, grepl("o__", Taxa)) %>% spread(Seq_ID, RA)

# after removing unknown counts, re-normalise relative abundances to sum to 1
species[,-1] <- lapply(species[,-1], function(x){ x/sum(x, na.rm=TRUE)})
genus[,-1] <- lapply(genus[,-1], function(x){ x/sum(x, na.rm=TRUE)})
family[,-1] <- lapply(family[,-1], function(x){ x/sum(x, na.rm=TRUE)})
phylum[,-1] <- lapply(phylum[,-1], function(x){ x/sum(x, na.rm=TRUE)})
kingdom[,-1] <- lapply(kingdom[,-1], function(x){ x/sum(x, na.rm=TRUE)})
class[,-1] <- lapply(class[,-1], function(x){ x/sum(x, na.rm=TRUE)})
order[,-1] <- lapply(order[,-1], function(x){ x/sum(x, na.rm=TRUE)})

# check samples now sum to 1
colSums(species[,-1])
colSums(genus[,-1])
colSums(family[,-1])
colSums(phylum[,-1])
colSums(kingdom[,-1])
colSums(class[,-1])
colSums(order[,-1])

# transpose data (format required for diversity metrics)
species_all_samples <- column_to_rownames(species, "Taxa") %>% t() %>% as.data.frame()
genus_all_samples <- column_to_rownames(genus, "Taxa") %>% t() %>% as.data.frame()
family_all_samples <- column_to_rownames(family, "Taxa") %>% t() %>% as.data.frame()
phylum_all_samples <- column_to_rownames(phylum, "Taxa") %>% t() %>% as.data.frame()
kingdom_all_samples <- column_to_rownames(kingdom, "Taxa") %>% t() %>% as.data.frame()
class_all_samples <- column_to_rownames(class, "Taxa") %>% t() %>% as.data.frame()
order_all_samples <- column_to_rownames(order, "Taxa") %>% t() %>% as.data.frame()

# Save taxonomy tables
save(species_all_samples, file = "data/processed/species_all_samples.Rdata")
save(genus_all_samples, file = "data/processed/genus_all_samples.Rdata")
save(family_all_samples, file = "data/processed/family_all_samples.Rdata")
save(phylum_all_samples, file = "data/processed/phylum_all_samples.Rdata")
save(kingdom_all_samples, file = "data/processed/kingdom_all_samples.Rdata")
save(class_all_samples, file = "data/processed/class_all_samples.Rdata")
save(order_all_samples, file = "data/processed/order_all_samples.Rdata")


# UNCERTAIN: taxonomy - each chapter  ------------------------------------------------

## baseline only 
base_taxonomy_table <- taxonomy %>%
  filter(Seq_ID %in% baseline_data$Seq_ID)

# Assuming the taxonomy information is in a column named 'Taxonomy'
taxonomy_counts <- base_taxonomy_table %>%
  mutate(
    Category = case_when(
      grepl("^k__", Taxa) ~ "Kingdom",
      grepl("^p__", Taxa) ~ "Phylum",
      grepl("^c__", Taxa) ~ "Class",
      grepl("^o__", Taxa) ~ "Order",
      grepl("^f__", Taxa) ~ "Family",
      grepl("^g__", Taxa) ~ "Genus",
      grepl("^s__", Taxa) ~ "Species",
      TRUE ~ "Unknown"
    )
  ) %>%
  group_by(Category) %>%
  summarize(Count = n())





# PART 2A: Metadata - stratifying --------------------------------------------------

# we are loading the main meta file that contains all sequence information of all the samples upto date OCT 2023 
# this file is created in the metadata_excel_sheet_combining.R script 
meta <- read_csv("data/processed/meta.csv")

# meta.csv file has an additional column called "Sheet" which we will remove
meta <- subset(meta, select = -Sheet)

# save as RData file for ease of use
save(meta, file = "data/processed/meta.Rdata") # contains all the metadata for m4efad (Oct 2023 update)

# adding a new column with the ARM of the study
meta_arm <- meta %>%
  mutate(Arm = paste(Condition, Age_months, sep = "_"))

save(meta_arm, file = "data/processed/meta_arm.Rdata") # contains metadata with an ARM column

# load the R data file as needed
load("data/processed/meta_arm.Rdata")

# creating a separate dataset for each of the testing facts
### 1: MAM v Healthy - chapter 1
### 2: 12 months v 15 months - chapter 2
### 3: Group 1 v Group 2 - chapter 1 & 2
### 4: Recovery v Non-recovered - chapter 2

## 1: MAM v Healthy
# baseline dataset for mam v healthy without anthropometrics
baseline_data <- meta %>% 
  filter(Condition %in% c("MAM", "Healthy") & Age_months == 12)

save(baseline_data, file = "data/processed/baseline_data.Rdata")

load("data/processed/baseline_data.Rdata")

# 2: 12 months v 15 months 

# we use the meta_arm dataset to filter as we need MAM 12mo + 15mo
refeed_data <- meta_arm %>% 
  filter(Arm %in% c("MAM_12", "MAM_15"))%>%
  arrange(Subject_ID) # at this stage, this data set contains all pairs both matched and unmatched

# identify the unmatched samples, i.e., ones without a second timepoint
unpaired_subjects <- refeed_data %>%
  group_by(Subject_ID) %>%
  filter(n_distinct(Age_months) == 1) # 9 samples have only one time point, unmatched

# remove the unpaired samples from the main dataset to only have paired samples
refeed_data <- refeed_data %>%
  anti_join(unpaired_subjects, by = "Subject_ID") %>% # anti_join removes the unpaired subjects
  mutate(Age_months = as.factor(Age_months)) # 296 samples have both time points

save(refeed_data, file = "data/processed/refeed_data.Rdata")
load("data/processed/refeed_data.Rdata")

# refeed_data contains the refeed dataset without anthro
refeed_data <- age12_age15 %>% 
  mutate(Age_months = as.factor(Age_months))

save(refeed_data, file = "data/refeed_data.Rdata")
load("data/refeed_data.Rdata")

# 3: Group 1 v Group 2
# to create a group specific dataset, we will load a new dataset that has this information 

LEAP_Anthro_Recovery_status_files_including_samples_batch4_5_1 <- read_excel("data/raw/LEAP Anthro Recovery status files including samples_batch4_5 1.xlsx", 
                                                                             col_types = c("skip", "text", "skip", 
                                                                                           "text", "skip", "text", "text", "text", 

                                                                                                                                                                                      "text", "text"))


# create a table for the dataset
group_data <- LEAP_Anthro_Recovery_status_files_including_samples_batch4_5_1

# Changed the column name on the recovery_meta file from SID to Subject_ID to match the datasets
colnames(group_data)[colnames(group_data) == "SID"] <- "Subject_ID"

# Changed the column name on the Feeds by Randomization to Feed_tpe
colnames(group_data)[colnames(group_data) == "Feeds by Randomization"] <- "Group"

# Combine the meta dataset with the group data
group_data <- inner_join(group_data, meta, by = "Subject_ID")

# filter out the 15 months samples
group_data <- group_data %>% filter(Age_months == "15")


save(group_data, file = "data/processed/group_data.Rdata") #this contains the data for refeed at 15 months where they began receiving either Gr1 or Gr1
load("data/processed/group_data.Rdata")


# 4: Recovered v Non-recovered
# the recovery dataset is also created using this same file - load again if needed

load("data/processed/meta.Rdata")


LEAP_Anthro_Recovery_status_files_including_samples_batch4_5_1 <- read_excel("data/raw/LEAP Anthro Recovery status files including samples_batch4_5 1.xlsx", 
                                                                             col_types = c("skip", "text", "skip", 
                                                                                           "text", "skip", "text", "text", "text", 
                                                                                           "text", "text"))


# create a table for the dataset
recovery_data <- LEAP_Anthro_Recovery_status_files_including_samples_batch4_5_1

# Changed the column name on the recovery_meta file from SID to Subject_ID to match the datasets
colnames(recovery_data)[colnames(recovery_data) == "SID"] <- "Subject_ID"

# Changed the column name on the Feeds by Randomization to Feed_tpe
colnames(recovery_data)[colnames(recovery_data) == "Feeds by Randomization"] <- "Group"

# Combine the meta dataset with the recovery data
recovery_data <- inner_join(recovery_data, meta, by = "Subject_ID")

# filter out the 15 months samples
recovery_data <- recovery_data %>% filter(Age_months == "15")

save(recovery_data, file = "data/processed/recovery_data.Rdata") # contains all the recovery metadata for m4efad 
load("data/processed/recovery_data.Rdata")

# 5: 3 yr data
yr3_data <- meta %>% 
  filter(Age_months == "36")

save(yr3_data, file = "data/processed/yr3_data.Rdata")

load("data/processed/yr3_data.Rdata")



# metadata and anthropometrics combined -----------------------------------


# combine all anthro  -----------------------------------------------------
# there are 3 separate files with anthropoemtric data taken at baseline
# we combine these files into one file and combine it to metadata
# we then filter the dataset as needed 
# these files are only suitable for the chp1 analyses ONLY
# as they contain one or two measurements of anthropometrics and NOTHING beyond that
# for other analyes beyond baseline, we will use the anthro_Healthy_2024 and anthro_MAM_2024 datasets
# or the combined version suitable for the relevant chapter


# Importing the excel files - the "X" is added in the front as R can't have a number as the object name to begin with

# 1 year HEALTHY

X1yr_Healthy_anthro <- read_excel("data/raw/1yr_Healthy_anthro.xlsx", 
                                  col_types = c("text", "date", "text", 
                                                "text", "numeric", "date", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "text"))
View(X1yr_Healthy_anthro)

# 1 year MAM

X1yr_MAM_anthro <- read_excel("data/raw/1yr_MAM_anthro.xlsx", 
                              col_types = c("text", "date", "text", 
                                            "text", "numeric", "date", "numeric", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "text"))
View(X1yr_MAM_anthro)

# 3 year MAM

X3yr_MAM_anthro <- read_excel("data/raw/3yr_MAM_anthro.xlsx", 
                              col_types = c("text", "date", "text", 
                                            "text", "numeric", "date", "numeric", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "text"))
View(X3yr_MAM_anthro)


# Combine by the rows

combined_anthro <- rbind(X1yr_MAM_anthro, X1yr_Healthy_anthro, X3yr_MAM_anthro)

# Combine the columns of combined_data and combined_anthro by Subject_ID

meta_anthro <- left_join(meta, combined_anthro, by = "Subject_ID")
save(meta_anthro, file = "data/processed/meta_anthro.Rdata")




# chp1: baseline anthro only ----------------------------------------------


# we will use the meta_anthro data file that contains all the anthropometrics data 
# we will filter for the Seq_IDs within baseline_data 
# it is important that the the number of rows are the same between baseline_species 
# and the new meta information data with weight 

# meta_anthro file currently contains all age groups 
# we will make some modifications to the table and select the columns of interest

meta_baseline <- meta_anthro %>% 
  filter(Age_months == "12") %>% 
  filter(An_Time == 1.00) %>% 
  select(Seq_ID, Sample_Name, Age_months, Condition, Subject_ID, 
         SEX, Weight, Length, MUAC, HC, WLZ_WHZ, Group, Seq_batch)


save(meta_baseline, file = "data/processed/meta_baseline.Rdata")



# PART 2B: Adding additional info to metadata -----------------------------

# This section focuses on adding additional and useful information to the main metadata file 
# for each specific dataset such as baseline, refeed, recovery and group metadata
# additional information added: Sex, Group, initial weight, length, muac, hc and WLZ_WHZ ratio

# load the data 

load("data/processed/baseline_data.Rdata")
load("data/processed/refeed_data.Rdata")
load("data/processed/recovery_data.Rdata")
load("data/processed/group_data.Rdata")

# load the raw files 
# the files that contain this information are separated into MAM and Healthy
# we will load them separately and then match to the appropriate Subject_ID

# contains Healthy information - used to add the values for baseline_data
anthro_Healthy_2024 <- read_excel("data/raw/anthro_Healthy_2024.xlsx", 
                                  col_types = c("text", "skip", "text", 
                                                "skip", "numeric", "skip", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "text")) # n = 600 
# first we will rename the columns as needed 
anthro_Healthy_2024 <- anthro_Healthy_2024 %>% 
  rename(Subject_ID = SID_LCC,
         Sex = SEX,
         base_Weight = Weight,
         base_Length = Length, 
         base_MUAC = MUAC,
         base_HC = HC,
         base_WLZ_WHZ = WLZ_WHZ)


# next, we only want the first measurement for baseline measurements
# so, we filter by An_Time where An_Time == 1 

anthro_Healthy_2024 <- anthro_Healthy_2024 %>% 
  filter(An_Time == 1.00) # n = 75

# contains all the MAM information 
anthro_MAM_2024 <- read_excel("data/raw/anthro_MAM_2024.xlsx", 
                              col_types = c("text", "skip", "text", 
                                            "text", "skip", "numeric", "skip", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "numeric", "text")) # n = 2167

  

# first we will rename the columns as needed 
anthro_MAM_2024 <- anthro_MAM_2024 %>% 
  rename(Subject_ID = SID_LCC,
         Sex = SEX,
         base_Weight = Weight,
         base_Length = Length, 
         base_MUAC = MUAC,
         base_HC = HC,
         base_WLZ_WHZ = WLZ_WHZ)


# next, we only want the first measurement for baseline measurements
# so, we filter by An_Time where An_Time == 1 
anthro_MAM_2024 <- anthro_MAM_2024 %>% 
  filter(An_Time == 1.00) # n = 159


# We will combine both the healthy and mam datasets to add onto the baseline data 
# other datasets will only need the anthro_MAM data as they do not have any healthy samples

# combine anthro_Healthy and anthro_MAM rows using rbind
all_anthro_2024 <- bind_rows(anthro_MAM_2024, anthro_Healthy_2024)

# Combine the anthro data to each dataset

# 1 - baseline_data
baseline_data <- baseline_data %>% 
  inner_join(all_anthro_2024, by = "Subject_ID")

# for the healthy samples, there is no Group value, therefore, add 0 
baseline_data <- baseline_data %>%
  mutate_all(~replace(., is.na(.), 0))

save(baseline_data, file ="data/processed/baseline_data.Rdata")


# 2 - refeed_data
refeed_data <- refeed_data %>% 
  inner_join(anthro_MAM_2024, by = "Subject_ID")

save(refeed_data, file ="data/processed/refeed_data.Rdata")


# 3 - recovery_data
recovery_data <- recovery_data %>% 
  inner_join(anthro_MAM_2024, by = "Subject_ID")

save(recovery_data, file ="data/processed/recovery_data.Rdata")


# 4 - group_data
group_data <- group_data %>% 
  inner_join(anthro_MAM_2024, by = "Subject_ID")

save(group_data, file ="data/processed/group_data.Rdata")


# PART 2C: Multiple Anthro measurements  ----------------------------------

# This part focuses on organising a dataset that contains multiple measurements
# of anthropometrics over a time period. For the purpopse of the PhD thesis
# we will filter the data to only contain upto and including 15months data. 

# the multiple measurement frequency is different for Healthy and MAM
# MAM were measured weekly and Healthy were measured at 3months after the enrolemnt measure 

# Therefore, as the analysis only focuses on up to 15months, we will filter out everything else

# reload the initial datasets 

# MAM
anthro_MAM <- read_excel("data/raw/anthro_MAM_2024.xlsx", 
                              col_types = c("text", "skip", "text", 
                                            "text", "skip", "numeric", "skip", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "numeric", "text")) 

# filter upto the last measurement for 15months
mam_anthro_measures <- anthro_MAM %>% 
  filter(An_Time <= 2.13) # 2.13 is the 15months mark for MAM

# rename as needed
mam_anthro_measures <- mam_anthro_measures %>% 
  rename(Subject_ID = SID_LCC,
         Sex = SEX)

save(mam_anthro_measures, file ="data/processed/mam_anthro_measures.Rdata")

# Healthy
anthro_Healthy <- read_excel("data/raw/anthro_Healthy_2024.xlsx", 
                                  col_types = c("text", "skip", "text", 
                                                "skip", "numeric", "skip", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "text")) 

# filter upto the last measurement for 15months
healthy_anthro_measures <- anthro_Healthy %>% 
  filter(An_Time <= 2.02) # 2.02 is the 15months mark for Healthy

# rename as needed
healthy_anthro_measures <- healthy_anthro_measures %>% 
  rename(Subject_ID = SID_LCC,
         Sex = SEX)

save(healthy_anthro_measures, file ="data/processed/healthy_anthro_measures.Rdata")

save(healthy_anthro_measures,mam_anthro_measures, file ="data/processed/anthro_measures.Rdata")


# PART 3: Taxa - Metadata  --------------------------------------------

# This section focuses on making datasets for each of the testing metadata sets
# we will use the orginal species_all_samples, genus_all_sample and phylum_all samples datasets
# as we created R data files, we will simply load them here


load("data/processed/species_all_samples.Rdata")
load("data/processed/genus_all_samples.Rdata")
load("data/processed/phylum_all_samples.Rdata")
load("data/processed/family_all_samples.Rdata")


# we turn these into objects to be used downstream

species_all_samples_edit <- species_all_samples %>% 
  rownames_to_column("Seq_ID")

genus_all_samples_edit <- genus_all_samples %>% 
  rownames_to_column("Seq_ID")

family_all_samples_edit <- family_all_samples %>% 
  rownames_to_column("Seq_ID")

phylum_all_samples_edit <- phylum_all_samples %>% 
  rownames_to_column("Seq_ID")

kingdom_all_samples_edit <- kingdom_all_samples %>% 
  rownames_to_column("Seq_ID")

order_all_samples_edit <- order_all_samples %>% 
  rownames_to_column("Seq_ID")

class_all_samples_edit <- class_all_samples %>% 
  rownames_to_column("Seq_ID")


# - baseline --------------------------------------------------------------

## SPECIES
# Prepare matching species data
baseline_species <- species_all_samples_edit %>% 
 semi_join(baseline_data, by ="Seq_ID")

save(baseline_species, file = "data/processed/baseline_species.Rdata") # save for future use


## GENUS
# Prepare matching species data
baseline_genus <- genus_all_samples_edit %>% 
  semi_join(baseline_data, by ="Seq_ID")

save(baseline_genus, file = "data/processed/baseline_genus.Rdata") # save for future use


## PHYLUM
# Prepare matching species data
baseline_phylum <- phylum_all_samples_edit %>% 
  semi_join(baseline_data, by ="Seq_ID")

save(baseline_phylum, file = "data/processed/baseline_phylum.Rdata") # save for future use



# - refeed --------------------------------------------------------------

## SPECIES
# Prepare matching species data
refeed_species <- species_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_species, file = "data/processed/refeed_species.Rdata") # save for future use


## GENUS
# Prepare matching species data
refeed_genus <- genus_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_genus, file = "data/processed/refeed_genus.Rdata") # save for future use


## PHYLUM
# Prepare matching species data
refeed_phylum <- phylum_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_phylum, file = "data/processed/refeed_phylum.Rdata") # save for future use


## KINGDOM
# Prepare matching species data
refeed_kingdom <- kingdom_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_kingdom, file = "data/processed/refeed_kingdom.Rdata") # save for future use

## ORDER
# Prepare matching species data
refeed_order <- order_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_order, file = "data/processed/refeed_order.Rdata") # save for future use

## CLASS
# Prepare matching species data
refeed_class <- class_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_class, file = "data/processed/refeed_class.Rdata") # save for future use

## FAMILY
# Prepare matching species data
refeed_family <- family_all_samples_edit %>% 
  semi_join(refeed_data, by ="Seq_ID")

save(refeed_family, file = "data/processed/refeed_family.Rdata") # save for future use


# - group --------------------------------------------------------------

## SPECIES
# Prepare matching species data
group_species <- species_all_samples_edit %>% 
  semi_join(group_data, by ="Seq_ID")

save(group_species, file = "data/processed/group_species.Rdata") # save for future use


## GENUS
# Prepare matching species data
group_genus <- genus_all_samples_edit %>% 
  semi_join(group_data, by ="Seq_ID")
  
save(group_genus, file = "data/processed/group_genus.Rdata") # save for future use


## PHYLUM
# Prepare matching species data
group_phylum <- phylum_all_samples_edit %>% 
  semi_join(group_data, by ="Seq_ID")

save(group_phylum, file = "data/processed/group_phylum.Rdata") # save for future use



# - recovery --------------------------------------------------------------

## SPECIES
# Prepare matching species data
recovery_species <- species_all_samples_edit %>% 
  semi_join(recovery_data, by ="Seq_ID")

save(recovery_species, file = "data/processed/recovery_species.Rdata") # save for future use


# conditions: mam v healthy -----------------------------------------------

# load the conditions specific metadata 

# mam
load("data/processed/mam_meta.Rdata")

mam_species <- species_all_samples_edit %>% 
  semi_join(mam_meta, by ="Seq_ID")

save(mam_species, file = "data/processed/mam_species.Rdata") # save for future use

# healthy
load("data/processed/healthy_meta.Rdata")

healthy_species <- species_all_samples_edit %>% 
  semi_join(healthy_meta, by ="Seq_ID")

save(healthy_species, file = "data/processed/healthy_species.Rdata") # save for future use

## GENUS
# Prepare matching species data
recovery_genus <- genus_all_samples_edit %>% 
  semi_join(recovery_data, by ="Seq_ID")

save(recovery_genus, file = "data/processed/recovery_genus.Rdata") # save for future use


## PHYLUM
# Prepare matching species data
recovery_phylum <- phylum_all_samples_edit %>% 
  semi_join(recovery_data, by ="Seq_ID")

save(recovery_phylum, file = "data/processed/recovery_phylum.Rdata") # save for future use


# PART 4: 10% Prevalence filter - Species  --------------------------------

# the prevalence filter is added to each individual dataset, ie., MAM v H, 12 v 15, Gr1 v Gr2, R v nR
# the total proportion will be different for each dataset as the total sample number changes 
# the species datasets have been made specific to each metadata set when preparing for Maaslin - refer to Maaslin speices

### 1: MAM v Healthy
load("data/processed/baseline_species.Rdata")

# we add a 10 % prevalence filter - this value can be changed by changing the filter line

# the data must be converted to the right format

# Check if baseline_species has row names, and if so, remove them
if (has_rownames(baseline_species)) {
  rownames(baseline_species) <- NULL
}

baseline_species <- baseline_species %>% 
  column_to_rownames("Seq_ID")

filtered_baseline_species <- baseline_species %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_baseline_species, file = "data/processed/filtered_baseline_species.Rdata") 

### 2: 12 v 15
load("data/refeed_species.Rdata")

# the data must be converted to the right format

# Check if baseline_species has row names, and if so, remove them
if (has_rownames(refeed_species)) {
  rownames(refeed_species) <- NULL
}

refeed_species <- refeed_species %>% 
  column_to_rownames("Seq_ID")

filtered_refeed_species <- refeed_species %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_refeed_species, file = "data/processed/filtered_refeed_species.Rdata") 

### 3: Gr1 v Gr2
load("data/processed/group_species.Rdata")

# the data must be converted to the right format

# Check if baseline_species has row names, and if so, remove them
if (has_rownames(group_species)) {
  rownames(group_species) <- NULL
}

group_species <- group_species %>% 
  column_to_rownames("Seq_ID")


filtered_group_species <- group_species %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_group_species, file = "data/processed/filtered_group_species.Rdata") 

### 4: R v nR
load("data/processed/recovery_species.Rdata")

# Check if baseline_species has row names, and if so, remove them
if (has_rownames(recovery_species)) {
  rownames(recovery_species) <- NULL
}

recovery_species <- recovery_species %>% 
  column_to_rownames("Seq_ID")


filtered_recovery_species <- recovery_species %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_recovery_species, file = "data/processed/filtered_recovery_species.Rdata") 

# PART 5: 10% Prevalence filter - Phylum ----------------------------------------------

# Add the filter for 10% prevalence filter for each dataset

### 1: MAM v Healthy
# load the feature for the filter
load("data/processed/baseline_phylum.Rdata")

# the data must be in the right format
if (has_rownames(baseline_phylum)) {
  rownames(baseline_phylum) <- NULL
}

baseline_phylum <- baseline_phylum %>% 
  column_to_rownames("Seq_ID")

filtered_baseline_phylum <- baseline_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_baseline_phylum, file = "data/processed/filtered_baseline_phylum.Rdata") 

### 2: 12 months v 15 months 
# load the feature for the filter
load("data/processed/refeed_phylum.Rdata")

# the data must be in the right format
if (has_rownames(refeed_phylum)) {
  rownames(refeed_phylum) <- NULL
}

refeed_phylum <- refeed_phylum %>% 
  column_to_rownames("Seq_ID")

filtered_refeed_phylum <- refeed_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_refeed_phylum, file = "data/processed/filtered_refeed_phylum.Rdata") 

### Group 1 v Group 2
# load the feature for the filter
load("data/group_phylum.Rdata")

# the data must be in the right format
if (has_rownames(group_phylum)) {
  rownames(group_phylum) <- NULL
}

group_phylum <- group_phylum %>% 
  column_to_rownames("Seq_ID")

filtered_group_phylum <- group_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_group_phylum, file = "data/processed/filtered_group_phylum.Rdata") 

### Recovered v Not recovered
# load the feature for the filter
load("data/processed/recovery_phylum.Rdata")

# the data must be in the right format
if (has_rownames(recovery_phylum)) {
  rownames(recovery_phylum) <- NULL
}

recovery_phylum <- recovery_phylum %>% 
  column_to_rownames("Seq_ID")

filtered_recovery_phylum <- recovery_phylum %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) # removing the last three rows that has count, sum and proportion

save(filtered_recovery_phylum, file = "data/processed/filtered_recovery_phylum.Rdata") 

# PART 6: Humann - Data Cleaning ------------------------------------------

pathways_sp_raw <- read_tsv("data/raw/m4efad_humann3_pathway_cpm_oct2023_stratified.tsv") #  pathways are stratified by species
pathways_raw <- read_tsv("data/raw/m4efad_humann3_pathway_cpm_oct2023_unstratified.tsv") # community pathway abundance

# DATA CLEANING
# Pathways

# wide format (puts the pathways as column names and the samples as rows)
pathways <- pathways_raw %>% 
  rename(Pathway = `# Pathway`) %>% 
  gather(Seq_ID, CPM, -Pathway) %>% 
  mutate(Seq_ID = str_remove(Seq_ID, "_Abundance")) %>% 
  spread(Pathway, CPM) %>% 
  column_to_rownames("Seq_ID") 

save(pathways, file = "data/processed/pathways.Rdata")


# long format, species-stratified (creates a table with a column for pathways, species, 
# sample and cpm of each pathway for each species)

pathways_sp_long <- pathways_sp_raw %>% 
  rename(Pathway_Species = `# Pathway`) %>% 
  gather(Seq_ID, CPM, -Pathway_Species) %>% 
  filter(CPM >0) %>% # keep only detectable rows
  separate(Pathway_Species, into = c("Pathway", "Species"), sep = "\\|") %>% 
  mutate(Seq_ID = str_remove(Seq_ID, "_Abundance"),
         Species = str_remove(Species, ".*\\."))

save(pathways_sp_long, file = "data/processed/pathways_sp_long.Rdata")

# long format (unstratified pathways made into a long format)
pathways_long <- pathways %>% 
  rownames_to_column("Seq_ID") %>% 
  gather(Pathway, CPM, -Seq_ID)

save(pathways_long, file = "data/processed/pathways_long.Rdata")

# Stratified pathways 

# load data
pathways_sp_raw <- read_tsv("data/raw/m4efad_humann3_pathway_cpm_oct2023_stratified.tsv") #  pathways are stratified by species

# wide format (puts the pathways as column names and the samples as rows)
# this is the format we generally use for Maaslin 
pathways_strat <- pathways_sp_raw %>% 
  rename(Pathway = `# Pathway`) %>% 
  gather(Seq_ID, CPM, -Pathway) %>% 
  mutate(Seq_ID = str_remove(Seq_ID, "_Abundance")) %>% 
  spread(Pathway, CPM) %>% 
  column_to_rownames("Seq_ID") 

save(pathways_strat, file = "data/processed/pathways_strat.Rdata")


# long format, species-stratified (creates a table with a column for pathways, species, 
# sample and cpm of each pathway for each species)

pathways_sp_long <- pathways_sp_raw %>% 
  rename(Pathway_Species = `# Pathway`) %>% 
  gather(Seq_ID, CPM, -Pathway_Species) %>% 
  filter(CPM >0) %>% # keep only detectable rows
  separate(Pathway_Species, into = c("Pathway", "Species"), sep = "\\|") %>% 
  mutate(Seq_ID = str_remove(Seq_ID, "_Abundance"),
         Species = str_remove(Species, ".*\\."))

save(pathways_sp_long, file = "data/pathways_sp_long.Rdata")



# Gene Tables 
# Note: Gene families table is a very large file
# to reduce its size, we will convert it to a long format and remove 0 quantification rows


genes <- read_tsv("data/raw/m4efad_humann3_genefam_cpm_oct2023_unstratified.tsv")

colnames(genes) <- gsub("_Abundance-RPKs", "", colnames(genes)) # remove the "_Abundance-RPKs" suffix from seq ID

colnames(genes) <- gsub("# Gene Family", "Gene_ID", colnames(genes)) # rename the "# Gene Family" column

genes_long <- melt(setDT(genes), id.vars = c("Gene_ID"), variable.name = "Seq_ID", value.name = "CPM") # convert to long format

setDT(genes_long) # convert to data table - more memory efficient 

genes_long <- genes_long[CPM >0] # remove non-detected rows

save(genes_long, file = "data/processed/genes_long.Rdata")




# PART 7: Controls - microbial and DNA ------------------------------------

# As part of the protocol, microbial community standard (extraction bias) and 
# DNA community standard (sequencing bias) were included with the samples

# here we separate these out for downstream analysis

# load the whole dataset file - this is the latest may 2024 file with all the sample information
load("data/processed/updated_meta.Rdata")

all_controls <- updated_meta %>% 
  filter(Condition == "control")

microbial_controls <- all_controls %>%
  filter(str_detect(Subject_ID, "^PC[1-5]$") | Subject_ID == "PC1_EH")

dna_controls <- all_controls %>%
  filter(str_detect(Subject_ID, "^D[1-5]$") | Subject_ID %in% c("D5_1", "D5_2"))

# save all the control related datasets within one Rdata object
save(all_controls, microbial_controls, dna_controls, file = "data/processed/controls.Rdata")

# this is the reference standard from the ZYMOBIOMICS Community Microbial and DNA standard 
reference_standard <- read_csv("data/raw/reference_standard.csv")

# the raw file contains unwanted rows, so we select them out
reference_standard <- reference_standard[1:10, ]

save(reference_standard, file = "data/processed/reference_standard.Rdata")

# PART 8: QC related datasets  --------------------------------------------

# The output from kneaddata contains all the age groups 
# we will not use the 3yr dataset in the course of the analyses, therefore 
# we will remove this

# load data 
# this file is the output from the kneaddata workflow summarised into a table 
kneaddata <- read_tsv("data/raw/m4efad_kneaddata_read_counts_april2024.tsv")

# rename the column "Sample" to match metadata
kneaddata_raw <- kneaddata %>% 
  rename(Seq_ID = Sample)


# - remove + create controls  ------------------------------------------------------

# load the controls Rdata file, this contains all the controls, we will remove them
# and add to a new data file
load("data/processed/controls.Rdata")

# controls removed dataset 
kneaddata_no_controls <- kneaddata_raw %>% 
  anti_join(all_controls, by = "Seq_ID")

# create separate controls only kneaddata
kneaddata_controls_only <- kneaddata_raw %>%
  inner_join(all_controls, by = "Seq_ID") 



# - remove 3yr samples  ---------------------------------------------------

# load the 3yr metadata 
load("data/processed/yr3_data.Rdata")

# remove these samples from the kneaddata
kneaddata_3yr_removed <- kneaddata_no_controls %>% 
 anti_join(yr3_data, by = "Seq_ID")



# combine metadata to samples  --------------------------------------------

kneaddata_metadata_combined <- kneaddata_3yr_removed %>% 
  inner_join(updated_meta, by = "Seq_ID")


# save the files  ---------------------------------------------------------

# raw = whole dataset incl 3yr samples + controls, 
# 3yr removed = 36 months samples + controls removed, 
# kneaddata_controls_only = kneaddata filtered for controls only,
# kneaddata_no_controls = controls removed but contains 3yr samples
# kneaddata_metadata_combined = 3yr removed samples + controls removed combined to metadata - this will be used mainly
save(kneaddata_raw, 
     kneaddata_3yr_removed, 
     kneaddata_controls_only, 
     kneaddata_no_controls, 
     kneaddata_metadata_combined, file = "data/processed/kneaddata.Rdata")



# Qc prior to sequencing - from Sequence provider -------------------------

# we receive qc reports from the sequence provider prior to library prep
# this gives us an indication of the DNA concentration in each sample and the quality
# of the DNA within each sample

# load the file with the data

# the data is stored in one excel files with multiple sheets for each batch

sheet1_qc_data_from_simon <- read_excel("data/qc_data_from_simon.xlsx", 
                                 sheet = "Sheet1", col_types = c("text", 
                                                                 "numeric", "numeric", "numeric", 
                                                                 "text"))


sheet2_qc_data_from_simon <- read_excel("data/qc_data_from_simon.xlsx", 
                                        sheet = "Sheet3", col_types = c("text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text"))


sheet3_qc_data_from_simon <- read_excel("data/qc_data_from_simon.xlsx", 
                                        sheet = "Sheet4", col_types = c("text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text"))


sheet4_qc_data_from_simon <- read_excel("data/qc_data_from_simon.xlsx", 
                                        sheet = "Sheet5", col_types = c("text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text"))


sheet5_qc_data_from_simon <- read_excel("data/qc_data_from_simon.xlsx", 
                                        sheet = "Sheet6", col_types = c("text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text"))


# Combine the datasets into one
qc_reports <- bind_rows(sheet1_qc_data_from_simon, 
                        sheet2_qc_data_from_simon, 
                        sheet3_qc_data_from_simon, 
                        sheet4_qc_data_from_simon, 
                        sheet5_qc_data_from_simon)

save(qc_reports, file = "data/processed/qc_reports.Rdata")



# DNA purity - lab extractions --------------------------------------------

# Note - script: lab_results_processing.R contains this component separately 

# # batch 1 ---------------------------------------------------------------

# NOTE: does not contain Seq_ID = C251 (DNA control)
# contains Seq_ID, Sample_Name and A260280 values 
QC_data_DNAonly_WL2022 <- read_excel("data/DNA-integrity-raw_files/QC_data_DNAonly_WL2022.xlsx", 
                                     sheet = "Sheet1")


b1 <- QC_data_DNAonly_WL2022 %>% 
  select(sequence_ID, sample_name, A260280, A260230) %>% 
  rename(Seq_ID = sequence_ID,
         Sample_Name = sample_name)

# Update Sample_Name for specified Seq_IDs
b1 <- b1 %>%
  left_join(updated_meta %>% select(Seq_ID, Sample_Name), by = "Seq_ID") %>%
  mutate(Sample_Name = if_else(Seq_ID %in% c("C247", "C248", "C249", "C250") & !is.na(Sample_Name.y), 
                               Sample_Name.y, Sample_Name.x)) %>%
  select(-Sample_Name.x, -Sample_Name.y)

# Add the row with Seq_ID = C251 from updated_meta
new_row <- updated_meta %>%
  filter(Seq_ID == "C251")

b1 <- bind_rows(b1, new_row)

b1 <- b1 %>% 
  select(Seq_ID, Sample_Name, A260280, A260230) %>% 
  mutate(across(everything(), ~replace_na(.x, 0)))



# # batch 2 ---------------------------------------------------------------

# two datasets need the sequencing IDs matched to the lab results 
# 1 - with all the DNA metrics and Sample name but no Seq_ID
batch3_DNA_all_new <- read_excel("data/DNA-integrity-raw_files/batch3_DNA_all_new.xlsx")

b2_A <- batch3_DNA_all_new %>%
  select(`Specimen ID`, A260280, A260230) %>%
  rename(Sample_Name = `Specimen ID`) %>%
  mutate(A260280 = as.numeric(A260280), 
         A260230 = as.numeric(A260230))

b2_A <- b2_A %>%
  mutate(Sample_Name = case_when(
    Sample_Name == "Negative control Standard" ~ "Negative control standard",
    Sample_Name == "Negative control Enhanced" ~ "Negative control enhanced",
    Sample_Name == "Micorbial Community Standard" ~ "Microbial Community Standard",
    TRUE ~ Sample_Name
  ))

# 2 - with seq IDS 
batch3_sequencing_data <- read_excel("data/DNA-integrity-raw_files/batch3_sequencing_data.xlsx")

b2_B <- batch3_sequencing_data %>% 
  select(Seq_ID, `Specimen ID`) %>%
  rename(Sample_Name = `Specimen ID`)

b2 <- full_join(b2_A, b2_B, by = "Sample_Name")

b2 <- b2 %>% 
  mutate(across(everything(), ~replace_na(.x, 0)))

# # batch 3  --------------------------------------------------------------

# NOTE: does not contain Seq_ID = J137 (DNA control)
batch4_DNA <- read_excel("data/DNA-integrity-raw_files/batch4_DNA.xlsx", 
                         sheet = "extraction_result_final")

b3 <- batch4_DNA %>% 
  select(Seq_ID, `Specimen ID`, A260280, A260230) %>%
  rename(Sample_Name = `Specimen ID`) 


# # batch 4 ---------------------------------------------------------------

# sheet with seq_ids
batch5_DNA <- read_excel("data/DNA-integrity-raw_files/batch5_DNA.xlsx", 
                         sheet = "sequencing")

b4_A <- batch5_DNA %>% 
  select(`Name on tube`, `Name on report`) %>% 
  rename(Seq_ID = `Name on tube`,
         Sample_Name = `Name on report`)

# sheet with SIDs to match the above Seq_IDs
batch5_DNA_2 <- read_excel("data/DNA-integrity-raw_files/batch5_DNA.xlsx", 
                           sheet = "all_results")

b4_B <- batch5_DNA_2 %>% 
  select(sample_name, `A260/280`, `A260/230`) %>% 
  rename(Sample_Name = sample_name,
         A260280 = `A260/280`,
         A260230 = `A260/230`)

b4_B <- b4_B %>%
  filter(!(Sample_Name == "negative control" & A260280 != 0))

b4_B <- b4_B %>%
  mutate(Sample_Name = case_when(
    Sample_Name == "positive control" ~ "ZymoBIOMICS Microbial Community Standard (P1)",
    Sample_Name == "negative control" ~ "Negative Control",
    TRUE ~ Sample_Name
  ))

b4 <- full_join(b4_A, b4_B, by = "Sample_Name")

b4 <- b4 %>% 
  mutate(across(everything(), ~replace_na(.x, 0)))

# # batch 5 ---------------------------------------------------------------

# contains only Seq_ID + A260280 ratios 
LEAP_Study_follow_up_samples_list_Shipment_6th_Batch_18_Dec_2023_MasterListResults <- read_excel("data/DNA-integrity-raw_files/LEAP Study_follow up samples list_Shipment 6th Batch_18 Dec 2023_MasterListResults.xlsx", 
                                                                                                 sheet = "ExtractionResults")

b5_A <- LEAP_Study_follow_up_samples_list_Shipment_6th_Batch_18_Dec_2023_MasterListResults %>% 
  select(`Tube ID`, `260/280`, `260/230`) %>% 
  rename(Seq_ID = `Tube ID`,
         A260280 = `260/280`,
         A260230 = `260/230`)


b5_A <- b5_A %>%
  mutate(Seq_ID = if_else(Seq_ID == "B1NTC", "NTC", Seq_ID))

b5_A <- b5_A %>%
  filter(Seq_ID != "B2NTC")

b5_A <- b5_A %>%
  filter(Seq_ID != "B3NTC")


# contains Sample_Name
batch_6_Leap_Study_Sample_Information_Hui_Hui_20240220 <- read_excel("data/DNA-integrity-raw_files/batch_6_Leap Study_Sample Information_Hui Hui_20240220.xlsx")

b5_B <- batch_6_Leap_Study_Sample_Information_Hui_Hui_20240220 %>% 
  select(`Name on tube`, `Name on report`) %>% 
  rename(Seq_ID = `Name on tube`,
         Sample_Name = `Name on report`)

b5_B <- b5_B %>%
  mutate(Seq_ID = if_else(Sample_Name == "Microbial Community DNA Std", "J674-1", Seq_ID))


b5 <- full_join(b5_A, b5_B, by = "Seq_ID")



# all batches -------------------------------------------------------------

all_batches <- bind_rows(b1, b2, b3, b4, b5)


## Load previously made qc report 
load("data/processed/qc_reports.Rdata")

## Load metadata file to match SID to Seq_ID
load("data/processed/updated_meta.Rdata")

# combine to updated_meta 
all_batches <- all_batches %>%
  left_join(updated_meta, by = "Seq_ID") %>%
  mutate(Sample_Name = coalesce(Sample_Name.y, Sample_Name.x)) %>%
  select(-Sample_Name.x, -Sample_Name.y)



# dna purity and quantity  ------------------------------------------------

dna_purity_quantity <- all_batches %>%
  left_join(qc_reports, by = "Seq_ID")

# save as Rdata 
save(dna_purity_quantity, file = "data/processed/dna_purity_quantity.Rdata")

# save as excel file
write_xlsx(dna_purity_quantity, "data/processed/dna_purity_quantity.xlsx")



# Datasets for Facet plots ------------------------------------------------

# this section creates datasets for facet plots using multiple variables

# use the original datasets for each without a filter 
load("data/processed/species.Rdata") # loading Metaphlan3 taxonomy data for species
load("data/processed/baseline_data.Rdata")
load("data/processed/refeed_data.Rdata")
load("data/processed/group_data.Rdata")
load("data/processed/recovery_data.Rdata")

# species RA table for the whole dataset - 461 observations (latest data from Oct 2023)
spp_whole_dataset <- species

## prepare metadata
# from each of the individual datasets, select the relevant columns

bd <- baseline_data %>% 
  select(Seq_ID, Subject_ID, Sample_Name, Age_months, Condition, Seq_batch)

rd <- refeed_data %>% 
  select(Seq_ID, Subject_ID, Sample_Name, Age_months, Condition, Seq_batch, Group)

gd <- group_data %>% 
  select(Seq_ID, Subject_ID, Sample_Name, Age_months, Condition, Seq_batch, Group)

recd <- recovery_data %>% 
  select(Seq_ID, Subject_ID, Sample_Name, Age_months, Condition, Seq_batch, Group.y, Recovery)



## clean up each dataset so they all match in the columns
# starting with the recovery dataset as it has the most columns


gd <- gd %>%
  mutate(Group = ifelse(Group == "Local RUSF (A)", 1, 
                        ifelse(Group == "ERUSF (B)", 2, 0)))

# add Group and Recovery columns for bd
bd <- bd %>%
  mutate(Recovery = ifelse(Condition == "Healthy", "Healthy control", "Pre-recovery"))

# add the correct Group value for the bd baseline dataset
bd <- left_join(bd, gd %>% select(Subject_ID, Group), by = "Subject_ID") %>%
  mutate(Group = coalesce(Group, 0))

# add the Recovery column for gd
gd <- left_join(gd, recd %>% select(Subject_ID, Recovery), by = "Subject_ID") 

# add Recovery column for rd
rd <- left_join(rd, recd %>% select(Subject_ID, Recovery), by = "Subject_ID") 

recd <- recd %>% 
  rename(Group = Group.y)

## combine all the four datasets

meta_whole_dataset <- rbind(bd, rd, gd, recd)
save(meta_whole_dataset, file = "data/processed/meta_whole_dataset.Rdata") # includes baseline healthy samples 

# in the whole dataset, create a column for health condition 
#(Baseline_healthy, Baseline_MAM, Recovered_MAM, Unrecovered_MAM, Baseline recovered and Baseline unrecovered)

meta_whole_dataset <- meta_whole_dataset %>%
  mutate(
    Health_condition = case_when(
      Age_months == 12 & Recovery == "Healthy control" ~ "Baseline healthy",
      Age_months == 12 & Recovery == "Pre-recovery" ~ "Baseline MAM",
      Age_months == 15 & Recovery == TRUE ~ "Infant Recovered",
      Age_months == 15 & Recovery == FALSE ~ "Infant Unrecovered",
      Age_months == 12 & Recovery == FALSE ~ "Baseline infant Unrecovered",
      Age_months == 12 & Recovery == TRUE ~ "Baseline infant Recovered"
    )
  )

# removing any NA rows with no data on recovery
meta_whole_dataset <- na.omit(meta_whole_dataset)

# filter out the healthy controls from the main dataset
MAM_meta_whole_dataset <- meta_whole_dataset %>% 
  filter(Condition == "MAM")

save(MAM_meta_whole_dataset, file = "data/processed/MAM_meta_whole_dataset.Rdata") # MAM only dataset- used in facets for refeed analysis


# Species Datasets for facet plots --------------------------------

load("data/processed/species.Rdata")
# species RA table for the whole dataset - 461 observations (latest data from Oct 2023)
spp_whole_dataset <- species

### Combine to the species dataset

# add a prevalence filter to the species dataset
filtered_spp_whole_dataset <- spp_whole_dataset %>%
  t() %>%
  as.data.frame() %>%
  mutate(count = rowSums(. > 0),  # Count only rows with values greater than zero
         Sum = ncol(.),
         proportion = count / Sum) %>% 
  filter(proportion >= 0.1) %>% 
  t() %>% #tranpose the axes back to original axes
  as.data.frame() %>% 
  rownames_to_column(var = "Seq_ID") %>% #add the Seq_ID column
  slice(1:(n() - 3)) 

# convert the rows in a column - if needed
# spp_whole_dataset <- spp_whole_dataset %>% 
  #column_to_rownames(var = "Seq_ID")

# join the spp dataset to the MAM whole dataset by Seq_ID

#wide format
# all including Healthy
wide_spp_whole <- filtered_spp_whole_dataset %>% 
  inner_join(meta_whole_dataset, by = "Seq_ID")

# MAM only
wide_spp_MAM_whole <- filtered_spp_whole_dataset %>% 
  inner_join(MAM_meta_whole_dataset, by = "Seq_ID")

# long format - all incl healthy controls
# long_spp_whole <- filtered_spp_whole_dataset %>% 
#  gather(Taxa, RA, -Seq_ID) %>%
#  group_by(Taxa) %>% 
#  inner_join(whole_data_arm) %>%
#  arrange(desc(Taxa))

#long format - MAM only
long_spp_MAM_whole <- filtered_spp_whole_dataset %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(MAM_meta_whole_dataset) %>%
  arrange(desc(Taxa)) %>%
  filter(Recovery != "Pre-recovery") # removing this as by default 12 months is pre recovery

save(long_spp_MAM_whole, file = "data/processed/long_spp_MAM_whole.Rdata") # used in facet plots 

load('data/processed/long_spp_MAM_whole.Rdata')


























