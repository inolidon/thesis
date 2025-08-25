## code summary:
# find which metaphlan4 database result has the highest representation of WAAFLE directed HGT species
# quantify the number of species prevalent in each sample, across the sample subsets
# investigate the correlation between number of HGT interactions and mean species relative abundance (RA)/prevalence
# compare HGT interactions of species by role classification
# compare mean RA/prevalence of species by role classification

## load libraries

library(tidyverse)
library(ggrepel)
library(ggbeeswarm)
library(patchwork)
library(dunn.test)

## set working directory

setwd("Z:/hgt_networks/clean_outputs/")

## load data

# species-specific HGT events (formatted WAAFLE results)
load("spp_dir_hgts_uc.RData")
load("spp_dir_hgts_gb.RData")

# full metadata for each trial, note: sample IDs must match metaphlan sample IDs (need to remove underscores from metaphlan samples IDs during processing)
load("formatted_metadata_uc.RData")
load("formatted_metadata_gb.RData")

# species network roles
load("joined_spp_roles.RData")

# formatted metadata, note: sample IDs must match metaphlan sample IDs (need to remove underscores from metaphlan samples IDs during processing)
load("sample_subsets.RData")

## load palettes

trial_palette <- c("FOCUS" = "#D55E00", "Gut Bugs" = "#0072B2")
role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")

## load functions

# get distinct species in metaphlan4 GTDB or SGB output
get_distinct_m4_spp <- function(formatted_metaphlan){
  ## formatted_metaphlan: path to metaphlan data, data contains the object 'species' where rownames = sample IDs, colnames = species (str)
  load(formatted_metaphlan)
  distinct_spp <- species %>%
    t() %>% # transpose data so rownames = species
    as.data.frame() %>% # convert from matrix to dataframe
    rownames_to_column("Species") %>%
    select(Species) %>%
    distinct() %>% # get distinct species
    mutate(Species = na_if(Species, "")) %>% # assign NA to missing species
    filter(!is.na(Species)) %>% # remove NA rows
    # note: both GTDB and SGB formatted species begin with 's__', SGB also have underscores within name
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
  distinct_spp
}

# get distinct species in WAAFLE directional HGT between species
get_all_dir_hgt_spp <- function(spp_hgts_formatted_samples){
  # spp_hgts_formatted_samples: directed HGT events between species, must contain columns 'CLADE_B' and 'CLADE_A' (df)
  # isolate donor species
  donor_spp <- spp_hgts_formatted_samples %>%
    select(CLADE_B) %>%
    rename(Species = CLADE_B)
  # isolate recipient species
  recip_spp <- spp_hgts_formatted_samples %>%
    select(CLADE_A) %>%
    rename(Species = CLADE_A)
  # join donor and recipient species
  distinct_spp <- donor_spp %>%
    bind_rows(recip_spp) %>%
    distinct() %>%
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
  distinct_spp
}

# find which WAAFLE directed HGT species are found in the metaphlan4 GTDB or SGB outputs
waafle_spp_in_m4 <- function(dir_hgt_spp, gtdb_spp, sgb_spp){
  ## dir_hgt_spp: distinct species involved in directed HGT from WAAFLE result (df)
  ## gtdb_spp: distinct species from metaphlan4 GTDB result (df)
  ## sgb_spp: distinct species from metaphlan4 SGB result (df)
  dir_hgt_spp %>%
    mutate(in_m4_gtdb = if_else(Species %in% gtdb_spp$Species, "Yes", "No")) %>% # check for waafle species in metaphlan4 gtdb species
    mutate(in_m4_sgb = if_else(Species %in% sgb_spp$Species, "Yes", "No")) # check for waafle species in metaphlan4 sgb species
}

# format the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs for plotting
format_waafle_m4_spp_ppn_data <- function(waafle_m4_spp_ppn, trial_name){
  ## waafle_m4_spp_ppn: percentage of WAAFLE species in m4 outputs by database used, contains columns 'in_m4_sgb' and 'in_m4_gtdb' (df)
  ## trial_name: name of trial (str)
  waafle_m4_spp_ppn %>%
    rename(SGB = in_m4_sgb,
           GTDB = in_m4_gtdb) %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rename(ppn_waafle_matches = V1,
           m4_database = rowname) %>% # MetaPhlAn database
    mutate(trial = trial_name)
}

# quantify the number of species prevalent in each sample, across the sample subsets
quantify_metaphlan_spp <- function(formatted_metaphlan, metadata){
  ## formatted_metaphlan: path to metaphlan data, where rownames = sample IDs, colnames = taxa (str)
  ## metadata: formatted metadata for sample subset, contains column 'Sample_ID' (df)
  load(formatted_metaphlan) # load formatted metaphlan taxonomy profiles for each trial
  species %>% # load taxonomy relative abundance profiles for each sample
    rownames_to_column("Sample_ID") %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(metadata %>% select(Sample_ID), by = "Sample_ID") %>% # join with subset sample data
    column_to_rownames("Sample_ID") %>% # no longer need sample ID column
    select(where(~ sum(.) > 0)) %>% # filter out columns (taxa) with sum(RA == 0) across all samples
    rowwise() %>% # for each sample..
    mutate(n_spp = sum(c_across(everything()) > 0)) %>% # count the number of species with RA > 0
    ungroup() %>%
    pull(n_spp)
}

# find the mean relative abundance of species with a certain role, across samples of a cohort
spp_roles_mean_RA <- function(metaphlan_RData, sample_subset, trial_cohort_name, joined_spp_roles){
  ## metaphlan_RData: path to processed metaphlan RData object for each trial containing 'metaphlan_long' (str)
  ## sample_subset: metadata for subset cohort samples (df)
  ## trial_cohort_name: trial cohort concatenated name (str)
  ## joined_spp_roles: joined species role data for all trial cohorts (df)
  load(metaphlan_RData) # load RData object
  spp_roles_mean_RA <- metaphlan_long %>%
    filter(grepl("s__", Taxa)) %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(sample_subset, by = "Sample_ID") %>% # join with subset metadata
    rename(Species = Taxa) %>%
    group_by(Group, Species) %>% # for each species in each group..
    summarise(mean_RA = mean(RA)) %>% # calculate mean RA of each species
    ungroup() %>%
    mutate(Trial_cohort = trial_cohort_name) %>%
    inner_join(joined_spp_roles, by = c("Species", "Trial_cohort")) %>% # join with species network role data
    mutate(Total_edges = Donor_count + Recipient_count) # total HGT interactions for each species
}

# find the mean sample prevalence of species with a certain role, across samples of a cohort
spp_roles_mean_prev <- function(metaphlan_RData, sample_subset, trial_cohort_name, joined_spp_roles){
  ## metaphlan_RData: path to processed metaphlan RData object for each trial containing 'metaphlan_long' (str)
  ## sample_subset: metadata for subset cohort samples (df)
  ## trial_cohort_name: trial cohort concatenated name (str)
  ## joined_spp_roles: joined species role data for all trial cohorts (df)
  n_samples <- nrow(sample_subset) # quantify number of samples in sample subset
  load(metaphlan_RData) # load RData object
  spp_roles_mean_prev <- metaphlan_long %>%
    filter(grepl("s__", Taxa)) %>%
    mutate(Sample_ID = str_remove(Sample_ID, "_")) %>% # remove underscores from metaphlan sample ID to match metadata
    right_join(sample_subset, by = "Sample_ID") %>% # join with subset metadata
    rename(Species = Taxa) %>%
    group_by(Group, Species) %>% # for each species in each group..
    summarise(Sample_count = sum(RA > 0, na.rm = TRUE)) %>% # quantify the number of samples with non-0 abundance
    ungroup() %>%
    mutate(mean_prev = Sample_count/n_samples) %>% # species mean sample prevalence
    select(-Sample_count) %>%
    mutate(Trial_cohort = trial_cohort_name) %>%
    inner_join(joined_spp_roles, by = c("Species", "Trial_cohort")) %>% # join with species network role data
    mutate(Total_edges = Donor_count + Recipient_count) # total HGT interactions for each species
}

## running

# find which metaphlan4 database has the most HGT species -----------------

# get distinct species in metaphlan4 GTDB or SGB output
distinct_m4_gtdb_spp_uc <- get_distinct_m4_spp("metaphlan4_gtdb_uc.Rdata") # m4 gtdb
distinct_m4_gtdb_spp_gb <- get_distinct_m4_spp("metaphlan4_gtdb_gb.Rdata") 
distinct_m4_sgb_spp_uc <- get_distinct_m4_spp("metaphlan4_sgb_uc.Rdata") # m4 sgb
distinct_m4_sgb_spp_gb <- get_distinct_m4_spp("metaphlan4_sgb_gb.Rdata") 

# subset HGT events for those that are directed and between species only
dir_spp_hgts_uc <- spp_dir_hgts_uc %>% filter(Directional == "Yes")
dir_spp_hgts_gb <- spp_dir_hgts_gb %>% filter(Directional == "Yes")

# get distinct species in WAAFLE directional HGT between species
distinct_dir_hgt_spp_uc <- get_all_dir_hgt_spp(dir_spp_hgts_uc)
distinct_dir_hgt_spp_gb <- get_all_dir_hgt_spp(dir_spp_hgts_gb)

# find the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
# note: ok to use full metaphlan profiles (not sample subsets) as only interested in proportion of WAAFLE representation here
waafle_m4_spp_matches_uc <- waafle_spp_in_m4(distinct_dir_hgt_spp_uc, gtdb_spp = distinct_m4_gtdb_spp_uc, sgb_spp = distinct_m4_sgb_spp_uc)
waafle_m4_spp_matches_gb <- waafle_spp_in_m4(distinct_dir_hgt_spp_gb, gtdb_spp = distinct_m4_gtdb_spp_gb, sgb_spp = distinct_m4_sgb_spp_gb)

# quantify the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
waafle_m4_spp_ppn_uc <- waafle_m4_spp_matches_uc %>% summarise(across(-Species, ~ mean(. == "Yes")))
waafle_m4_spp_ppn_gb <- waafle_m4_spp_matches_gb %>% summarise(across(-Species, ~ mean(. == "Yes")))

# format the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs for plotting
waafle_m4_spp_ppn_uc_formatted <- format_waafle_m4_spp_ppn_data(waafle_m4_spp_ppn_uc, "FOCUS")
waafle_m4_spp_ppn_gb_formatted <- format_waafle_m4_spp_ppn_data(waafle_m4_spp_ppn_gb, "Gut Bugs")

# plot the proportion of exact matching species from WAAFLE output in metaphlan4 GTDB and SGB outputs
waafle_m4_spp_ppn_plot <- waafle_m4_spp_ppn_uc_formatted %>% 
  bind_rows(waafle_m4_spp_ppn_gb_formatted) %>% # bind data for all trials
  ggplot(aes(x = m4_database, y = ppn_waafle_matches)) +
  geom_bar(stat = "identity") +
  facet_grid(~trial) +
  xlab("MetaPhlAn database") + ylab("Proportion of WAAFLE HGT species represented") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "../clean_figures/waafle_m4_spp_ppn.svg", 
#       plot = waafle_m4_spp_ppn_plot, units = "cm", width=12, height=12)

# note: WAAFLE HGT species are better represented in SGB database results for both trials, use metaphlan4 SGB database results from here on

# quantify the number of species prevalent in each sample, across the sample subsets
uc_n_species <- quantify_metaphlan_spp("metaphlan4_sgb_uc.Rdata", metadata_uc) # FOCUS
mean(uc_n_species) # 80.6
sd(uc_n_species) # 41.4

gb_n_species <- quantify_metaphlan_spp("metaphlan4_sgb_gb.Rdata", metadata_gb) # Gut Bugs
mean(gb_n_species) # 199
sd(gb_n_species) # 62.6


# investigate mean RA of different network roles --------------------------

# find the mean relative abundance of species with a certain role, across samples of each cohort
spp_roles_mean_RA_uc_donor <- spp_roles_mean_RA("metaphlan4_sgb_uc.Rdata", donor_samples_uc, "FOCUS Donor", joined_spp_roles)
spp_roles_mean_RA_uc_pre_fmt <- spp_roles_mean_RA("metaphlan4_sgb_uc.Rdata", pre_fmt_samples_uc, "FOCUS Pre-FMT", joined_spp_roles)
spp_roles_mean_RA_uc_post_fmt <- spp_roles_mean_RA("metaphlan4_sgb_uc.Rdata", post_fmt_samples_uc, "FOCUS Post-FMT", joined_spp_roles)
spp_roles_mean_RA_uc_pre_placebo <- spp_roles_mean_RA("metaphlan4_sgb_uc.Rdata", pre_placebo_samples_uc, "FOCUS Pre-placebo", joined_spp_roles)
spp_roles_mean_RA_uc_post_placebo <- spp_roles_mean_RA("metaphlan4_sgb_uc.Rdata", post_placebo_samples_uc, "FOCUS Post-placebo", joined_spp_roles)

spp_roles_mean_RA_gb_donor <- spp_roles_mean_RA("metaphlan4_sgb_gb.Rdata", donor_samples_gb, "Gut Bugs Donor", joined_spp_roles)
spp_roles_mean_RA_gb_pre_fmt <- spp_roles_mean_RA("metaphlan4_sgb_gb.Rdata", pre_fmt_samples_gb, "Gut Bugs Pre-FMT", joined_spp_roles)
spp_roles_mean_RA_gb_post_fmt <- spp_roles_mean_RA("metaphlan4_sgb_gb.Rdata", post_fmt_samples_gb, "Gut Bugs Post-FMT", joined_spp_roles)
spp_roles_mean_RA_gb_pre_placebo <- spp_roles_mean_RA("metaphlan4_sgb_gb.Rdata", pre_placebo_samples_gb, "Gut Bugs Pre-placebo", joined_spp_roles)
spp_roles_mean_RA_gb_post_placebo <- spp_roles_mean_RA("metaphlan4_sgb_gb.Rdata", post_placebo_samples_gb, "Gut Bugs Post-placebo", joined_spp_roles)

# join data for all cohorts
joined_spp_roles_mean_RA <- bind_rows(spp_roles_mean_RA_uc_donor, spp_roles_mean_RA_uc_pre_fmt, spp_roles_mean_RA_uc_post_fmt, 
                                      spp_roles_mean_RA_uc_pre_placebo, spp_roles_mean_RA_uc_post_placebo, spp_roles_mean_RA_gb_donor,
                                      spp_roles_mean_RA_gb_pre_fmt, spp_roles_mean_RA_gb_post_fmt, spp_roles_mean_RA_gb_pre_placebo,
                                      spp_roles_mean_RA_gb_post_placebo)

joined_spp_roles_mean_RA$Type <- factor(joined_spp_roles_mean_RA$Type, levels = c("Source", "Conduit", "Sink"))

# get mean RA of species with each role across all cohorts
joined_spp_roles_mean_RA_all <- joined_spp_roles_mean_RA %>%
  select(Species, Type, mean_RA, Total_edges) %>%
  group_by(Species, Type) %>%
  summarise(mean_RA = mean(mean_RA), # calculate mean RA of each species in each role across all cohorts
            Total_edges = sum(Total_edges)) %>% 
  ungroup()

# plot correlation between number of HGT interactions and mean species relative abundance
spp_RA_edges_corr <- joined_spp_roles_mean_RA_all %>%
  ggplot(aes(x = mean_RA, y = Total_edges, fill = Type), colour = "black") +
  geom_label_repel(data = joined_spp_roles_mean_RA_all %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale" | Species == "Collinsella aerofaciens"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label conserved sinks
                   nudge_x = -0.2,
                   nudge_y = 0.35,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  geom_point(shape = 21, size = 3.5) +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300)) + # log y axis to separate points more
  xlab("Mean relative abundance") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "../clean_figures/spp_RA_edges_corr.svg", 
#       plot = spp_RA_edges_corr, units = "cm", width=10, height=10)


# investigate mean prevalence of different network roles ------------------

# find the mean prevalence of species with a certain role, across samples of each cohort
spp_roles_mean_prev_uc_donor <- spp_roles_mean_prev("metaphlan4_sgb_uc.Rdata", donor_samples_uc, "FOCUS Donor", joined_spp_roles)
spp_roles_mean_prev_uc_pre_fmt <- spp_roles_mean_prev("metaphlan4_sgb_uc.Rdata", pre_fmt_samples_uc, "FOCUS Pre-FMT", joined_spp_roles)
spp_roles_mean_prev_uc_post_fmt <- spp_roles_mean_prev("metaphlan4_sgb_uc.Rdata", post_fmt_samples_uc, "FOCUS Post-FMT", joined_spp_roles)
spp_roles_mean_prev_uc_pre_placebo <- spp_roles_mean_prev("metaphlan4_sgb_uc.Rdata", pre_placebo_samples_uc, "FOCUS Pre-placebo", joined_spp_roles)
spp_roles_mean_prev_uc_post_placebo <- spp_roles_mean_prev("metaphlan4_sgb_uc.Rdata", post_placebo_samples_uc, "FOCUS Post-placebo", joined_spp_roles)

spp_roles_mean_prev_gb_donor <- spp_roles_mean_prev("metaphlan4_sgb_gb.Rdata", donor_samples_gb, "Gut Bugs Donor", joined_spp_roles)
spp_roles_mean_prev_gb_pre_fmt <- spp_roles_mean_prev("metaphlan4_sgb_gb.Rdata", pre_fmt_samples_gb, "Gut Bugs Pre-FMT", joined_spp_roles)
spp_roles_mean_prev_gb_post_fmt <- spp_roles_mean_prev("metaphlan4_sgb_gb.Rdata", post_fmt_samples_gb, "Gut Bugs Post-FMT", joined_spp_roles)
spp_roles_mean_prev_gb_pre_placebo <- spp_roles_mean_prev("metaphlan4_sgb_gb.Rdata", pre_placebo_samples_gb, "Gut Bugs Pre-placebo", joined_spp_roles)
spp_roles_mean_prev_gb_post_placebo <- spp_roles_mean_prev("metaphlan4_sgb_gb.Rdata", post_placebo_samples_gb, "Gut Bugs Post-placebo", joined_spp_roles)

# join data for all cohorts
joined_spp_roles_mean_prev <- bind_rows(spp_roles_mean_prev_uc_donor, spp_roles_mean_prev_uc_pre_fmt, spp_roles_mean_prev_uc_post_fmt, 
                                        spp_roles_mean_prev_uc_pre_placebo, spp_roles_mean_prev_uc_post_placebo, spp_roles_mean_prev_gb_donor,
                                        spp_roles_mean_prev_gb_pre_fmt, spp_roles_mean_prev_gb_post_fmt, spp_roles_mean_prev_gb_pre_placebo,
                                        spp_roles_mean_prev_gb_post_placebo)

joined_spp_roles_mean_prev$Type <- factor(joined_spp_roles_mean_prev$Type, levels = c("Source", "Conduit", "Sink"))

# get mean prevalence of species with each role across all cohorts
joined_spp_roles_mean_prev_all <- joined_spp_roles_mean_prev %>%
  select(Species, Type, mean_prev, Total_edges) %>%
  group_by(Species, Type) %>%
  summarise(mean_prev = mean(mean_prev), # calculate mean prevalence of each species in each role across all cohorts
            Total_edges = sum(Total_edges)) %>% # calculate total interactions of each species in each role across all cohorts
  ungroup()

# plot correlation between number of HGT interactions and mean species prevalence
spp_prev_edges_corr <- joined_spp_roles_mean_prev_all %>%
  ggplot(aes(x = mean_prev, y = Total_edges), colour = "black") +
  geom_label_repel(data = joined_spp_roles_mean_prev_all %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale" | Species == "Collinsella aerofaciens"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label conserved sinks
                   nudge_x = -0.2,
                   nudge_y = 0.25,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  geom_point(aes(fill = Type), shape = 21, size = 3.5) +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300)) + # log y axis to separate points more
  xlab("Mean prevalence") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "../clean_figures/spp_prev_edges_corr.svg", 
#       plot = spp_prev_edges_corr, units = "cm", width=10, height=10)

# remake figure just for legend
spp_edges_corr_legend <- joined_spp_roles_mean_prev_all %>%
  ggplot(aes(x = mean_prev, y = Total_edges), colour = "black") +
  geom_point(aes(fill = Type), shape = 21, size = 3.5) +
  geom_label_repel(data = joined_spp_roles_mean_prev_all %>% 
                     mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
                            Species = str_replace_all(Species, "_", " ")) %>% # remove underscores from species names
                     filter(Species == "Faecalibacterium prausnitzii" | Species == "Eubacterium rectale"), 
                   aes(label = Species, fill = Type), min.segment.length = 0, # label outliers
                   nudge_x = -0.2,
                   nudge_y = 0.25,
                   size = 2.5,
                   show.legend = FALSE) + # still show all points as point shape and fill in legend
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  scale_y_log10() + # log y axis to separate points more
  xlab("Mean species prevalence") + ylab("Total HGT interactions") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "../clean_figures/spp_edges_corr_legend.svg", 
#       plot = spp_edges_corr_legend, units = "cm", width=10, height=10)


# investigate association between role and HGT interactions ---------------

joined_spp_roles_HGT_int_all <- joined_spp_roles_mean_RA_all %>%
  select(Species, Type, Total_edges) # not interested in RA or prevalence, just HGT interactions for each species in each role

spp_roles_HGT_int_all <- joined_spp_roles_HGT_int_all %>%
  ggplot(aes(x = Type, y = Total_edges, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  scale_y_log10(breaks = c(3, 10, 30, 100, 300),
                expand = expansion(mult = c(0.05, 0.15))) + # log y axis, add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab("Species HGT network role") + ylab("Total HGT interactions") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
ggsave(file = "../clean_figures/spp_roles_HGT_int_all.svg", 
       plot = spp_roles_HGT_int_all, units = "cm", width=12, height=12)

# assess statistical significance of association between HGT interactions and network role

# normality testing
joined_spp_roles_HGT_int_all %>% filter(Type == "Source") %>% pull(Total_edges) %>% shapiro.test() # p = 1.385e-09

# statistically compare HGT interactions of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(Total_edges ~ Type, data = joined_spp_roles_HGT_int_all) # p = 0.0006396 ***

# post-hoc testing to see which roles differ
dunn_test_int <- dunn.test(joined_spp_roles_HGT_int_all$Total_edges, joined_spp_roles_HGT_int_all$Type, method = "bonferroni")
role_comparisons_int <- dunn_test_int$comparisons
role_pvals_int <- dunn_test_int$P.adj

sig_role_comparisons_int <- role_comparisons_int[role_pvals_int < 0.05]
sig_role_pvals_int <- role_pvals_int[role_pvals_int < 0.05]
data.frame(Comparison = sig_role_comparisons_int, P.adj = sig_role_pvals_int)
# Comparison        P.adj
# Conduit - Sink 0.000521459 ***
# Conduit - Source 0.034266151 *


# investigate association between role and mean RA/prev -------------------

# plot mean RA and prevalence of species in each network role across all cohorts
# note: plots need to be made separately due to only one requiring log-transformation of y axis
spp_roles_mean_RA_all <- joined_spp_roles_mean_RA_all %>%
  rename(Value = mean_RA) %>%
  mutate(metric = "Mean relative abundance") %>%
  #filter(Value != 0) %>% # remove species with mean RA of 0 as this cannot be log transformed
  mutate(Value = if_else(Value == 0, Value+1e-6, Value)) %>% # add a pseudo count to mean RA of 0 as this cannot be log transformed
  ggplot(aes(x = metric, y = Value, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(metric~., scales = "free", space = "free") +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.15))) + # log y axis, add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(), # remove x axis labels
        axis.ticks.x = element_blank(), # remove x axis ticks
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

spp_roles_mean_prev_all <- joined_spp_roles_mean_prev_all %>%
  rename(Value = mean_prev) %>%
  mutate(metric = "Mean prevalence") %>%
  ggplot(aes(x = metric, y = Value, fill = Type)) +
  geom_quasirandom(shape = 21, dodge.width = 0.75, size = 2) +
  geom_boxplot(alpha = 0.85, outlier.colour = NA, position = position_dodge2(preserve = "single")) +
  facet_grid(metric~., scales = "free", space = "free") +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), 
                     expand = expansion(mult = c(0.05, 0.15))) + # add 15% extra space above points to add significance bars if needed
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab(NULL) + 
  theme_bw() +
  theme(legend.position = "none",
        axis.text.x = element_blank(), # remove x axis labels
        axis.ticks.x = element_blank(), # remove x axis ticks
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))

# join plots together
spp_roles_RA_prev_all <- spp_roles_mean_RA_all / spp_roles_mean_prev_all
#ggsave(file = "../clean_figures/spp_roles_RA_prev_all.svg", 
#       plot = spp_roles_RA_prev_all, units = "cm", width=11, height=20)

# assess statistical significance of association between species mean RA and network role

# normality testing
joined_spp_roles_mean_RA_all %>% filter(Type == "Source") %>% pull(mean_RA) %>% shapiro.test() # p = 9.15e-09

# statistically compare mean RA of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(mean_RA ~ Type, data = joined_spp_roles_mean_RA_all) # p = 7.025e-06 ***

# post-hoc testing to see which roles differ
dunn_test_RA <- dunn.test(joined_spp_roles_mean_RA_all$mean_RA, joined_spp_roles_mean_RA_all$Type, method = "bonferroni")
role_comparisons_RA <- dunn_test_RA$comparisons
role_pvals_RA <- dunn_test_RA$P.adj

sig_role_comparisons_RA <- role_comparisons_RA[role_pvals_RA < 0.05]
sig_role_pvals_RA <- role_pvals_RA[role_pvals_RA < 0.05]
data.frame(Comparison = sig_role_comparisons_RA, P.adj = sig_role_pvals_RA)
# Comparison        P.adj
# Conduit - Sink 1.208653e-05 ***
# Sink - Source 2.699521e-05 ***

# assess statistical significance of association between species mean prevalence and network role

# normality testing
joined_spp_roles_mean_prev_all %>% filter(Type == "Source") %>% pull(mean_prev) %>% shapiro.test() # p = 0.001865

# statistically compare mean prevalence of species across all trial cohorts by role classification, using non parametric tests
kruskal.test(mean_prev ~ Type, data = joined_spp_roles_mean_prev_all) # p = 3.234e-06 ***

# post-hoc testing to see which roles differ
dunn_test_prev <- dunn.test(joined_spp_roles_mean_prev_all$mean_prev, joined_spp_roles_mean_prev_all$Type, method = "bonferroni")
role_comparisons_prev <- dunn_test_prev$comparisons
role_pvals_prev <- dunn_test_prev$P.adj

sig_role_comparisons_prev <- role_comparisons_prev[role_pvals_prev < 0.05]
sig_role_pvals_prev <- role_pvals_prev[role_pvals_prev < 0.05]
data.frame(Comparison = sig_role_comparisons_prev, P.adj = sig_role_pvals_prev)
# Comparison        P.adj
# Conduit - Sink 1.349431e-06 ***
# Sink - Source 1.608037e-04 ***
