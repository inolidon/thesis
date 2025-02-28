## code summary:
# for each group (MAM v Healthy), plot the proportion of all HGT events, and those between species, with direction
# for each group, plot the proportion of species-specific HGT events occurring within each taxonomic level
# for the whole m4efad study - this compare against other studies

## load libraries

library(tidyverse)

## set working directory

# setwd("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/data/anna-hgt")

## load functions

# define directional HGT events
define_dir_hgts <- function(joined_hgt_df, metadata, trial_name){
  ## joined_hgt_df: single dataframe containing all contigs with HGT events detected by waafle, across samples (df)
  ## metadata: subset of metadata containing all samples for investigation (df)
  ## trial_name: trial name (str)
  joined_hgt_df %>%
    inner_join(metadata, by = "Sample_ID") %>% # join with subset metadata by extracted sample ID
    mutate(Directional = if_else(DIRECTION == "B>A", "Yes", "No"), # identify directional HGTs (designated A?B or B>A)
           Type = "All HGTs", # define HGT type
           Trial = trial_name) # define trial name (study name)
}

# subset directional HGT events between species
subset_spp_dir_hgts <- function(dir_hgts, trial_name){
  ## dir_hgts: output of define_dir_hgts() (df)
  ## trial_name: trial name (str)
  dir_hgts %>%
    filter(grepl("^s__", CLADE_A) & grepl("^s__", CLADE_B)) %>% # get only events between species
    mutate(Type = "Species-specific HGTs", # define HGT type
           Trial = trial_name) # define trial name
}

## load data

# metadata on local PC
# load("data/processed/clean_outputs/metadata_m4efad.RData")

# on github repo
load("hgt_analyses/data/metadata_m4efad.RData")

# waafle
# load("data/processed/hgt_df_m4efad.RData")

# on github repo
load("hgt_analyses/data/hgt_df_m4efad.RData")


## load palettes

trial_palette <- c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")
trial_palette2 <- c("M4EFaD" = "darkkhaki")

## running

# investigate the proportion of HGT events with direction -----------------

# define directional HGT events
# whole dataset
dir_hgts_m4efad <- define_dir_hgts(hgt_df_m4efad, formatted_metadata_m4efad, trial_name = "M4EFaD")

# subsets/groups considered "trials" - note: only using the 12 months baseline samples 4/Feb/2025 
dir_hgts_mam12 <- define_dir_hgts(hgt_df_m4efad, MAM_12_samples, trial_name = "MAM_12")
dir_hgts_healthy12 <- define_dir_hgts(hgt_df_m4efad, Healthy_12_samples, trial_name = "Healthy_12")

# subset directional HGT events between species
# whole dataset
spp_dir_hgts_m4efad <- subset_spp_dir_hgts(dir_hgts_m4efad, trial_name = "M4EFaD")

# subsets/groups considered "trials" - note: only using 12 months baseline samples
spp_dir_hgts_mam12 <- subset_spp_dir_hgts(dir_hgts_mam12, trial_name = "MAM_12")
spp_dir_hgts_healthy12 <- subset_spp_dir_hgts(dir_hgts_healthy12, trial_name = "Healthy_12")


# whole m4efad data -------------------------------------------------------

# plot the proportion of all HGT events, and those between species, with direction
ppn_hgt_directional_plot2 <- dir_hgts_m4efad %>% 
  bind_rows(spp_dir_hgts_m4efad) %>% # bind data for all trials
  select(Type, Trial, Directional) %>%
  group_by(Type, Trial, Directional) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = Trial, y = n, fill = Directional)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = n), position = position_fill(vjust=0.5), colour="white") +
  scale_fill_manual(values = c("#999999", "#4ABE56"), name = "Directed") +
  facet_grid(~Type) +
  xlab("Trial") + ylab("Proportion of HGT events") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg format
# ggsave(file = "plots/chp3/clean_plots/m4efad_ppn_hgt_directional.svg", 
      # plot = ppn_hgt_directional_plot2, units = "cm", width=12, height=12)
# pdf format
# ggsave(file = "plots/chp3/clean_plots/m4efad_ppn_hgt_directional.pdf", 
      # plot = ppn_hgt_directional_plot2, units = "cm", width=12, height=12)


# subset groups: 12 - MAM v Healthy ---------------------------------------

# plot the proportion of all HGT events, and those between species, with direction
ppn_hgt_directional_plot <- dir_hgts_mam12 %>% 
  bind_rows(dir_hgts_healthy12, spp_dir_hgts_mam12, spp_dir_hgts_healthy12) %>% # bind data for all trials
  select(Type, Trial, Directional) %>%
  group_by(Type, Trial, Directional) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = Trial, y = n, fill = Directional)) +
  geom_bar(position = "fill", stat = "identity") +
  geom_text(aes(label = n), position = position_fill(vjust=0.5), colour="white") +
  scale_fill_manual(values = c("#999999", "#4ABE56"), name = "Directed") +
  facet_grid(~Type) +
  xlab("Group") + ylab("Proportion of HGT events") + 
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg format
# ggsave(file = "plots/chp3/clean_plots/ppn_hgt_directional.svg", 
#      plot = ppn_hgt_directional_plot, units = "cm", width=12, height=12)
# pdf format
# ggsave(file = "plots/chp3/clean_plots/ppn_hgt_directional.pdf", 
#       plot = ppn_hgt_directional_plot, units = "cm", width=12, height=12)

# quantify HGT events

# all m4efad
nrow(dir_hgts_m4efad) # m4efad all HGTs: 48276
nrow(dir_hgts_m4efad %>% filter(Directional == "Yes"))/nrow(dir_hgts_m4efad)*100 # 12.66882 % directed

nrow(spp_dir_hgts_m4efad) # m4efad species-specific HGTs: 32765
nrow(spp_dir_hgts_m4efad %>% filter(Directional == "Yes"))/nrow(spp_dir_hgts_m4efad)*100 # 11.81749 % directed

# MAM_12
nrow(dir_hgts_mam12) # mam 12 HGTs: 9255
nrow(dir_hgts_mam12 %>% filter(Directional == "Yes"))/nrow(dir_hgts_mam12)*100 # 14.01405 % directed

nrow(spp_dir_hgts_mam12) # mam 12 species-specific HGTs: 6167
nrow(spp_dir_hgts_mam12 %>% filter(Directional == "Yes"))/nrow(spp_dir_hgts_mam12)*100 # 13.79925 % directed


# Healthy_12
nrow(dir_hgts_healthy12) # healthy 12 HGTs: 4851
nrow(dir_hgts_healthy12 %>% filter(Directional == "Yes"))/nrow(dir_hgts_healthy12)*100 # 11.11111 % directed

nrow(spp_dir_hgts_healthy12) # healthy 12 HGTs: 2885
nrow(spp_dir_hgts_healthy12 %>% filter(Directional == "Yes"))/nrow(spp_dir_hgts_healthy12)*100 # 11.75043 % directed


# investigate the proportion of intra-taxon HGT events --------------------

# whole m4efad data  ------------------------------------------------------

# all HGT events from m4efad 
all_hgt_events <- dir_hgts_m4efad 

# bind rows of taxonomy a and taxonomy b information to capture all taxa involved in HGT (sample subset)
hgt_taxa_a <- all_hgt_events %>%
  select(CLADE_A, TAXONOMY_A) %>%
  rename(Clade = CLADE_A,
         Taxonomy = TAXONOMY_A)

hgt_taxa_b <- all_hgt_events %>%
  select(CLADE_B, TAXONOMY_B) %>%
  rename(Clade = CLADE_B,
         Taxonomy = TAXONOMY_B)

# create taxonomy lookup table for all species with HGT interactions (sample subset)
all_hgt_taxa_distinct <- hgt_taxa_a %>%
  bind_rows(hgt_taxa_b) %>%
  distinct() %>% # only need distinct taxa to make up lookup table
  mutate(Lowest_taxonomic_level = case_when( # identify lowest taxonomic level - use for filtering lookup table later
    str_starts(Clade, "s__") ~ "Species",
    str_starts(Clade, "g__") ~ "Genus",
    str_starts(Clade, "f__") ~ "Family",
    str_starts(Clade, "o__") ~ "Order",
    str_starts(Clade, "c__") ~ "Class",
    str_starts(Clade, "p__") ~ "Phylum",
    str_starts(Clade, "k__") ~ "Kingdom")) %>%
  # extract kingdom
  mutate(Kingdom_tmp = str_extract(Taxonomy, "(?<=k__).*"), # k__ is before match, .* match everything
         Kingdom = if_else(grepl("\\|", Kingdom_tmp), str_extract(Kingdom_tmp, ".*(?=\\|p__)"), Kingdom_tmp)) %>% # if there are more levels after, extract this part
  select(-Kingdom_tmp) %>% # remove temporary column
  # extract phylum if available
  mutate(Phylum_tmp = str_extract(Taxonomy, "(?<=p__).*"),
         Phylum = if_else(grepl("\\|", Phylum_tmp), str_extract(Phylum_tmp, ".*(?=\\|c__)"), Phylum_tmp)) %>%
  select(-Phylum_tmp) %>%
  # extract class if available
  mutate(Class_tmp = str_extract(Taxonomy, "(?<=c__).*"),
         Class = if_else(grepl("\\|", Class_tmp), str_extract(Class_tmp, ".*(?=\\|o__)"), Class_tmp)) %>%
  select(-Class_tmp) %>%
  # extract order if available
  mutate(Order_tmp = str_extract(Taxonomy, "(?<=o__).*"),
         Order = if_else(grepl("\\|", Order_tmp), str_extract(Order_tmp, ".*(?=\\|f__)"), Order_tmp)) %>%
  select(-Order_tmp) %>%
  # extract family if available
  mutate(Family_tmp = str_extract(Taxonomy, "(?<=f__).*"),
         Family = if_else(grepl("\\|", Family_tmp), str_extract(Family_tmp, ".*(?=\\|g__)"), Family_tmp)) %>%
  select(-Family_tmp) %>%
  # extract genus if available
  mutate(Genus_tmp = str_extract(Taxonomy, "(?<=g__).*"),
         Genus = if_else(grepl("\\|", Genus_tmp), str_extract(Genus_tmp, ".*(?=\\|s__)"), Genus_tmp)) %>%
  select(-Genus_tmp) %>%
  # extract species if available
  mutate(Species_tmp = str_extract(Taxonomy, "(?<=s__).*"),
         Species = str_replace_all(Species_tmp, "_", " ")) %>% # remove underscores from species names
  select(-Species_tmp) %>%
  
  select(-Taxonomy) # no longer need taxonomy column

# find the proportion of species-specific HGT events occurring within each taxonomic level
prop_within_taxa_hgts <- spp_dir_hgts_m4efad %>%
  filter(Directional == "Yes") %>% # filter for directional HGT events
  mutate(Donor_species = str_remove(CLADE_B, "s__"), # remove prefix from species names
         Recipient_species = str_remove(CLADE_A, "s__")) %>%
  mutate(Donor_species = str_replace_all(Donor_species, "_", " "), # remove underscores from species names
         Recipient_species = str_replace_all(Recipient_species, "_", " ")) %>%
  select(Donor_species, Recipient_species, Trial, Timepoint_general) %>%
  left_join(all_hgt_taxa_distinct %>% select(-c(Clade, Lowest_taxonomic_level)), by = c("Donor_species" = "Species")) %>% # join source with taxonomy data
  rename(Source_kingdom = Kingdom,
         Source_phylum = Phylum,
         Source_class = Class,
         Source_order = Order,
         Source_family = Family,
         Source_genus = Genus,
         Source_species = Donor_species) %>%
  left_join(all_hgt_taxa_distinct %>% select(-c(Clade, Lowest_taxonomic_level)), by = c("Recipient_species" = "Species")) %>% # join sink with taxonomy data
  rename(Sink_kingdom = Kingdom,
         Sink_phylum = Phylum,
         Sink_class = Class,
         Sink_order = Order,
         Sink_family = Family,
         Sink_genus = Genus,
         Sink_species = Recipient_species) %>%
  group_by(Trial, Timepoint_general) %>% # group by trial, cohort
  # calculate proportion of HGT events in each cohort that occur within each taxonomic level
  summarise(Kingdom = mean(Source_kingdom == Sink_kingdom),
            Phylum = mean(Source_phylum == Sink_phylum),
            Class = mean(Source_class == Sink_class),
            Order = mean(Source_order == Sink_order),
            Family = mean(Source_family == Sink_family),
            Genus = mean(Source_genus == Sink_genus)) %>%
  ungroup() %>%
  pivot_longer(cols = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), # convert data to long format for plotting
               names_to = "Taxonomic_level", # Taxonomic level
               values_to = "Proportion_HGTs") %>% # proportion of HGT events occurring within taxonomic group
  mutate(Percentage_HGTs = Proportion_HGTs*100) # convert to percentage

prop_within_taxa_hgts$Timepoint_general <- factor(prop_within_taxa_hgts$Timepoint_general, levels = c("MAM_12", "MAM_15", "MAM_24", "Healthy_12", "Healthy_24"))

# quantify average rate of intra-taxon HGT events
prop_within_taxa_hgts %>%
  group_by(Taxonomic_level) %>%
  summarise(mean = mean(Percentage_HGTs))

# plot the proportion of directed species-specific HGT events occurring within each taxonomic level
prop_within_taxa_hgts_plot <- prop_within_taxa_hgts %>%
  ggplot(aes(x = factor(Taxonomic_level, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")), y = Percentage_HGTs, fill = Trial)) +
  #geom_line(aes(group = interaction(Trial, Timepoint_general), colour = Trial), show.legend = FALSE) +
  geom_point(aes(group = interaction(Trial, Timepoint_general),  shape = Timepoint_general), size = 3) +
  #scale_colour_manual(name = "FMT trial", values = trial_palette) +
  scale_fill_manual(name = "M4EFaD", values = trial_palette2) +
  scale_shape_manual(name = "Groups", values = c(21, 22, 23, 24, 25)) +
  scale_y_continuous(name = "HGT events occurring within taxonomic group (%)", limits = c(20, 100)) +
  xlab(NULL) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
#ggsave(file = "../clean_figures/prop_within_taxa_hgts.svg", 
#       plot = prop_within_taxa_hgts_plot, units = "cm", width=14, height=12)

# get donor species in species-specific HGT events
dir_hgt_donor_spp <- spp_dir_hgts_m4efad %>%
  #bind_rows(spp_dir_hgts_gb) %>% # bind data for all trials
  filter(Directional == "Yes") %>% # filter for directional HGT events
  mutate(Donor_species = str_remove(CLADE_B, "s__"), # remove prefix from species names
         Donor_species = str_replace_all(Donor_species, "_", " ")) %>% # remove underscores from species names
  select(Species = Donor_species)

# get recipient species in species-specific HGT events
dir_hgt_recipient_spp <- spp_dir_hgts_m4efad %>%
 # bind_rows(spp_dir_hgts_gb) %>% # bind data for all trials
  filter(Directional == "Yes") %>% # filter for directional HGT events
  mutate(Recipient_species = str_remove(CLADE_A, "s__"), # remove prefix from species names
         Recipient_species = str_replace_all(Recipient_species, "_", " ")) %>% # remove underscores from species names
  select(Species = Recipient_species)

# join and get distinct species
distinct_dir_hgt_spp <- bind_rows(dir_hgt_donor_spp, dir_hgt_recipient_spp) %>% distinct() # 298 distinct species involved in HGT across all samples

# count number of distinct taxa involved in directed species-specific HGT events at each taxonomic level
n_distinct_dir_hgt_taxa <- distinct_dir_hgt_spp %>%
  left_join(all_hgt_taxa_distinct, by = "Species") %>%
  select(-c(Species, Clade, Lowest_taxonomic_level)) %>%
  summarise(across(everything(), n_distinct)) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("Taxonomic_level") %>%
  rename(n_distinct_taxa = V1)

# plot the number of distinct taxa involved in directed species-specific HGT events at each taxonomic level
hgt_taxa_count_plot <- n_distinct_dir_hgt_taxa %>%
  mutate(Group = "a") %>% # grouping variable for line plot
  ggplot(aes(x = factor(Taxonomic_level, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")), y = n_distinct_taxa)) +
  geom_line(aes(group = Group)) +
  scale_y_reverse(name = "Number of distinct taxa", limits = c(100, 0), position = "right") + # reverse y axis
  xlab(NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
ggsave(file = "plots/chp3//clean_plots/hgt_taxa_count_plot.svg", 
     plot = hgt_taxa_count_plot, units = "cm", width=14, height=12)

# note: ggplot doesn't support plots with two independent y axes, so we make above plots with same dimensions and overlay afterwards

# remake proportion plot just to get legend
prop_within_taxa_hgts_legend <- prop_within_taxa_hgts %>%
  ggplot(aes(x = factor(Taxonomic_level, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")), y = Percentage_HGTs, fill = Trial)) +
  #geom_line(aes(group = interaction(Trial, Timepoint_general), colour = Trial), show.legend = FALSE) +
  geom_point(aes(group = interaction(Trial, Timepoint_general),  shape = Timepoint_general), size = 3) +
  #scale_colour_manual(name = "FMT trial", values = trial_palette) +
  scale_fill_manual(name = "M4EFaD", values = trial_palette2) +
  scale_shape_manual(name = "Groups", values = c(21, 22, 23, 24, 25)) +
  scale_y_continuous(name = "HGT events occurring within taxonomic group (%)", limits = c(20, 100)) +
  xlab(NULL) +
  guides(fill = guide_legend("Trial", override.aes = list(shape = 21, color = "black")),
         shape = guide_legend("Group")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
ggsave(file = "plots/chp3/clean_plots/prop_within_taxa_hgts_legend.svg", 
   plot = prop_within_taxa_hgts_legend, units = "cm", width=18, height=12)



# save data files ---------------------------------------------------------
# Note: These files are already in the hgt_analyses/data folder in the github repo

#save(spp_dir_hgts_m4efad, file = "data/processed/clean_outputs/spp_dir_hgts_m4efad.RData")
#write.table(all_hgt_taxa_distinct, "data/processed/clean_outputs/hgt_taxa_lookup.txt", sep = "\t", quote = FALSE)





