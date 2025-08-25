## code summary:
# compare the impact of  intervention on HGT network membership
# investigate the restoration of 'healthy' cohort species network roles, post-intervention
# compare the restoration of species network role to 'healthy' cohort role, in FMT vs placebo recipients

## load libraries

library(tidyverse)
library(ggalluvial)

## set working directory

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad/")

## load data

# species network roles
load("data/processed/clean_outputs/joined_spp_roles.RData")

## load palettes

role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")


## load functions

# investigate the restoration of 'healthy' cohort species network roles, post-intervention
spp_role_restoration <- function(joined_spp_roles, donor_cohort, pre_intervention_cohort, post_intervention_cohort){
  ## joined_spp_roles: joined species role data for all trial cohorts (df)
  ## donor_cohort: name of trial donor cohort (str)
  ## pre_intervention_cohort: name of trial pre-intervention cohort (str)
  ## post_intervention_cohort: name of trial post-intervention cohort (str)
  spp_role_restoration <- joined_spp_roles %>%
    filter(Trial_cohort == donor_cohort) %>%
    select(Species, Type) %>%
    rename(Donor_role = Type) %>% # species network role in donor cohort
    inner_join(joined_spp_roles %>% filter(Trial_cohort == pre_intervention_cohort) %>% select(Species, Type), by = "Species") %>%
    rename(Pre_intervention_role = Type) %>% # species network role in pre-intervention cohort
    inner_join(joined_spp_roles %>% filter(Trial_cohort == post_intervention_cohort) %>% select(Species, Type), by = "Species") %>%
    rename(Post_intervention_role = Type) %>% # species network role in post-intervention cohort
    mutate(Role_restored = ifelse(Donor_role == Post_intervention_role & 
                                    Donor_role != Pre_intervention_role & Post_intervention_role != Pre_intervention_role, Donor_role, "No"))
  # fix factor levels
  spp_role_restoration$Donor_role <- factor(spp_role_restoration$Donor_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Pre_intervention_role <- factor(spp_role_restoration$Pre_intervention_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Post_intervention_role <- factor(spp_role_restoration$Post_intervention_role, levels = c("Source", "Conduit", "Sink"))
  
  spp_role_restoration # return data
}

# plot the restoration of 'healthy' cohort species network roles, post-intervention 
role_restoration_alluvial <- function(role_restoration_data, pre_intervention_name, post_intervention_name){
  ## role_restoration_data: output of spp_role_restoration() (df)
  ## pre_intervention_name: e.g. "Pre-FMT" or "Pre-placebo" (str)
  ## post_intervention_name: e.g. "Post-FMT" or "Post-placebo" (str)
  role_restoration_data %>%
    ggplot(aes(axis1 = factor(Donor_role), axis2 = factor(Pre_intervention_role), axis3 = factor(Post_intervention_role), y = 1)) +
    geom_alluvium(aes(fill = Role_restored), alpha = 0.6) + 
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(name = "Donor ID", drop = FALSE, values = role_palette) +
    scale_x_continuous(breaks=c(1, 2, 3), 
                       labels=c("Donor", pre_intervention_name, post_intervention_name)) +
    ylab("Number of shared species") +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          legend.background = element_rect(color="black", size=0.4))
}


# for the m4efad trial - the healthy cohort only contains 2 timepoints, 12 months and 24 months. 
# there are now only 2 axes for the plot and the function logic is renewed to match the 2 axes. 

# Function to analyze and visualize species role restoration in the healthy cohort
spp_role_restoration_healthy <- function(joined_spp_roles, donor_cohort, post_intervention_cohort){
  ## joined_spp_roles: joined species role data for all trial cohorts (df)
  ## donor_cohort: name of trial donor cohort (str)
  ## post_intervention_cohort: name of trial post-intervention cohort (str)
  
  spp_role_restoration <- joined_spp_roles %>%
    filter(Trial_cohort == donor_cohort) %>%
    select(Species, Type) %>%
    rename(Donor_role = Type) %>%  # Species network role in donor cohort
    inner_join(joined_spp_roles %>% filter(Trial_cohort == post_intervention_cohort) %>% select(Species, Type), by = "Species") %>%
    rename(Post_intervention_role = Type) %>%  # Species network role in post-intervention cohort
    mutate(Role_restored = ifelse(Donor_role == Post_intervention_role, Donor_role, "No"))  # Check restoration
  
  # Convert roles to factors
  spp_role_restoration$Donor_role <- factor(spp_role_restoration$Donor_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Post_intervention_role <- factor(spp_role_restoration$Post_intervention_role, levels = c("Source", "Conduit", "Sink"))
  
  return(spp_role_restoration)
}

# Function to create an alluvial plot for role restoration (Donor → Post-Intervention)
role_restoration_alluvial_healthy <- function(role_restoration_data, post_intervention_name){
  ## role_restoration_data: output of spp_role_restoration_healthy() (df)
  ## post_intervention_name: e.g. "Post-FMT" or "Post-placebo" (str)
  
  role_restoration_data %>%
    ggplot(aes(axis1 = factor(Donor_role), axis2 = factor(Post_intervention_role), y = 1)) +
    geom_alluvium(aes(fill = Role_restored), alpha = 0.6) +  
    geom_stratum() +  
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +  
    scale_fill_manual(name = "Donor ID", drop = FALSE, values = role_palette) +  
    scale_x_continuous(breaks = c(1, 2), labels = c("Donor", post_intervention_name)) +  
    ylab("Number of shared species") +  
    theme_minimal() +  
    theme(legend.position = "none",  
          panel.grid.major = element_blank(),  
          panel.grid.minor = element_blank(),  
          panel.border = element_blank(),  
          legend.background = element_rect(color = "black", size = 0.4))
}



## running

# investigate restoration of network roles after ----------------------

# investigate the restoration of 'healthy' cohort species network roles, post-intervention
# for the m4efad trial, the equivalent to UC and GB placebo would be investigating role restoration between 12 v 24 healthy cohort 

# mam - baseline to post intervention 
mam_role_restoration <- spp_role_restoration(joined_spp_roles, donor_cohort = "M4EFaD MAM_12", 
                                                pre_intervention_cohort = "M4EFaD MAM_15", post_intervention_cohort = "M4EFaD MAM_24")


# healthy post time - 12 to 24 months
healthy_restoration <- spp_role_restoration(joined_spp_roles, donor_cohort = "M4EFaD Healthy_12", 
                                                    pre_intervention_cohort = "M4EFaD Healthy_12", post_intervention_cohort = "M4EFaD Healthy_24")

## function specific for healthy with two axes
healthy_restoration_only <- spp_role_restoration_healthy(joined_spp_roles, donor_cohort = "M4EFaD Healthy_12", 
                                            post_intervention_cohort = "M4EFaD Healthy_24")



# plot the restoration of 'mam' and 'healthy' cohort species network roles, post-intervention (mam) and post time (healthy)

## MAM
mam_role_restoration_plot <- role_restoration_alluvial(mam_role_restoration, 
                                                       pre_intervention_name = "MAM pre-intervention", 
                                                       post_intervention_name = "MAM post-intervention")
ggsave(file = "plots/chp3/clean_plots/mam_role_restoration_plot.svg", plot = mam_role_restoration_plot, units = "cm", width=14, height=8) 
ggsave(file = "plots/chp3/clean_plots/mam_role_restoration_plot.pdf", plot = mam_role_restoration_plot, units = "cm", width=14, height=8) 

## Healthy 
healthy_role_restoration_plot <- role_restoration_alluvial(healthy_restoration, 
                                                              pre_intervention_name = "Healthy_12", post_intervention_name = "Healthy_24")
ggsave(file = "plots/chp3/clean_plots/healthy_role_restoration_plot.svg", plot = healthy_role_restoration_plot, units = "cm", width=14, height=8) 
ggsave(file = "plots/chp3/clean_plots/healthy_role_restoration_plot.pdf", plot = healthy_role_restoration_plot, units = "cm", width=14, height=8) 


## healthy using 2 axes only 
healthy_role_restoration_only_plot <- role_restoration_alluvial_healthy(healthy_restoration_only, 
                                                           donor_role = "Healthy_12", post_intervention_name = "Healthy_24")

healthy_role_restoration_only_plot <- role_restoration_alluvial_healthy(healthy_restoration_only, "Healthy_24")

ggsave(file = "plots/chp3/clean_plots/healthy_only_role_restoration_plot.svg", plot = healthy_role_restoration_only_plot, units = "cm", width=14, height=8) 
ggsave(file = "plots/chp3/clean_plots/healthy_only_role_restoration_plot.pdf", plot = healthy_role_restoration_only_plot, units = "cm", width=14, height=8) 


# for each trial, statistically compare the restoration of species network role to 'healthy' cohort role, in FMT vs placebo recipients
# calculate proportions
mean(fmt_role_restoration_uc$Role_restored != "No") # 0.118
mean(placebo_role_restoration_uc$Role_restored != "No") # 0.139

mean(fmt_role_restoration_gb$Role_restored != "No") # 0.206
mean(placebo_role_restoration_gb$Role_restored != "No") # 0.174

# compare proportions with Fisher's exact test due to small sample sizes
n_restored_fmt_uc <- sum(fmt_role_restoration_uc$Role_restored != "No")
n_total_fmt_uc <- nrow(fmt_role_restoration_uc)
n_restored_placebo_uc <- sum(placebo_role_restoration_uc$Role_restored != "No")
n_total_placebo_uc <- nrow(placebo_role_restoration_uc)

fisher.test(matrix(c(n_restored_fmt_uc, n_total_fmt_uc - n_restored_fmt_uc,
                     n_restored_placebo_uc, n_total_placebo_uc - n_restored_placebo_uc), nrow = 2)) # p = 1

n_restored_fmt_gb <- sum(fmt_role_restoration_gb$Role_restored != "No")
n_total_fmt_gb <- nrow(fmt_role_restoration_gb)
n_restored_placebo_gb <- sum(placebo_role_restoration_gb$Role_restored != "No")
n_total_placebo_gb <- nrow(placebo_role_restoration_gb)

fisher.test(matrix(c(n_restored_fmt_gb, n_total_fmt_gb - n_restored_fmt_gb,
                     n_restored_placebo_gb, n_total_placebo_gb - n_restored_placebo_gb), nrow = 2)) # p = 0.5962



# TESTING -----------------------------------------------------------------

## Load palettes
role_palette <- c("Source-Sink Change" = "#DD4444",  # Highlight species changing between Source and Sink
                  "Other" = "grey")  # Everything else stays grey

## Function to analyze species role transitions across cohorts
analyze_role_transitions <- function(joined_spp_roles, donor_cohort, pre_intervention_cohort, post_intervention_cohort){
  role_transitions <- joined_spp_roles %>%
    filter(Trial_cohort == donor_cohort) %>%
    select(Species, Type) %>%
    rename(Donor_role = Type) %>%
    inner_join(joined_spp_roles %>% filter(Trial_cohort == pre_intervention_cohort) %>%
                 select(Species, Type), by = "Species") %>%
    rename(Pre_intervention_role = Type) %>%
    inner_join(joined_spp_roles %>% filter(Trial_cohort == post_intervention_cohort) %>%
                 select(Species, Type), by = "Species") %>%
    rename(Post_intervention_role = Type) %>%
    mutate(Role_change_highlight = ifelse(
      (Donor_role == "Source" & Post_intervention_role == "Sink") |
        (Donor_role == "Sink" & Post_intervention_role == "Source"),
      "Source-Sink Change", "Other"  # Only highlight Source <-> Sink transitions
    ))
  
  # Fix factor levels
  role_transitions$Donor_role <- factor(role_transitions$Donor_role, levels = c("Source", "Conduit", "Sink"))
  role_transitions$Pre_intervention_role <- factor(role_transitions$Pre_intervention_role, levels = c("Source", "Conduit", "Sink"))
  role_transitions$Post_intervention_role <- factor(role_transitions$Post_intervention_role, levels = c("Source", "Conduit", "Sink"))
  role_transitions$Role_change_highlight <- factor(role_transitions$Role_change_highlight, levels = c("Source-Sink Change", "Other"))
  
  role_transitions  # Return data
}

## Function to plot species role transitions as an alluvial diagram
plot_role_transitions <- function(role_transitions, pre_intervention_name, post_intervention_name){
  role_transitions %>%
    ggplot(aes(axis1 = factor(Donor_role), 
               axis2 = factor(Pre_intervention_role), 
               axis3 = factor(Post_intervention_role), 
               y = 1)) +
    geom_alluvium(aes(fill = Role_change_highlight), alpha = 0.6) +  # Color only specific changes
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(name = "Role Change", drop = FALSE, values = role_palette) +  # Updated legend
    scale_x_continuous(breaks=c(1, 2, 3), 
                       labels=c("Donor", pre_intervention_name, post_intervention_name)) +
    ylab("Number of shared species") +
    theme(legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.border = element_blank(),
          legend.background = element_rect(color="black", size=0.4))
}


### HEALTHY 

# Function to analyze species role restoration in the healthy cohort
analyze_role_restoration_healthy <- function(joined_spp_roles, donor_cohort, post_intervention_cohort){
  spp_role_restoration <- joined_spp_roles %>%
    filter(Trial_cohort == donor_cohort) %>%
    select(Species, Type) %>%
    rename(Donor_role = Type) %>%
    inner_join(joined_spp_roles %>% filter(Trial_cohort == post_intervention_cohort) %>%
                 select(Species, Type), by = "Species") %>%
    rename(Post_intervention_role = Type) %>%
    mutate(Role_restored = case_when(
      Donor_role == Post_intervention_role ~ "No Change",
      Donor_role == "Conduit" | Post_intervention_role == "Conduit" ~ "No Change",
      Donor_role == "Sink" & Post_intervention_role == "Source" ~ "Role Switch",
      Donor_role == "Source" & Post_intervention_role == "Sink" ~ "Role Switch",
      TRUE ~ "No Change"
    ))
  
  # Convert roles to factors
  spp_role_restoration$Donor_role <- factor(spp_role_restoration$Donor_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Post_intervention_role <- factor(spp_role_restoration$Post_intervention_role, levels = c("Source", "Conduit", "Sink"))
  spp_role_restoration$Role_restored <- factor(spp_role_restoration$Role_restored, levels = c("No Change", "Role Switch"))
  
  return(spp_role_restoration)
}



plot_role_restoration_healthy <- function(role_restoration_data, post_intervention_name){
  role_restoration_data %>%
    ggplot(aes(axis1 = factor(Donor_role), 
               axis2 = factor(Post_intervention_role), 
               y = 1)) +
    geom_alluvium(aes(fill = Role_restored), alpha = 0.6) +  
    geom_stratum() +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_fill_manual(name = "Role Change", 
                      values = c("No Change" = "grey", "Role Switch" = "#DD4444"), 
                      drop = FALSE) +
    scale_x_continuous(breaks = c(1, 2), 
                       labels = c("Donor", post_intervention_name)) +
    ylab("Number of shared species") +
    theme_minimal() +
    theme(legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.background = element_rect(color = "black", size = 0.4))
}





## RUNNING


# MAM - Baseline to Post-Intervention  
mam_role_transitions <- analyze_role_transitions(joined_spp_roles, 
                                                 donor_cohort = "M4EFaD MAM_12", 
                                                 pre_intervention_cohort = "M4EFaD MAM_15", 
                                                 post_intervention_cohort = "M4EFaD MAM_24")

# Plot the role transitions for MAM cohort  
mam_role_transitions_plot <- plot_role_transitions(mam_role_transitions,  
                                                   pre_intervention_name = "MAM pre-intervention",  
                                                   post_intervention_name = "MAM post-intervention")



# Analyze role restoration for the healthy cohort
healthy_role_restoration <- analyze_role_restoration_healthy(joined_spp_roles, 
                                                             donor_cohort = "M4EFaD Healthy_12", 
                                                             post_intervention_cohort = "M4EFaD Healthy_24")

# Generate the alluvial plot
healthy_role_restoration_plot <- plot_role_restoration_healthy(healthy_role_restoration, "Healthy_24")



### SINK SOURCE CHANGE 

sink_source_species <- joined_spp_roles %>%
  filter(Type %in% c("Sink", "Source")) %>%      # Only keep rows where role is Sink or Source
  group_by(Species) %>%
  filter(n_distinct(Type) == 2) %>%               # Keep species that have both Sink and Source
  ungroup() %>%
  arrange(Species, Cohort) %>%                    # Arrange by species and cohort
  select(Species, Type, Cohort)                   # Select the desired columns


# Create the heatmap
ggplot(sink_source_species, aes(x = Cohort, y = Species, fill = Type)) +
  geom_tile(color = "white") +  # White borders for better separation
  scale_fill_manual(values = c("Sink" = "#EE8866", "Source" = "#44BB99")) +  # Set custom colors for Sink and Source
  labs(title = "Species Switching between Sink and Source Roles Across Cohorts",
       x = "Cohort",
       y = "Species",
       fill = "Role") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




#### CALCULATING THE ROLE SWITCH PERCENTAGES 

library(dplyr)
library(ggplot2)

species_switch <- joined_spp_roles %>%
  mutate(Group = case_when(
    grepl("Healthy", Cohort) ~ "Healthy",
    grepl("MAM", Cohort) ~ "MAM",
    TRUE ~ NA_character_
  )) %>%
  filter(Type %in% c("Sink", "Source")) %>%  # Only consider Sink and Source roles
  group_by(Group, Species) %>%
  summarise(n_roles = n_distinct(Type), .groups = "drop") %>%
  mutate(Switch = if_else(n_roles == 2, "Switched", "No Switch"))

# Step 2: Summarise percentages by group
switch_summary <- species_switch %>%
  group_by(Group, Switch) %>%
  summarise(n_species = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(total_species = sum(n_species),
         percent = (n_species / total_species) * 100) %>%
  ungroup()

print(switch_summary)

# Step 3: Plot the percentages as a bar graph
ggplot(switch_summary, aes(x = Group, y = percent, fill = Switch)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage of Species with Sink ↔ Source Role Switches",
       x = "Group (MAM vs Healthy)",
       y = "Percentage of Species",
       fill = "Role Switch") +
  scale_fill_manual(values = c("Switched" = "#DD4444", "No Switch" = "grey")) +
  theme_minimal()



## TESTING - how many unique subjects had the Coprococcus species as a donor or reciever

# Filter for rows where s__Coprococcus_comes appears in CLADE_A or CLADE_B
coprococcus_counts <- dir_spp_hgts_m4efad %>%
  filter(CLADE_A == "s__Coprococcus_comes" | CLADE_B == "s__Coprococcus_comes") %>%
  summarise(unique_participants = n_distinct(Participant_ID))

# Print the result
print(coprococcus_counts)

coprococcus_counts_in_samples <- dir_spp_hgts_m4efad %>%
  filter(CLADE_A == "s__Coprococcus_comes" | CLADE_B == "s__Coprococcus_comes") %>%
  summarise(unique_samples = n_distinct(Sample_ID))

coprococcus_counts_in_contig <- dir_spp_hgts_m4efad %>%
  filter(CLADE_A == "s__Coprococcus_comes" | CLADE_B == "s__Coprococcus_comes") %>%
  summarise(unique_contig = n_distinct(CONTIG_NAME))

# n = 58


# Identify species that have different Type in different Cohorts
species_different_roles <- joined_spp_roles %>%
  group_by(Species) %>%
  filter(n_distinct(Type) > 1) %>%
  ungroup()

# Print the result
print(species_different_roles)


# Identify species that have different Type in MAM vs Healthy cohorts
species_different_roles_mam_healthy <- joined_spp_roles %>%
  mutate(Group = ifelse(grepl("MAM", Cohort), "MAM", "Healthy")) %>% # Categorize cohorts into MAM and Healthy
  group_by(Species, Group) %>%
  summarise(distinct_roles = n_distinct(Type), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = distinct_roles, values_fill = 0) %>%
  filter(MAM > 0 & Healthy > 0) %>% # Select species that appear in both groups
  inner_join(joined_spp_roles, by = "Species") %>%
  select(Species, Cohort, Type)

# Print the result
print(species_different_roles_mam_healthy)



# Identify species that have different "Type" between any MAM and any Healthy cohort
species_different_roles_mam_healthy <- joined_spp_roles %>%
  mutate(Group = ifelse(grepl("MAM", Cohort), "MAM", "Healthy")) %>% # Classify into MAM or Healthy
  group_by(Species, Group) %>%
  summarise(unique_types = list(unique(Type)), .groups = "drop") %>% # Get unique types per group
  pivot_wider(names_from = Group, values_from = unique_types) %>% # Spread MAM and Healthy
  filter(!map2_lgl(MAM, Healthy, ~setequal(.x, .y))) %>% # Keep species where MAM and Healthy roles are different
  select(Species) %>%
  inner_join(joined_spp_roles, by = "Species") %>% # Rejoin with the dataset to retain details
  arrange(Species)

# Print the result
print(species_different_roles_mam_healthy)




# Species w/ different roles within each Cohort  ------------------------------

library(dplyr)
library(tidyr)
library(ggplot2)

# Identify species that switch between Sink and Source within each Cohort
species_sink_source_within_cohort <- joined_spp_roles %>%
  group_by(Cohort, Species) %>%
  filter("Sink" %in% Type & "Source" %in% Type) %>% # Ensure species have both Sink and Source roles
  ungroup() %>%
  arrange(Cohort, Species) # Arrange by Cohort and Species


library(ggplot2)
library(dplyr)

# Calculate summary statistics for dr_ratio within each Cohort
dr_ratio_trend <- joined_spp_roles %>%
  group_by(Cohort) %>%
  summarise(
    mean_dr_ratio = mean(dr_ratio, na.rm = TRUE),
    median_dr_ratio = median(dr_ratio, na.rm = TRUE),
    sd_dr_ratio = sd(dr_ratio, na.rm = TRUE),
    count = n()
  )

# Print the summary statistics
print(dr_ratio_trend)

# Plot the distribution of dr_ratio within each Cohort
ggplot(joined_spp_roles, aes(x = Cohort, y = dr_ratio)) +
  geom_boxplot(aes(fill = Cohort), alpha = 0.7) + 
  geom_jitter(aes(color = Cohort), width = 0.2, alpha = 0.5) + 
  theme_minimal() +
  labs(
    title = "DR Ratio Trend Across Cohorts",
    x = "Cohort",
    y = "Donor-Receiver Ratio"
  ) +
  theme(legend.position = "none")



## How many species maintained their role within each Cohort and how many switched to Conduit?

library(dplyr)

# Process the dataset
role_proportions <- joined_spp_roles %>%
  mutate(Group = sub("_.*", "", Cohort)) %>%  # Remove suffix to get Group
  group_by(Group, Species) %>%
  summarise(
    n_roles = n_distinct(Type),  # Count distinct roles per species within the group
    .groups = "drop"
  ) %>%
  mutate(
    Role_Change = ifelse(n_roles > 1, "Changed Role", "Maintained Role")  # Identify species that changed roles
  ) %>%
  group_by(Group, Role_Change) %>%
  summarise(
    Species_Count = n(),  # Count species in each category
    .groups = "drop"
  ) %>%
  group_by(Group) %>%
  mutate(
    Total_Species = sum(Species_Count),
    Proportion = Species_Count / Total_Species * 100  # Calculate proportion
  ) %>%
  ungroup()

# Print the result
print(role_proportions)

# Plot the results
ggplot(role_proportions, aes(x = Group, y = Proportion, fill = Role_Change)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Proportion of Species Maintaining or Changing Roles",
       x = "Group",
       y = "Proportion (%)",
       fill = "Role Status") +
  scale_fill_manual(values = c("Changed Role" = "#DD4444", "Maintained Role" = "#44BB99")) +
  theme_minimal()

# Extract species that maintained their roles
maintained_roles_species <- joined_spp_roles %>%
  mutate(Group = sub("_.*", "", Cohort)) %>%
  group_by(Group, Species) %>%
  summarise(n_roles = n_distinct(Type), .groups = "drop") %>%
  filter(n_roles == 1) %>%  # Species that have the same role across all instances
  select(Group, Species)

# Extract species that changed roles
changed_roles_species <- joined_spp_roles %>%
  mutate(Group = sub("_.*", "", Cohort)) %>%
  group_by(Group, Species) %>%
  summarise(n_roles = n_distinct(Type), .groups = "drop") %>%
  filter(n_roles > 1) %>%  # Species that switched roles
  select(Group, Species)

# Print results
print("Species that maintained roles:")
print(maintained_roles_species)

print("Species that changed roles:")
print(changed_roles_species)



# Run prop test  ----------------------------------------------------------

library(dplyr)
library(broom)

# Prepare data for the proportion test
role_proportion_test_data <- joined_spp_roles %>%
  mutate(Group = sub("_.*", "", Cohort)) %>%  # Extract Group name (MAM or Healthy)
  group_by(Group, Species) %>%
  summarise(n_roles = n_distinct(Type), .groups = "drop") %>%
  mutate(Role_Change = ifelse(n_roles > 1, "Changed Role", "Maintained Role")) %>%
  group_by(Group, Role_Change) %>%
  summarise(Species_Count = n(), .groups = "drop") %>%
  group_by(Group) %>%
  mutate(Total_Species = sum(Species_Count)) %>%
  ungroup()

# Create a summary table
role_summary <- role_proportion_test_data %>%
  select(Group, Role_Change, Species_Count, Total_Species) %>%
  pivot_wider(names_from = Role_Change, values_from = Species_Count, values_fill = 0)

print(role_summary)  # View data for testing


# Run the proportion test
prop_test_result <- prop.test(
  x = c(role_summary$`Changed Role`[role_summary$Group == "MAM"],
        role_summary$`Changed Role`[role_summary$Group == "Healthy"]),
  n = c(role_summary$Total_Species[role_summary$Group == "MAM"],
        role_summary$Total_Species[role_summary$Group == "Healthy"])
)

# Print the test result
print(prop_test_result)
tidy(prop_test_result)  # Cleaned-up results using broom









