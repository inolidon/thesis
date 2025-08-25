## code summary:
# compare species pairs common and unique between Groups (MAM v Healthy) 
# compare species pairs across groups and within timepoints for each group
# identify core shared species pairs 
# identify most prevalent species pair interactions 
# compare the change in prevalence for the shared species pairs between groups


# load packages -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggvenn)
library(grid)
library(rstatix)
library(ggbeeswarm)

# load color palettes  ----------------------------------------------------

# Palette for groups
group_palette <- c(
  "MAM_12" = "darkgoldenrod1",
  "MAM_24" = "cornsilk3",
  "Healthy_12" = "cyan4",
  "Healthy_24" = "chocolate"
)

categories_palette <- c(
  "MAM_12" = "darkgoldenrod1",
  "MAM_24" = "cornsilk3",
  "Healthy_12" = "cyan4",
  "Healthy_24" = "chocolate",
  "Super transferome" = "deeppink"
)

# Palette for timepoints
timepoints_palette <- c(
  "12" = "darkgoldenrod1",
  "15" = "darkmagenta",
  "24" = "cornsilk3"
)



# load functions ----------------------------------------------------------

# Function: compare_interaction_pairs
# Description: 
#   Given a dataset with columns 'CLADE_A', 'CLADE_B', and 'Group', this function:
#   - Filters the data for two specified groups.
#   - Creates a new column "Interaction" (formatted as "CLADE_B -> CLADE_A") that preserves direction.
#   - Identifies the unique interaction pairs for each group.
#   - Computes the common interactions (present in both groups) and the interactions unique to each group.
#   - Calculates the proportion of shared interactions relative to the union of all interactions.
#   - Returns a list containing a summary data frame and a pairwise classification data frame.
compare_interaction_pairs <- function(data, group1, group2) {
  
  # Filter data for the two groups of interest and create the 'Interaction' column
  data_filtered <- data %>%
    filter(Group %in% c(group1, group2)) %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A))
  
  # Get unique interactions for each group
  interactions_group1 <- data_filtered %>%
    filter(Group == group1) %>%
    distinct(Interaction) %>%
    pull(Interaction)
  
  interactions_group2 <- data_filtered %>%
    filter(Group == group2) %>%
    distinct(Interaction) %>%
    pull(Interaction)
  
  # Identify common interactions and those unique to each group
  common_interactions <- intersect(interactions_group1, interactions_group2)
  unique_group1 <- setdiff(interactions_group1, interactions_group2)
  unique_group2 <- setdiff(interactions_group2, interactions_group1)
  
  # Get the union of all interactions
  union_interactions <- union(interactions_group1, interactions_group2)
  
  # Calculate proportions relative to the union
  prop_common <- length(common_interactions) / length(union_interactions)
  prop_unique_group1 <- length(unique_group1) / length(union_interactions)
  prop_unique_group2 <- length(unique_group2) / length(union_interactions)
  
  # Create summary data frame
  summary_df <- data.frame(
    Interaction_Type = c("Common", paste0("Unique_", group1), paste0("Unique_", group2)),
    Count = c(length(common_interactions), length(unique_group1), length(unique_group2)),
    Proportion = c(prop_common, prop_unique_group1, prop_unique_group2)
  )
  
  # Create pairwise classification data frame
  pairwise_df <- data.frame(Interaction = union_interactions, stringsAsFactors = FALSE) %>%
    mutate(Classification = case_when(
      Interaction %in% common_interactions ~ "Common",
      Interaction %in% unique_group1 ~ paste0("Unique_", group1),
      Interaction %in% unique_group2 ~ paste0("Unique_", group2),
      TRUE ~ NA_character_
    ))
  
  return(list(summary = summary_df, pairwise = pairwise_df))
}


# Function: compare_interaction_pairs_timepoint
#
# Description:
#   Given a dataset with columns 'CLADE_A', 'CLADE_B', 'Timepoint', and 'Group',
#   this function filters the data for a specified timepoint and for two groups.
#   It creates a new "Interaction" column (formatted as "CLADE_B -> CLADE_A") to
#   preserve direction, then extracts the unique interaction pairs for each group,
#   identifies common interactions (present in both groups) and those unique to each group,
#   and calculates the proportion of shared interactions relative to the union of all interactions.
#
# Arguments:
#   data: A data frame with at least the columns: CLADE_A, CLADE_B, Timepoint, and Group.
#   timepoint: A character (or numeric) value specifying the timepoint to filter on (e.g., "12" or "24").
#   group1: A character string for the first group (e.g., "MAM").
#   group2: A character string for the second group (e.g., "Healthy").
#
# Returns:
#   A list with two elements:
#     - summary: A data frame summarizing the counts and proportions of interactions classified as
#                "Common", "Unique_<group1>", and "Unique_<group2>".
#     - pairwise: A data frame listing each unique interaction (as "CLADE_B -> CLADE_A") and its classification.
compare_interaction_pairs_timepoint <- function(data, timepoint, group1, group2) {
  
  # Filter data for the specified timepoint and for the two groups of interest
  data_filtered <- data %>%
    filter(Timepoint == timepoint, Group %in% c(group1, group2)) %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A))
  
  # Extract unique interaction pairs for each group
  interactions_group1 <- data_filtered %>%
    filter(Group == group1) %>%
    distinct(Interaction) %>%
    pull(Interaction)
  
  interactions_group2 <- data_filtered %>%
    filter(Group == group2) %>%
    distinct(Interaction) %>%
    pull(Interaction)
  
  # Identify common and unique interaction pairs
  common_interactions <- intersect(interactions_group1, interactions_group2)
  unique_group1 <- setdiff(interactions_group1, interactions_group2)
  unique_group2 <- setdiff(interactions_group2, interactions_group1)
  
  # Create the union of all interactions
  union_interactions <- union(interactions_group1, interactions_group2)
  
  # Calculate proportions relative to the union
  prop_common <- if(length(union_interactions) > 0) length(common_interactions) / length(union_interactions) else NA
  prop_unique_group1 <- if(length(union_interactions) > 0) length(unique_group1) / length(union_interactions) else NA
  prop_unique_group2 <- if(length(union_interactions) > 0) length(unique_group2) / length(union_interactions) else NA
  
  # Create summary data frame
  summary_df <- data.frame(
    Interaction_Type = c("Common", paste0("Unique_", group1), paste0("Unique_", group2)),
    Count = c(length(common_interactions), length(unique_group1), length(unique_group2)),
    Proportion = c(prop_common, prop_unique_group1, prop_unique_group2)
  )
  
  # Create pairwise classification data frame
  pairwise_df <- data.frame(Interaction = union_interactions, stringsAsFactors = FALSE) %>%
    mutate(Classification = case_when(
      Interaction %in% common_interactions ~ "Common",
      Interaction %in% unique_group1 ~ paste0("Unique_", group1),
      Interaction %in% unique_group2 ~ paste0("Unique_", group2),
      TRUE ~ NA_character_
    ))
  
  return(list(summary = summary_df, pairwise = pairwise_df))
}


# Function: compare_species_pairs_four_groups
# Description:
#   Given a dataset with columns:
#     - CLADE_A (recipient species)
#     - CLADE_B (donor species)
#     - Timepoint_general (e.g., "MAM_12", "MAM_24", "Healthy_12", "Healthy_24")
#   This function:
#     1. Filters the data for the specified groups.
#     2. Creates an "Interaction" column (formatted as "CLADE_B -> CLADE_A") that preserves direction.
#     3. For each unique interaction, determines in which groups it occurs.
#     4. Classifies each interaction as "Common" (if present in all groups), 
#        "Unique_[Group]" (if only in one group), or "Shared_[...]" if in two or three groups.
#     5. Returns a list with:
#         - pairwise: A data frame with each unique interaction, the groups in which it is present, and its classification.
#         - summary: A summary table of counts and proportions for each classification.
#
# Args:
#   data: A data frame with at least the columns "CLADE_A", "CLADE_B", and "Timepoint_general".
#   groups: A character vector of group names to compare (default: c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24")).
#
# Returns:
#   A list with two elements:
#     - pairwise: Detailed data frame for each unique interaction.
#     - summary: Summary table with counts and proportions for each classification.

compare_species_pairs_four_groups <- function(data, 
                                              groups = c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24")) {
  # Step 1: Filter data for the specified groups and create an Interaction column preserving direction
  data_filtered <- data %>%
    filter(Timepoint_general %in% groups) %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A))
  
  # Step 2: For each unique interaction, determine the groups (Timepoint_general) in which it appears
  interaction_groups <- data_filtered %>%
    group_by(Interaction) %>%
    summarise(Groups_Present = list(sort(unique(Timepoint_general))), .groups = "drop")
  
  # Step 3: Classify each interaction using rowwise evaluation
  interaction_groups <- interaction_groups %>%
    rowwise() %>%
    mutate(Classification = {
      g <- Groups_Present  # a vector of groups where this interaction is present
      if(setequal(g, sort(groups))) {
        "Common"
      } else if(length(g) == 1) {
        paste0("Unique_", g)
      } else if(length(g) == 2) {
        paste0("Shared_", paste(g, collapse = "_"))
      } else if(length(g) == 3) {
        paste0("Shared_three_", paste(g, collapse = "_"))
      } else {
        NA_character_
      }
    }) %>%
    ungroup()
  
  # Step 4: Create a summary table: count the number of interactions for each classification 
  summary_df <- interaction_groups %>%
    group_by(Classification) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(Proportion = Count / sum(Count))
  
  # Step 5: Create a pairwise data frame that splits the Interaction column into Donor and Recipient
  pairwise_df <- interaction_groups %>%
    separate(Interaction, into = c("Donor", "Recipient"), sep = " -> ")
  
  return(list(summary = summary_df, pairwise = pairwise_df))
}






# Split the Interaction column into Donor and Recipient 
# applicable for comparisons between two sets of data (not to be used in four group timepoint comparison)
split_interactions <- function(pairwise_df) {
  pairwise_df %>%
    separate(Interaction, into = c("Donor", "Recipient"), sep = " -> ")
}

# Function: plot_interaction_venn_custom
# Description:
#   Given a data frame with columns "Donor", "Recipient", and "Classification",
#   this function creates an "Interaction" column (formatted as "Donor -> Recipient"),
#   constructs unique interaction sets for two groups ("MAM" and "Healthy"),
#   and then plots a Venn diagram comparing these sets using custom colors and aesthetics.
#
# Dataset Requirements:
#   - Donor: Character column representing the donor species.
#   - Recipient: Character column representing the recipient species.
#   - Classification: A factor or character column with values like "Common", "Unique_MAM", "Unique_Healthy".
#
# Returns:
#   A ggplot object displaying the Venn diagram.
plot_interaction_venn_custom <- function(data) {
  
  # Create the Interaction column preserving direction
  data <- data %>%
    mutate(Interaction = paste(Donor, "->", Recipient))
  
  # Build the sets for each group:
  venn_data <- list(
    MAM = data %>% 
      filter(Classification %in% c("Common", "Unique_MAM")) %>% 
      pull(Interaction) %>% 
      unique(),
    Healthy = data %>% 
      filter(Classification %in% c("Common", "Unique_Healthy")) %>% 
      pull(Interaction) %>% 
      unique()
  )
  
  # Plot the Venn diagram using ggvenn with custom aesthetics
  plot <- ggvenn(venn_data,
                 fill_color = c("darkgoldenrod1", "cyan4"),
                 stroke_size = 1,
                 set_name_size = 5,
                 text_size = 4,
                 fill_alpha = 0.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  return(plot)
}


# Function: plot_interaction_venn_four
# Description:
#   Creates a Venn diagram comparing unique directed interaction pairs across
#   four groups (e.g., "MAM_12", "MAM_24", "Healthy_12", "Healthy_24").
#
# Dataset Requirements:
#   - CLADE_A: Recipient species.
#   - CLADE_B: Donor species.
#   - Timepoint_general: Group identifier (e.g., "MAM_12", "MAM_24", "Healthy_12", "Healthy_24").
#
# Args:
#   data: A data frame with the required columns.
#   groups: A character vector specifying the groups to compare.
#           Default is c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24").
#
# Returns:
#   A ggplot object displaying the Venn diagram.
plot_interaction_venn_four <- function(data, 
                                       groups = c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24")) {
  
  # Create the Interaction column preserving direction: "CLADE_B -> CLADE_A"
  data <- data %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A))
  
  # Build a list of unique interactions for each specified group using the Timepoint_general column
  venn_list <- lapply(groups, function(g) {
    data %>%
      filter(Timepoint_general == g) %>%
      pull(Interaction) %>%
      unique()
  })
  names(venn_list) <- groups
  
  # Plot the Venn diagram using ggvenn with custom aesthetics, using default colors
  plot <- ggvenn(venn_list,
                 stroke_size = 1,
                 set_name_size = 5,
                 text_size = 4,
                 fill_alpha = 0.5) +
    theme(plot.title = element_text(hjust = 0.5, size = 14))
  
  # Draw the Venn diagram on a new grid page
  grid.newpage()
  grid.draw(plot)
  
  return(plot)
}


## Super transferome - Prevalence calculation 
# use the main dataframe that contains the CLADE_A and CLADE_B columns with Sample_ID 

calculate_prevalence_dir_spp <- function(data) {
  data %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A)) %>%  # Combine CLADE_B and CLADE_A into a unique interaction
    group_by(Interaction) %>%
    summarize(
      total_sample_ids = n_distinct(Sample_ID),  # Count unique Sample_IDs where the interaction occurs
      total_pairs = n(),                         # Total occurrences of the interaction pair
      prevalence = total_sample_ids / n_distinct(data$Sample_ID) * 100  # Calculate prevalence as a percentage
    ) %>%
    arrange(desc(prevalence))  # Arrange by highest prevalence
}


## Prevalence calculation for each Group (eg: timepoint, pre post treatement etc.)

calculate_prevalence_timepoint <- function(data, timepoint) {
  # Filter the data for the specified timepoint (e.g., "MAM_12", "MAM_24", "Healthy_12", "Healthy_24")
  data_filtered <- data %>% 
    filter(Timepoint_general == timepoint)
  
  # Calculate the total number of unique samples at that timepoint
  total_samples <- n_distinct(data_filtered$Sample_ID)
  
  # Compute prevalence for each interaction pair within this timepoint
  prevalence_table <- data_filtered %>%
    mutate(Interaction = paste(CLADE_B, "->", CLADE_A)) %>%  # create the Interaction column
    group_by(Interaction) %>%
    summarize(
      total_sample_ids = n_distinct(Sample_ID),  # number of unique samples containing the interaction
      total_pairs = n(),                         # total occurrences of the interaction
      prevalence = total_sample_ids / total_samples * 100,  # prevalence as a percentage
      .groups = "drop"
    ) %>%
    arrange(desc(prevalence))
  
  return(prevalence_table)
}


# Function: plot_species_pair_prevalence
# Requirements: data must contain Donor, Recipient, prevalence, and Category columns.
# The function removes "s__" from species names, creates a Species_Pair column,
# and plots a dot plot of prevalence using a custom color palette.
plot_species_pair_prevalence <- function(data, palette = categories_palette) {
  data <- data %>%
    mutate(
      Donor         = gsub("^s__", "", Donor),
      Recipient     = gsub("^s__", "", Recipient),
      Species_Pair  = paste(Donor, Recipient, sep = " - ")
    ) %>%
    arrange(Species_Pair)
  
  ggplot(data, aes(x = prevalence, y = reorder(Species_Pair, prevalence), color = Category)) +
    geom_point(size = 3) +
    scale_color_manual(values = palette) +
    labs(
      title = "Prevalence of Species Pairs by Category",
      x     = "Prevalence (%)",
      y     = "Species Pair"
    ) +
    theme_minimal() +
    theme(
      plot.title           = element_text(hjust = 0.5, size = 14),
      axis.text.y          = element_text(size = 9),
      axis.title.x         = element_text(size = 9),
      axis.title.y         = element_text(size = 9),
      # grid lines
      panel.grid.major.x   = element_line(color = "grey90"),
      panel.grid.minor     = element_blank(),
      panel.grid.major.y   = element_line(color = "grey90"),
      # draw axis lines
      axis.line.x          = element_line(color = "black"),
      axis.line.y          = element_line(color = "black")
    )
}



# Function: plot_species_abundance
# Requirements: The input data should have columns:
#   - Sample_ID: Unique identifier for each sample (rows in species_long)
#   - Taxa: Species names (e.g., "s__Coprococcus_comes")
#   - RA: Relative abundance of each species in each sample.
# The metadata (formatted_metadata_m4efad) should contain a Sample_ID column and a Timepoint_general column.
#
# Args:
#   species_name: A character string specifying the species of interest.
#   data: The long-format species abundance dataset (default: species_long).
#   metadata: The formatted metadata dataset (default: formatted_metadata_m4efad).
#   palette: A named vector of colors for each Timepoint_general category (default: categories_palette).
#
# Returns:
#   A ggplot object showing the boxplot of RA by Timepoint_general for the specified species.
plot_species_abundance <- function(species_name,
                                   data = species_long,
                                   metadata = formatted_metadata_m4efad,
                                   palette = categories_palette) {
  
  # Merge species abundance data with metadata by Sample_ID
  merged_data <- data %>%
    inner_join(metadata, by = "Sample_ID")
  
  # Filter for the species of interest
  species_data <- merged_data %>%
    filter(Taxa == species_name)
  
  # Create the plot with contaminant_conditions_plot aesthetics
  p <- ggplot(species_data, aes(x = factor(Timepoint_general), y = RA, fill = factor(Timepoint_general))) +
    geom_quasirandom(
      dodge.width = 0.75,
      shape       = 16,
      size        = 1,
      alpha       = 0.7,
      aes(color = factor(Timepoint_general))
    ) +
    geom_boxplot(
      outlier.colour = NA,
      alpha          = 0.5,
      color          = "black",
      linewidth      = 0.3,
      position       = position_dodge(width = 0.75)
    ) +
    scale_fill_manual(values = palette) +
    scale_color_manual(values = palette) +
    labs(
      title = paste("Relative Abundance of", species_name, "by Timepoint"),
      x     = "",
      y     = "Relative Abundance (RA)"
    ) +
    theme_bw() +
    theme(
      panel.border    = element_blank(),
      panel.grid      = element_blank(),
      axis.text       = element_text(),
      axis.line       = element_line(),
      axis.title      = element_text(size = 9),
      axis.title.y    = element_text(size = 9),
      axis.title.x    = element_text(size = 9),
      legend.position = "none"
    )
  
  return(p)
}


# Function: run_species_stats
# Description:
#   This function merges the species abundance data with metadata,
#   filters for the specified species, and runs a Kruskal-Wallis test 
#   to assess differences in relative abundance (RA) across Timepoint_general.
#   It then performs pairwise Dunn's tests with BH correction.
#
# Requirements:
#   - The species_long dataset must contain: Sample_ID, Taxa, and RA.
#   - The formatted_metadata_m4efad dataset must contain: Sample_ID and Timepoint_general.
#
# Args:
#   species_name: A character string with the species name of interest.
#   data: The long-format species abundance dataset (default: species_long).
#   metadata: The metadata dataset (default: formatted_metadata_m4efad).
#
# Returns:
#   A list with two elements:
#     - kw: The result of the global Kruskal-Wallis test.
#     - pairwise: The results of the pairwise Dunn's tests.
run_species_stats <- function(species_name, data = species_long, metadata = formatted_metadata_m4efad) {
  
  # Merge species abundance data with metadata by Sample_ID
  merged_data <- data %>%
    inner_join(metadata, by = "Sample_ID") %>%
    filter(Taxa == species_name)
  
  # Run the global Kruskal-Wallis test for RA across Timepoint_general
  kw_result <- merged_data %>%
    kruskal_test(RA ~ Timepoint_general)
  
  # Run pairwise Dunn's test for RA across Timepoint_general with BH adjustment
  pairwise_result <- merged_data %>%
    dunn_test(RA ~ Timepoint_general, p.adjust.method = "BH")
  
  return(list(kw = kw_result, pairwise = pairwise_result))
}

# load data ---------------------------------------------------------------

# dataset with directional species level HGTs for the whole dataset
load("data/processed/clean_outputs/dir_spp_hgts_m4efad.RData")

# contains the metaphlan 3 taxpnomy table post processing (col names and placemet)
load("data/processed/species_all_samples.Rdata")

# contains metadata with correct column names 
# Note: once loaded dataset is called: formatted_metadata_m4efad
load("data/processed/clean_outputs/metadata_m4efad.Rdata")


# running -----------------------------------------------------------------


# Investigate the overall species pairs unique and shared between  --------

# Suppose your dataset is called 'dir_spp_hgts_m4efad'
# and it has a column 'Group' that may contain values like "MAM" or "Healthy" (or others).
# To compare interactions between "MAM" and "Healthy", run:
results <- compare_interaction_pairs(dir_spp_hgts_m4efad, group1 = "MAM", group2 = "Healthy")
print(results$summary) # provides proportions of unique and common interactions between groups
print(results$pairwise) # provides species pairs identities 

# Create separate dataset with species pairs identities - to be used downstream 
species_pairs_m4efad <- split_interactions(results$pairwise)

# plot Venn diagram
# Assuming your dataset is called species_pairs_m4efad:
species_pairs_m4efad_plot <- plot_interaction_venn_custom(species_pairs_m4efad)

# svg format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_plot.svg", 
       plot = species_pairs_m4efad_plot, units = "cm", width=10, height=8)
# pdf format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_plot.pdf", 
       plot = species_pairs_m4efad_plot, units = "cm", width=10, height=8)

# Investigate the species pairs unique and shared between timepoints of groups  --------

# Timepoint comparison 
# Compare interaction pairs at timepoint "12" between MAM and Healthy groups:
result_12 <- compare_interaction_pairs_timepoint(dir_spp_hgts_m4efad, timepoint = "12", group1 = "MAM", group2 = "Healthy")
print(result_12$summary)
print(result_12$pairwise)

# Create separate dataset with species pairs identities - to be used downstream 
species_pairs_m4efad_12 <- split_interactions(result_12$pairwise)

# Compare interaction pairs at timepoint "24" between MAM and Healthy groups:
result_24 <- compare_interaction_pairs_timepoint(dir_spp_hgts_m4efad, timepoint = "24", group1 = "MAM", group2 = "Healthy")
print(result_24$summary)
print(result_24$pairwise)

# Create separate dataset with species pairs identities - to be used downstream 
species_pairs_m4efad_24 <- split_interactions(result_24$pairwise)

# plot Venn diagram - 12 months 
# Assuming your dataset is called species_pairs_m4efad:
species_pairs_m4efad_12_plot <- plot_interaction_venn_custom(species_pairs_m4efad_12)

# svg format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_12_plott.svg", 
       plot = species_pairs_m4efad_12_plot, units = "cm", width=10, height=8)
# pdf format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_12_plot.pdf", 
       plot = species_pairs_m4efad_12_plot, units = "cm", width=10, height=8)



# plot Venn diagram - 24 months
# Assuming your dataset is called species_pairs_m4efad:
species_pairs_m4efad_24_plot <- plot_interaction_venn_custom(species_pairs_m4efad_24)

# svg format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_24_plot.svg", 
       plot = species_pairs_m4efad_24_plot, units = "cm", width=10, height=8)
# pdf format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_24_plot.pdf", 
       plot = species_pairs_m4efad_24_plot, units = "cm", width=10, height=8)


# Identify the core species pair that are shared between all timepoints of all groups  --------

## Comparing species pairs interaction between four timepoints of two groups 

# Assuming your dataset is called dir_spp_hgts_m4efad:
result_four_groups <- compare_species_pairs_four_groups(dir_spp_hgts_m4efad)
print(result_four_groups$summary)
print(result_four_groups$pairwise)

# Convert the outputs to data frames (if needed)
summary_four_groups <- as.data.frame(result_four_groups$summary)
pairwise_four_groups <- as.data.frame(result_four_groups$pairwise)


## plot Venn diagram - all four timepoints
# Assuming your dataset is called dir_spp_hgts_m4efad and it contains the required columns:
# colors are assigned by ggvenn default 
species_pairs_m4efad_four_plot <- plot_interaction_venn_four(dir_spp_hgts_m4efad)

# svg format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_four_plot.svg", 
       plot = species_pairs_m4efad_four_plot, units = "cm", width=17, height=17)
# pdf format 
ggsave(file = "plots/chp3/clean_plots/species_pairs_m4efad_four_plot.pdf", 
       plot = species_pairs_m4efad_four_plot, units = "cm", width=15, height=15)



# Identify core shared species pairs  -------------------------------------

core_shared_interactions <- pairwise_four_groups %>%
  filter(Classification == "Common")



# Investigate the prevalence of species pairs  ----------------------------

# For each Group (timepoint)
# For MAM_12:
prevalence_MAM_12 <- calculate_prevalence_timepoint(dir_spp_hgts_m4efad, "MAM_12")
print(prevalence_MAM_12)

# For MAM_24:
prevalence_MAM_24 <- calculate_prevalence_timepoint(dir_spp_hgts_m4efad, "MAM_24")
print(prevalence_MAM_24)

# For Healthy_12:
prevalence_Healthy_12 <- calculate_prevalence_timepoint(dir_spp_hgts_m4efad, "Healthy_12")
print(prevalence_Healthy_12)

# For Healthy_24:
prevalence_Healthy_24 <- calculate_prevalence_timepoint(dir_spp_hgts_m4efad, "Healthy_24")
print(prevalence_Healthy_24)


# Combine and modify datasets

# Combine datasets with the appropriate Category labels
combined_prevalence <- bind_rows(
  prevalence_MAM_12 %>% mutate(Category = "MAM_12"),
  prevalence_MAM_24 %>% mutate(Category = "MAM_24"),
  prevalence_Healthy_12 %>% mutate(Category = "Healthy_12"),
  prevalence_Healthy_24 %>% mutate(Category = "Healthy_24")
)

# Process the combined dataset: split the Interaction column into Donor and Recipient columns
combined_prevalence <- combined_prevalence %>%
  separate(Interaction, into = c("Donor", "Recipient"), sep = "->", remove = TRUE) %>%
  mutate(
    Donor = trimws(Donor),
    Recipient = trimws(Recipient)
  )

# Filter the dataset for species pairs with prevalence >= 10% (standard for MB studies)
combined_prev_10_percent <- combined_prevalence %>%
  filter(prevalence >= 10)

# Core shared species pairs with their prevalence in each Group 
core_shared_prevalence <- combined_prevalence %>%
  semi_join(core_shared_interactions, by = c("Donor", "Recipient"))

# Filter the core shared species pairs to only include species pairs with a prevalence of 10%
core_shared_prev_10_percent <- core_shared_prevalence %>% 
  filter(prevalence >= 10)

# modify column entries for plotting
core_shared_prev_10_percent <- core_shared_prev_10_percent %>%
  mutate(
    Donor     = Donor     %>%
      str_remove("^s__") %>%
      str_replace_all("_", " "),
    Recipient = Recipient %>%
      str_remove("^s__") %>%
      str_replace_all("_", " ")
  )

# for viewing only 
grouped_species_pairs <- core_shared_prev_10_percent %>%
  arrange(Donor, Recipient)

# plot the differences in species pairs with >10% prevalence across categories 
# Assuming 'core_shared_prev_10_percent' is your dataset with columns: Donor, Recipient, prevalence, and Category.
prev_plot <- plot_species_pair_prevalence(core_shared_prev_10_percent)

# svg format 
# ggsave(file = "plots/chp3/clean_plots/prev_plot.svg", 
#       plot = prev_plot, units = "cm", width=20, height=10)
# pdf format 
ggsave(file = "plots/chp3/clean_plots/prev_plot.pdf", 
      plot = prev_plot, units = "cm", width=20, height=10)




# investigate the abundance of the core prevalent species pairs -----------

# Modify the species table
# convert species table to long format 
species_long <- species_all_samples %>%
  rownames_to_column(var = "Sample_ID") %>%  # Convert rownames to a column called "Sample_ID"
  pivot_longer(cols = -Sample_ID,           # Pivot all columns except Sample_ID
               names_to = "Taxa",           # New column "Taxa" will contain species names (old column names)
               values_to = "RA")            # New column "RA" will contain the relative abundance values


# remove all RA = 0 rows 
species_long <- species_long %>%
  filter(RA != 0)

# remove the 15 months timepoint from the metadata
formatted_metadata_m4efad <- formatted_metadata_m4efad %>%
  filter(Timepoint_general != "MAM_15")

# plot the abundance of selected species across groups 
# species that are the most prevalent and are part of the conserved species across groups have been selected

# List of species to plot
species_list <- c("s__Coprococcus_comes", "s__Collinsella_aerofaciens", 
                  "s__Prevotella_stercorea", "s__Prevotella_copri", "s__Bifidobacterium_longum", 
                  "s__Faecalibacterium_prausnitzii")

# Generate and print plots for each species
plots <- lapply(species_list, plot_species_abundance)
plots  # This will display a list of ggplot objects, one for each species

# Assuming `plots` is a list of ggplot objects returned by plot_species_abundance()
# and `species_list` is defined as:
# species_list <- c(
#   "s__Coprococcus_comes", "s__Collinsella_aerofaciens",
#   "s__Prevotella_stercorea", "s__Prevotella_copri",
#   "s__Bifidobacterium_longum", "s__Faecalibacterium_prausnitzii"
# )

# Create output directory if it doesn't exist
out_dir <- "plots/chp3/clean_plots/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Loop over each species and save its plot
for (i in seq_along(species_list)) {
  species <- species_list[i]
  p <- plots[[i]]
  
  # Clean up species name for filename:
  #   remove "s__", replace underscores with hyphens
  species_clean <- species %>%
    gsub("^s__", "", .) %>%
    gsub("_", "-", .)
  
  # Construct a descriptive filename
  file_name <- file.path(
    out_dir,
    paste0(species_clean, "_abundance.pdf")
  )
  
  # Save the plot
  ggsave(
    filename = file_name,
    plot     = p,
    units    = "cm",
    width    = 5,
    height   = 6
  )
}


# Run statistics for s__Coprococcus_comes
stats_coprococcus <- run_species_stats("s__Coprococcus_comes")
print(stats_coprococcus$kw)
print(stats_coprococcus$pairwise)

# Run statistics for s__Collinsella_aerofaciens
stats_collinsella <- run_species_stats("s__Collinsella_aerofaciens")
print(stats_collinsella$kw)
print(stats_collinsella$pairwise)

# Run statistics for s__Prevotella_stercorea
stats_prevotella_stercorea <- run_species_stats("s__Prevotella_stercorea")
print(stats_prevotella_stercorea$kw)
print(stats_prevotella_stercorea$pairwise)

# Run statistics for s__Prevotella_copri
stats_prevotella_copri <- run_species_stats("s__Prevotella_copri")
print(stats_prevotella_copri$kw)
print(stats_prevotella_copri$pairwise)

# Run statistics for s__Bifidobacterium_longum
stats_Bifidobacterium_longum <- run_species_stats("s__Bifidobacterium_longum")
print(stats_Bifidobacterium_longum$kw)
print(stats_Bifidobacterium_longum$pairwise)

# Run statistics for s__Faecalibacterium_prausnitzii
stats_Faecalibacterium_prausnitzii <- run_species_stats("s__Faecalibacterium_prausnitzii")
print(stats_Faecalibacterium_prausnitzii$kw)
print(stats_Faecalibacterium_prausnitzii$pairwise)

# save data ---------------------------------------------------------------

# contains all species pairs unique and shared between Groups
save(species_pairs_m4efad, file = "data/processed/clean_outputs/species_pairs_m4efad.RData")
# contains 12 months species pairs unique and shared between Groups
save(species_pairs_m4efad_12, file = "data/processed/clean_outputs/species_pairs_m4efad_12.RData")
# contains 24 months species pairs unique and shared between Groups
save(species_pairs_m4efad_24, file = "data/processed/clean_outputs/species_pairs_m4efad_24.RData")
# contains prevalence of species pairs for the super transferome and each timepoint
save(combined_prevalence, file = "data/processed/clean_outputs/combined_prevalence.RData")



