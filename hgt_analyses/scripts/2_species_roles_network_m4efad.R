
## code summary:
# create a super HGT network with all pooled samples from all trial cohorts (all groups of m4efad)
# calculate donor-receiver ratio for all nodes in 10k random networks and look at the distribution 
# select the donor-receiver ratio threshold (to define sources/sinks) according to the tails of the distribution
# define species as sources, conduits or sinks based on their role in HGT interactions across cohort samples (subset metadata) 
# calculate the proportion of species with each role across cohorts
# investigate the conservation of network roles across cohorts

## NOTE: all large data files are stored in the 'save data files' section at the end of this script

# load packages -----------------------------------------------------------
library(tidyverse)
library(igraph)
library(patchwork)
library(ComplexUpset)
library(broom)


# load data ---------------------------------------------------------------

# HGT events between species for sample subset (individual donors, 
# recipients with pre- and post-intervention), contains column with Directional == "Yes" / "No"
load("data/processed/clean_outputs/spp_dir_hgts_m4efad.RData")

# taxonomy lookup table
all_hgt_taxa_distinct <- read.table("data/processed/clean_outputs/hgt_taxa_lookup.txt", sep = "\t")
unique_phyla <- unique(all_hgt_taxa_distinct$Phylum)
na_phylum_rows <- all_hgt_taxa_distinct %>% filter(is.na(Phylum))

# Convert unique phyla into a dataframe
phyla_list_m4fad <- data.frame(Phylum = unique_phyla)

# formatted metadata, note: sample IDs must match extracted WAAFLE sample IDs
load("data/processed/clean_outputs/sample_subsets.RData")




# load palettes -----------------------------------------------------------

role_palette <- c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")

phylum_palette <- c(
  "Firmicutes" = "#009E73", 
  "Bacteroidetes" = "#56B4E9", 
  "Actinobacteria" = "#E69F00", 
  "Proteobacteria" = "#CC79A7", 
  "Verrucomicrobia" = "#F0E442", 
  "Microsporidia" = "#9999CC", 
  "Spirochaetes" = "#B65D40",
  "Ascomycota" = "olivedrab3", 
  "Fusobacteria" = "maroon2", 
  "Synergistetes" = "wheat2", 
  "Apicomplexa" = "thistle2", 
  "Chloroflexi" = "steelblue4", 
  "Tenericutes" = "mediumorchid",
  "Chlamydiae" = "dodgerblue3", 
  "Eukaryota_noname" = "gray50", 
  "Candidatus_Saccharibacteria" = "darkorange2",
  "NA" = "black"
)




# load functions  ---------------------------------------------------------

# subset directed HGTs between species by sample
get_hgts_per_sample <- function(dir_spp_hgts, metadata){
  ## dir_spp_hgts: dataframe containing WAAFLE HGT events for each sample (df)
  ## metadata: subset of sample metadata, must contain column 'Sample_ID' (df)
  dir_spp_hgts %>%
    select(CLADE_B, CLADE_A, SYNTENY, ANNOTATIONS.UNIREF50, ANNOTATIONS.UNIREF90, Sample_ID) %>% # simplify dataframe
    inner_join(metadata, by = c("Sample_ID")) %>% # join with metadata by extracted sample ID
    select(Participant_ID, Sample_ID, Group, CLADE_B, CLADE_A, SYNTENY, ANNOTATIONS.UNIREF50, ANNOTATIONS.UNIREF90) %>% # simplify dataframe
    rename(Donor = CLADE_B, # renaming as events are directional
           Recipient = CLADE_A)
}

# count the total number of edges (HGT interactions) of each species across the sample subset and filter low counts 
get_spp_edge_count <- function(hgts_per_sample){
  # Step 1: Count the number of edges (HGT interactions) for each species
  edge_counts <- hgts_per_sample %>%
    group_by(Species = Donor) %>%
    summarise(Donor_count = n(), .groups = "drop") %>%
    full_join(
      hgts_per_sample %>%
        group_by(Species = Recipient) %>%
        summarise(Recipient_count = n(), .groups = "drop"),
      by = "Species"
    ) %>%
    # Replace NAs with 0
    mutate(Donor_count = ifelse(is.na(Donor_count), 0, Donor_count),
           Recipient_count = ifelse(is.na(Recipient_count), 0, Recipient_count))
  
  # Step 2: Calculate thresholds using median and MAD (k = 1)
  thresholds <- edge_counts %>%
    summarise(
      median_donor = median(Donor_count, na.rm = TRUE),
      mad_donor = mad(Donor_count, na.rm = TRUE),
      threshold_donor = median_donor + 1 * mad_donor,
      median_recipient = median(Recipient_count, na.rm = TRUE),
      mad_recipient = mad(Recipient_count, na.rm = TRUE),
      threshold_recipient = median_recipient + 1 * mad_recipient
    )
  
  donor_thresh <- thresholds$threshold_donor
  recipient_thresh <- thresholds$threshold_recipient
  
  # Step 3: Filter out rows with low edge counts (i.e. counts below thresholds)
  filtered_edge_counts <- edge_counts %>%
    filter(Donor_count >= donor_thresh | Recipient_count >= recipient_thresh)
  
  return(filtered_edge_counts)
}

# get network node layout for plotting
network_node_coordinates <- function(igraph_object, species_phylum_table){
  ## igraph_object: igraph object for the HGT network (igraph)
  ## species_phylum_table: lookup table containing phylum classification for each species (df)
  set.seed(1234) # set seed
  node_layout <- layout_with_fr(igraph_object) # Ruchterman-Reingold layout
  node_coords <- data.frame(x = node_layout[, 1], # node coordinates
                            y = node_layout[, 2],
                            Species = V(igraph_object)$name)
  node_coords %>% # return node coordinates
    left_join(species_phylum_table %>% select(Species, Phylum), by = "Species") # get phylum of each species
}

# get network edge layout for plotting
network_edge_coordinates <- function(igraph_object, node_coords){
  ## igraph_object: igraph object for the HGT network (igraph)
  ## node_coords: output from network_node_coordinates() (df)
  set.seed(1234) # set seed
  edge_coords <- as.data.frame(as_edgelist(igraph_object)) # get edge list
  colnames(edge_coords) <- c("from", "to") # rename columns
  # add coordinates to edges dataframe
  edge_coords$x1 <- node_coords$x[match(edge_coords$from, node_coords$Species)]
  edge_coords$y1 <- node_coords$y[match(edge_coords$from, node_coords$Species)]
  edge_coords$x2 <- node_coords$x[match(edge_coords$to, node_coords$Species)]
  edge_coords$y2 <- node_coords$y[match(edge_coords$to, node_coords$Species)]
  edge_coords # return edge coordinates
}

# make a HGT network plot from a list of directed HGTs between species for each sample in a particular cohort 
make_hgt_network_plot <- function(dir_spp_hgts, species_phylum_table, plot_title) {
  ## dir_spp_hgts: directed HGT events for each sample in a cohort (df)
  ## species_phylum_table: lookup table containing phylum classification for each species (df)
  ## plot_title: title for the network plot (str)
  
  set.seed(1234) # for reproducibility
  
  # Step 1: Preprocess data - select and clean donor/recipient names
  dir_spp_hgts <- dir_spp_hgts %>%
    select(Donor, Recipient) %>%
    mutate(
      Donor = str_remove(Donor, "s__"),
      Donor = str_replace_all(Donor, "_", " "),
      Recipient = str_remove(Recipient, "s__"),
      Recipient = str_replace_all(Recipient, "_", " ")
    )
  
  # Step 2: Calculate edge counts per species (as Donor and Recipient)
  edge_counts <- dir_spp_hgts %>%
    group_by(Species = Donor) %>%
    summarise(Donor_count = n(), .groups = "drop") %>%
    full_join(
      dir_spp_hgts %>%
        group_by(Species = Recipient) %>%
        summarise(Recipient_count = n(), .groups = "drop"),
      by = "Species"
    ) %>%
    mutate(
      Donor_count = ifelse(is.na(Donor_count), 0, Donor_count),
      Recipient_count = ifelse(is.na(Recipient_count), 0, Recipient_count)
    )
  
  # Step 3: Calculate thresholds using median and MAD (with k = 1)
  thresholds <- edge_counts %>%
    summarise(
      median_donor = median(Donor_count, na.rm = TRUE),
      mad_donor = mad(Donor_count, na.rm = TRUE),
      threshold_donor = median_donor + 1 * mad_donor,
      median_recipient = median(Recipient_count, na.rm = TRUE),
      mad_recipient = mad(Recipient_count, na.rm = TRUE),
      threshold_recipient = median_recipient + 1 * mad_recipient
    )
  
  donor_thresh <- thresholds$threshold_donor
  recipient_thresh <- thresholds$threshold_recipient
  
  # Step 4: Filter species based on threshold criteria
  filtered_species <- edge_counts %>%
    filter(Donor_count >= donor_thresh | Recipient_count >= recipient_thresh)
  
  # Filter the HGT events to only include edges where either donor or recipient is in filtered_species
  dir_spp_hgts <- dir_spp_hgts %>%
    filter(Donor %in% filtered_species$Species | Recipient %in% filtered_species$Species)
  
  # Step 5: Create the directed network graph from the filtered edges
  directed_network <- graph_from_data_frame(d = dir_spp_hgts[, c(1, 2)], directed = TRUE)
  
  # Weight network edges by the number of HGT interactions
  E(directed_network)$weight <- count_multiple(directed_network)
  
  # Step 6: Get node and edge layouts for plotting (using pre-defined helper functions)
  node_coords <- network_node_coordinates(directed_network, species_phylum_table)
  edge_coords <- network_edge_coordinates(directed_network, node_coords)
  
  # Step 7: Plot the HGT network to visualize interactions between taxa
  ggplot() +
    geom_segment(data = edge_coords,
                 aes(x = x1, y = y1, xend = x2, yend = y2),
                 arrow = arrow(type = "closed", length = unit(1.5, "mm")),
                 colour = "grey") +
    geom_point(data = node_coords, aes(x = x, y = y, fill = Phylum),
               size = 2, shape = 21) +
    scale_fill_manual(values = phylum_palette, name = "Phylum") +
    ggtitle(plot_title) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 12),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.background = element_rect(color = "black", size = 0.4))
}

# randomise 10,000 HGT networks and extract edge lists
randomise_hgt_networks <- function(igraph_object){
  ## igraph_object: igraph object for the super HGT network created from all pooled samples (igraph)
  n_nodes <- gorder(igraph_object) # calculate the number of nodes in the true network
  n_edges <- gsize(igraph_object) # calculate the number of edges in the true network
  set.seed(1234) # set the seed for randomisation
  results <- tibble() # initialise an empty tibble
  for (i in 1:10000) { # loop 10,000 times
    random_graph <- erdos.renyi.game(n = n_nodes, p.or.m = n_edges, type = "gnm", directed = TRUE) # create a random network based on true network
    edges <- as_tibble(as_edgelist(random_graph, names = TRUE)) # extract the source-target edge list from the random network
    colnames(edges) <- c("Donor", "Recipient") # rename edge list columns
    edges <- edges %>%
      mutate(loop_counter = i) # add a loop counter
    results <- bind_rows(results, edges) # bind the results output
  }
  as.data.frame(results) # return overall result as dataframe
}

# calculate the donor-receiver ratio for all nodes in 10,000 random networks, adding pseudo count of 1 to all donor/recipient edge counts
randomised_dr_ratio <- function(randomised_network_edges){
  ## randomised_network_edges: output from randomise_hgt_networks() (df)
  donor_count <- randomised_network_edges %>% # calculate randomised out-edges for each node
    group_by(loop_counter, Species = Donor) %>%
    summarise(Donor_count = n()) %>%
    ungroup()
  recipient_count <- randomised_network_edges %>% # calculate randomised in-edges for each node
    group_by(loop_counter, Species = Recipient) %>%
    summarise(Recipient_count = n()) %>%
    ungroup()
  joined_count <- full_join(donor_count, recipient_count, by = c("Species", "loop_counter")) # join results
  # after the join, nodes with 0 donor count or 0 recipient count in a given loop will be represented by NA
  joined_count %>% # calculate donor-receiver ratio for each node in each randomised network
    replace(is.na(.), 0) %>% # replace NAs with 0
    mutate(Donor_count_pseudo = Donor_count+1, # add pseudo count of 1 to all out edges
           Recipient_count_pseudo = Recipient_count+1) %>% # add a pseudo count of 1 to all in edges
    mutate(dr_ratio = Donor_count_pseudo / Recipient_count_pseudo) # calculate donor-receiver ratio for all nodes
}

# define species network roles based on donor-receiver ratio thresholds
define_spp_roles <- function(spp_edge_count, lower_threshold, upper_threshold, type, trial_name, cohort_name){
  ## spp_edge_count: output from get_spp_edge_count() contains observed in and out edge counts for each species (df)
  ## lower_threshold: lower threshold for species HGT donor-receiver ratio to define sinks vs conduits (num)
  ## upper_threshold: upper threshold for species HGT donor-receiver ratio to define sources vs conduits (num)
  ## trial_name: trial name (chr)
  ## cohort_name: cohort name (chr)
  spp_edge_count_pseudo <- spp_edge_count %>% # NAs are already replaced with 0
    mutate(Donor_count_pseudo = Donor_count+1, # add pseudo count of 1 to all out edges
           Recipient_count_pseudo = Recipient_count+1) %>% # add a pseudo count of 1 to all in edges
    mutate(dr_ratio = Donor_count_pseudo / Recipient_count_pseudo)
  # define species network roles based on DR ratio thresholds defined from simulated networks
  spp_type <- spp_edge_count_pseudo %>%
    mutate(Type = case_when(dr_ratio <= lower_threshold ~ "Sink", # DR ratio less than or equal to lower threshold = sink
                            dr_ratio >= upper_threshold ~ "Source", # DR ratio greater than or equal to upper threshold = source
                            dr_ratio > lower_threshold & dr_ratio < upper_threshold ~ "Conduit")) %>% # DR ratio in between = conduit
    mutate(Trial = trial_name, # define trial name
           Cohort = cohort_name) # define cohort name
  spp_type # return species role types
}

# calculate network metrics 
calculate_network_metrics_randomized <- function(dir_spp_hgts, n_iter = 1000) {
  # Preprocess data
  processed_data <- dir_spp_hgts %>%
    select(Donor, Recipient) %>%
    mutate(Donor = str_remove(Donor, "s__"),
           Donor = str_replace_all(Donor, "_", " "),
           Recipient = str_remove(Recipient, "s__"),
           Recipient = str_replace_all(Recipient, "_", " "))
  
  # Create the observed directed network
  obs_network <- graph_from_data_frame(d = processed_data, directed = TRUE)
  
  # Calculate number of nodes and edges
  n_nodes <- gorder(obs_network)
  n_edges <- gsize(obs_network)
  
  # Compute observed metrics
  obs_undirected <- as.undirected(obs_network, mode = "collapse")
  obs_community <- cluster_louvain(obs_undirected)
  observed_modularity <- modularity(obs_community)
  observed_clustering <- transitivity(obs_network, type = "global")
  
  # Initialize vectors to store metrics from randomized networks
  modularity_values <- numeric(n_iter)
  clustering_values <- numeric(n_iter)
  
  set.seed(1234)  # Ensure reproducibility
  
  # Generate random networks with same number of nodes and edges, compute metrics for each
  for (i in 1:n_iter) {
    random_net <- erdos.renyi.game(n = n_nodes, p.or.m = n_edges, type = "gnm", directed = TRUE)
    undirected_random <- as.undirected(random_net, mode = "collapse")
    community_rand <- cluster_louvain(undirected_random)
    
    modularity_values[i] <- modularity(community_rand)
    clustering_values[i] <- transitivity(random_net, type = "global")
  }
  
  # Calculate the average and standard deviation of the metrics across random networks
  avg_random_modularity <- mean(modularity_values)
  sd_random_modularity <- sd(modularity_values)
  avg_random_clustering <- mean(clustering_values)
  sd_random_clustering <- sd(clustering_values)
  
  # Return a list containing observed metrics and the random distributions summary
  return(list(
    observed_modularity = observed_modularity,
    observed_clustering = observed_clustering,
    avg_random_modularity = avg_random_modularity,
    sd_random_modularity = sd_random_modularity,
    avg_random_clustering = avg_random_clustering,
    sd_random_clustering = sd_random_clustering,
    random_modularity_distribution = modularity_values,
    random_clustering_distribution = clustering_values
  ))
}

# Permutation function
run_permutation_test <- function(metrics_result, alternative = "greater") {
  # metrics_result should be a list with:
  # - observed_modularity, observed_clustering
  # - random_modularity_distribution, random_clustering_distribution
  
  n_iter <- length(metrics_result$random_modularity_distribution)
  
  if (alternative == "greater") {
    p_modularity <- mean(metrics_result$random_modularity_distribution >= metrics_result$observed_modularity)
    p_clustering <- mean(metrics_result$random_clustering_distribution >= metrics_result$observed_clustering)
  } else if (alternative == "less") {
    p_modularity <- mean(metrics_result$random_modularity_distribution <= metrics_result$observed_modularity)
    p_clustering <- mean(metrics_result$random_clustering_distribution <= metrics_result$observed_clustering)
  } else if (alternative == "two.sided") {
    mod_mean <- mean(metrics_result$random_modularity_distribution)
    clust_mean <- mean(metrics_result$random_clustering_distribution)
    p_modularity <- mean(abs(metrics_result$random_modularity_distribution - mod_mean) >= 
                           abs(metrics_result$observed_modularity - mod_mean))
    p_clustering <- mean(abs(metrics_result$random_clustering_distribution - clust_mean) >= 
                           abs(metrics_result$observed_clustering - clust_mean))
  } else {
    stop("Invalid alternative hypothesis. Choose 'greater', 'less', or 'two.sided'.")
  }
  
  # If no random network produced an extreme value, set the p-value to 1/n_iter
  if(p_modularity == 0) p_modularity <- 1 / n_iter
  if(p_clustering == 0) p_clustering <- 1 / n_iter
  
  return(list(p_modularity = p_modularity, p_clustering = p_clustering))
}

# run permutation test for network metrics where i = 1000
run_permutation_tests_for_all <- function(metrics_list, alternative = "greater", correction_method = "BH") {
  # Run the permutation test for each dataset in the list
  results <- lapply(metrics_list, function(metrics_result) {
    run_permutation_test(metrics_result, alternative = alternative)
  })
  
  # Create a data frame with the raw p-values
  results_df <- data.frame(
    Dataset = names(results),
    p_modularity = sapply(results, function(x) x$p_modularity),
    p_clustering = sapply(results, function(x) x$p_clustering),
    row.names = NULL
  )
  
  # Apply multiple test correction on each set of p-values
  results_df$p_modularity_adj <- p.adjust(results_df$p_modularity, method = correction_method)
  results_df$p_clustering_adj <- p.adjust(results_df$p_clustering, method = correction_method)
  
  return(results_df)
}

# run proportion test for independent variables 
# Helper function to run proportion tests for a pair of cohorts
run_prop_test <- function(df, cohort1, cohort2) {
  # Filter data for the two cohorts
  data_filtered <- df %>% 
    filter(Cohort %in% c(cohort1, cohort2))
  
  # Calculate total counts per cohort
  totals <- data_filtered %>%
    group_by(Cohort) %>%
    summarise(total = sum(Role_count), .groups = "drop")
  
  # Merge the totals back into the filtered data
  data_filtered <- data_filtered %>%
    left_join(totals, by = "Cohort") %>%
    rename(cohort_total = total)
  
  # For each role (Type), run a proportion test comparing the two cohorts
  prop_tests <- data_filtered %>%
    group_by(Type) %>%
    summarise(
      count_cohort1 = Role_count[Cohort == cohort1],
      count_cohort2 = Role_count[Cohort == cohort2],
      total_cohort1 = cohort_total[Cohort == cohort1],
      total_cohort2 = cohort_total[Cohort == cohort2],
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      test = list(prop.test(c(count_cohort1, count_cohort2), 
                            c(total_cohort1, total_cohort2))),
      tidy_test = list(broom::tidy(test))
    ) %>%
    unnest(tidy_test)
  
  return(prop_tests)
}


# make a presence absence matrix with colnames = trial cohort, and one column with the species, for each network role type
format_upset_data <- function(joined_spp_roles, role_type){
  ## joined_spp_roles: joined species role data for all trial cohorts, with minimum 'Species', 'Type' and 'Trial_cohort' columns (df)
  ## role_type: network species role (str)
  joined_spp_roles %>%
    filter(Type == role_type) %>%
    select(Species, Trial_cohort) %>%
    mutate(Present = 1L) %>%
    pivot_wider(names_from = Trial_cohort, values_from = Present, values_fill = 0L)
}


# make complex upset plot for species with each network role
make_upset_plot <- function(
    upset_data,
    trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4,
    trial_1_colour, trial_2_colour,
    species_role_name
) {
  # upset_data: data frame prepared for ComplexUpset
  # trial_cohort_1–trial_cohort_4: names of the four cohort columns
  # trial_1_colour, trial_2_colour: hex codes for the two trial colours
  # species_role_name: label for the set‐size y‐axis (if you re‐enable it)
  
  ComplexUpset::upset(
    upset_data,
    rev(c(trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4)),
    
    queries = list(
      upset_query(set = trial_cohort_1, fill = trial_1_colour),
      upset_query(set = trial_cohort_2, fill = trial_1_colour),
      upset_query(set = trial_cohort_3, fill = trial_1_colour),
      upset_query(set = trial_cohort_4, fill = trial_2_colour)
    ),
    
    height_ratio = 0.8,
    width_ratio  = 0.2,
    name         = NULL,
    
    base_annotations = list(
      ' ' = (
        intersection_size(
          counts = FALSE,
          bar_number_threshold = 1,
          width = 0.8
        ) +
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
          theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = 'black')
          )
      )
    ),
    
    stripes    = 'white',
    matrix     = intersection_matrix(
      geom = geom_point(shape = 'circle filled', size = 2, stroke = 0.45)
    ),
    set_sizes  = FALSE,
    sort_sets        = FALSE,
    sort_intersections = 'descending'
  )
}

# find conserved species network roles across all trial cohorts
conserved_network_roles <- function(upset_data){
  ## upset_data: output from format_upset_data() (df)
  upset_data %>%
    rowwise() %>%
    filter(sum(c_across(2:5)) == 4) %>% # find species with role in all 5 trial cohorts
    ungroup() %>%
    select(Species) %>%
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
}


compare_species_roles2 <- function(data, cohort1, cohort2) {
  # Filter data for the two cohorts and only for the roles of interest
  data_filtered <- data %>%
    filter(Cohort %in% c(cohort1, cohort2)) %>%
    filter(Type %in% c("Source", "Sink", "Conduit"))
  
  # For each cohort and each role, get the list of unique species
  role_lists <- data_filtered %>%
    group_by(Cohort, Type) %>%
    summarise(Species_List = list(unique(Species)), .groups = "drop") %>%
    pivot_wider(names_from = Cohort, values_from = Species_List)
  
  # For each role, calculate shared and unique species, and the proportion of shared species
  results <- role_lists %>%
    rowwise() %>%
    mutate(
      common = length(intersect(get(cohort1), get(cohort2))),
      unique_cohort1 = length(setdiff(get(cohort1), get(cohort2))),
      unique_cohort2 = length(setdiff(get(cohort2), get(cohort1))),
      total_union = common + unique_cohort1 + unique_cohort2,
      prop_common = ifelse(total_union > 0, common / total_union, NA_real_)
    ) %>%
    ungroup() %>%
    mutate(Comparison = paste(cohort1, "vs", cohort2)) %>%
    select(Comparison, Type, common, unique_cohort1, unique_cohort2, total_union, prop_common)
  
  return(results)
}




# investigate species behaviour in a super HGT network --------------------

# join data from both trials
super_dir_spp_hgts <- spp_dir_hgts_m4efad %>%
  filter(Directional == "Yes") %>% # must filter for directional events
  select(CLADE_B, CLADE_A) %>% # only need donor and receiver species of each HGT interaction
  rename(Donor = CLADE_B, # renaming as events are directional
         Recipient = CLADE_A) %>%
  mutate(Donor = str_remove(Donor, "s__"), # remove prefix from species names
         Donor = str_replace_all(Donor, "_", " ")) %>% # remove underscores from species names
  mutate(Recipient = str_remove(Recipient, "s__"), # remove prefix from species names
         Recipient = str_replace_all(Recipient, "_", " ")) # remove underscores from species names

# count number of in and out edges for each species
super_spp_edge_count <- get_spp_edge_count(super_dir_spp_hgts)

# create a directed graph from the donor species - receiver species edge list
set.seed(1234) # set seed for this section, also set in functions to get network node and edge coordinates
super_directed_network <- graph_from_data_frame(d = super_dir_spp_hgts[, c(1, 2)], directed = TRUE) # col1 = donor species, col2 = recipient species
E(super_directed_network)$weight <- count_multiple(super_directed_network) # weight network by number of HGT interactions

super_network_n_nodes <- gorder(super_directed_network) # number of nodes = 442
super_network_n_edges <- gsize(super_directed_network) # number of edges = 3872

# get the network node layout for plotting
super_node_coords <- network_node_coordinates(super_directed_network, all_hgt_taxa_distinct)

# get the network edge layout for plotting
super_edge_coords <- network_edge_coordinates(super_directed_network, super_node_coords)

# plot the super HGT network to visualise interactions between taxa
super_hgt_network <- ggplot() +
  geom_segment(data = super_edge_coords, 
               aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(type = "closed", length = unit(2.5, "mm")),
               colour = "grey") + 
  geom_point(data = super_node_coords, aes(x = x, y = y, fill = Phylum), size = 3, shape = 21) +
  scale_fill_manual(values = phylum_palette, name = "Phylum") + # phylum mapped to node outline
  theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg format 
# ggsave(file = "plots/chp3/clean_plots/super_hgt_network.svg", 
#       plot = super_hgt_network, units = "cm", width=14, height=10)
# pdf format 
# ggsave(file = "plots/chp3/clean_plots/super_hgt_network.pdf", 
#       plot = super_hgt_network, units = "cm", width=14, height=10)

# randomise 10,000 HGT networks based on the super network and extract edge lists (note: takes close to 45 mins)
randomised_network_edges <- randomise_hgt_networks(super_directed_network)

# calculate the donor-receiver ratio for all nodes in 10,000 random networks
randomised_dr_ratios <- randomised_dr_ratio(randomised_network_edges)
randomised_dr_ratios$dr_ratio_log <- log10(randomised_dr_ratios$dr_ratio) # log donor-receiver ratios

# determine thresholds for source/sink definitions based on 2SD either side of mean log10 donor-receiver ratio (normal distribution)
log_mean <- mean(randomised_dr_ratios$dr_ratio_log) # 9.544999e-06
log_sd <- sd(randomised_dr_ratios$dr_ratio_log) # 0.1985073
log_lower_threshold <- log_mean - 2 * log_sd # -0.397005
log_upper_threshold <- log_mean + 2 * log_sd # 0.3970241
lower_threshold <- 10^log_lower_threshold # 0.4008621
upper_threshold <- 10^log_upper_threshold # 2.494733

# plot donor-receiver ratios for all nodes in 10,000 random networks (log10 transformed)
role_shading <- data.frame(xmin = c(-Inf, log_lower_threshold, log_upper_threshold), # create shaded areas for each role
                           xmax = c(log_lower_threshold, log_upper_threshold, Inf),
                           fill = c("Sink", "Conduit", "Source"))

randomised_dr_ratios_hist_log <- ggplot(randomised_dr_ratios, aes(x = dr_ratio_log)) +
  geom_rect(data = role_shading,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
            alpha = 0.31, inherit.aes = FALSE) +
  geom_histogram(binwidth = 0.1, alpha = 0.85, fill = "#36454F", colour = "black") +
  geom_vline(xintercept = log_lower_threshold, color = "#EE8866", linetype = "dashed", size = 0.5) +
  geom_vline(xintercept = log_upper_threshold, color = "#44BB99", linetype = "dashed", size = 0.5) +
  xlab("Donor-receiver ratio (log base 10)") + ylab("Frequency") +
  scale_fill_manual(values = role_palette) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg format 
# ggsave(file = "plots/chp3/clean_plots/randomised_dr_ratios_log.svg", 
#       plot = randomised_dr_ratios_hist_log, units = "cm", width=12, height=10)
# pdf format
ggsave(file = "plots/chp3/clean_plots/randomised_dr_ratios_log.pdf", 
      plot = randomised_dr_ratios_hist_log, units = "cm", width=10, height=10)


# define species network roles using above thresholds
super_spp_roles <- define_spp_roles(super_spp_edge_count, lower_threshold, upper_threshold, trial_name = NA, cohort_name = NA) %>%
  select(-c(Trial, Cohort))

# quantify percentage of species in super HGT network with each network role
nrow(super_spp_roles %>% filter(Type == "Source"))/nrow(super_spp_roles)*100 # 41.70616% = sources
nrow(super_spp_roles %>% filter(Type == "Conduit"))/nrow(super_spp_roles)*100 # 35.07109% = conduits
nrow(super_spp_roles %>% filter(Type == "Sink"))/nrow(super_spp_roles)*100 # 23.22275% = sinks


# investigate species behaviour in cohort-specific HGT networks -----------

# subset HGT events for those that are directed and between species only
dir_spp_hgts_m4efad <- spp_dir_hgts_m4efad %>% filter(Directional == "Yes")
# dir_spp_hgts_healthy12 <- spp_dir_hgts_healthy12 %>% filter(Directional == "Yes")

# get HGT events for each sample, in each cohort
hgts_per_sample_mam12 <- get_hgts_per_sample(dir_spp_hgts_m4efad, MAM_12_samples) # switch out specific metadata
hgts_per_sample_mam24 <- get_hgts_per_sample(dir_spp_hgts_m4efad, MAM_24_samples)
hgts_per_sample_healthy12 <- get_hgts_per_sample(dir_spp_hgts_m4efad, Healthy_12_samples)
hgts_per_sample_healthy24 <- get_hgts_per_sample(dir_spp_hgts_m4efad, Healthy_24_samples)

# make individual HGT networks for each trial cohort
mam12_hgt_network <- make_hgt_network_plot(hgts_per_sample_mam12, all_hgt_taxa_distinct, "MAM_12")
mam24_hgt_network <- make_hgt_network_plot(hgts_per_sample_mam24, all_hgt_taxa_distinct, "MAM_24")
healthy12_hgt_network <- make_hgt_network_plot(hgts_per_sample_healthy12, all_hgt_taxa_distinct, "Healthy_12")
healthy24_hgt_network <- make_hgt_network_plot(hgts_per_sample_healthy24, all_hgt_taxa_distinct, "Healthy_24")


# join plots together
hgt_networks <- wrap_plots(mam12_hgt_network, mam24_hgt_network, healthy12_hgt_network, healthy24_hgt_network,
                           ncol = 2, nrow = 2)
# svg
# ggsave(file = "plots/chp3/clean_plots/hgt_networks.svg", 
#       plot = hgt_networks, units = "cm", width=20, height=12)
# pdf 
ggsave(file = "plots/chp3/clean_plots/hgt_networks.pdf", 
       plot = hgt_networks, units = "cm", width=10, height=10)


# Calculate network metrics 
# Apply the randomized metrics function to each dataset
mam12_res <- calculate_network_metrics_randomized(hgts_per_sample_mam12, n_iter = 1000)
mam15_res <- calculate_network_metrics_randomized(hgts_per_sample_mam15, n_iter = 1000)
mam24_res <- calculate_network_metrics_randomized(hgts_per_sample_mam24, n_iter = 1000)
healthy12_res <- calculate_network_metrics_randomized(hgts_per_sample_healthy12, n_iter = 1000)
healthy24_res <- calculate_network_metrics_randomized(hgts_per_sample_healthy24, n_iter = 1000)

# Combine the results into a comparison table
metrics_comparison <- data.frame(
  Dataset = c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24"),
  Observed_Modularity = c(mam12_res$observed_modularity,
                          mam24_res$observed_modularity,
                          healthy12_res$observed_modularity,
                          healthy24_res$observed_modularity),
  Observed_Clustering = c(mam12_res$observed_clustering,
                          mam24_res$observed_clustering,
                          healthy12_res$observed_clustering,
                          healthy24_res$observed_clustering),
  Avg_Random_Modularity = c(mam12_res$avg_random_modularity,
                            mam24_res$avg_random_modularity,
                            healthy12_res$avg_random_modularity,
                            healthy24_res$avg_random_modularity),
  SD_Random_Modularity = c(mam12_res$sd_random_modularity,
                           mam24_res$sd_random_modularity,
                           healthy12_res$sd_random_modularity,
                           healthy24_res$sd_random_modularity),
  Avg_Random_Clustering = c(mam12_res$avg_random_clustering,
                            mam24_res$avg_random_clustering,
                            healthy12_res$avg_random_clustering,
                            healthy24_res$avg_random_clustering),
  SD_Random_Clustering = c(mam12_res$sd_random_clustering,
                           mam24_res$sd_random_clustering,
                           healthy12_res$sd_random_clustering,
                           healthy24_res$sd_random_clustering)
)


# format table (optional)
metrics_table <- metrics_comparison %>%
  pivot_longer(cols = -Dataset, names_to = "Metric", values_to = "Value") %>%
  pivot_wider(names_from = Dataset, values_from = Value)

# calculate the p value for observed versus actual network metrics for each dataset 
# Assuming you have computed network metrics for each dataset:
# mam12_res, mam24_res, healthy12_res, healthy24_res
metrics_list <- list(
  MAM_12 = mam12_res,
  MAM_24 = mam24_res,
  Healthy_12 = healthy12_res,
  Healthy_24 = healthy24_res
)

# Run permutation tests with multiple test correction (default is BH)
permutation_results <- run_permutation_tests_for_all(metrics_list, alternative = "greater", correction_method = "BH")


# Shannon diversity - 12 v 24  --------------------------------------------

# This section looks at the alpha diveristy between MAM v Healthy  
load("data/processed/hgt/sp_only_counts_shannon_and_richness.Rdata")

shannon_12_24 <- sp_only_counts_shannon_richness %>%
  filter(Age_months %in% c(12, 24))

shannon_overall_all_summary <- shannon_12_24 %>%
  group_by(Age_months) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    .groups = "drop"
  )

shannon_all_summary <- shannon_12_24 %>%
  group_by(Age_months, Condition) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    .groups = "drop"
  )


# Modify the plot with updated aesthetics and color scheme
shannon_12_24_plot <- ggplot(shannon_12_24, aes(x = Age_months, y = Shannon, fill = interaction(Condition, Age_months))) +
  geom_quasirandom(
    dodge.width = 0.75, shape = 16, size = 1, alpha = 0.7,
    aes(color = interaction(Condition, Age_months))
  ) +
  geom_boxplot(
    outlier.colour = NA, alpha = 0.5, color = "black",
    linewidth = 0.3, position = position_dodge(width = 0.75)
  ) +
  facet_wrap(~Condition) +
  scale_fill_manual(values = c(
    "MAM.12" = "darkgoldenrod1", "MAM.24" = "darkkhaki",
    "Healthy.12" = "cyan4", "Healthy.24" = "chocolate"
  )) +
  scale_color_manual(values = c(
    "MAM.12" = "darkgoldenrod1", "MAM.24" = "darkkhaki",
    "Healthy.12" = "cyan4", "Healthy.24" = "chocolate"
  )) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("Timepoint") +
  ylab("Shannon's diversity index")

# pdf
ggsave(file = "plots/chp3/clean_plots/shannon_12_24_plot.pdf", 
       plot = shannon_12_24_plot, units = "cm", width=6, height=8)


# Perform Wilcoxon tests for both MAM and Healthy between 12 and 24 months
wilcox_shannon_overall_12_24 <- wilcox.test(Shannon ~ Age_months, data = shannon_12_24)
wilcox_shannon_mam <- wilcox.test(Shannon ~ Age_months, data = shannon_12_24 %>% filter(Condition == "MAM"))
wilcox_shannon_healthy <- wilcox.test(Shannon ~ Age_months, data = shannon_12_24 %>% filter(Condition == "Healthy"))





# Defining species roles in each cohort -----------------------------------

# count number of in and out edges for each species across all samples in cohort
spp_edge_count_mam12 <- get_spp_edge_count(hgts_per_sample_mam12)
spp_edge_count_mam24 <- get_spp_edge_count(hgts_per_sample_mam24)
spp_edge_count_healthy12 <- get_spp_edge_count(hgts_per_sample_healthy12)
spp_edge_count_healthy24 <- get_spp_edge_count(hgts_per_sample_healthy24)

# define species network roles across all samples in cohort
spp_roles_mam12 <- define_spp_roles(spp_edge_count_mam12, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "MAM_12")
spp_roles_mam24 <- define_spp_roles(spp_edge_count_mam24, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "MAM_24")
spp_roles_healthy12 <- define_spp_roles(spp_edge_count_healthy12, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "Healthy_12")
spp_roles_healthy24 <- define_spp_roles(spp_edge_count_healthy24, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "Healthy_24")

# join data for all trial cohorts
joined_spp_roles <- bind_rows(spp_roles_mam12, spp_roles_mam24, spp_roles_healthy12, spp_roles_healthy24) %>%
  mutate(Trial_cohort = paste0(Trial, " ", Cohort)) # concatenate trial and cohort names

# calculate the frequency of each species network role in each trial cohort
spp_role_freq <- joined_spp_roles %>%
  select(Type, Trial, Cohort) %>%
  group_by(Trial, Cohort, Type) %>%
  summarise(Role_count = n()) %>% # count number of species in each trial / cohort with each role
  ungroup() %>%
  group_by(Trial, Cohort) %>%
  mutate(Role_freq = Role_count/sum(Role_count)) %>% # get the frequency of each role in each trial / cohort
  ungroup()

# reorder factor orders for plotting (reversed due to coord_flip)
spp_role_freq$Type <- factor(spp_role_freq$Type, levels = c("Source", "Sink", "Conduit"))
spp_role_freq$Trial <- factor(spp_role_freq$Trial, levels = c("All", "M4EFaD"))
spp_role_freq$Cohort <- factor(spp_role_freq$Cohort, levels = c("MAM_12", "MAM_24", "Healthy_12", "Healthy_24", "Super transferome"))

# plot the frequency of each species network role in each trial cohort
spp_role_freq_plot <- spp_role_freq %>%
  ggplot(aes(x = Cohort, y = Role_freq, fill = Type)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(name = "Species HGT\nnetwork role", values = role_palette) +
  xlab(NULL) + ylab("Proportion of species") + 
  coord_flip() + # make horizontal stacked bar
  facet_grid(Trial~., scales = "free", space = "free") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg 
# ggsave(file = "plots/chp3/clean_plots/spp_role_freq.svg", 
#       plot = spp_role_freq_plot, units = "cm", width=12, height=10)

# pdf 
ggsave(file = "plots/chp3/clean_plots/spp_role_freq.pdf", 
    plot = spp_role_freq_plot, units = "cm", width=12, height=10)


# Run the proportion tests for MAM_12 vs Healthy_12
results_MAM12_Healthy12 <- run_prop_test(spp_role_freq, "MAM_12", "Healthy_12")
# Run the proportion tests for MAM_24 vs Healthy_24
results_MAM24_Healthy24 <- run_prop_test(spp_role_freq, "MAM_24", "Healthy_24")


# investigate overlap of HGT network roles across cohorts -----------------

# format upset plot data for each species network role
source_upset_data <- format_upset_data(joined_spp_roles, "Source")
conduit_upset_data <- format_upset_data(joined_spp_roles, "Conduit")
sink_upset_data <- format_upset_data(joined_spp_roles, "Sink")

# plot overlap in sources between cohorts
source_spp_overlap <- make_upset_plot(source_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                      "darkgoldenrod1", "cyan4", "Source species") 
#svg
# ggsave(file = "plots/chp3/clean_plots/source_spp_overlap.svg", 
#       plot = source_spp_overlap, units = "cm", width=15, height=10)

#pdf
ggsave(file = "plots/chp3/clean_plots/source_spp_overlap.pdf", 
     plot = source_spp_overlap, units = "cm", width=10, height=9)

# plot overlap in conduits between cohorts
conduit_spp_overlap <- make_upset_plot(conduit_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                       "darkgoldenrod1", "cyan4", "Conduit species")
# svg
# ggsave(file = "plots/chp3/clean_plots/conduit_spp_overlap.svg", 
#       plot = conduit_spp_overlap, units = "cm", width=15, height=10)
# pdf
ggsave(file = "plots/chp3/clean_plots/conduit_spp_overlap.pdf", 
     plot = conduit_spp_overlap, units = "cm", width=10, height=9)

# plot overlap in sinks between cohorts
sink_spp_overlap <- make_upset_plot(sink_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                    "darkgoldenrod1", "cyan4", "Sink species")
# svg 
# ggsave(file = "plots/chp3/clean_plots/sink_spp_overlap.svg", 
#       plot = sink_spp_overlap, units = "cm", width=15, height=10)

# pdf 
ggsave(file = "plots/chp3/clean_plots/sink_spp_overlap.pdf", 
     plot = sink_spp_overlap, units = "cm", width=10, height=9)

# find conserved species network roles across all trial cohorts
conserved_source_spp <- conserved_network_roles(source_upset_data) # 6 species
conserved_conduit_spp <- conserved_network_roles(conduit_upset_data) # 2 species
conserved_sink_spp <- conserved_network_roles(sink_upset_data) # 8 species


# investigate species that had all three roles  ---------------------------
## Species with 3 roles 

three_roles_species <- joined_spp_roles %>%
  select(Species, Type, Cohort) %>%
  group_by(Species) %>%
  mutate(n_distinct_roles = n_distinct(Type)) %>%
  filter(n_distinct_roles == 3) %>%
  arrange(Species) %>%
  ungroup()


# how many unique subjects had the Coprococcus species as a donor or reciever

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


# investigate species w/ different roles within each Cohort  ------------------------------

# Identify species that switch between Sink and Source within each Cohort
species_sink_source_within_cohort <- joined_spp_roles %>%
  group_by(Cohort, Species) %>%
  filter("Sink" %in% Type & "Source" %in% Type) %>% # Ensure species have both Sink and Source roles
  ungroup() %>%
  arrange(Cohort, Species) # Arrange by Cohort and Species

# no species switch from sink to source or vice versa within Cohorts, ie. MAM or healthy


## Q: How many species maintained their role within each Cohort and how many switched to Conduit?

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
role_switch_conduit_plot <- ggplot(role_proportions, aes(x = Group, y = Proportion, fill = Role_Change)) +
  geom_bar(
    stat     = "identity",
    position = "stack",
    linewidth = 0.3,
    alpha     = 0.7
  ) +
  scale_fill_manual(
    name   = "Role Status",
    values = c(
      "Changed Role"    = "steelblue",
      "Maintained Role" = "grey"
    ),
    labels = c(
      "Changed Role"    = "Conduit switch",
      "Maintained Role" = "Role unchanged"
    )
  ) +
  theme_bw() +
  theme(
    panel.border    = element_blank(),
    panel.grid      = element_blank(),
    axis.text       = element_text(),
    axis.line       = element_line(),
    axis.title      = element_text(size = 9),
    axis.title.y    = element_text(size = 9),
    legend.position = "right",
    legend.title    = element_text(size = 9),
    legend.text     = element_text(size = 8)
  ) +
  xlab("") +
  ylab("Proportion (%)")

# pdf 
ggsave(file = "plots/chp3/clean_plots/role_switch_conduit_plot.pdf", 
       plot = role_switch_conduit_plot, units = "cm", width=8, height=9)

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


## What proportions changed from sink to conduit, source to conduit and conduit to sink, conduit to source?
# 1) Pivot out the 12 mo and 24 mo roles for each species & group
transitions <- joined_spp_roles %>%
  mutate(Group = sub("_.*", "", Cohort)) %>% 
  filter(Cohort %in% c("MAM_12","MAM_24","Healthy_12","Healthy_24")) %>%
  select(Group, Species, Cohort, Type) %>%
  pivot_wider(
    names_from  = Cohort,
    values_from = Type
  ) %>%
  # Create a single transition string
  mutate(
    transition = case_when(
      !is.na(MAM_12) & !is.na(MAM_24) & Group == "MAM" ~ paste(MAM_12, "→", MAM_24),
      !is.na(Healthy_12) & !is.na(Healthy_24) & Group == "Healthy" ~ paste(Healthy_12, "→", Healthy_24),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(transition))

# 2) Count and compute proportions for each transition in each group
transition_props <- transitions %>%
  group_by(Group, transition) %>%
  summarise(
    n_species = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    total_changed = sum(n_species),             # total species that changed role
    prop_changed  = n_species / total_changed * 100
  ) %>%
  ungroup()

print(transition_props)




# Run proportion test  ----------------------------------------------------------

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



# investigate the identity of the species that changed roles --------------

# Compare MAM_12 vs Healthy_12
result_12 <- compare_species_roles2(joined_spp_roles, "MAM_12", "Healthy_12")
# Compare MAM_24 vs Healthy_24
result_24 <- compare_species_roles2(joined_spp_roles, "MAM_24", "Healthy_24")

# Combine results into one table if desired
final_comparison2 <- bind_rows(result_12, result_24)
print(final_comparison2)


prop_tests_by_role <- final_comparison2 %>%
  filter(Comparison %in% c("MAM_12 vs Healthy_12", "MAM_24 vs Healthy_24")) %>%
  group_by(Type) %>%
  summarise(
    x12 = common[Comparison == "MAM_12 vs Healthy_12"],
    n12 = total_union[Comparison == "MAM_12 vs Healthy_12"],
    x24 = common[Comparison == "MAM_24 vs Healthy_24"],
    n24 = total_union[Comparison == "MAM_24 vs Healthy_24"],
    p_value = prop.test(c(x12, x24), c(n12, n24))$p.value,
    .groups = "drop"
  ) %>%
  ungroup() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH")
  )

print(prop_tests_by_role)


# reshape for plotting
plot_df <- final_comparison2 %>%
  select(Type, Comparison, prop_common) %>%
  # split "MAM_12 vs Healthy_12" into Cohort1="MAM_12", Cohort2="Healthy_12"
  separate(Comparison, into = c("Cohort1", "Cohort2"), sep = " vs ") %>%
  # derive a simple Timepoint label from Cohort1
  mutate(Timepoint = ifelse(grepl("12", Cohort1), "12", "24"))

# Now you can plot as a connected‐dot / slope graph:
comparison_roles_12_24_plot <- ggplot(plot_df, aes(x = Timepoint, y = prop_common, group = Type, color = Type)) +
  geom_line(size = 0.7) +
  geom_point(size = 3) +
  facet_wrap(~ Type) +
  scale_color_manual(values = role_palette) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    x     = "Age (months)",
    y     = "Shared Species (%)",
    color = "Role"
  ) +
  theme_bw() +
  theme(
    # retain the grey facet strip background
    strip.text        = element_text(face = "bold"),
    # thin border around each facet to separate roles
    panel.border      = element_rect(color = "grey80", fill = NA, size = 0.3),
    panel.spacing     = unit(0.5, "lines"),
    # keep only horizontal grid lines
    panel.grid.major.y = element_line(color = "grey80"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    # draw x and y axis lines
    axis.line         = element_line(color = "black"),
    axis.title        = element_text(size = 9),
    axis.text         = element_text(size = 10),
    legend.position   = "none"
  ) +
  xlab("") +
  ylab("Shared Species (%)")

# pdf 
ggsave(file = "plots/chp3/clean_plots/comparison_roles_12_24_plot.pdf", 
       plot = comparison_roles_12_24_plot, units = "cm", width=8, height=9)



# save data files ---------------------------------------------------------

save(super_spp_roles, file = "data/processed/clean_outputs/super_spp_roles.RData")
save(log_lower_threshold, log_upper_threshold, lower_threshold, upper_threshold, file = "data/processed/clean_outputs/spp_role_thresholds.RData")
save(hgts_per_sample_mam12, hgts_per_sample_mam15, hgts_per_sample_mam24, hgts_per_sample_healthy12, hgts_per_sample_healthy24,
     file = "data/processed/clean_outputs/dir_hgts_per_sample.RData")
save(joined_spp_roles, file = "data/processed/clean_outputs/joined_spp_roles.RData")
save(randomised_network_edges, file = "data/processed/clean_outputs/randomised_network_edges.RData")




