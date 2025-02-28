## code summary:
# create a super HGT network with all pooled samples from all trial cohorts (all groups of m4efad)
# calculate donor-receiver ratio for all nodes in 10k random networks and look at the distribution 
# select the donor-receiver ratio threshold (to define sources/sinks) according to the tails of the distribution
# define species as sources, conduits or sinks based on their role in HGT interactions across cohort samples (subset metadata) 
# calculate the proportion of species with each role across cohorts
# investigate the conservation of network roles across cohorts

## load libraries

library(tidyverse)
library(igraph)
library(patchwork)
library(ComplexUpset)
library(broom)


## set working directory

setwd("~/Library/CloudStorage/OneDrive-TheUniversityofAuckland/Desktop/IWD - PhD/r_projects/m4efad")

## load data

# HGT events between species for sample subset (individual donors, recipients with pre- and post-intervention), contains column with Directional == "Yes" / "No"
# on local pc 
#load("data/processed/clean_outputs/spp_dir_hgts_m4efad.RData")
# load("spp_dir_hgts_gb.RData")

# on github repo 
load("hgt_analyses/data/spp_dir_hgts_m4efad.RData")

# formatted metadata, note: sample IDs must match extracted WAAFLE sample IDs
# load("data/processed/clean_outputs/sample_subsets.RData")

# on github repo 
load("hgt_analyses/data/sample_subsets.RData")

# taxonomy lookup table
# on local pc
# all_hgt_taxa_distinct <- read.table("data/processed/clean_outputs/hgt_taxa_lookup.txt", sep = "\t")

# on github repo 
all_hgt_taxa_distinct <- read.table("hgt_analyses/data/hgt_taxa_lookup.txt", sep = "\t")

unique_phyla <- unique(all_hgt_taxa_distinct$Phylum)
na_phylum_rows <- all_hgt_taxa_distinct %>% filter(is.na(Phylum))

# Convert unique phyla into a dataframe
phyla_list_m4fad <- data.frame(Phylum = unique_phyla)


## load palettes

#phylum_palette <- c("Firmicutes" = "#009E73", "Bacteroidetes" = "#56B4E9", "Actinobacteria" = "#E69F00", "Proteobacteria" = "#CC79A7", 
#              "Verrucomicrobia" = "#F0E442", "Microsporidia" = "#9999CC", "Spirochaetes" = "#B65D40")
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

## load functions

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
make_upset_plot <- function(upset_data, trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4, trial_cohort_5,
                            trial_1_colour, trial_2_colour, species_role_name){
  # note: code assumes 2 trials with 5 cohorts each
  ## upset_data: output of format_upset_data() (df)
  ## trial_cohort_1 - trial_cohort_5: cohorts from trial 1 (str)
  ## trial_cohort_6 - trial_cohort_10: cohorts from trial 1 (str)
  ## trial_1_colour: hexcode for trial 1 (str)
  ## trial_2_colour hexcode for trial 2 (str)
  ComplexUpset::upset(upset_data, # call package as function conflicts with UpSetR package if that is loaded
                      rev(c(trial_cohort_1, trial_cohort_2, trial_cohort_3, trial_cohort_4, trial_cohort_5)),  # reverses order in plot
                      queries=list( # trial cohorts in same order as above vector
                        upset_query(set=trial_cohort_1, fill=trial_1_colour),
                        upset_query(set=trial_cohort_2, fill=trial_1_colour),
                        upset_query(set=trial_cohort_3, fill=trial_1_colour),
                        upset_query(set=trial_cohort_4, fill=trial_2_colour),
                        upset_query(set=trial_cohort_5, fill=trial_2_colour)),
                      height_ratio = 0.8, width_ratio = 0.2,
                      name = NULL, # remove interactions label
                      base_annotations=list(' '=(intersection_size(
                        counts=FALSE, # to remove numbers on top of bars
                        bar_number_threshold=1,  # show all numbers on top of bars
                        width=0.8) + # reduce width of the bars
                          scale_y_continuous(expand=expansion(mult=c(0, 0.05))) + # add some space on the top of the bars
                          theme(panel.grid.major=element_blank(), # hide grid lines
                                panel.grid.minor=element_blank(), 
                                axis.line=element_line(colour='black')))), # show axis lines
                      stripes='white', # disable stripe colouring
                      matrix=intersection_matrix(geom=geom_point(shape='circle filled', size=2, stroke=0.45)),
                      set_sizes = FALSE, # to remove set size plot on left
                      #set_sizes=(upset_set_size(geom=geom_bar(width=0.4)) +
                      #             theme(axis.title = NULL,
                      #                   axis.title.x = NULL,
                      #                   axis.line.x=element_line(colour='black'),
                      #                   axis.ticks.x=element_line()) +
                      #             ylab(species_role_name)),
                      sort_sets=FALSE, sort_intersections='descending')
}

# find conserved species network roles across all trial cohorts
conserved_network_roles <- function(upset_data){
  ## upset_data: output from format_upset_data() (df)
  upset_data %>%
    rowwise() %>%
    filter(sum(c_across(2:6)) == 5) %>% # find species with role in all 5 trial cohorts
    ungroup() %>%
    select(Species) %>%
    mutate(Species = str_remove(Species, "s__"), # remove prefix from species names
           Species = str_replace_all(Species, "_", " ")) # remove underscores from species names
}


# Running -----------------------------------------------------------------

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
# ggsave(file = "plots/chp3/clean_plots/randomised_dr_ratios_log.pdf", 
#       plot = randomised_dr_ratios_hist_log, units = "cm", width=12, height=10)


# define species network roles using above thresholds
super_spp_roles <- define_spp_roles(low_count_removed_super_spp_edge_count, lower_threshold, upper_threshold, trial_name = NA, cohort_name = NA) %>%
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
hgts_per_sample_mam15 <- get_hgts_per_sample(dir_spp_hgts_m4efad, MAM_15_samples)
hgts_per_sample_mam24 <- get_hgts_per_sample(dir_spp_hgts_m4efad, MAM_24_samples)
hgts_per_sample_healthy12 <- get_hgts_per_sample(dir_spp_hgts_m4efad, Healthy_12_samples)
hgts_per_sample_healthy24 <- get_hgts_per_sample(dir_spp_hgts_m4efad, Healthy_24_samples)

# make individual HGT networks for each trial cohort
mam12_hgt_network <- make_hgt_network_plot(hgts_per_sample_mam12, all_hgt_taxa_distinct, "MAM_12")
mam15_hgt_network <- make_hgt_network_plot(hgts_per_sample_mam15, all_hgt_taxa_distinct, "MAM_15")
mam24_hgt_network <- make_hgt_network_plot(hgts_per_sample_mam24, all_hgt_taxa_distinct, "MAM_24")
healthy12_hgt_network <- make_hgt_network_plot(hgts_per_sample_healthy12, all_hgt_taxa_distinct, "Healthy_12")
healthy24_hgt_network <- make_hgt_network_plot(hgts_per_sample_healthy24, all_hgt_taxa_distinct, "Healthy_24")


# join plots together
hgt_networks <- wrap_plots(mam12_hgt_network, mam15_hgt_network, mam24_hgt_network, healthy12_hgt_network, healthy24_hgt_network,
                           ncol = 3, nrow = 2)
# svg
# ggsave(file = "plots/chp3/clean_plots/hgt_networks.svg", 
#       plot = hgt_networks, units = "cm", width=20, height=12)
# pdf 
# ggsave(file = "plots/chp3/clean_plots/hgt_networks.pdf", 
#       plot = hgt_networks, units = "cm", width=20, height=12)


# Calculate network metrics 
# Apply the randomized metrics function to each dataset
mam12_res <- calculate_network_metrics_randomized(hgts_per_sample_mam12, n_iter = 1000)
mam15_res <- calculate_network_metrics_randomized(hgts_per_sample_mam15, n_iter = 1000)
mam24_res <- calculate_network_metrics_randomized(hgts_per_sample_mam24, n_iter = 1000)
healthy12_res <- calculate_network_metrics_randomized(hgts_per_sample_healthy12, n_iter = 1000)
healthy24_res <- calculate_network_metrics_randomized(hgts_per_sample_healthy24, n_iter = 1000)

# Combine the results into a comparison table
metrics_comparison <- data.frame(
  Dataset = c("MAM_12", "MAM_15", "MAM_24", "Healthy_12", "Healthy_24"),
  Observed_Modularity = c(mam12_res$observed_modularity,
                          mam15_res$observed_modularity,
                          mam24_res$observed_modularity,
                          healthy12_res$observed_modularity,
                          healthy24_res$observed_modularity),
  Observed_Clustering = c(mam12_res$observed_clustering,
                          mam15_res$observed_clustering,
                          mam24_res$observed_clustering,
                          healthy12_res$observed_clustering,
                          healthy24_res$observed_clustering),
  Avg_Random_Modularity = c(mam12_res$avg_random_modularity,
                            mam15_res$avg_random_modularity,
                            mam24_res$avg_random_modularity,
                            healthy12_res$avg_random_modularity,
                            healthy24_res$avg_random_modularity),
  SD_Random_Modularity = c(mam12_res$sd_random_modularity,
                           mam15_res$sd_random_modularity,
                           mam24_res$sd_random_modularity,
                           healthy12_res$sd_random_modularity,
                           healthy24_res$sd_random_modularity),
  Avg_Random_Clustering = c(mam12_res$avg_random_clustering,
                            mam15_res$avg_random_clustering,
                            mam24_res$avg_random_clustering,
                            healthy12_res$avg_random_clustering,
                            healthy24_res$avg_random_clustering),
  SD_Random_Clustering = c(mam12_res$sd_random_clustering,
                           mam15_res$sd_random_clustering,
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
# mam12_res, mam15_res, mam24_res, healthy12_res, healthy24_res
metrics_list <- list(
  MAM_12 = mam12_res,
  MAM_15 = mam15_res,
  MAM_24 = mam24_res,
  Healthy_12 = healthy12_res,
  Healthy_24 = healthy24_res
)

# Run permutation tests with multiple test correction (default is BH)
permutation_results <- run_permutation_tests_for_all(metrics_list, alternative = "greater", correction_method = "BH")

# count number of in and out edges for each species across all samples in cohort
spp_edge_count_mam12 <- get_spp_edge_count(hgts_per_sample_mam12)
spp_edge_count_mam15 <- get_spp_edge_count(hgts_per_sample_mam15)
spp_edge_count_mam24 <- get_spp_edge_count(hgts_per_sample_mam24)
spp_edge_count_healthy12 <- get_spp_edge_count(hgts_per_sample_healthy12)
spp_edge_count_healthy24 <- get_spp_edge_count(hgts_per_sample_healthy24)

# define species network roles across all samples in cohort
spp_roles_mam12 <- define_spp_roles(spp_edge_count_mam12, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "MAM_12")
spp_roles_mam15 <- define_spp_roles(spp_edge_count_mam15, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "MAM_15")
spp_roles_mam24 <- define_spp_roles(spp_edge_count_mam24, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "MAM_24")
spp_roles_healthy12 <- define_spp_roles(spp_edge_count_healthy12, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "Healthy_12")
spp_roles_healthy24 <- define_spp_roles(spp_edge_count_healthy24, lower_threshold, upper_threshold, trial_name = "M4EFaD", cohort_name = "Healthy_24")

# join data for all trial cohorts
joined_spp_roles <- bind_rows(spp_roles_mam12, spp_roles_mam15, spp_roles_mam24, spp_roles_healthy12, spp_roles_healthy24) %>%
  mutate(Trial_cohort = paste0(Trial, " ", Cohort)) # concatenate trial and cohort names

# calculate the frequency of each species network role in each trial cohort
spp_role_freq <- joined_spp_roles %>%
  select(Type, Trial, Cohort) %>%
  bind_rows(super_spp_roles %>% select(Type) %>% mutate(Trial = "All", Cohort = "Super transferome")) %>% # also include super transferome data
  group_by(Trial, Cohort, Type) %>%
  summarise(Role_count = n()) %>% # count number of species in each trial / cohort with each role
  ungroup() %>%
  group_by(Trial, Cohort) %>%
  mutate(Role_freq = Role_count/sum(Role_count)) %>% # get the frequency of each role in each trial / cohort
  ungroup()

# reorder factor orders for plotting (reversed due to coord_flip)
spp_role_freq$Type <- factor(spp_role_freq$Type, levels = c("Source", "Sink", "Conduit"))
spp_role_freq$Trial <- factor(spp_role_freq$Trial, levels = c("All", "M4EFaD"))
spp_role_freq$Cohort <- factor(spp_role_freq$Cohort, levels = c("MAM_12", "MAM_15", "MAM_24", "Healthy_12", "Healthy_24", "Super transferome"))

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
# ggsave(file = "plots/chp3/clean_plots/spp_role_freq.pdf", 
#       plot = spp_role_freq_plot, units = "cm", width=12, height=10)


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
source_spp_overlap <- make_upset_plot(source_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_15", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                      "darkgoldenrod1", "cyan4", "Source species") 
#svg
# ggsave(file = "plots/chp3/clean_plots/source_spp_overlap.svg", 
#       plot = source_spp_overlap, units = "cm", width=15, height=10)

#pdf
# ggsave(file = "plots/chp3/clean_plots/source_spp_overlap.pdf", 
#       plot = source_spp_overlap, units = "cm", width=15, height=10)

# plot overlap in conduits between cohorts
conduit_spp_overlap <- make_upset_plot(conduit_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_15", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                       "darkgoldenrod1", "cyan4", "Conduit species")
# svg
# ggsave(file = "plots/chp3/clean_plots/conduit_spp_overlap.svg", 
#       plot = conduit_spp_overlap, units = "cm", width=15, height=10)
# pdf
# ggsave(file = "plots/chp3/clean_plots/conduit_spp_overlap.pdf", 
#       plot = conduit_spp_overlap, units = "cm", width=15, height=10)

# plot overlap in sinks between cohorts
sink_spp_overlap <- make_upset_plot(sink_upset_data, "M4EFaD MAM_12", "M4EFaD MAM_15", "M4EFaD MAM_24", "M4EFaD Healthy_12", "M4EFaD Healthy_24",
                                    "darkgoldenrod1", "cyan4", "Sink species")
# svg 
# ggsave(file = "plots/chp3/clean_plots/sink_spp_overlap.svg", 
#       plot = sink_spp_overlap, units = "cm", width=15, height=10)

# pdf 
# ggsave(file = "plots/chp3/clean_plots/sink_spp_overlap.pdf", 
#       plot = sink_spp_overlap, units = "cm", width=15, height=10)

# find conserved species network roles across all trial cohorts
conserved_source_spp <- conserved_network_roles(source_upset_data) # 6 species
conserved_conduit_spp <- conserved_network_roles(conduit_upset_data) # 2 species
conserved_sink_spp <- conserved_network_roles(sink_upset_data) # 7 species

# find proportion of species with changed network roles across cohorts (not just missing from HGT in all cohort samples)
n_spp_roles <- joined_spp_roles %>% 
  select(Species, Type) %>%
  group_by(Species) %>%
  mutate(n_cohorts = if_else(n() == 1, "Only 1 cohort", "Across >1 cohort")) %>% # identify species with HGT interactions in at least 2 cohorts
  mutate(n_distinct_roles = n_distinct(Type)) %>% # count number of roles for each species across all cohorts
  mutate(Conduit_switch = if_else(n_distinct_roles == 2 & any(Type == "Conduit"), 
                                  "Conduit\nswitch", NA)) %>% # for species with 2 roles, was it a conduit switch?
  ungroup() %>%
  select(-Type) %>%
  distinct() %>%
  group_by(n_cohorts, n_distinct_roles, Conduit_switch) %>%
  summarise(Species_count = n()) %>% # summarise number of distinct species with each role count
  ungroup() %>%
  mutate(Species_prop = Species_count/sum(Species_count)) # proportion of distinct species with each role count

# reorder factor orders for plotting
n_spp_roles$n_cohorts <- factor(n_spp_roles$n_cohorts, levels = c("Only 1 cohort", "Across >1 cohort"))

# plot proportion of species with different network role counts across cohorts
n_spp_roles_plot <- n_spp_roles %>%
  ggplot(aes(x = n_distinct_roles, y = Species_prop, fill = Conduit_switch)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = "#99DDFF") +
  geom_text(aes(label = Conduit_switch), position = position_stack(vjust=0.5), na.rm = TRUE, colour="black") +
  facet_grid(~n_cohorts, scales = "free", space = "free") +
  scale_x_continuous(breaks = 1:3) + 
  xlab("Number of HGT network roles") + ylab("Proportion of species") + 
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), #formats background
        panel.grid.minor = element_blank(),
        legend.background = element_rect(color="black", size=0.4))
# svg
# ggsave(file = "plots/chp3/clean_plots/n_spp_roles.svg", 
#       plot = n_spp_roles_plot, units = "cm", width=12, height=12)

# pdf
# ggsave(file = "plots/chp3/clean_plots/n_spp_roles.pdf", 
#       plot = n_spp_roles_plot, units = "cm", width=12, height=12)


## Species with 3 roles 

three_roles_species <- joined_spp_roles %>%
  select(Species, Type) %>%
  group_by(Species) %>%
  mutate(n_distinct_roles = n_distinct(Type)) %>%
  filter(n_distinct_roles == 3) %>%
  ungroup()


three_roles_species <- joined_spp_roles %>%
  select(Species, Type, Cohort) %>%
  group_by(Species) %>%
  mutate(n_distinct_roles = n_distinct(Type)) %>%
  filter(n_distinct_roles == 3) %>%
  arrange(Species) %>%
  ungroup()


# save data files ---------------------------------------------------------

save(super_spp_roles, file = "data/processed/clean_outputs/super_spp_roles.RData")
save(log_lower_threshold, log_upper_threshold, lower_threshold, upper_threshold, file = "data/processed/clean_outputs/spp_role_thresholds.RData")
save(hgts_per_sample_mam12, hgts_per_sample_mam15, hgts_per_sample_mam24, hgts_per_sample_healthy12, hgts_per_sample_healthy24,
     file = "data/processed/clean_outputs/dir_hgts_per_sample.RData")
save(joined_spp_roles, file = "data/processed/clean_outputs/joined_spp_roles.RData")
save(randomised_network_edges, file = "data/processed/clean_outputs/randomised_network_edges.RData")



# TEST CODES - Do Not Run  ------------------------------------------------

## Test 1 - without Giant Component Filter 

library(igraph)

# Function to calculate modularity and clustering coefficient
calculate_network_metrics <- function(dir_spp_hgts) {
  # Preprocess data
  dir_spp_hgts <- dir_spp_hgts %>%
    select(Donor, Recipient) %>%
    mutate(Donor = str_remove(Donor, "s__"),
           Donor = str_replace_all(Donor, "_", " "),
           Recipient = str_remove(Recipient, "s__"),
           Recipient = str_replace_all(Recipient, "_", " "))
  
  # Create directed network
  network <- graph_from_data_frame(d = dir_spp_hgts, directed = TRUE)
  
  # Convert to undirected for modularity calculation
  undirected_network <- as.undirected(network, mode = "collapse")
  
  # Calculate modularity using community detection
  community <- cluster_louvain(undirected_network)
  modularity_value <- modularity(community)
  
  # Calculate clustering coefficient
  clustering_coeff <- transitivity(network, type = "global")
  
  return(list(modularity = modularity_value, clustering_coefficient = clustering_coeff))
}

# Apply function to each dataset
mam12_metrics <- calculate_network_metrics(hgts_per_sample_mam12)
mam15_metrics <- calculate_network_metrics(hgts_per_sample_mam15)
mam24_metrics <- calculate_network_metrics(hgts_per_sample_mam24)
healthy12_metrics <- calculate_network_metrics(hgts_per_sample_healthy12)
healthy24_metrics <- calculate_network_metrics(hgts_per_sample_healthy24)

# Print results
list(MAM_12 = mam12_metrics, 
     MAM_15 = mam15_metrics, 
     MAM_24 = mam24_metrics, 
     Healthy_12 = healthy12_metrics, 
     Healthy_24 = healthy24_metrics)


network_metrics <- list(
  MAM_12 = mam12_metrics, 
  MAM_15 = mam15_metrics, 
  MAM_24 = mam24_metrics, 
  Healthy_12 = healthy12_metrics, 
  Healthy_24 = healthy24_metrics
)

# Convert to a data frame
network_metrics_df <- do.call(rbind, lapply(network_metrics, as.data.frame))
network_metrics_df <- data.frame(Network = rownames(network_metrics_df), network_metrics_df, row.names = NULL)

# Print the dataset
network_metrics_df


# Compute the standard deviation for each metric
sd_modularity <- sd(network_metrics_df$modularity)
sd_clustering <- sd(network_metrics_df$clustering_coefficient)

# Combine the results into a new data frame
std_metrics_df <- data.frame(
  Metric = c("Modularity", "Clustering Coefficient"),
  Standard_Deviation = c(sd_modularity, sd_clustering)
)

# Print the standard deviation results
std_metrics_df




## Test 2 - Giant Component Filter Added 

library(igraph)
library(dplyr)
library(stringr)

# Function to calculate modularity and clustering coefficient with giant component filter
calculate_network_metrics <- function(dir_spp_hgts) {
  # Preprocess data
  dir_spp_hgts <- dir_spp_hgts %>%
    select(Donor, Recipient) %>%
    mutate(Donor = str_remove(Donor, "s__"),
           Donor = str_replace_all(Donor, "_", " "),
           Recipient = str_remove(Recipient, "s__"),
           Recipient = str_replace_all(Recipient, "_", " "))
  
  # Create directed network
  network <- graph_from_data_frame(d = dir_spp_hgts, directed = TRUE)
  
  # Convert to undirected for modularity calculation
  undirected_network <- as.undirected(network, mode = "collapse")
  
  # Filter to keep only the giant component
  giant_component <- induced_subgraph(undirected_network, which(clusters(undirected_network)$membership == which.max(table(clusters(undirected_network)$membership))))
  
  # Calculate modularity using community detection on the giant component
  community <- cluster_louvain(giant_component)
  modularity_value <- modularity(community)
  
  # Calculate clustering coefficient for the giant component
  clustering_coeff <- transitivity(giant_component, type = "global")
  
  return(list(modularity = modularity_value, clustering_coefficient = clustering_coeff))
}

# Apply function to each dataset
mam12_metrics <- calculate_network_metrics(hgts_per_sample_mam12)
mam15_metrics <- calculate_network_metrics(hgts_per_sample_mam15)
mam24_metrics <- calculate_network_metrics(hgts_per_sample_mam24)
healthy12_metrics <- calculate_network_metrics(hgts_per_sample_healthy12)
healthy24_metrics <- calculate_network_metrics(hgts_per_sample_healthy24)

# Print results
list(MAM_12 = mam12_metrics, 
     MAM_15 = mam15_metrics, 
     MAM_24 = mam24_metrics, 
     Healthy_12 = healthy12_metrics, 
     Healthy_24 = healthy24_metrics)


network_metrics <- list(
  MAM_12 = mam12_metrics, 
  MAM_15 = mam15_metrics, 
  MAM_24 = mam24_metrics, 
  Healthy_12 = healthy12_metrics, 
  Healthy_24 = healthy24_metrics
)

# Convert to a data frame
network_metrics_df <- do.call(rbind, lapply(network_metrics, as.data.frame))
network_metrics_df <- data.frame(Network = rownames(network_metrics_df), network_metrics_df, row.names = NULL)

# Print the dataset
network_metrics_df


# Compute the standard deviation for each metric
sd_modularity <- sd(network_metrics_df$modularity)
sd_clustering <- sd(network_metrics_df$clustering_coefficient)

# Combine the results into a new data frame
std_metrics_df <- data.frame(
  Metric = c("Modularity", "Clustering Coefficient"),
  Standard_Deviation = c(sd_modularity, sd_clustering)
)

# Print the standard deviation results
std_metrics_df


## Testing three role species 

library(ggplot2)
library(dplyr)

# Prepare data: Ensure species are ordered for better visualization
heatmap_data <- three_roles_species %>%
  select(Species, Cohort, Type) %>%
  arrange(Species, Cohort)

# Create heatmap
ggplot(heatmap_data, aes(x = Cohort, y = Species, fill = Type)) +
  geom_tile(color = "white") +  # White borders for better separation
  scale_fill_manual(values = c("Source" = "#44BB99",  "Conduit" = "#99DDFF", "Sink" = "#EE8866")) +  # Customize colors
  labs(title = "Species with 3 Roles Across Cohorts",
       x = "Cohort",
       y = "Species",
       fill = "Role") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



#### Pulling out the three role species counts 

three_roles_species <- joined_spp_roles %>%
  
  
  
  sp_three_roles <- joined_spp_roles %>%
  inner_join(three_roles_species %>% select(Species) %>% distinct(), by = "Species")

sp_three_roles <- sp_three_roles %>%
  ungroup() %>%
  arrange(Species)


# Instances of true donor or reciever counts = 0 
sp_three_roles_filtered <- sp_three_roles %>%
  filter(Donor_count == 0 | Recipient_count == 0)


# All times when "Conduit" is assigned 

sp_three_roles_filtered_conduit <- sp_three_roles %>%
  filter(Type == "Conduit")


## Checking min_max_score, avg_max_score for these species 

# how confident are we that the species with atleast one count is real vs random?
# we check using waafle's scoring method


# extact the orginal scores for each of the 7 species 

score_check_three_roles_sp <- hgt_sp_level %>%
  filter(CLADE_A %in% sp_three_roles$Species | CLADE_B %in% sp_three_roles$Species)

score_check_three_roles_sp <- score_check_three_roles_sp %>% 
  arrange(CLADE_A) %>% 
  filter(DIRECTION != "A?B")


## PLotting the distribution of contig length and scores 

library(ggplot2)

# Plot the distribution of CONTIG_LENGTH
plot_contig_length <- ggplot(score_check_three_roles_sp, aes(x = CONTIG_LENGTH)) +
  geom_histogram(fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of CONTIG_LENGTH",
       x = "CONTIG_LENGTH",
       y = "Frequency") +
  theme_minimal()

# Plot the distribution of MIN_MAX_SCORE
plot_min_max_score <- ggplot(score_check_three_roles_sp, aes(x = MIN_MAX_SCORE)) +
  geom_histogram(fill = "forestgreen", color = "black", alpha = 0.7) +
  labs(title = "Distribution of MIN_MAX_SCORE",
       x = "MIN_MAX_SCORE",
       y = "Frequency") +
  theme_minimal()

# Plot the distribution of AVG_MAX_SCORE
plot_avg_max_score <- ggplot(score_check_three_roles_sp, aes(x = AVG_MAX_SCORE)) +
  geom_histogram(fill = "firebrick", color = "black", alpha = 0.7) +
  labs(title = "Distribution of AVG_MAX_SCORE",
       x = "AVG_MAX_SCORE",
       y = "Frequency") +
  theme_minimal()

# To display the plots, simply print them:
plot_contig_length
plot_min_max_score
plot_avg_max_score

# Calculate the mean and standard deviation for each metric and store in a new dataset
metrics_summary <- score_check_three_roles_sp %>%
  summarise(
    mean_contig_length = mean(CONTIG_LENGTH, na.rm = TRUE),
    sd_contig_length   = sd(CONTIG_LENGTH, na.rm = TRUE),
    mean_min_max_score = mean(MIN_MAX_SCORE, na.rm = TRUE),
    sd_min_max_score   = sd(MIN_MAX_SCORE, na.rm = TRUE),
    mean_avg_max_score = mean(AVG_MAX_SCORE, na.rm = TRUE),
    sd_avg_max_score   = sd(AVG_MAX_SCORE, na.rm = TRUE)
  )

# Print the resulting dataset
metrics_summary


#### CHECKING LCA


library(dplyr)
library(tidyr)

lca_check <- score_check_three_roles_sp %>%
  # Pivot the CLADE_A and CLADE_B columns into a long format
  pivot_longer(cols = c(CLADE_A, CLADE_B), names_to = "clade", values_to = "Species_match") %>%
  # Filter for rows where the species match is one of the species in sp_three_roles
  filter(Species_match %in% sp_three_roles$Species) %>%
  # Select the species and the corresponding LCA, renaming Species_match to Species
  select(Species = Species_match, LCA) %>%
  # Keep only distinct combinations so that multiple occurrences with the same LCA are not repeated
  distinct()

# View the resulting dataset
lca_check



## Pulling out ancestral transfers 
filtered_k_bacteria <- score_check_three_roles_sp %>%
  filter(grepl("^k__Bacteria", LCA))



# Low HGT count removal  --------------------------------------------------

# use super_spp_edge_count dataset which contains a list of species with their D and R
# counts at species level transfers, n = 442 unique species for the whole m4efad trial 

library(dplyr)
library(ggplot2)

# Calculate summary statistics for Donor_count and Recipient_count using super_spp_edge_count
summary_stats <- super_spp_edge_count %>%
  summarise(
    mean_donor = mean(Donor_count, na.rm = TRUE),
    sd_donor   = sd(Donor_count, na.rm = TRUE),
    mean_recipient = mean(Recipient_count, na.rm = TRUE),
    sd_recipient   = sd(Recipient_count, na.rm = TRUE)
  )

summary_stats <- super_spp_edge_count %>%
  summarise(
    mean_donor = mean(Donor_count, na.rm = TRUE),
    sd_donor   = sd(Donor_count, na.rm = TRUE),
    median_donor = median(Donor_count, na.rm = TRUE),
    mean_recipient = mean(Recipient_count, na.rm = TRUE),
    sd_recipient   = sd(Recipient_count, na.rm = TRUE),
    median_recipient = median(Recipient_count, na.rm = TRUE)
  )



# Print summary statistics
print(summary_stats)

# Plot distribution for Donor_count
plot_donor <- ggplot(super_spp_edge_count, aes(x = Donor_count)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Donor_count", x = "Donor_count", y = "Frequency") +
  theme_minimal()

# Plot distribution for Recipient_count
plot_recipient <- ggplot(super_spp_edge_count, aes(x = Recipient_count)) +
  geom_histogram(binwidth = 1, fill = "firebrick", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Recipient_count", x = "Recipient_count", y = "Frequency") +
  theme_minimal()

# Display the plots
plot_donor
plot_recipient




ggsave(file = "plots/chp3//plot_donor.svg", 
       plot = plot_donor, units = "cm", width=10, height=10)
ggsave(file = "plots/chp3//plot_donor.pdf", 
       plot = plot_donor, units = "cm", width=10, height=10)


ggsave(file = "plots/chp3/plot_recipient.svg", 
       plot = plot_recipient, units = "cm", width=10, height=10)

ggsave(file = "plots/chp3/plot_recipient.pdf", 
       plot = plot_recipient, units = "cm", width=10, height=10)



# MEAN-based filtering  ---------------------------------------------------


## Using the mean count, we filter for species greater than the mean 

filtered_super_spp_edge_count <- super_spp_edge_count %>%
  filter(Donor_count > 8.76 | Recipient_count > 8.76)

# we plot the distribution of the new filtered species set that justifies above 

# Plot distribution for Donor_count in the filtered dataset
plot_donor_filtered <- ggplot(filtered_super_spp_edge_count, aes(x = Donor_count)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Donor_count (Filtered)", 
       x = "Donor_count", 
       y = "Frequency") +
  theme_minimal()

# Plot distribution for Recipient_count in the filtered dataset
plot_recipient_filtered <- ggplot(filtered_super_spp_edge_count, aes(x = Recipient_count)) +
  geom_histogram(binwidth = 1, fill = "firebrick", color = "black", alpha = 0.7) +
  labs(title = "Distribution of Recipient_count (Filtered)", 
       x = "Recipient_count", 
       y = "Frequency") +
  theme_minimal()

# Display the plots
plot_donor_filtered
plot_recipient_filtered

# we get the summary stats 
summary_stats_filtered <- filtered_super_spp_edge_count %>%
  summarise(
    mean_donor = mean(Donor_count, na.rm = TRUE),
    sd_donor = sd(Donor_count, na.rm = TRUE),
    min_donor = min(Donor_count, na.rm = TRUE),
    max_donor = max(Donor_count, na.rm = TRUE),
    range_donor = max_donor - min_donor,
    
    mean_recipient = mean(Recipient_count, na.rm = TRUE),
    sd_recipient = sd(Recipient_count, na.rm = TRUE),
    min_recipient = min(Recipient_count, na.rm = TRUE),
    max_recipient = max(Recipient_count, na.rm = TRUE),
    range_recipient = max_recipient - min_recipient
  )

summary_stats_filtered

# save the files 
ggsave(file = "plots/chp3/plot_donor_filtered.svg", 
       plot = plot_donor_filtered, units = "cm", width=10, height=10)
ggsave(file = "plots/chp3/plot_donor_filtered.pdf", 
       plot = plot_donor_filtered, units = "cm", width=10, height=10)


ggsave(file = "plots/chp3/plot_recipient_filtered.svg", 
       plot = plot_recipient_filtered, units = "cm", width=10, height=10)

ggsave(file = "plots/chp3/plot_recipient_filtered.pdf", 
       plot = plot_recipient_filtered, units = "cm", width=10, height=10)




# MEDIAN-based filtering --------------------------------------------------

# when filtering the dataset using mean, we lose a large proprtion of the data
# as the mean is relatively large and the variability is very high
# therefore, we use a more robust method that is sensitive to large variations 
# and is not affected by extreme values
# Median Absolute Deviation 


mad_stats <- super_spp_edge_count %>%
  summarise(
    mad_donor = mad(Donor_count, na.rm = TRUE),
    mad_recipient = mad(Recipient_count, na.rm = TRUE)
  )

mad_stats


# after deterning the MAD for each of the donor counts and recipient counts
# we will calculate the threshold 

# Threshold = Median + ( k × MAD )
# k = 1, this is to capture the maximum amount of information


thresholds <- super_spp_edge_count %>%
  summarise(
    median_donor = median(Donor_count, na.rm = TRUE),
    mad_donor = mad(Donor_count, na.rm = TRUE),
    threshold_donor = median_donor + 1 * mad_donor,
    
    median_recipient = median(Recipient_count, na.rm = TRUE),
    mad_recipient = mad(Recipient_count, na.rm = TRUE),
    threshold_recipient = median_recipient + 1 * mad_recipient
  )

thresholds


# donor count threshold = 4.9652
# recipient count threshold = 2.4826

## Filter the dataset to remove the low counts

# Extract threshold values
donor_thresh <- thresholds$threshold_donor
recipient_thresh <- thresholds$threshold_recipient

# Filter dataset to only include rows where both Donor_count and Recipient_count 
# are greater than or equal to their respective thresholds

low_count_removed_super_spp_edge_count <- super_spp_edge_count %>%
  filter(Donor_count >= donor_thresh | Recipient_count >= recipient_thresh)

# Plot the distribution of the new dataset 

# Convert the dataset into long format for plotting
long_counts <- low_count_removed_super_spp_edge_count %>%
  select(Donor_count, Recipient_count) %>%
  pivot_longer(cols = c(Donor_count, Recipient_count),
               names_to = "Count_Type",
               values_to = "Count")

# Plot the distributions, facetted by Count_Type
low_counts_removed_distribution <- ggplot(long_counts, aes(x = Count)) +
  geom_histogram(fill = "steelblue", color = "black", binwidth = 1, alpha = 0.7) +
  facet_wrap(~ Count_Type, scales = "free_x") +
  labs(title = "Distribution of Donor_count and Recipient_count",
       x = "Count",
       y = "Frequency") +
  theme_minimal()

# save the files 
ggsave(file = "plots/chp3/clean_plots/low_counts_removed_distribution.svg", 
       plot = low_counts_removed_distribution, units = "cm", width=10, height=10)
ggsave(file = "plots/chp3/clean_plots/low_counts_removed_distribution.pdf", 
       plot = low_counts_removed_distribution, units = "cm", width=10, height=10)


# Calculate summary statistics as a single row
metrics_summary <- low_count_removed_super_spp_edge_count %>%
  summarise(
    Mean_Donor = mean(Donor_count, na.rm = TRUE),
    SD_Donor = sd(Donor_count, na.rm = TRUE),
    Median_Donor = median(Donor_count, na.rm = TRUE),
    Range_Donor = max(Donor_count, na.rm = TRUE) - min(Donor_count, na.rm = TRUE),
    Mean_Recipient = mean(Recipient_count, na.rm = TRUE),
    SD_Recipient = sd(Recipient_count, na.rm = TRUE),
    Median_Recipient = median(Recipient_count, na.rm = TRUE),
    Range_Recipient = max(Recipient_count, na.rm = TRUE) - min(Recipient_count, na.rm = TRUE)
  )

# Transform the summary into a table with rows as metrics and columns for Donor and Recipient
summary_table <- metrics_summary %>%
  pivot_longer(everything(), names_to = "Metric_Group", values_to = "Value") %>%
  separate(Metric_Group, into = c("Metric", "Group"), sep = "_") %>%
  pivot_wider(names_from = Group, values_from = Value) %>%
  rename(Donor = Donor, Recipient = Recipient)


## save as an Rdata file 




















