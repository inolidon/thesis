HGT Analyses – README

Important Notice
These scripts are a work in progress and not the final versions. Updates and refinements are ongoing.

Script 1: 1_format_metadata_m4efad.R

Purpose:

This script is only required if you need to modify your own metadata to replicate the study.
For reproduction of the m4efad analyses, the preprocessed output files are already available in hgt_analyses/data.
Functionality:

Standardizes metadata column names to align with other studies (Gut Bugs and Ulcerative Colitis).
Generates sample subsets based on nutritional status (MAM, Healthy) at different timepoints (12, 15, and 24 months).
Output Files (saved in hgt_analyses/data):

metadata_m4efad.RData – Processed metadata file with standardized column names.
hgt_df_m4efad.RData – Reformatted HGT dataframe, aligning column names (e.g., Seq_ID → Participant_ID).
sample_subsets.RData – Subset of samples categorized by nutritional status and timepoints.
Script 2: 2_hgt_patterns_m4efad.R

Purpose:

Selects species-level HGT interactions.
Analyzes the proportion of total HGT events compared to species-level HGT events.
Visualizes HGT patterns with directional relationships.
Usage:

This script relies on preprocessed data files in hgt_analyses/data and can be executed directly.
Plot saving via ggsave() is currently commented out—modify as needed.
Output Files:

Processed data outputs are stored in hgt_analyses/data for use in the next script.
Script 3: 3_species_roles_m4efad.R

Purpose:

Assigns species roles in HGT networks.
Constructs and analyzes network structures.
Examines species role behavior within study cohorts.
Key Output:

joined_spp_roles.RData – A processed dataset containing species role classifications.
This file is required for subsequent analyses and is stored in hgt_analyses/data.
Reproducibility Notes

The dataset and preprocessed outputs have been provided for ease of replication.
Any modifications to metadata should be made by running 1_format_metadata_m4efad.R before proceeding with downstream scripts.
Ensure that all required dependencies and R packages are installed before execution.
For any questions regarding script usage or data processing, please refer to the associated documentation or contact the study authors.
