HGT Analyses – README

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


