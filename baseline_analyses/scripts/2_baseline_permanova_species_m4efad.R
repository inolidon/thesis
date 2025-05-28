## code summary 
## data can be directly loaded by using the files in the "save datasets" section 
# Load baseline species data and associated metadata.

# Calculate Bray-Curtis dissimilarity matrix from baseline species table.
# Reshape dissimilarity matrix for pairwise comparisons and NMDS ordination.

# Perform NMDS for visualization and plot Bray-Curtis distances colored by Condition.

# Run PERMANOVA:
# (1) Basic test for effect of Condition.
# (2) Multivariable model including Condition + demographic + clinical + technical covariates.
# (3) Individual PERMANOVAs to assess effect of each covariate separately on species composition.

# Adjust p-values for multiple testing (FDR correction) on individual PERMANOVA results.

# Save dissimilarity matrices and NMDS output.




# load packages  ----------------------------------------------------------

library(tidyverse) # required for data wrangling and ggplot
library(vegan) # required for ordination
library(ggbeeswarm) # required for geom_quasirandom() which adds jittered points to plot
library(Maaslin2) # required for differential pathway analysis
library(smplot2) # fits line to scatter plots and adds pearson correlation
library(scales) # required for label_comma()
library(data.table) # require for gene table cleaning
library(dplyr)
library(stringr)
library(readxl)
library(reshape2)
library(igraph)
library(pheatmap)
library(stats)
library(broom)




# load data  --------------------------------------------------------------

# baseline metadata 
load("data/processed/baseline_data.Rdata")
# baseline species - raw 
load("data/processed/baseline_species.Rdata")
# prevalence 10% filtered species 
load("data/processed/filtered_baseline_species.Rdata")
load("data/processed/all_meta_baseline.Rdata")
# contains all permanova values 
load("data/processed/baseline_permanova.Rdata") 




# beta diversity  ---------------------------------------------------------

# dataset needs to be in this format for the matrix to be created
baseline_species <- baseline_species %>% 
  column_to_rownames("Seq_ID")

# Calculate Bray-Curtis Dissimilarity Index (vegan package is used)
BC_dis_mam_healthy <- vegdist(baseline_species, method = "bray") 

# Data re-structuring 
BC_dis_long_mam_healthy <- BC_dis_mam_healthy %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column("Sample1") %>%
  gather(Sample2, BC_dis, -Sample1) %>%
  mutate(ID = apply(cbind(as.character(Sample1), as.character(Sample2)), 1, function(x) {str_c(sort(x), collapse = ":")})) %>%
  distinct(ID, .keep_all = T) %>% # removes duplicates
  filter(!(Sample1 == Sample2)) %>% # removes rows where samples are compared to themselves
  select(-ID)


# For visualisation of the Brays-Curtis dissimilarity index - NMDS
# Perform Ordination using Nonmetric Multidimensional Scaling (nMDS)
set.seed(280) # set a seed to ensure reproducible results (applicable to methods involving randomness)


# sp ~ condition ----------------------------------------------------------

# Running mds on the Bray-Curtis distance
mds_mam_healthy <- metaMDS(BC_dis_mam_healthy)

mds_data_mam_healthy <- mds_mam_healthy$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Seq_ID") %>% 
  inner_join(baseline_data)



# BC plot -----------------------------------------------------------------

BC_plot <- ggplot(mds_data_mam_healthy, aes(x = MDS1, y = MDS2)) + 
  geom_point(
    aes(fill = Condition, color = Condition),
    shape = 21,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  scale_color_manual(values = c("MAM" = "darkgoldenrod1", "Healthy" = "cyan4")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.text = element_text(),
    axis.line = element_line(),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")

# pdf
ggsave(file = "plots/chp1/clean_plots/BC_plot.pdf", 
       plot = BC_plot, units = "cm", width=9, height=9)

# run Permanova -----------------------------------------------------------

# Run PERMANOVA using the precomputed Bray-Curtis distance matrix
condition_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, Condition) 

permanova_result <- adonis2(
  BC_dis_mam_healthy ~ Condition,
  data = condition_data,
  permutations = 999,
  method = "bray"
)


# multi-permanova ---------------------------------------------------------

# Prepare metadata with selected covariates
metadata_multi <- all_meta_baseline %>%
  select(Seq_ID, Subject_ID, Condition, Sex,Place_of_birth, Breastfeeding_duration_months,
         Delivery_Mode, Weight, Length, MUAC, HC, WLZ_WHZ, Seq_batch)

# Ensure the row order matches the distance matrix
metadata_multi <- metadata_multi %>%
  filter(Seq_ID %in% rownames(baseline_species)) %>%
  arrange(match(Seq_ID, rownames(baseline_species)))

# Run PERMANOVA with the new covariates
permanova_multi <- adonis2(
  BC_dis_mam_healthy ~ Condition + Sex + Weight + WLZ_WHZ + Length + MUAC + HC + Place_of_birth + 
    Breastfeeding_duration_months + Delivery_Mode + Seq_batch,
  data = metadata_multi,
  permutations = 999,
  method = "bray"
)

# View result
print(permanova_multi)

# p value = 0.165
# R2 = 0.05313



# individual variable-permanova -------------------------------------------

# sp ~ WLZ_WHZ ratio  -----------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?
wlz_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, base_WLZ_WHZ) 

wlz_data2 <- wlz_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")



permanova_sp_WLZ <- adonis2(baseline_species ~ base_WLZ_WHZ, 
                            data = wlz_data2,
                            by = "margin", 
                            permutations = 999,
                            method = "bray")
# SIGNIFICANT DIFFERENCE - the WLZ_WHZ ratio significantly explains 0.9% of the variation
# in the species composition of the baseline dataset 

# p = 0.043
# r2 = 0.00751582


#  sp ~ weight   ----------------------------------------------- 

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?
weight_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, base_Weight) 

weight_data2 <- weight_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_weight <- adonis2(baseline_species ~ base_Weight, 
                               data = weight_data2,
                               by = "margin", 
                               permutations = 999,
                               method = "bray")

# p = 0.096	
# r = 0.006530109

#  sp ~ sex   -----------------------------------------------

sex_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, Sex) 

sex_data2 <- sex_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_sex <- adonis2(baseline_species ~ Sex, 
                            data = sex_data,
                            by = "margin", 
                            permutations = 999,
                            method = "bray")
# p = 0.741
# r = 0.003106412


#  sp ~ length   -----------------------------------------------

length_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, base_Length) 

length_data2 <- length_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_length <- adonis2(baseline_species ~ base_Length, 
                               data = length_data2,
                               by = "margin", 
                               permutations = 999,
                               method = "bray")


# p = 0.353
# r = 0.004699863

#  sp ~ muac   -----------------------------------------------

muac_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, base_MUAC) 

muac_data2 <- muac_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_muac <- adonis2(baseline_species ~ base_MUAC, 
                             data = muac_data2,
                             by = "margin", 
                             permutations = 999,
                             method = "bray")

# p = 0.189
# r = 0.005649538

#  sp ~ hc   -----------------------------------------------

hc_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, base_HC) 

hc_data2 <- hc_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_hc <- adonis2(baseline_species ~ base_HC, 
                           data = hc_data,
                           by = "margin", 
                           permutations = 999,
                           method = "bray")


# p = 0.212
# r = 0.005474936


# sp ~ Place of birth  -----------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?

pob_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID) 

pob_data2 <- pob_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")

permanova_sp_POB <- adonis2(baseline_species ~ Place_of_birth, 
                            data = pob_data2,
                            by = "margin", 
                            permutations = 999,
                            method = "bray")

# p = 0.838
# r2 = 0.002870808


# sp ~ Delivery mode -----------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?
dmode_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID) 

dmode_data2 <- dmode_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")


permanova_sp_DMode <- adonis2(baseline_species ~ Delivery_Mode, 
                              data = dmode_data2,
                              by = "margin", 
                              permutations = 999,
                              method = "bray")

# p = 0.146
# r2 = 0.005814564


# sp ~ Seq_batch -----------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?
seqb_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID, Seq_batch) 



permanova_sp_seq_batch <- adonis2(baseline_species ~ Seq_batch, 
                                  data = seqb_data,
                                  by = "margin", 
                                  permutations = 999,
                                  method = "bray")

# p = 0.422
# r2 = 0.00430593


# sp ~ Breast_feeding -----------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between MAM v Healthy group?
bf_data <- baseline_data %>% 
  select(Seq_ID, Subject_ID) 

bf_data2 <- bf_data %>% 
  left_join(all_meta_baseline, by = "Subject_ID")


permanova_sp_BFD <- adonis2(baseline_species ~ Breastfeeding_duration_months, 
                            data = bf_data2,
                            by = "margin", 
                            permutations = 999,
                            method = "bray")

# p = 0.418
# r2 = 0.004316728



# Adjust p values -----------------------------------------------

# The dataset baseline_permanova contains all the p values and r2 values for above permanovas 
# Assuming 'baseline_permanova' dataset is loaded with columns Variable, pvalue, and r2

# Adjust the p values 
baseline_permanova <- baseline_permanova %>%
  mutate(
    p_fdr = p.adjust(pvalue, method = "fdr")
  )




































# save data ---------------------------------------------------------------

# Save Bray Curtis values and long table and mds data
save(BC_dis_mam_healthy, file = "data/processed/BC_dis_mam_healthy.Rdata")
save(BC_dis_long_mam_healthy, file = "data/processed/BC_dis_long_mam_healthy.Rdata")
save(mds_data_mam_healthy, file = "data/processed/mds_data_mam_healthy.Rdata") 








































































