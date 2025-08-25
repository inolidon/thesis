# code summary
# the infants that recovered have varying trajectories 
# we stratify them into Fast and Slow and analyse MB metrics 




# Load Packages  ----------------------------------------------------------
library(tidyverse)
library(ggbeeswarm)
library(rstatix)
library(DescTools)
library(lme4)
library(lmerTest)
library(circlize)
library(stringr)
library(dplyr)
library(readr)
library(purrr)
library(vegan)
library(lme4)
library(broom.mixed)
library(knitr)
library(car)
library(AMR)


# Load data ---------------------------------------------------------------

# anthro cleaned files - the process of making these file are in the load and clean data section below
load("data/processed/anthro_2Q_MAM.Rdata")
load("data/processed/anthro_4Q_H.Rdata")
# dataset used in recovery strat plots
load("data/processed/anthro_recov.Rdata")

# metadata
load("data/processed/recovery_data.Rdata")

# alpha div metrics for Recovery speed category
load("data/processed/alpha_div_R_category.Rdata")

# recovery species dataset
load("data/processed/recovery_species.Rdata")
# genus data for P/B ratio 
load("data/processed/genus_all_samples.Rdata")
# phylum for F/B ratio 
load("data/processed/phylum_all_samples.Rdata")

# Z SCORE  ----------------------------------------------------------------

# 1: Load & Clean data ------------------------------------------------------

## The data we will use:
# anthro_2024 - for all MAM data 
# anthro_H - for all Healthy data

### This dataset is the 31 January updated dataset for anthro
anthro_MAM_2024 <- read_excel("data/raw/anthro_MAM_2024.xlsx", 
                              col_types = c("text", "date", "text", 
                                            "text", "numeric", "numeric", "date", 
                                            "numeric", "numeric", "numeric", 
                                            "numeric", "numeric", "text"))


# updated anthro dataset - USE THIS 
anthro_2024 <- anthro_MAM_2024 %>% 
  select(SID_LCC, WLZ_WHZ, An_Time, SEX, Group, Weight, Length, MUAC, HC) %>% 
  rename(Subject_ID = SID_LCC) %>% 
  filter(WLZ_WHZ != 99.99) %>% # remove all the rows that don't have multiple measures
  arrange(Subject_ID, An_Time) %>%
  group_by(Subject_ID) %>%
  mutate(measure_no = row_number()) # this is to have a new row with a number for each measure

## remove any Subject_IDs tha only have one Measure 
anthro_2024 <- anthro_2024 %>%
  group_by(Subject_ID) %>%
  filter(n_distinct(measure_no) > 1)


## Healthy controls
anthro_Healthy_2024 <- read_excel("data/raw/anthro_Healthy_2024.xlsx", 
                                  col_types = c("text", "date", "text", 
                                                "skip", "numeric", "skip", "numeric", 
                                                "numeric", "numeric", "numeric", 
                                                "numeric", "skip"))

anthro_H <- anthro_Healthy_2024 %>% 
  select(SID_LCC, WLZ_WHZ, An_Time, SEX, Weight, Length, MUAC, HC) %>% 
  rename(Subject_ID = SID_LCC) %>% 
  filter(WLZ_WHZ != 99.99) %>% # remove all the rows that don't have multiple measures
  arrange(Subject_ID, An_Time) %>%
  group_by(Subject_ID) %>%
  mutate(measure_no = row_number())

# remove any rows that have only one measure
anthro_H <- anthro_H %>%
  group_by(Subject_ID) %>%
  filter(n_distinct(measure_no) > 1)


# 2: Prepare data ---------------------------------------------------------

#### ADD A NEW WEEK COLUMN

# The data in anthro_2024 and anthro_H needs a column called Weeks

## first we will rename the original datasets to relevant names 

## anthro_2024 --> anthro_2Q_MAM : 2 Q means the dataset contains measures upto the 2nd Quater (An_TIme = 3.02), MAM = MAM only
## anthro_H --> anthro_4Q_H : H = Healthy 

anthro_2Q_MAM <- anthro_2024

anthro_4Q_H <- anthro_H

# MAM dataset first

# convert the An_Time column to a character, as it is originally numeric
anthro_2Q_MAM$An_Time <- as.character(anthro_2Q_MAM$An_Time)

# we will split the An_Time into columns to calculate in weeks.
anthro_2Q_MAM <- anthro_2Q_MAM %>%
  separate(An_Time, into = c("set", "num_w"), sep = "\\.", remove = FALSE)

# for set = 1, the num_w values comes up as NA.
# we need to change this to 0 to represent the first measurement at baseline
# for An_Time = 2.10, the set = 2, and num_w = 1
# we need to convert 1 to 10 to represent 10 weeks 

anthro_2Q_MAM$num_w <- ifelse(is.na(anthro_2Q_MAM$num_w ), "0",  # Replace NA with "0"
                              ifelse(anthro_2Q_MAM$num_w  == "1", "10",  # Replace "1" with "10"
                                     anthro_2Q_MAM$num_w ))  # Otherwise, keep the original value


# When set = 3, the measuremtn is not weekly, it is 3 monthly.
# so, when 3.01, 3 means 3 monthly measure and 01 means first
# therefore, 3.01 is 3 months or 12 weeks after the 2.13 or week 13
# this means 13 + 12weeks = 25 weeks or week 25
# and, for 3.02, thats another 3 months added to the previous measure
# so, 25 + 12 weeks = 37 weeks or week 37

anthro_2Q_MAM <- anthro_2Q_MAM %>%
  mutate(num_w = ifelse(set == "3" & num_w == "01", "25", #mutate enables modification of the column
                        ifelse(set == "3" & num_w == "02", "37", num_w)))


# we will convert the num_w column to numeric 
anthro_2Q_MAM$num_w <- as.numeric(anthro_2Q_MAM$num_w)


# HEALTHY dataset second

# convert the An_Time column to a character, as it is originally numeric
anthro_4Q_H$An_Time <- as.character(anthro_4Q_H$An_Time)

# we will split the An_Time into columns to calculate in weeks.
anthro_4Q_H <- anthro_4Q_H %>%
  separate(An_Time, into = c("set", "num_w"), sep = "\\.", remove = FALSE)

# the healthy dataset has different timing for anthro measures
# 1.00 is enrolment and 2.00 is one month later first measure, so its week 4
# 2.01, is the second measure another MONTH later for 4 + 4 weeks = 8, week 8 (monthly measures only for 3 months so 2.02)
# 3.01 is the first measure in 3 months, so last measure 2.02 is 16 weeks, 3.01 will be 16 + 12 weeks (3 mo) = 28 weeks, week 28

anthro_4Q_H <- anthro_4Q_H %>%
  mutate(num_w = ifelse(set == "1" & is.na(num_w), "0",
                        ifelse(set == "2" & is.na(num_w), "4",
                               ifelse(set == "2" & num_w == "01", "8",
                                      ifelse(set == "2" & num_w == "02", "12",
                                             ifelse(set == "2" & num_w == "03", "16",
                                                    ifelse(set == "3" & num_w == "01", "28",
                                                           ifelse(set == "3" & num_w == "02", "40",
                                                                  ifelse(set == "3" & num_w == "03", "52", 
                                                                         ifelse(set == "3" & num_w == "04", "64", num_w))))))))))


# we will convert the num_w column to numeric 
anthro_4Q_H$num_w <- as.numeric(anthro_4Q_H$num_w)

# save the files 
save(anthro_2Q_MAM, file = "data/processed/anthro_2Q_MAM.Rdata")
save(anthro_4Q_H, file = "data/processed/anthro_4Q_H.Rdata")

# 3: Q - Trajectory Plot  -----------------------------------------------------

# MAM
# plot with trajectory for each Subject_ID
ggplot(anthro_2Q_MAM, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_MAM_Zscores.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot

ggplot(anthro_2Q_MAM, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_MAM_Zscores.pdf", h = 20, w = 20) 


# HEALTHY
# plot with trajectory for each Subject_ID
ggplot(anthro_4Q_H, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 16, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 16
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_4Q_H_Zscores.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot

ggplot(anthro_4Q_H, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 16, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 16
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_4Q_H_Zscores.pdf", h = 20, w = 20) 


# 4: Q MAM - Recovery  ----------------------------------------------------------

# Next, combine the data with recovery state data to identify the infants that 
# relapsed post being identified as "Recovered" or "Unrecovered"


# select the columns of interest and save in a new dataset 
recov_state <- recovery_data %>% 
  select(Subject_ID, Recovery, Group.x, Sample_Name, Seq_ID)

recov_state < recov_state %>% 
  rename(Group = Group.x)


# combine recov to anthro_2Q_MAM (only MAM infants have the stratification of R v nR)
anthro_recov <- anthro_2Q_MAM %>% 
  inner_join(recov_state, by = "Subject_ID")

save(anthro_recov, file ="data/processed/anthro_recov.Rdata")

# filter for R only 
R_anthro_rec <- anthro_recov %>% 
  filter(Recovery == "TRUE")

# filter for nR only 
nR_anthro_rec <- anthro_recov %>% 
  filter(Recovery == "FALSE")


# plot --------------------------------------------------------------------

# R only 

# Point plot with trajectory for each Subject_ID
ggplot(R_anthro_rec, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Recovered - Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_R_Zscores.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(R_anthro_rec, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Recovered - Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_R_Zscores.pdf", h = 20, w = 20) 


# nR

# Point plot with trajectory for each Subject_ID
ggplot(nR_anthro_rec, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Not recovered - Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_nR_Zscores.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(nR_anthro_rec, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Not recovered - Trajectory of Z scores of WLZ/WHZ ratio",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_nR_Zscores.pdf", h = 20, w = 20) 



# 5: WK13 -----------------------------------------------------------------


# Trajectory  -------------------------------------------------------------

## WK 13 ALL 

# Point plot with trajectory for each Subject_ID
ggplot(wk13_MAM, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio upto Week 13",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_nR_Zscores.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(wk13_MAM, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio upto Week 13",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_nR_Zscores.pdf", h = 20, w = 20) 


## WK 13 - R

# Line plot for the same data - neater plot
ggplot(wk13_R, aes(x = num_w, y = WLZ_WHZ, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Trajectory of Z scores of WLZ/WHZ ratio upto Week 13",
    x = "Number of Weeks",
    y = "WLZ_WHZ"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_wk13_R_Zscores.pdf", h = 20, w = 20) 



# WEIGHT ------------------------------------------------------------------

### All MAM (Weight - y axis is log transformed)

# Point plot with trajectory for each Subject_ID
ggplot(anthro_2Q_MAM, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "All MAM infants - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  theme_minimal() +
  scale_y_log10() +  # Log-transform y-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_MAM_weight.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(anthro_2Q_MAM, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "All MAM infants - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  scale_y_log10() +  # Log-transform x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_MAM_weight.pdf", h = 20, w = 20) 


#### R only 

# Point plot with trajectory for each Subject_ID
ggplot(R_anthro_rec, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Recovered ONLY - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  theme_minimal() +
  scale_y_log10() +  # Log-transform y-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_R_weight.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(R_anthro_rec, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Recovered ONLY - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  scale_y_log10() +  # Log-transform x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_R_weight.pdf", h = 20, w = 20) 


#### nR

# Point plot with trajectory for each Subject_ID
ggplot(nR_anthro_rec, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_point(position = position_dodge(width = 0.8)) +
  geom_line(position = position_dodge(width = 0.8)) +
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Not recovered ONLY - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  theme_minimal() +
  scale_y_log10() +  # Log-transform y-axis
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")  # Facet by Subject_ID with free y-axis scales

ggsave("plots/v3/draft/point_2Q_nR_weight.pdf", h = 20, w = 20) 

# Line plot for the same data - neater plot
ggplot(nR_anthro_rec, aes(x = num_w, y = Weight, group = Subject_ID)) +
  geom_line() +   # Keep only geom_line for a line chart
  #geom_hline(yintercept = -1, linetype = "dashed", color = "red") +  # Add intercept line
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +  # Add vertical line at num_w = 13
  labs(
    title = "Not recovered ONLY - Trajectory of Weight (kg)",
    x = "Number of Weeks",
    y = "Weight (kg) (log10 scale)"
  ) +
  scale_y_log10() +  # Log-transform x-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +  # Rotate x-axis labels if needed
  facet_wrap(~Subject_ID, scales = "free_y")

ggsave("plots/v3/draft/line_2Q_nR_weight.pdf", h = 20, w = 20) 



# *Histograms ----------------------------------------------------

# all MAM data
# toggle between WLZ_WHZ, weight and length
ggplot(anthro_2Q_MAM, aes(x = Length)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +  # Adjust binwidth and colors as needed
  labs(
    title = "Distribution of Length (cm) - all MAM",
    x = "Length (cm)",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("plots/v3/draft/anthro_2Q_MAM_Length_histogram.pdf", h = 6, w = 6) 


# Healthy
# toggle between WLZ_WHZ, weight and length
ggplot(anthro_4Q_H, aes(x = Length)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +  # Adjust binwidth and colors as needed
  labs(
    title = "Distribution of Length (cm) - Healthy",
    x = "Length (cm)",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("plots/v3/draft/anthro_4Q_H_Length_histogram.pdf", h = 6, w = 6) 

# R only
# toggle between WLZ_WHZ, weight and length
ggplot(R_anthro_rec, aes(x = Weight)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +  # Adjust binwidth and colors as needed
  labs(
    title = "Distribution of Weight (kg) - Recovered",
    x = "Weight",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("plots/v3/draft/anthro_2Q_R_Weight_histogram.pdf", h = 6, w = 6) 


# nR only
# toggle between WLZ_WHZ, weight and length
ggplot(nR_anthro_rec, aes(x = WLZ_WHZ)) +
  geom_histogram(binwidth = 0.5, fill = "skyblue", color = "black") +  # Adjust binwidth and colors as needed
  labs(
    title = "Distribution of WLZ_WHZ ratio - Not recovered",
    x = "WLZ_WHZ",
    y = "Frequency"
  ) +
  theme_minimal()
ggsave("plots/v3/draft/anthro_2Q_nR_WLZ_WHZ_histogram.pdf", h = 6, w = 6) 


# RECOVERY RATE (wks) - 13wk ----------------------------------------

# We are going to stratify the data into Fast and Slow recovery
# on the standard or enhanced refeed


# 1: Data prep ------------------------------------------------------------

# First, we msut remove all measurements after the 13wk 
# new dataset that contains all the MAM measures upto week 13 only (An_Time = 13)

wk13_MAM <- anthro_recov %>%
  filter(An_Time <= 2.13)

# next, we will filter the R only 
wk13_R <- wk13_MAM %>% 
  filter(Recovery == "TRUE")

# nR dataset 
wk13_nR <- wk13_MAM %>% 
  filter(Recovery == "FALSE")


# 2: Rate calculation -----------------------------------------------------

R_rate_wks <- wk13_R %>%
  group_by(Subject_ID) %>%
  summarise(Weeks = sum(!is.na(num_w))) %>% 
  mutate(category = ifelse(Weeks > 8, "Fast", "Slow")) %>%
  arrange(desc(Weeks))

# what is the average number of weeks to recovery in the recovered data?
# 8.17 wks +/- 2.88
R_rate_wks %>% summarize(mean_num_weeks = mean(Weeks, na.rm = TRUE))
R_rate_wks %>% summarize(sd_num_weeks = sd(Weeks, na.rm = TRUE))


# 3: Plot  ----------------------------------------------------------------

# Plot growth rates for comparison
plot_growth_rate <- ggplot(R_rate_wks, aes(x = reorder(Subject_ID, Weeks), y = Weeks, fill = category)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 8, linetype = "dashed", color = "red") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_fill_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  labs(
    x = "Subject ID",
    y = "Weeks to Recovery"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    plot.title = element_blank(),
    legend.position = "none"
  )

# pdf
ggsave(file = "plots/chp2/clean_plots/plot_growth_rate.pdf", 
       plot = plot_growth_rate, units = "cm", width=15, height=8)

# Plot histogram of growth rates
plot_growth_hist <- ggplot(R_rate_wks, aes(x = Weeks)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black", linewidth = 0.3) +
  labs(
    x = "Weeks to Recovery",
    y = "Frequency"
  ) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text.x = element_text(size = 9),
    axis.text.y = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    plot.title = element_blank()
  )

# pdf
ggsave(file = "plots/chp2/clean_plots/plot_growth_hist.pdf", 
       plot = plot_growth_hist, units = "cm", width=6, height=8)


# MB METRICS --------------------------------------------------------------

# prep data ---------------------------------------------------------------

# R_rate_wks dataset does not have Seq_ID column, so we must join it with a dataset that has recovery data
# recovery_data has all the info

R_rate_wks <- R_rate_wks %>% 
  inner_join(recovery_data, by = "Subject_ID")

R_rate_wks %>% 
  count(category)

# 1: Alpha div ---------------------------------------------------------------

# Calculate alpha diversity metrics (Species richness, Evenness, and Shannon's diversity index)

recovery_species <-recovery_species %>% 
  column_to_rownames ("Seq_ID")

alpha_div_R_category <- tibble(Seq_ID = rownames(recovery_species),
                               Richness = specnumber(recovery_species), # how many species 
                               Evenness = diversity(recovery_species, index = "shannon")/log(Richness),
                               Shannon = diversity(recovery_species, index = "shannon")) %>% 
  inner_join(R_rate_wks)

save(alpha_div_R_category, file = "data/processed/alpha_div_R_category.Rdata")


# Reshape to long format
alpha_div_long_category <- alpha_div_R_category %>%
  select(Seq_ID, category, Shannon, Richness, Evenness) %>%
  pivot_longer(
    cols = c(Shannon, Richness, Evenness),
    names_to = "Alpha_metric",
    values_to = "Value"
  )

# Set order for consistent facet layout
alpha_div_long_category$Alpha_metric <- factor(
  alpha_div_long_category$Alpha_metric,
  levels = c("Shannon", "Richness", "Evenness")
)

# Plot
alpha_div_plot_category <- ggplot(alpha_div_long_category, aes(x = category, y = Value, fill = category, color = category)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 21,
    size = 1,
    alpha = 0.7
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(0.75)
  ) +
  facet_wrap(~ Alpha_metric, scales = "free_y") +
  scale_fill_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  scale_color_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("Alpha Diversity Metric")

ggsave(file = "plots/chp2/clean_plots/alpha_div_plot_category.pdf",
       plot = alpha_div_plot_category, units = "cm", width = 8, height = 6)


#  Stats  -----------------------------------------------------------------

# Perform the Wilcoxon rank-sum test 
wilcox_alpha_div_R_category <- wilcox.test(Shannon ~ category, data = alpha_div_R_category)
# p-value = 0.8932

## Spp Richness 

# Perform the Wilcoxon rank-sum test 
wilcox_alpha_div_R_category_richness <- wilcox.test(Richness ~ category, data = alpha_div_R_category)
# p = 0.3087

## Spp Evenness

# Perform the Wilcoxon rank-sum test 
wilcox_alpha_div_R_category_evenness <- wilcox.test(Evenness ~ category, data = alpha_div_R_category)
# p = 0.4385


# generate descriptive stats for alpha diversity 
descript_stats_Rwks <- alpha_div_R_category %>% 
  group_by(category) %>% 
  summarise(n = n(),
            mean = mean(Shannon),
            sd = sd(Shannon),
            median = median(Shannon),
            min = min(Shannon),
            max = max(Shannon),
            range = max-min)

# summary stats
alpha_div_summary_R_category <- alpha_div_R_category %>%
  group_by(category) %>%
  summarise(
    mean_shannon = mean(Shannon, na.rm = TRUE),
    sd_shannon = sd(Shannon, na.rm = TRUE),
    mean_richness = mean(Richness, na.rm = TRUE),
    sd_richness = sd(Richness, na.rm = TRUE),
    mean_evenness = mean(Evenness, na.rm = TRUE),
    sd_evenness = sd(Evenness, na.rm = TRUE),
    .groups = "drop"
  )


# 2: Beta div -------------------------------------------------------------

# Calculate Bray-Curtis Dissimilarity Index (vegan package is used)
BC_dis_recovery <- vegdist(recovery_species, method = "bray") 

# Data re-structuring 
BC_dis_long_recovery <- BC_dis_recovery %>% 
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

# Running mds on the Bray-Curtis distance
mds_recovery <- metaMDS(BC_dis_recovery)

mds_data_R_category <- mds_recovery$points %>% 
  as.data.frame() %>% 
  rownames_to_column("Seq_ID") %>% 
  inner_join(R_rate_wks)

# Save Bray Curtis values and long table and mds data
save(BC_dis_recovery, file = "data/BC_dis_recovery.Rdata")
save(BC_dis_long_recovery, file = "data/BC_dis_long_recovery.Rdata")
save(mds_data_recovery, file = "data/mds_data_recovery.Rdata") 


# PLOT Brays Curtis
plot_mds_R_category <- ggplot(mds_data_R_category, aes(x = MDS1, y = MDS2)) + 
  geom_point(
    aes(fill = category, color = category),
    shape = 21,
    size = 2,
    alpha = 0.7
  ) +
  scale_fill_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  scale_color_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none",
    aspect.ratio = 1
  ) +
  xlab("NMDS1") +
  ylab("NMDS2")

ggsave(file = "plots/chp2/clean_plots/plot_mds_R_category.pdf", 
       plot = plot_mds_R_category, units = "cm", width=6, height=6)

# Permanova ---------------------------------------------------------------

#species against condition (y ~ x), Is there a significant difference in the species between Recovered v Not-recovered group?
load("data/recovery_species.Rdata")
# convert R_rate_wks column Seq_ID into a row


R_rate_wks <- R_rate_wks %>% 
  column_to_rownames("Seq_ID")

R_species <- filter(species, rownames(species) %in% rownames(R_rate_wks))

permanova_rfastvslow <- adonis2(R_species ~ category, 
                                data = R_rate_wks,
                                by = "margin", 
                                permutations = 999,
                                method = "bray")

# p = 0.891
# R2 = 0.008871782

# 3: P/B ratio ------------------------------------------------------------

# Make sure that the R_rate_wks dataset has a column called "Seq_ID" - in the right format
R_rate_wks <- R_rate_wks %>% 
  rownames_to_column("Seq_ID")

# Prepare the P/B ratio dataset for R category (Fast vs Slow)
PBratio_R_category <- genus_all_samples %>%
  select(g__Prevotella, g__Bacteroides) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    g__Prevotella = g__Prevotella + 1e-6,
    g__Bacteroides = g__Bacteroides + 1e-6,
    PBratio = g__Prevotella / g__Bacteroides
  ) %>%
  inner_join(R_rate_wks, by = "Seq_ID")

# Plot
PB_plot_R_category <- ggplot(PBratio_R_category, aes(x = category, y = PBratio, fill = category, color = category)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 21,
    size = 1,
    alpha = 0.7
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(
    limits = c(1e-6, 200),
    breaks = c(1e-6, 1e-4, 1e-2, 1, 100)
  ) +
  scale_fill_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  scale_color_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("P/B Ratio (log10-transformed)")

# Save
ggsave("plots/chp2/clean_plots/PB_plot_R_category.pdf", 
       plot = PB_plot_R_category, units = "cm", width = 4, height = 6)


# Statistic check
wilcox_result_R_category_p_b <- wilcox.test(PBratio_R_category$PBratio ~ PBratio_R_category$category)
# p-value = 0.4601


# 4: Top Phyla ------------------------------------------------------------
# top phyla tbale
top_phyla_R_category <- phylum %>%
  rownames_to_column("Seq_ID") %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(Gr1_2) %>%
  group_by(Taxa, Group) %>%
  summarise(MeanRA = mean(RA) * 100,
            SD = sd(RA) * 100) %>%
  arrange(desc(MeanRA)) %>%
  group_by(Group) %>%
  slice(1:5)

# structure the data for the plot
top_phyla_R_category <- filtered_recovery_phyla %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(R_rate_wks) 

# Violin plot with mean and log-scaled y-axis
ggplot(top_phyla_R_category, aes(x = Taxa, y = RA, fill = category)) +
  geom_violin(alpha = 0.5) +
  geom_quasirandom(dodge.width = 0.9, shape = 21, size = 2) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), alpha = 0.5, outlier.colour = NA) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 3,
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c("skyblue2", "goldenrod")) +
  labs(
    title = "Comparison of Phylum Relative Abundance between Recovery rates",
    x = "Phyla",
    y = "Relative Abundance (log10 scale)"
  ) +
  scale_y_log10() +  # Add log scale to the y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("plots/v3/draft/R_category_phylaRA.pdf", h = 6, w = 9)

# phylum data table to run stat
Gr1_Gr2_phylum_data <- phylum %>%
  rownames_to_column("Seq_ID") %>%
  gather(Taxa, RA, -Seq_ID) %>%
  inner_join(Gr1_2) %>%
  mutate(RA = RA * 100,
         Taxa_plot = ifelse(RA > 1, Taxa, "Other (< 1%)"))

top_phyla_R_category %>%
  #filter(Taxa %in% c("p__Bacteroidetes", "p__Firmicutes", "p__Actinobacteria", "p__Proteobacteria", "p__Spirochaetes", "p__Verrucomicrobia", "p__Viruses_unclassified")) %>%
  group_by(Taxa) %>%
  summarise(p_value = wilcox.test(RA ~ category)$p.value)



# 5: Species RA -----------------------------------------------------------
## NOT USEFUL PLOT 

# change the rows to column called Seq_ID 
R_rate_wks <- R_rate_wks %>% rownames_to_column("Seq_ID")

# structure the data for the plot
spp_R_category <- filtered_recovery_species %>% 
  gather(Taxa, RA, -Seq_ID) %>%
  group_by(Taxa) %>% 
  inner_join(R_rate_wks) 

# Violin plot with mean and log-scaled y-axis
ggplot(spp_R_category, aes(x = Taxa, y = RA, fill = category)) +
  geom_boxplot(position = position_dodge(width = 0.75)) +
  geom_point(aes(group = category), position = position_dodge(width = 0.75), size = 0.001) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 18,
    size = 0.001,
    position = position_dodge(width = 0.75)
  ) +
  scale_fill_manual(values = c("skyblue2", "goldenrod")) +
  labs(
    title = "Comparison of Species Relative Abundance between Recovery rates",
    x = "Phyla",
    y = "Relative Abundance (log10 scale)"
  ) +
  scale_y_log10() +  # Add log scale to the y-axis
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave("plots/v2/R_v_nR_speciesRA.pdf", h = 20, w = 40)



# 6: Maaslin - Spp --------------------------------------------------------
# Prepare metadata
R_rate_wks <- R_rate_wks %>% 
  rownames_to_column("Seq_ID")

R_cat_data <- R_rate_wks %>% 
  column_to_rownames("Seq_ID")

# Prepare matching species data
R_species <- filter(species, rownames(species) %in% rownames(R_rate_wks))

# there is an additonal data point in the recovery_data set, therefore, we remove it to match
R_cat_data <- R_cat_data %>% rownames_to_column("Seq_ID")
R_species <- R_species %>% rownames_to_column("Seq_ID")
R_cat_data <- semi_join(R_cat_data, R_species, by = "Seq_ID")

# convert the column Seq_ID back to a row for data entry into Maaslin
R_cat_data <- R_cat_data %>%  column_to_rownames("Seq_ID")
R_species <- R_species %>%  column_to_rownames("Seq_ID")

# Rownames need to be in the same order for both metadata and phyla data
R_cat_data <- R_cat_data[sort(rownames(R_cat_data)),]

R_species <- R_species[sort(rownames(R_species)),]
save(recovery_species, file = "data/recovery_species.Rdata")

all(rownames(R_cat_data) == rownames(R_species)) # check they match up

# Run maaslin2
Maaslin2(input_data = R_species, 
         input_metadata = R_cat_data,
         output = "data/maaslin2/R_category", # specify directory to save output files
         fixed_effects = c("category"), # testing variable (can add multiple variables)
         #random_effects = "Feed_type", # have to add given we have multiple samples from the same individual (i.e. non-independent)
         transform = "log", # log transform pathway counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.1, # species needs to be found in atleast 10%%% of samples to be included
         cores = 4) # speeds things up



# Read in output files
all_results_species_R <- read_tsv("data/maaslin2/R_category/all_results.tsv")
sig_results_species_R <- read_tsv("data/maaslin2/R_category/significant_results.tsv")

n_distinct(all_results_species_R$feature) # 148 species included in testing
n_distinct(sig_results_species_R$feature) # 0 species found to be significantly different 

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(all_results_species_R, coef >0)) # 67 species increased
n_distinct(filter(all_results_species_R, coef <0)) # 81 species decreased

# top 10 species that increased in Group 2
head(arrange(all_results_species_recovery, desc(coef)), 10) 

# how many species are reduced in Group 2
n_distinct(filter(all_results_species_recovery, coef <0)) # 71 species are reduced 

# PLOT

all_results_species_R %>% ggplot(aes(x = coef, y = feature, fill = factor(coef > 0))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("salmon2", "yellowgreen")) +
  labs(x = "Effect size and direction of change in abundance",
       y = "Species") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, vjust=0.5, hjust=1, size = 8),
        axis.text.y = element_text(size = 8)) +
  guides(fill = FALSE)  # Remove the legend

ggsave("plots/final/maaslin_sig_species_Gr1_Gr2.jpg", h = 2, w = 6)

# Ordered x axis

ordered_plot <- all_results_species_recovery %>%
  ggplot(aes(x = coef, y = reorder(feature, -coef), fill = factor(coef > 0))) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("salmon2", "yellowgreen")) +
  labs(x = "Effect size and direction of change in abundance", y = "Species") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, vjust = 1, hjust = 1, size = 8),
    axis.text.y = element_text(size = 6)
  ) +
  guides(fill = FALSE)  # Remove the legend

ggsave("plots/v2/maaslin_species_recoveryNOTsig_ordered.jpg", h = 8, w = 6)
ggsave("plots/v2/maaslin_species_recoveryNOTsig_ordered.pdf", h = 8, w = 6)




# 7: Maaslin - Paths ------------------------------------------------------

# Prepare metadata
R_cat_data <- R_rate_wks %>% 
  column_to_rownames("Seq_ID")

# Prepare matching pathways data
R_pathways <- filter(pathways, rownames(pathways) %in% rownames(R_cat_data))

# there is an additonal data point in the group_data set, therefore, we remove it to match
R_cat_data <- R_cat_data %>% rownames_to_column("Seq_ID")
R_pathways <- R_pathways %>% rownames_to_column("Seq_ID")
R_cat_data <- semi_join(R_cat_data, R_pathways, by = "Seq_ID")

# convert the column Seq_ID back to a row for data entry into Maaslin
R_cat_data <-R_cat_data %>%  column_to_rownames("Seq_ID")
R_pathways <- R_pathways %>%  column_to_rownames("Seq_ID")


# Rownames need to be in the same order for both metadata and phyla data
R_cat_data <- R_cat_data[sort(rownames(R_cat_data)),]

R_pathways <- R_pathways[sort(rownames(R_pathways)),]

all(rownames(R_cat_data) == rownames(R_pathways)) # check they match up

# Run maaslin2
Maaslin2(input_data = R_pathways, 
         input_metadata = R_cat_data,
         output = "data/maaslin2/pathways_R", # specify directory to save output files
         fixed_effects = "category", # testing variable (can add multiple variables)
         # random_effects = "Subject_ID", # have to add given we have multiple samples from the same individual (i.e. non-independent)
         transform = "log", # log transform pathway counts
         normalization = "none",
         min_abundance = 0, # no minimum abundance required
         min_prevalence = 0.05, # pathway needs to be found in atleast 5% of samples to be included
         cores = 4) # speeds things up


# Read in output files
all_results_R_pathways <- read_tsv("data/maaslin2/pathways_R/all_results.tsv")
sig_results_R_pathways <- read_tsv("data/maaslin2/pathways_R/significant_results.tsv")

n_distinct(all_results_R_pathways$feature) # 460 pathways included in testing
n_distinct(sig_results_R_pathways$feature) # 0 pathways found to be significantly different between MAM and Healthy

# we can look at the coefficient value to work out whether phyla increased/decreased in relative abundance
n_distinct(filter(all_results_refeed_group_pathways, coef >0)) # 265 pathways increased
n_distinct(filter(all_results_refeed_group_pathways, coef <0)) # 196 pathways decreased

# top 10 pathways that increased following refeed
head(arrange(all_results_refeed_group_pathways, desc(coef)), 10) 

# top 10 pathways that decreased
head(arrange(all_results_refeed_group_pathways, (coef)), 20) 



# PLOT

all_results_refeed_group_pathways %>% ggplot(aes(x = feature, y = coef, fill = factor(coef > 0))) +
  geom_bar(stat = "identity") +
  labs(title = "Community-level differential pathway abundance in Refeed Group 2 comapred to Group 1",
       x = "Community-level pathways",
       y = "Effect size and direction of change in abundance") +
  scale_fill_manual(values = c("salmon2", "yellowgreen")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1, size = 8)) +
  guides(fill = FALSE)  # Remove the legend

ggsave("plots/final/differential_pathways_Gr1_2.jpg", h = 10, w = 25)



# 8: F/B ratio ------------------------------------------------------------
# F/B Ratio calculation for Recovery Category
FBratio_R_category <- phylum_all_samples %>%
  select(p__Firmicutes, p__Bacteroidetes) %>%
  rownames_to_column("Seq_ID") %>%
  mutate(
    p__Firmicutes = p__Firmicutes + 1e-6,
    p__Bacteroidetes = p__Bacteroidetes + 1e-6,
    FBratio = p__Firmicutes / p__Bacteroidetes
  ) %>%
  inner_join(R_rate_wks, by = "Seq_ID")

# Plot
FB_plot_R_category <- ggplot(FBratio_R_category, aes(x = category, y = FBratio, fill = category, color = category)) +
  geom_quasirandom(
    dodge.width = 0.75,
    shape = 21,
    size = 1,
    alpha = 0.7
  ) +
  geom_boxplot(
    outlier.colour = NA,
    alpha = 0.5,
    color = "black",
    linewidth = 0.3,
    position = position_dodge(0.75)
  ) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.7) +
  scale_y_log10(
    limits = c(1e-6, 200),
    breaks = c(1e-6, 1e-4, 1e-2, 1, 100)
  ) +
  scale_fill_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  scale_color_manual(values = c("Fast" = "skyblue2", "Slow" = "goldenrod2")) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.line = element_line(),
    axis.text = element_text(size = 9),
    axis.title = element_text(size = 9),
    axis.title.y = element_text(size = 9),
    legend.position = "none"
  ) +
  xlab("") +
  ylab("F/B Ratio (log10-transformed)")

# Save
ggsave(
  file = "plots/chp2/clean_plots/FB_plot_R_category.pdf",
  plot = FB_plot_R_category,
  units = "cm",
  width = 4,
  height = 6
)

# Wilcoxon test
wilcox_result_R_category <- wilcox.test(FBratio ~ category, data = FBratio_R_category)
# p-value = 0.5517

# Summary
FB_summary_R_category <- FBratio_R_category %>%
  group_by(category) %>%
  summarise(
    mean_firmicutes = mean(p__Firmicutes),
    mean_bacteroidetes = mean(p__Bacteroidetes),
    F_B_ratio = mean_firmicutes / mean_bacteroidetes,
    .groups = "drop"
  )

# category mean_firmicutes mean_bacteroidetes F_B_ratio
# Fast               0.189              0.388     0.488
# Slow               0.223              0.392     0.567

