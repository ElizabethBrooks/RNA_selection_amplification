#!/usr/bin/env Rscript

# R script to create family data tables for plotting
# usage: 01_ST2_family_tables.R

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/ST2_family_table"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# read in cluster family sequence data
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv")

# read in sequence family count data
#seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification/family_identities_atLeast90.csv")
r8_seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/r8_S8_L001_family_identities_max_atLeast90.csv")

# set the round num
roundNum <- 8

# list of cluster IDs
cluster_list <- r8_peaks$cluster_ID

# color blind safe plotting palette
#safe_colors <- c(carto_pal(name="Safe"), "#000000")
safe_colors <- carto_pal(name="Safe")[1:length(cluster_list)]

# create cluster data frame
cluster_data <- data.frame(
  cluster_ID = cluster_list,
  family_ID = rep(NA, length(cluster_list)),
  seq_counts = rep(NA, length(cluster_list)),
  read_counts = rep(NA, length(cluster_list)),
  read_abun = rep(NA, length(cluster_list)),
  sequence = r8_peaks$sequence,
  cluster_color = safe_colors
)

# add read and sequence counts
for (cluster_num in min(cluster_data$cluster_ID):max(cluster_data$cluster_ID)) {
  # add read counts
  cluster_data$read_counts[cluster_num+1] <-  sum(r8_seqs_family[r8_seqs_family$peak_cluster_ID == cluster_num,"read_counts"])
  # add sequence counts
  cluster_data$seq_counts[cluster_num+1] <-  nrow(r8_seqs_family[r8_seqs_family$peak_cluster_ID == cluster_num,])
}

# add abundances
cluster_data$read_abun <- 100*cluster_data$read_counts/quality[roundNum]

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# add family numbers
cluster_data$family_ID <- seq(from = 1, to = length(cluster_list), by = 1)
  
# export data
write.csv(cluster_data, file = paste(out_dir, "/r8_family_count_data.csv", sep = ""), row.names = FALSE, quote = FALSE)
