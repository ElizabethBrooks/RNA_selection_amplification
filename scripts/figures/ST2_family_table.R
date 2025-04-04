#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/ST2_family_table"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# read in cluster family sequence data
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized/r8_S8_L001_cluster_peaks_table.csv")
#seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification/family_identities_atLeast90.csv")
seqs_family <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification/family_identities.csv")

# read in sequence count data
counts_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# subset to the round 8 data
r8_seqs_family <- seqs_family[seqs_family$run_name == "8",]
r8_counts_input <- counts_input[counts_input$run_name == "8" & counts_input$counts_run_name == "r8_S8_L001",]

# list of unique seqeunces
unique_seqs <- r8_counts_input[unique(r8_counts_input$sequence),"sequence"]

# loop over each unique sequence
for (seq in 1:length(unique_seqs)) {
  
}
# To-do: keep all counts and add sequence identities
# merge the counts with the identities
r8_family_data <- merge(r8_counts_input, r8_seqs_family, by = "sequence_ID", all.x = TRUE)

# set the round num
roundNum <- 8

# list of cluster IDs
#cluster_list <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
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
for (cluster_num in 0:max(cluster_data$cluster_ID)) {
  # To-do: make sure r8 sequence family data is sorted by increasing cluster ID
  # add read counts
  cluster_data$read_counts[cluster_num] <-  sum(r8_seqs_family[r8_seqs_family$peak_cluster_ID == cluster_num,])
  # add sequence counts
  cluster_data$seq_counts[cluster_num] <-  nrow(r8_seqs_family[r8_seqs_family$peak_cluster_ID == cluster_num,])
}

# add abundances
cluster_data_out$read_abun <- 100*cluster_data_out$read_counts/quality[roundNum]

# sort cluster data
cluster_data_out <- cluster_data_out[order(cluster_data_out$read_abun, decreasing = TRUE),]  

# add family numbers
cluster_data_out$fam_num <- seq(from = 1, to = length(cluster_list), by = 1)
  
# export data
write.csv(cluster_data_out, file = paste(out_dir, "/r8_family_count_data.csv", sep = ""), row.names = FALSE, quote = FALSE)
