#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/S1_cluster_sequence_identity"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in cluster identity data
r8_seqs_identity <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09_identified/07a_clustered/r8_S8_L001_formatted_above9_cluster_sequences_identity_table.csv")

# list of cluster IDs
cluster_list <- unique(r8_seqs_identity$cluster_ID)

# calculate cluster IDs and abundances
r8_high_quality <- 889374
cluster_read_counts <- c(40999, 256601, 34117, 82079, 14439, 26345, 6555, 12539, 15921, 2572, 4227, 8543, 4043)
cluster_seq_counts <- c(108, 1114, 98, 116, 269, 85, 33, 49, 75, 153, 275, 286, 271)
cluster_read_abun <- 100*cluster_read_counts/r8_high_quality

# create cluster data frame
cluster_data <- data.frame(
  cluster_ID = cluster_list,
  read_counts = cluster_read_counts,
  seq_counts = cluster_seq_counts,
  read_abun = cluster_read_abun
)

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$read_abun, decreasing = TRUE),]  

# add family numbers
cluster_data$fam_num <- seq(from = 1, to = 13, by = 1)

# setup data frame length
data_length <- length(cluster_list)

# data frame of cluster abundances
cluster_identities <- data.frame(
  fam_ID = rep(NA, data_length),
  cluster_ID = rep(NA, data_length), 
  avg_ID = rep(NA, data_length),
  highest_ID = rep(NA, data_length),
  lowest_ID = rep(NA, data_length)
)

# loop over each cluster
for (cluster_num in 0:max(cluster_list)) {
  # setup index
  index <- cluster_num+1
  # update family number for publishing
  cluster_out <- cluster_data[cluster_data$cluster_ID == cluster_num, "fam_num"]
  # add family number
  cluster_identities$fam_ID[index] <- cluster_out
  # add cluster ID
  cluster_identities$cluster_ID[index] <- cluster_num
  # add avg, lowest, and highest IDs 
  cluster_identities$avg_ID[index] <- mean(r8_seqs_identity[r8_seqs_identity$cluster_ID == cluster_num, "percent_ID"])
  cluster_identities$highest_ID[index] <- max(r8_seqs_identity[r8_seqs_identity$cluster_ID == cluster_num, "percent_ID"])
  cluster_identities$lowest_ID[index] <- min(r8_seqs_identity[r8_seqs_identity$cluster_ID == cluster_num, "percent_ID"])
}

# line plots with round 8 clustering identity info
r8_peaks_identity_plot <- ggplot(data=cluster_identities, aes(x = cluster_ID)) +
  geom_line(aes(y = avg_ID, color = "Avg"), size = 1) +
  #geom_line(aes(y = highest_ID, color = "Highest")) +
  geom_line(aes(y = lowest_ID, color = "Lowest"), size = 1) +
  geom_point(aes(y = avg_ID, color = "Avg")) +
  #geom_point(aes(y = highest_ID, color = "Highest")) +
  geom_point(aes(y = lowest_ID, color = "Lowest")) +
  #scale_color_manual(name = "Statistic", values = c("Avg" = safe_colors[5], "Highest" = safe_colors[9], "Lowest" = safe_colors[4])) +
  scale_color_manual(name = "Statistic", values = c("Avg" = safe_colors[5], "Lowest" = safe_colors[4])) +
  theme_classic() +
  scale_x_continuous(labels = as.character(cluster_identities$cluster_ID), breaks = cluster_identities$cluster_ID) +
  ylab("Percent Identity to Peak") +
  xlab("Cluster Number")
# save the plot
exportFile <- paste(out_dir, "/r8_sequence_identity.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(r8_peaks_identity_plot)
dev.off()

# export plotting data
write.csv(cluster_identities, file = paste(out_dir, "/r8_sequence_identity.csv", sep = ""), row.names = FALSE, quote = FALSE)
