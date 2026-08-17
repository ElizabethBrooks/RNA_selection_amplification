#!/usr/bin/env Rscript

# R script to create family base conservation plots
# usage: 04_F4A_family_base_conservation.R

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(stringr)

# set outputs directory
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2_r8"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2_unique"
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2_r8_unique"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in family data
fam_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/ST2_family_table/r8_family_count_data.csv")

# read in cluster family data
cluster_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max_atLeast90.csv")

# subset the counts data
cluster_data <- cluster_data[cluster_data$run_name == 8 & cluster_data$counts_run_name == "r8_S8_L001",]
cluster_data <- cluster_data[!duplicated(cluster_data$sequence),]

# loop over each cluster
for (cluster_num in 0:9) {
  # set the family IDs
  cluster_data[cluster_data$peak_cluster_ID == cluster_num, "family_ID"] <- fam_data[fam_data$cluster_ID == cluster_num, "family_ID"]
}

# re-name cluster ID column
names(cluster_data)[names(cluster_data) == 'peak_cluster_ID'] <- 'cluster_ID'

# output sequence data with family IDs
write.csv(cluster_data, file = paste(out_dir, "/family_sequences.csv", sep = ""), row.names = FALSE, quote = FALSE)

# list of cluster IDs
cluster_list <- unique(cluster_data$cluster_ID)

# initialize data frame for base counts
base_counts <- data.frame(
  cluster_ID = rep(NA, 40*4),
  fam_ID = rep(NA, 40*4),
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4),
  read_conservation = rep(NA, 40*4)
)
base_counts_out <- data.frame()

# initialize data frame for sequences
seqs_matrix <- data.frame()

# loop over each cluster
for (cluster_num in min(cluster_list):max(cluster_list)) {
#for (cluster_num in 1:1) {
  #cluster_num <- 1
  # convert list of sequences into a matrix
  seqs_matrix <- do.call(rbind, type.convert(strsplit(cluster_data[cluster_data$cluster_ID == cluster_num, "sequence"], ""), as.is = TRUE))
  # set family number
  fam_num <- head(cluster_data[cluster_data$cluster_ID == cluster_num, "family_ID"], 1)
  # loop over each base
  for (base_num in 1:40) {
    # update indicies
    index <- ((base_num-1)*4)+1
    index_max <- index+3
    # add cluster ID
    base_counts$cluster_ID[index:index_max] <- rep(cluster_num, 4)
    # add family ID
    base_counts$fam_ID[index:index_max] <- rep(fam_num, 4)
    # add base ID
    base_counts$base_ID[index:index_max] <- rep(base_num, 4)
    # add percent conservation of each base character
    base_counts$conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
    base_counts$conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
    base_counts$conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
    base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
    # to-do: update read conservation
    # add percent read conservation of each base character
    base_counts$read_conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
    base_counts$read_conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
    base_counts$read_conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
    base_counts$read_conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
    # add each base character
    base_counts$base[index] <- "A"
    base_counts$base[index+1] <- "C"
    base_counts$base[index+2] <- "G"
    base_counts$base[index_max] <- "U"
  }
  # change zeros to NAs for plotting
  base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
  # change the bases to factors for plotting
  #base_counts$base <- factor(base_counts$base, levels=c("U", "G", "C", "A"))
  # set family plot title
  run_title <- paste("Family", fam_num, "Nucleotide Conservation")
  # create heatmap of base conservation
  base_counts_plot <- ggplot(data = base_counts, aes(reorder(as.character(base_ID), base_ID), base, fill= conservation_na)) + 
    theme_classic(base_size = 18) +
    geom_tile(colour = "black") +
    # left P2
    annotate("rect", xmin = c(5.5), xmax = c(12.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1.5) +
    # right P2
    annotate("rect", xmin = c(25.5), xmax = c(32.5), ymin = c(0.5), ymax = c(4.5), 
             colour = safe_colors[2], fill = "transparent", linewidth = 1.5) +
    # overhang compliment
    annotate("rect", xmin = c(15.5), xmax = c(23.5), ymin = c(0.5), ymax = c(4.5),
             #xmin = c(15.5, 16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5), xmax = c(16.5, 17.5, 18.5, 19.5, 20.5, 21.5, 22.5, 23.5, 24.5, 25.5), 
             #ymin = c(3.5, 2.5, 1.5, 2.5, 2.5, 0.5, 0.5, 3.5, 2.5, 1.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5, 1.5, 1.5, 4.5, 3.5, 2.5), 
             #xmin = c(15.5, 16.5, 17.5, 18.5, 19.5), xmax = c(16.5, 17.5, 18.5, 19.5, 20.5), 
             #ymin = c(3.5, 2.5, 1.5, 2.5, 2.5), ymax = c(4.5, 3.5, 2.5, 3.5, 3.5), 
             colour = safe_colors[1], fill = "transparent", linewidth = 1.5) +
    ylab("Nucleotide") +
    xlab("Nucleotide Position") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 14), 
          axis.text.y = element_text(size = 14),
          axis.title = element_text(size = 18)) +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5, size = 18)) +
    scale_fill_gradient2(name = "Percent\nAbundance",
                         low = "#F0E442",
                         mid = safe_colors[7],
                         high = safe_colors[5],
                         midpoint = max(base_counts$conservation)/2,
                         na.value = "white") +
    coord_fixed() +
    geom_text(aes(label = round(conservation_na, digits = 2)), size = 4, colour =ifelse(base_counts$conservation_na < 50, "black", "white"))
  # save the plot
  exportFile <- paste(out_dir, "/family", fam_num, "_base_conservation.png", sep = "")
  png(exportFile, units="in", width=20, height=3, res=300)
  print(base_counts_plot)
  dev.off()
  # add current family data to outputs
  base_counts_out <- rbind(base_counts_out, base_counts)
}

# export plotting data
write.csv(base_counts_out, file = paste(out_dir, "/family_base_conservation.csv", sep = ""), row.names = FALSE, quote = FALSE)
