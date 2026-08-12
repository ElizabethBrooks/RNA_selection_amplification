#!/usr/bin/env Rscript

# R script to create family base conservation plots
# usage: 06_family1_co_occuring_bases.R

# turn of scientific notation
options(scipen=9999)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(stringr)
library(dplyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/family1_co_occuring_bases/alternates_b20"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# set working directory
setwd(out_dir)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# read in family data
fam_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/ST2_family_table/r8_family_count_data.csv")

# read in cluster family data
cluster_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max_atLeast90.csv")

# subset the counts data to round 8
cluster_data <- cluster_data[cluster_data$run_name == 8 & cluster_data$counts_run_name == "r8_S8_L001",]

# subset the counts data to family 1 (cluster 1)
cluster_data <- cluster_data[cluster_data$peak_cluster_ID == 1,]

# remove duplicates
#cluster_data <- cluster_data[!duplicated(cluster_data$sequence),]

# move sequence column to front
cluster_data <- cluster_data %>% relocate(sequence) 

###
## most common base (base 20 = G) alternate U
###

# subset to counts data
cluster_data_conserved <- cluster_data[substr(cluster_data$sequence, 20, 20) == "T", ]

# initialize data frame for base counts
base_counts <- data.frame(
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4)
)

# initialize data frame for sequences
seqs_matrix <- data.frame()

# convert list of sequences into a matrix
seqs_matrix <- do.call(rbind, type.convert(strsplit(cluster_data_conserved[, "sequence"], ""), as.is = TRUE))

# loop over each base
for (base_num in 1:40) {
  # update indicies
  index <- ((base_num-1)*4)+1
  index_max <- index+3
  # add base ID
  base_counts$base_ID[index:index_max] <- rep(base_num, 4)
  # add percent conservation of each base character
  base_counts$conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
  base_counts$conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
  base_counts$conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
  base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
  # add each base character
  base_counts$base[index] <- "A"
  base_counts$base[index+1] <- "C"
  base_counts$base[index+2] <- "G"
  base_counts$base[index_max] <- "U"
}

# change zeros to NAs for plotting
base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
# set family plot title
run_title <- "Family 1 Nucleotide Conservation"
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
exportFile <-"family_1_base_conservation_b20_U.png"
png(exportFile, units="in", width=20, height=3, res=300)
print(base_counts_plot)
dev.off()

# clean up sequences data frame
cluster_data_conserved_out <- select(cluster_data_conserved, -c(sequence_ID, counts, counts_run_name))

# write sequences to file
write.csv(cluster_data_conserved_out, "family_1_sequences_b20_U.csv", quote = FALSE, row.names = FALSE)

###
## most common base (base 20 = G) alternate C
###

# subset to counts data
cluster_data_conserved <- cluster_data[substr(cluster_data$sequence, 20, 20) == "C", ]

# initialize data frame for base counts
base_counts <- data.frame(
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4)
)

# initialize data frame for sequences
seqs_matrix <- data.frame()

# convert list of sequences into a matrix
seqs_matrix <- do.call(rbind, type.convert(strsplit(cluster_data_conserved[, "sequence"], ""), as.is = TRUE))

# loop over each base
for (base_num in 1:40) {
  # update indicies
  index <- ((base_num-1)*4)+1
  index_max <- index+3
  # add base ID
  base_counts$base_ID[index:index_max] <- rep(base_num, 4)
  # add percent conservation of each base character
  base_counts$conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
  base_counts$conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
  base_counts$conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
  base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
  # add each base character
  base_counts$base[index] <- "A"
  base_counts$base[index+1] <- "C"
  base_counts$base[index+2] <- "G"
  base_counts$base[index_max] <- "U"
}

# change zeros to NAs for plotting
base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
# set family plot title
run_title <- "Family 1 Nucleotide Conservation"
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
exportFile <-"family_1_base_conservation_b20_C.png"
png(exportFile, units="in", width=20, height=3, res=300)
print(base_counts_plot)
dev.off()

# clean up sequences data frame
cluster_data_conserved_out <- select(cluster_data_conserved, -c(sequence_ID, counts, counts_run_name))

# write sequences to file
write.csv(cluster_data_conserved_out, "family_1_sequences_b20_C.csv", quote = FALSE, row.names = FALSE)

###
## most common base (base 20 = G) alternate A
###

# subset to counts data
cluster_data_conserved <- cluster_data[substr(cluster_data$sequence, 20, 20) == "A", ]

# initialize data frame for base counts
base_counts <- data.frame(
  base_ID = rep(NA, 40*4),
  base = rep(NA, 40*4),
  conservation = rep(NA, 40*4)
)

# initialize data frame for sequences
seqs_matrix <- data.frame()

# convert list of sequences into a matrix
seqs_matrix <- do.call(rbind, type.convert(strsplit(cluster_data_conserved[, "sequence"], ""), as.is = TRUE))

# loop over each base
for (base_num in 1:40) {
  # update indicies
  index <- ((base_num-1)*4)+1
  index_max <- index+3
  # add base ID
  base_counts$base_ID[index:index_max] <- rep(base_num, 4)
  # add percent conservation of each base character
  base_counts$conservation[index] <- 100*sum(str_count(seqs_matrix[,base_num], "A"))/nrow(seqs_matrix)
  base_counts$conservation[index+1] <- 100*sum(str_count(seqs_matrix[,base_num], "C"))/nrow(seqs_matrix)
  base_counts$conservation[index+2] <- 100*sum(str_count(seqs_matrix[,base_num], "G"))/nrow(seqs_matrix)
  base_counts$conservation[index_max] <- 100*sum(str_count(seqs_matrix[,base_num], "T"))/nrow(seqs_matrix)
  # add each base character
  base_counts$base[index] <- "A"
  base_counts$base[index+1] <- "C"
  base_counts$base[index+2] <- "G"
  base_counts$base[index_max] <- "U"
}

# change zeros to NAs for plotting
base_counts$conservation_na <- ifelse(base_counts$conservation == 0, NA, base_counts$conservation)
# set family plot title
run_title <- "Family 1 Nucleotide Conservation"
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
exportFile <-"family_1_base_conservation_b20_A.png"
png(exportFile, units="in", width=20, height=3, res=300)
print(base_counts_plot)
dev.off()

# clean up sequences data frame
cluster_data_conserved_out <- select(cluster_data_conserved, -c(sequence_ID, counts, counts_run_name))

# write sequences to file
write.csv(cluster_data_conserved_out, "family_1_sequences_b20_A.csv", quote = FALSE, row.names = FALSE)
