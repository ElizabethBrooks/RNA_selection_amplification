#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# supress cowplot package messages
suppressMessages( require(cowplot) )

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/plots/04_top10_sequence_counts"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# round numbers
rounds <- c(1, 2, 3, 4, 5, 6, 7, 8)

# % diversity per round
#diversity_doped <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20, 97.30, 92.40, 86.43)
diversity <- c(99.67, 99.66, 99.65, 99.62, 98.44, 54.61, 15.86, 12.20)

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in sequence count data
seqs_counts <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_top10/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))

# reverse complement the sequences
#seqs_counts$sequence <- rev(chartr("ATGC","TACG",seqs_counts$sequence))

# add the fraction abundances
seqs_counts$frac_abundance <- NA
for (run_num in 1:8) {
  seqs_counts[seqs_counts$counts_run_name == run_num, "frac_abundance"] <- seqs_counts[seqs_counts$counts_run_name == run_num, "counts"]/quality[run_num]
}

# add log values
seqs_counts$log_counts <- log(seqs_counts$counts)
#seqs_counts$log_frac_abundance <- log(seqs_counts$frac_abundance)

# add percent abundance
seqs_counts$perc_abundance <- 100*seqs_counts$frac_abundance

# set infinite and NA values equal to zero
is.na(seqs_counts)<-sapply(seqs_counts, is.infinite)

# re-order the data for plotting
seqs_counts <- seqs_counts[order(seqs_counts$log_counts, decreasing=TRUE),]

# setup midpoint values for plotting
seqs_counts_noNA <- seqs_counts
seqs_counts_noNA[is.na(seqs_counts_noNA)] <- 0
mid_log_counts <- max(seqs_counts_noNA$log_counts)/2
#mid_log_frac_abundance <- log(min(seqs_counts_noNA[seqs_counts_noNA$frac_abundance != 0,"frac_abundance"]))/2
mid_perc_abundance <- max(seqs_counts$perc_abundance)/2

# initialize plot list
counts_heatmap_list <- vector('list', 8)
percent_heatmap_list <- vector('list', 8)

# loop over each round and create heatmaps
for (run_num in 1:8) {
  # subset seq data
  seqs_counts_subset <- seqs_counts[seqs_counts$run_name == run_num,]
  # remove run tags
  seqs_counts_subset[c("run", "ID")] <- do.call(rbind, strsplit(seqs_counts_subset$sequence_ID, "_"))
  # set round plot title
  run_title <- paste("Round", run_num)
  # create heatmap of log counts
  counts_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(ID), log_counts), fill= log_counts)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Sequence ID") +
    xlab("Round Number") +
    scale_fill_gradient2(name = "Log Counts",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = mid_log_counts,
                         na.value = "white")
  # save the plot
  exportFile <- paste(out_dir, "/r", run_num, "_top10_sequence_log_counts.png", sep = "")
  png(exportFile, units="in", width=5, height=5, res=300)
  print(counts_heatmap_subset)
  dev.off()
  # add the log counts heatmap to the plot list
  counts_heatmap_list[[run_num]] <- counts_heatmap_subset
  # create heatmap of percent abundances
  percent_heatmap_subset <- ggplot(data = seqs_counts_subset, aes(counts_run_name, reorder(as.character(ID), perc_abundance), fill= perc_abundance)) + 
    theme_bw() +
    geom_tile(colour = "black") +
    ggtitle(run_title) +
    theme(plot.title = element_text(hjust = 0.5)) +
    ylab("Sequence ID") +
    xlab("Round Number") +
    scale_fill_gradient2(name = "PA",
                         low = safe_colors[3],
                         mid = safe_colors[4],
                         high = safe_colors[5],
                         midpoint = mid_perc_abundance,
                         na.value = "white")
  # save the plot
  exportFile <- paste(out_dir, "/r", run_num, "_top10_sequence_percent_abundances.png", sep = "")
  png(exportFile, units="in", width=5, height=5, res=300)
  print(percent_heatmap_subset)
  dev.off()
  # add the heatmap to the plot list
  percent_heatmap_list[[run_num]] <- percent_heatmap_subset
}

# combine heatmaps with the log counts for each of the top 10 sequences per round
grid_log_counts_plot <- plot_grid(plotlist=counts_heatmap_list,  ncol = 4, align = 'v')
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_log_counts.png", sep = "")
png(exportFile, units="in", width=20, height=10, res=300)
print(grid_log_counts_plot)
dev.off()

# combine heatmaps with the percent abundances for each of the top 10 sequences per round
grid_perc_abundance_plot <- plot_grid(plotlist=percent_heatmap_list,  ncol = 4, align = 'v')
# save the plot
exportFile <- paste(out_dir, "/top10_sequence_percent_abundances.png", sep = "")
png(exportFile, units="in", width=20, height=10, res=300)
print(grid_perc_abundance_plot)
dev.off()

# export plotting data
write.csv(seqs_counts, file = paste(out_dir, "/top10_sequence_counts.csv", sep = ""), row.names = FALSE, quote = FALSE)
