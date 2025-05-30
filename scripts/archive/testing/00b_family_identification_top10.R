#!/usr/bin/env Rscript

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
#library(plyr)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/00b_family_identification_top10"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)

# read in cluster family sequence data
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized/07a_clustered/r8_S8_L001_formatted_above9_cluster_peaks_table.csv")

# read in sequence count data
seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11_quantified_top10/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
  
# retrieve data from run 8 and remove duplicate sequence data
#seqs_input <- seqs_input[seqs_input$counts_run_name == "1",]
seqs_input <- seqs_input[!duplicated(seqs_input$sequence_ID),]

# setup data frame length
data_length <- nrow(seqs_input) * nrow(r8_peaks)

# data frame of sequence data
seqs_data <- data.frame(
  run_name = rep(NA, data_length),
  sequence_ID = rep(NA, data_length),
  sequence = rep(NA, data_length),
  peak_cluster_ID = rep(NA, data_length),
  peak_identity = rep(NA, data_length)
)

# set the sequence length
seqLength <- 40

# loop over each sequence and compare with the peak to note if it has >= 90% similarity
for (seq_num in 0:(data_length-1)) {
  for (peak_num in 1:nrow(r8_peaks)) {
    # compare the current sequence with the current peak
    numMatch <- mapply(
      function(x, y) {
        len <- length(x)
        sum(x[1:len] == y[1:len])
      }, 
      strsplit(seqs_input$sequence[seq_num+1], ''), 
      strsplit(r8_peaks$sequence[peak_num], '')
    )
    # determine identity percent
    identity_perc <- 100*numMatch/seqLength
    # set the index
    index <- ((seq_num)*nrow(r8_peaks)) + peak_num
    # record data
    seqs_data$run_name[index] <- seqs_input$run_name[seq_num+1]
    seqs_data$sequence_ID[index] <- seqs_input$sequence_ID[seq_num+1]
    seqs_data$sequence[index] <- seqs_input$sequence[seq_num+1]
    seqs_data$peak_cluster_ID[index] <- r8_peaks$cluster_ID[peak_num]
    seqs_data$peak_identity[index] <- identity_perc
  }
}

# export data
write.csv(seqs_data, file = paste(out_dir, "/family_identities_top10.csv", sep = ""), row.names = FALSE, quote = FALSE)

# check how many sequences have at >= 90% identity to each peak
for (cluster_num in 0:(nrow(r8_peaks)-1)) {
  print(cluster_num)
  print(nrow(seqs_data[seqs_data[seqs_data$peak_cluster_ID == cluster_num,]$peak_identity >= 90,]))
}
