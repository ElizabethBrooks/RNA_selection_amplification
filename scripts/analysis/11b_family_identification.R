#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

args = commandArgs(trailingOnly=TRUE)

# R script to identify sequence families

# install missing packages
#install.packages("dplyr")

# libraries
library(dplyr)

# turn of scientific notation
options(scipen=10000)

# set the sequence length
seqLength <- 40

# set the input round number
#round_num <- "1"
round_num <- args[1]
round_name <- paste("r", round_num, "_S", round_num, "_L001", sep = "")

# create outputs directory
#out_dir <- "/scratch365/ebrooks5/RNA_evolution/outputs/11a_family_identification_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification_above2"
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11a_family_identification_all"
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in cluster family sequence data
#peaksFile <- "/scratch365/ebrooks5/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"
#peaksFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv"
peaksFile <- args[3]
r8_peaks <- read.csv(peaksFile)

# read in sequence count data for the specified round
#seqsFile <- "/scratch365/ebrooks5/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv"
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09b_quantified_above2/counts_plot_table.csv"
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/counts_plot_table.csv"
seqsFile <- args[4]
seqs_input <- read.csv(seqsFile, colClasses=c("run_name"="character", "counts_run_name"="character"))

# subset to the round data
seqs_input <- seqs_input[seqs_input$run_name == round_num & seqs_input$counts_run_name == round_name,]

# remove NAs
#seqs_input <- na.omit(seqs_input)

# reverse complement the sequences
#seqs_input$sequence <- rev(chartr("ATGC","TACG",seqs_input$sequence))

# setup data frame length
#data_length <- nrow(seqs_input) * nrow(r8_peaks)
data_length <- nrow(seqs_input)

# data frame of sequence data
seqs_data <- data.frame()

# loop over each sequence and compare with the peak sequences
for (seq_num in 0:(data_length-1)) {
  # loop over each peak sequence
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
    # record data
    curr_data <- seqs_input[seq_num+1,]
    curr_data$peak_cluster_ID <- r8_peaks$cluster_ID[peak_num]
    curr_data$peak_identity <- identity_perc
    seqs_data <- rbind(seqs_data, curr_data)
  }
}

# export data
write.csv(seqs_data, file = paste(out_dir, "/", round_name, "_family_identities.csv", sep = ""), row.names = FALSE, quote = FALSE)

# data frame to track how many sequences have >= 90% identity to each peak
#identity_check <- data.frame(
#  cluster = rep(NA, nrow(r8_peaks)),
#  identity_count = rep(NA, nrow(r8_peaks))
#)

# check how many sequences have >= 90% identity to each peak
#for (cluster_num in 0:(nrow(r8_peaks)-1)) {
#  clusterIndex <- cluster_num+1
#  identity_check$cluster[clusterIndex] <- cluster_num
#  clusterSize <- 0
#  clusterSize <- clusterSize + nrow(seqs_data[seqs_data[seqs_data$peak_cluster_ID == cluster_num,]$peak_identity >= 90,])
#  identity_check$identity_count[clusterIndex] <- clusterSize
#}
#identity_check

# keep sequences that have >= 90% identity to any peak
seqs_out <- seqs_data[seqs_data$peak_identity >= 90,]

# export data
write.csv(seqs_out, file = paste(out_dir, "/", round_name, "_family_identities_atLeast90.csv", sep = ""), row.names = FALSE, quote = FALSE)

# subset to the round 8 data
#r8_seqs_family <- seqs_data[seqs_data$run_name == "8" & seqs_data$counts_run_name == "r8_S8_L001",]

# set of unique sequences
#unique_seqs <- seqs_data[!duplicated(seqs_data$sequence),]
unique_seqs <- seqs_data %>% distinct(sequence)

# data frame for sequence counts and identities
family_data <- data.frame()

# loop over each unique sequence
for (seq in 1:nrow(unique_seqs)) {
  # retrieve current sequence data
  curr_seq <- unique_seqs$sequence[seq]
  curr_seq_data <- seqs_data[seqs_data$sequence == curr_seq,]
  # retrieve max identity data
  max_id_data <- curr_seq_data[which.max(curr_seq_data$peak_identity),]
  # add data for current sequence to the outputs
  family_data <- rbind(family_data, max_id_data)
}

# export data
write.csv(family_data, file = paste(out_dir, "/", round_name, "_family_identities_max.csv", sep = ""), row.names = FALSE, quote = FALSE)

# keep sequences that have >= 90% identity to any peak
family_data_out <- family_data[family_data$peak_identity >= 90,]

# export data
write.csv(family_data_out, file = paste(out_dir, "/", round_name, "_family_identities_max_atLeast90.csv", sep = ""), row.names = FALSE, quote = FALSE)
