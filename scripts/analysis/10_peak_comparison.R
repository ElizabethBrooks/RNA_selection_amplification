#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

# R script to create analysis plots for the RNA evolution project

# turn of scientific notation
options(scipen=10000)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/10b_peak_comparison"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# read in cluster family sequence data
#r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/08_summarized_1500/r8_S8_L001_cluster_peaks_table.csv")
r8_peaks <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/07_clustered_cd_hit/r8_S8_L001_formatted_above2.clstr_peak_data.txt")

# set the sequence length
seqLength <- 40

# setup data frame length
data_length <- nrow(r8_peaks) * nrow(r8_peaks)

# data frame of sequence data
seqs_data <- data.frame(
  run_name = rep(NA, data_length),
  sequence_ID = rep(NA, data_length),
  sequence = rep(NA, data_length),
  peak_cluster_ID = rep(NA, data_length),
  peak_identity = rep(NA, data_length)
)

# loop over each sequence and compare with the peak to note if it has >= 95% similarity
for (seq_num in 1:data_length) {
  # loop over each peak sequence
  for (peak_num in 1:nrow(r8_peaks)) {
    # compare the current sequence with the current peak
    numMatch <- mapply(
      function(x, y) {
        len <- length(x)
        sum(x[1:len] == y[1:len])
      }, 
      strsplit(r8_peaks$sequence[seq_num], ''), 
      strsplit(r8_peaks$sequence[peak_num], '')
    )
    # determine identity percent
    identity_perc <- 100*numMatch/seqLength
    # set the index
    index <- ((seq_num-1)*nrow(r8_peaks)) + peak_num
    # record data
    seqs_data$run_name[index] <- r8_peaks$run_name[seq_num]
    seqs_data$sequence_ID[index] <- r8_peaks$sequence_ID[seq_num]
    seqs_data$sequence[index] <- r8_peaks$sequence[seq_num]
    seqs_data$peak_cluster_ID[index] <- r8_peaks$cluster_ID[peak_num]
    seqs_data$peak_identity[index] <- identity_perc
  }
}

# export data
write.csv(seqs_data, file = paste(out_dir, "/peak_identities.csv", sep = ""), row.names = FALSE, quote = FALSE)

# check how many sequences have >= 90% identity to each peak
for (cluster_num in 0:(nrow(r8_peaks)-1)) {
  print(cluster_num)
  print(nrow(seqs_data[seqs_data[seqs_data$peak_cluster_ID == cluster_num,]$peak_identity >= 90,]))
}

# keep sequences that have >= 90% identity to any peak
seqs_out <- seqs_data[seqs_data$peak_identity >= 90,]

# export data
write.csv(seqs_out, file = paste(out_dir, "/peak_identities_atLeast90.csv", sep = ""), row.names = FALSE, quote = FALSE)
