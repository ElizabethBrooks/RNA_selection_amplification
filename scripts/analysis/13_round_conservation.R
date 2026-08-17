#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

args = commandArgs(trailingOnly=TRUE)

# R script to identify concerved overhang sequences

# set the library paths for the CRC servers
.libPaths("/afs/crc.nd.edu/user/e/ebrooks5/R/x86_64-pc-linux-gnu-library/4.4")

# turn of scientific notation
options(scipen=10000)

# import libraries
#library(ggplot2)
library(scales)
#library(rcartocolor)
library(stringr)

# set the input round number
#round_num <- "8"
round_num <- args[1]
round_name <- paste("r", round_num, "_S", round_num, "_L001", sep = "")

# set outputs directory
#out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/figures/F5_overhang_conservation_all"
out_dir <- args[2]
dir.create(out_dir, showWarnings = FALSE)

# read in sequence data
#seqs_input <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified/counts_plot_table_noDoped.csv", colClasses=c("run_name"="character", "counts_run_name"="character"))
#seqsFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/09a_quantified_all/r8_S8_L001_in_r8_S8_L001_counts_plot_table.csv"
seqsFile <- args[3]
seqs_input <- read.csv(seqsFile, colClasses=c("run_name"="character", "counts_run_name"="character"))

# subset to the round data
seqs_input <- seqs_input[seqs_input$run_name == round_num & seqs_input$counts_run_name == round_name,]

# pad the sequences for the sliding window
seqs_input$sequence <- paste0("XXXXXX", seqs_input$sequence, "XXXXXX")

# color blind safe plotting palette
#safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# A pairs with U and G pairs with C, but G also pairs with U
# G (C, U)
# A (U)
# U (A, G)
# C (G)

# overhang sequence: AUUCCGCA
# overhang wobble positions: 2, 3, and 6
# complement wobble positions: 3, 6, and 7
# 3
# 3, 6
# 3, 7
# 3, 6, 7
# 6
# 6, 7
# 7

# set the overhang sequence
overhang <- rev(c("ATTCCGCA"))

# store the overhang in an array
overhang <- unlist(strsplit(overhang, ""))

# reverse the overhang for comparing
rev_overhang <- rev(overhang)

# set expected overhang (ATTCCGCA) reverse complement (TGCGGAAT) sequence
#complement_seq <- c("T","G","C","G","G","A","A","T","G","C")
#complement_seq <- c("TGCGGAATGC")
#complement_seq <- c("T","G","C","G","G")
complement_seq <- c("TGCGGAAT")

# store the expected complement in an array
complement_seq <- unlist(strsplit(complement_seq, ""))

# wobble positions (TG[C,T]GG[A,G][A,G]T)
# TG[T]GGAAT
# TG[T]GG[G]AT
# TG[T]GGA[G]T
# TG[T]GG[G][G]T
# TGCGG[G]AT
# TGCGG[G][G]T
# TGCGGA[G]T

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#filtered_reads <- c(18, 19, 26, 27, 1585, 10626, 7230, 6315)

# list of round IDs
round_list <- unique(seqs_input$run_name)

# convert list of sequences into a matrix
seq_matrix <- do.call(rbind, type.convert(strsplit(seqs_input$sequence, ""), as.is = TRUE))

# trim the sequence to keep the bases of the overhang
#seq_matrix <- seqs_input_matrix[,16:20]

# set seqeunce and overhang complement lengths
seq_length <- 6 + 40 + 6
complement_length <- 8

# set minimum length and identity
min_length <- 2
min_subset_length <- 3

# set length for window sliding
sliding_length <- seq_length - (complement_length-1)

# set overhang data length
seq_data_length <- nrow(seqs_input)

# initialize data frame for base counts
complement_data <- data.frame(
  run_name = rep(NA, seq_data_length),
  sequence_ID = rep(NA, seq_data_length),
  sequence = rep(NA, seq_data_length),
  counts = rep(NA, seq_data_length),
  counts_run_name = rep(NA, seq_data_length),
  complement = rep(NA, seq_data_length),
  identity = rep(0, seq_data_length),
  identity_subset = rep(0, seq_data_length),
  tag = rep(NA, seq_data_length),
  tag_subset = rep(NA, seq_data_length),
  gap = rep(NA, seq_data_length),
  wobble = rep(NA, seq_data_length),
  location = rep(NA, seq_data_length),
  all_locations = rep(NA, seq_data_length),
  all_identities = rep(NA, seq_data_length),
  all_tags = rep(NA, seq_data_length)
)

# initialize loop variable
sequence_ID_last <- "0_0"
wobble_flag <- NA

# loop over each sequence
for (seq_num in 1:seq_data_length) {
  # add run ID
  complement_data$run_name[seq_num] <- seqs_input$run_name[seq_num]
  # add sequence ID
  complement_data$sequence_ID[seq_num] <- seqs_input$sequence_ID[seq_num]
  # add sequence
  complement_data$sequence[seq_num] <- str_remove_all(paste(seq_matrix[seq_num,], collapse=""), "X")
  # add counts
  complement_data$counts[seq_num] <- seqs_input$counts[seq_num]
  # add counts run name
  complement_data$counts_run_name[seq_num] <- seqs_input$counts_run_name[seq_num]
  # check if the current sequence has already been analyzed
  # if so, set the current sequence data to be the same as the previous and stop parsing
  if (seqs_input$sequence_ID[seq_num] == sequence_ID_last) {
    # set the window sequence
    complement_data$complement[seq_num] <- complement_data$complement[seq_num-1]
    # set the percent identity
    complement_data$identity[seq_num] <- complement_data$identity[seq_num-1]
    # set the longest subset window identity
    complement_data$identity_subset[seq_num] <- complement_data$identity_subset[seq_num-1]
    # set the tag
    complement_data$tag[seq_num] <- complement_data$tag[seq_num-1]
    # set the tag subset
    complement_data$tag_subset[seq_num] <- complement_data$tag_subset[seq_num-1]
    # set window gap flag
    complement_data$gap[seq_num] <- complement_data$gap[seq_num-1]
    # set the wobble flag
    complement_data$wobble[seq_num] <- complement_data$wobble[seq_num-1]
    # set the location
    complement_data$location[seq_num] <- complement_data$location[seq_num-1]
    # set all complementary locations
    complement_data$all_locations[seq_num] <- complement_data$all_locations[seq_num-1]
    # set all complementary identities
    complement_data$all_identities[seq_num] <- complement_data$all_identities[seq_num-1]
    # set all complementary tags
    complement_data$all_tags[seq_num] <- complement_data$all_tags[seq_num-1]
    # jump to the end of the loop and stop parsing the current sequence
    next
  }
  # set last sequence_ID
  sequence_ID_last <- seqs_input$sequence_ID[seq_num]
  # loop over each base of the sequence
  for (base_index in 1:sliding_length) {
    # set end index
    end_index <- base_index + (complement_length-1)
    # get the next 7 bases to create the 8bp sliding window
    slide_window <- seq_matrix[seq_num,base_index:end_index]
    # reset flags
    wobble_flag <- "no"
    gap_flag <- "no"
    # initialize subset length variable and mismatch flag
    subset_length <- 0
    subset_length_list <- c(0,0)
    subset_longest <- 0
    subset_second_longest <- 0
    mismatch_flag <- 0
    # loop over each consecutive base of the window
    for (window_index in 1:complement_length) {
      # check if the current base is a match
      # A pairs with U and G pairs with C, but G also pairs with U
      if (rev_overhang[window_index] == "A" & slide_window[window_index] == "T") { # A (U)
        # increment subset length
        subset_length <- subset_length+1
      }else if (rev_overhang[window_index] == "C" & slide_window[window_index] == "G"){ # C (G)
        # increment subset length
        subset_length <- subset_length+1
      }else if (rev_overhang[window_index] == "G" & slide_window[window_index] == "C"){ # G (C)
        # increment subset length
        subset_length <- subset_length+1
      }else if (rev_overhang[window_index] == "G" & slide_window[window_index] == "T"){ # G (U), wobble
        # increment subset length
        subset_length <- subset_length+1
        # flag wobble
        wobble_flag <- "yes"
      }else if (rev_overhang[window_index] == "T" & slide_window[window_index] == "A"){ # U (A)
        # increment subset length
        subset_length <- subset_length+1
      }else if (rev_overhang[window_index] == "T" & slide_window[window_index] == "G"){ # U (G), wobble
        # increment subset length
        subset_length <- subset_length+1
        # flag wobble
        wobble_flag <- "yes"
      }else{ # mismatch
        # flag mismatch
        mismatch_flag <- 1
      }
      # check if end of window or mismatch
      if (window_index == complement_length || mismatch_flag == 1){ 
        # store current subset length
        subset_length_list <- c(subset_length_list, subset_length)
        # reset subset length
        subset_length <- 0
        # reset flag mismatch
        mismatch_flag <- 0
      }
    }
    # sort and retrieve the longest and second longest subset lengths
    subset_lengths <- subset_length_list[subset_length_list >= min_length]
    subset_longest <- ifelse(is.na(subset_lengths[1]), 0, subset_lengths[1])
    subset_second_longest <- ifelse(is.na(subset_lengths[2]), 0, subset_lengths[2])
    # check if the total length is not at least the minimum length
    if (subset_longest < min_length) {
      # jump to the end of the loop and stop parsing the current window
      next
    }
    # check if the first and second longest consecutive matches are not at least 3bp
    if (subset_longest >= min_subset_length & subset_second_longest >= min_subset_length) {
      # set the total length
      subset_total_length <- subset_longest + subset_second_longest
    } else {
      # set the total length
      subset_total_length <- subset_longest
    }
    # set longest window subset identity
    subset_total_identity <- 100*(subset_total_length)/complement_length
    subset_longest_identity <- 100*subset_longest/complement_length
    # check if there is a gap
    if (subset_longest_identity < subset_total_identity) {
      # set the gap flag
      gap_flag <- "yes"
      # set the total and location tags
      total_tag <- paste(subset_longest, subset_second_longest, sep = "_")
    } else {
      # set the total and location tags
      total_tag <- subset_longest
    }
    # set the location tag
    loc_tag <- paste(base_index, end_index, sep = "-")
    # check if the identity of the longest consecutive subset is larger than the previous largest window identity
    if (subset_total_identity >= complement_data$identity[seq_num] & gap_flag == "no") {
      # store the current window sequence as the complement
      complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
      # add percent identity to expected overhang complement
      complement_data$identity[seq_num] <- subset_total_identity
      # add subset percent identity to expected overhang complement
      complement_data$identity_subset[seq_num] <- subset_longest_identity
      # set the tag
      complement_data$tag[seq_num] <- total_tag
      # set the tag subset
      complement_data$tag_subset[seq_num] <- subset_longest
      # flag that the current window does not have a gap
      complement_data$gap[seq_num] <-  gap_flag
      # set the wobble flag
      complement_data$wobble[seq_num] <- wobble_flag
      # set the location
      complement_data$location[seq_num] <- paste(base_index, end_index, sep = "-")
      # update all complementary locations
      complement_data$all_locations[seq_num] <- paste(complement_data$all_locations[seq_num], loc_tag, sep = ";")
      # update all complementary identities
      complement_data$all_identities[seq_num] <- paste(complement_data$all_identities[seq_num], subset_total_identity, sep = ";")
      # set all complementary tags
      complement_data$all_tags[seq_num] <- paste(complement_data$all_tags[seq_num], total_tag, sep = ";")
    } else if (subset_total_identity > complement_data$identity[seq_num]) {
      # store the current window sequence as the complement
      complement_data$complement[seq_num] <- paste(seq_matrix[seq_num,base_index:end_index], collapse="")
      # add percent identity to expected overhang complement
      complement_data$identity[seq_num] <- subset_total_identity
      # add subset percent identity to expected overhang complement
      complement_data$identity_subset[seq_num] <- subset_longest_identity
      # set the tag
      complement_data$tag[seq_num] <- total_tag
      # set the tag subset
      complement_data$tag_subset[seq_num] <- subset_longest
      # flag that the current window does not have a gap
      complement_data$gap[seq_num] <-  gap_flag
      # set the wobble flag
      complement_data$wobble[seq_num] <- wobble_flag
      # set the location
      complement_data$location[seq_num] <- paste(base_index, end_index, sep = "-")
      # update all complementary locations
      complement_data$all_locations[seq_num] <- paste(complement_data$all_locations[seq_num], loc_tag, sep = ";")
      # update all complementary identities
      complement_data$all_identities[seq_num] <- paste(complement_data$all_identities[seq_num], subset_total_identity, sep = ";")
      # set all complementary tags
      complement_data$all_tags[seq_num] <- paste(complement_data$all_tags[seq_num], total_tag, sep = ";")
    }
  }
}

# remove initializing NAs
complement_data$all_locations <- gsub("NA;", "", complement_data$all_locations)
complement_data$all_identities <- gsub("NA;", "", complement_data$all_identities)
complement_data$all_tags <- gsub("NA;", "", complement_data$all_tags)

# export data
write.csv(complement_data, file = paste(out_dir, "/", round_name, "_overhang_data_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
