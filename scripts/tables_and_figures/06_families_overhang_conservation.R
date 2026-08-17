#!/usr/bin/env Rscript

# created by: Elizabeth Brooks

# R script to create overhang conservation plots
# usage: 07_family_overhang_conservation.R

# turn of scientific notation
options(scipen=10000)

# import libraries
library(ggplot2)
library(rcartocolor)
library(stringr)
library(dplyr)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), palette.colors(palette = "Okabe-Ito"))

# set the number of rounds
num_rounds <- 8

# numbers of high quality reads
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique_reads <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)
#above_two_reads <- c(18, 19, 26, 27, 1585, 10626, 7230, 6315)

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/families_overhang_conservation_above2"
dir.create(out_dir, showWarnings = FALSE)

# read in identity data
identityFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/13b_overhang_conservation_above2/overhang_data_wobble.csv"
complement_data_rounds <- read.csv(identityFile)

# read in family data
familyFile <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F4A_family_base_conservation_above2_r8_unique/family_sequences.csv"
family_data <- read.csv(familyFile)

# add complement data
complement_data_family <- merge(complement_data_rounds, family_data, by = "sequence")

# remove unnecessary columns
complement_data_family <- subset(complement_data_family, select = -c(run_name.y, sequence_ID.y, counts.y, counts_run_name.y))

# re-name columns
names(complement_data_family)[names(complement_data_family) %in% c("run_name.x", "sequence_ID.x", "counts.x", "counts_run_name.x")] <- c("run_name", "sequence_ID", "counts", "counts_run_name")

# replace NAs with zeros
complement_data_family[is.na(complement_data_family)] <- 0
#complement_data_rounds[is.na(complement_data_rounds)] <- 0

# subset to round 8 data
complement_data_family <- complement_data_family[complement_data_family$run_name == 8,]

# copy data frame for later plotting
complement_data <- complement_data_family
#complement_data <- complement_data_rounds

# replace 0's through 4's with <5
complement_data[complement_data$tag == 0, "tag"] <- "<5"
complement_data[complement_data$tag == 2, "tag"] <- "<5"
complement_data[complement_data$tag == 3, "tag"] <- "<5"
complement_data[complement_data$tag == 4, "tag"] <- "<5"

# vectors of bins (total, consecutive, gaped)
#tag_bins <- unique(complement_data$tag)
tag_bins <- c("<5", "5", "3_3", "6", "4_3", "3_4", "7", "8")

# set data lengths
data_length <- length(tag_bins)

# list of family IDs
fam_list <- unique(complement_data$family_ID)

# initialize data frame for identity bin counts
complement_counts_sorted <- data.frame(
  family_ID = rep(NA, data_length),
  tag = rep(NA, data_length),
  counts = rep(NA, data_length),
  counts_unique = rep(NA, data_length),
  frac_abundance = rep(NA, data_length),
  frac_abundance_unique = rep(NA, data_length)
)
complement_counts_out <- data.frame()

# loop over each family
for (fam_num in min(fam_list):max(fam_list)) {
  # loop over tag bins
  for (bin_index in 1:data_length) {
    # retrieve current tag
    cur_tag <- tag_bins[bin_index]
    # add family ID
    complement_counts_sorted$family_ID[bin_index] <- fam_num
    # add tag
    complement_counts_sorted$tag[bin_index] <- cur_tag
    # set the number of unique sequences
    seq_data_length <- nrow(complement_data[complement_data$family_ID == fam_num,])
    # add overhang complement counts
    complement_counts_sorted$counts[bin_index] <- sum(complement_data[complement_data$tag == cur_tag & complement_data$family_ID == fam_num, "counts"])
    # add overhang complement unique counts
    complement_counts_sorted$counts_unique[bin_index] <- nrow(complement_data[complement_data$tag == cur_tag & complement_data$family_ID == fam_num,])
    # add fraction abundance
    complement_counts_sorted$frac_abundance_unique[bin_index] <- complement_counts_sorted$counts_unique[bin_index]/seq_data_length
    complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/quality[fam_num]
    #complement_counts_sorted$frac_abundance[bin_index] <- complement_counts_sorted$counts[bin_index]/above_two_reads[fam_num]
  }
  # add current run data
  complement_counts_out <- rbind(complement_counts_out, complement_counts_sorted)
}

# add percent counts
complement_counts_out$perc_abundance_unique <- 100*complement_counts_out$frac_abundance_unique
complement_counts_out$perc_abundance <- 100*complement_counts_out$frac_abundance

# export data
write.csv(complement_counts_out, file = paste(out_dir, "/", "families_overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)

# get lists of tags
#tag_list <- unique(complement_counts_out$tag)
tag_list <- c("8", "7", "4_3", "3_4", "6", "3_3", "5", "<5")

# add mapping table
identity_mappings <- data.frame(
  tag = tag_list,
  colors = safe_colors[1:8]
)

# add placeholder columns
complement_counts_out$colors <- NA

# loop over each possible tag label
for (label_num in 1:nrow(identity_mappings)) {
  # set tag colors for plotting
  ifelse(
    complement_counts_out$tag == identity_mappings$tag[label_num], 
    complement_counts_out[complement_counts_out$tag == identity_mappings$tag[label_num], "colors"] <- identity_mappings$colors[label_num], 
    NA
  )
}

# setup data for plotting
complement_counts_out$tag <- as.factor(complement_counts_out$tag)
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 3,]
#complement_counts_out <- complement_counts_out[complement_counts_out$tag != 4,]

# sort the data for plotting
complement_counts_total <-  complement_counts_out %>% arrange(factor(tag, levels = tag_list))
complement_counts_total$tag <- factor(complement_counts_total$tag, levels = c("8", "7", "4_3", "3_4", "6", "3_3", "5", "<5"))

# set order of family IDs for plotting
desired_order <- seq(1, 10)

## plots using sequencing read counts

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=tag, y=perc_abundance_unique, x=as.character(family_ID))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  scale_fill_manual(breaks = unique(complement_counts_total$tag), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  scale_x_discrete(limits = desired_order) +
  labs(fill = "Complementarity") +
  ylab("Proportion") +
  xlab("Family")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(base_counts_plot)
dev.off()

# create bar plot of total overhang identity percent
base_counts_plot <- ggplot(complement_counts_total, aes(fill=tag, y=perc_abundance_unique, x=as.character(family_ID))) + 
  geom_bar(position="stack", stat="identity") +
  theme_classic(base_size = 16) +
  theme(panel.background = element_rect(fill = 'black')) +
  scale_fill_manual(breaks = unique(complement_counts_total$tag), values = unique(complement_counts_total$colors), labels = unique(complement_counts_total$tag)) +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper")) +#, x = guide_axis(cap = "upper")) +
  geom_text(aes(label=paste0(sprintf("%1.1f", perc_abundance_unique),"%")),
            position=position_stack(vjust=0.5), size = 3, color = "white") +
  scale_x_discrete(limits = desired_order) +
  labs(fill = "Complementarity") +
  ylab("Proportion") +
  xlab("Family")
# save the plot
exportFile <- paste(out_dir, "/overhang_percent_abundance_unique_total_chart_perc.png", sep = "")
png(exportFile, units="in", width=6, height=4, res=300)
print(base_counts_plot)
dev.off()

# export data
write.csv(complement_counts_sorted, file = paste(out_dir, "/overhang_conservation_wobble.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(complement_counts_total, file = paste(out_dir, "/overhang_conservation_wobble_total.csv", sep = ""), row.names = FALSE, quote = FALSE)
