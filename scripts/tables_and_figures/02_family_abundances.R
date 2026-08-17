#!/usr/bin/env Rscript

# R script to create family abundance plots
# usage: 02_F2A_F2B_F2C_family_abundances.R

# turn of scientific notation
options(scipen=10000)

# inport libraries
library(ggplot2)
library(scales)
library(rcartocolor)
library(dplyr)
#library(plyr)

# color blind safe plotting palette
safe_colors <- c(carto_pal(name="Safe"), "#000000", "#D55E00", "#0072B2", "#F0E442", "#E69F00", "#009E73", "#999999")

# set outputs directory
out_dir <- "/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/F2A_F2B_F2C_family_abundances_above2"

# create outputs directory
dir.create(out_dir, showWarnings = FALSE)

# read in ligation rates
ligation_rates <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/ligation_rates_table.csv")

# numbers of high quality reads
#quality_doped <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374, 865509, 807849, 1143871)
quality <- c(1039660, 1067585, 1033048, 866423, 981844, 916485, 582260, 889374)
#unique <- c(1036229, 1063996, 1029483, 863123, 966495, 500507, 92366, 108529)

# read in family data
fam_data <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/tables_and_figures/ST2_family_table/r8_family_count_data.csv")

# read in sequences that have at least 90% identity to any peak
seqs_identities <- read.csv("/Users/bamflappy/PfrenderLab/RNA_evolution/outputs/11b_family_identification_above2/family_identities_max_atLeast90.csv")

# data exploration
r8_id <- seqs_identities[seqs_identities$run_name == 8,]
r8_id_similar <- r8_id[r8_id$peak_identity >= 90,]
r8_id_dissimilar <- r8_id[r8_id$peak_identity < 90,]
write.csv(r8_id_similar, file = paste(out_dir, "/round8_family_sequence_data.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(r8_id_dissimilar, file = paste(out_dir, "/round8_orphan_sequence_data.csv", sep = ""), row.names = FALSE, quote = FALSE)

# list of cluster and run IDs
cluster_list <- unique(seqs_identities$peak_cluster_ID)
run_list <- unique(seqs_identities$run_name)

# set the data length
data_length <- length(cluster_list)*length(run_list)
  
# create cluster data frame
cluster_data <- data.frame(
  run_name = rep(NA, data_length),
  cluster_ID = rep(NA, data_length),
  family_ID = rep(NA, data_length),
  seq_counts = rep(NA, data_length),
  read_counts = rep(NA, data_length),
  read_abun = rep(NA, data_length),
  sequence = rep(NA, data_length),
  cluster_color = rep(NA, data_length)
)

# initialize the index
index <- 1

# loop over each run
for (run_num in min(run_list):max(run_list)) {
  # loop over each cluster
  for (cluster_num in min(cluster_list):max(cluster_list)) {
    # set the current run, cluster, and family IDs
    cluster_data$run_name[index] <- run_num
    cluster_data$cluster_ID[index] <- cluster_num
    cluster_data$family_ID[index] <- fam_data[fam_data$cluster_ID == cluster_num, "family_ID"]
    # add sequence counts
    cluster_data$seq_counts[index] <-  nrow(seqs_identities[seqs_identities$peak_cluster_ID == cluster_num & seqs_identities$run_name == run_num,])
    # add read counts
    cluster_data$read_counts[index] <-  sum(seqs_identities[seqs_identities$peak_cluster_ID == cluster_num & seqs_identities$run_name == run_num,"read_counts"])
    # add read abundances
    cluster_data$read_abun[index] <- 100*cluster_data$read_counts[index]/quality[run_num]
    # add sequence
    cluster_data$sequence[index] <- fam_data[fam_data$cluster_ID == cluster_num, "sequence"]
    # add cluster color
    cluster_data$cluster_color[index] <- fam_data[fam_data$cluster_ID == cluster_num, "cluster_color"]
    # increment index
    index <- index+1
  }
}

# sort cluster data
cluster_data <- cluster_data[order(cluster_data$family_ID, decreasing = FALSE),]  

# setup ligation rates and diversity data for plotting
ligation_rates_data <- cbind(ligation_rates[ligation_rates$stat == "rate",], ligation_rates[ligation_rates$stat == "error", "value"])
ligation_rates_data <- ligation_rates_data[ , !(names(ligation_rates_data) %in% c("stat", "stat_color"))]
colnames(ligation_rates_data) <- c("run_name","rate","error")
diversity_data <- ligation_rates[ligation_rates$stat == "diversity",]
ligation_driversity_data <- merge(ligation_rates_data, diversity_data)
colnames(ligation_driversity_data) <- c("run_name","rate","error","diversity","stat","stat_color")
ligation_driversity_data <- ligation_driversity_data[,c("run_name","rate","error","diversity")]
#coeff <- (max(ligation_driversity_data$rate)+max(ligation_driversity_data$error))/max(ligation_driversity_data$diversity)
coeff <- 0.05/max(ligation_driversity_data$diversity)

# setup axis title
axis_title <- bquote(italic('k')[obs](h^-1))

# combined line plot of ligation rates with diversity
ligation_rates_diversity <- ggplot(ligation_driversity_data, aes(run_name)) +
  geom_line(aes(y = diversity, linetype = "Diversity"), size = 1.25, color = safe_colors[15]) +
  geom_point(aes(y = diversity), size = 2.75, color = safe_colors[15]) +
  geom_line(aes(y = rate/coeff, linetype = "Activity"), size = 1.25, color = safe_colors[2]) + 
  geom_point(aes(y = rate/coeff), size = 2.75, color = safe_colors[2]) +
  geom_errorbar(aes(ymin=(rate-error)/coeff, ymax=(rate+error)/coeff), width=.2,
                position=position_dodge(0.05), color = safe_colors[2]) +
  theme_classic(base_size = 16) +
  guides(y = guide_axis(cap = "upper")) +
  scale_y_continuous(
    name = "Percent Unique", breaks=seq(0, 100, 20),# labels = function(x) paste0(x, "%"),
    sec.axis = sec_axis(~.*coeff, name=axis_title, guide = guide_axis(cap = "upper"))
  ) +
  scale_x_continuous("Round", labels = as.character(ligation_driversity_data$run_name), breaks = ligation_driversity_data$run_name) +
  xlab("Round") +
  scale_linetype_manual(values=c("solid", "solid"), name = "Statistic")
# save the plot
exportFile <- paste(out_dir, "/sequence_diversity_ligation_rates.png", sep = "")
png(exportFile, units="in", width=6, height=4, res=300)
print(ligation_rates_diversity)
dev.off()

# line plot with percent abundance per round for each of the families
cluster_abundances_plot <- ggplot(data=cluster_data, aes(x=as.character(run_name), y=read_abun, group=family_ID, color=cluster_color))+
  geom_line(size = 1.25) +
  geom_point(size = 2.75) +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Family", labels = cluster_data$family_ID, breaks = cluster_data$cluster_color, guide = "legend") +
  scale_y_continuous(limits=c(0, 40), breaks=seq(0, 40, 5)) +#, labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  theme(axis.line = element_line()) +
  ylab("Percent Abundance") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/family_percent_abundances.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(cluster_abundances_plot)
dev.off()

# retrieve peak counts and abundances
peaks_data <- seqs_identities[seqs_identities$sequence %in% fam_data$sequence,] 

# loop over each run
for (run_num in 1:8) {
  # add read abundances
  peaks_data[peaks_data$run_name == run_num, "read_abun"] <- 100*peaks_data[peaks_data$run_name == run_num, "counts"]/quality[run_num]
}
#peaks_data$read_abun <- as.numeric(peaks_data$read_abun)
peaks_data <- peaks_data[,c("run_name","read_counts","read_abun","sequence")]

# rename columns
colnames(peaks_data) <- c("run_name","peak_counts","peak_abun","sequence")

# add peaks counts and abundances
peaks_data <- merge(fam_data, peaks_data)

# sort cluster data
peaks_data <- peaks_data[order(peaks_data$family_ID, decreasing = FALSE),]  

# create table of peaks data
peaks_table <- peaks_data[peaks_data$sequence %in% fam_data$sequence & peaks_data$run_name == 8,] 

# calculate average identity within families
peaks_table$avg_identity <- NA
peaks_table$min_identity <- NA
peaks_table$max_identity <- NA
for (clust_num in 0:9) {
  peaks_table$avg_identity[clust_num+1] <- sum(seqs_identities[seqs_identities$peak_cluster_ID == clust_num & seqs_identities$run_name == 8, "peak_identity"])/length(seqs_identities[seqs_identities$peak_cluster_ID == clust_num & seqs_identities$run_name == 8, "peak_identity"])
  peaks_table$min_identity[clust_num+1] <- min(seqs_identities[seqs_identities$peak_cluster_ID == clust_num & seqs_identities$run_name == 8, "peak_identity"])
  peaks_table$max_identity[clust_num+1] <- max(seqs_identities[seqs_identities$peak_cluster_ID == clust_num & seqs_identities$run_name == 8, "peak_identity"])
}

# create table with identity data
identity_table <- data.frame(
  run_num = rep(NA, 80),
  cluster_ID = rep(NA, 80),
  family_ID = rep(NA, 80),
  family_color = rep(NA, 80),
  avg_identity = rep(NA, 80),
  min_identity = rep(NA, 80),
  max_identity = rep(NA, 80)
)

# calculate average identity within families for each round
index <- 1
for (clust in 0:9) {
  for (run in 1:8) {
    # add the data
    identity_table$run_num[index] <- run
    identity_table$cluster_ID[index] <- clust
    identity_table$family_ID[index] <- peaks_table[peaks_table$cluster_ID == clust, "family_ID"]
    identity_table$family_color[index] <- peaks_table[peaks_table$cluster_ID == clust, "cluster_color"]
    identity_table$avg_identity[index] <- sum(seqs_identities[seqs_identities$peak_cluster_ID == clust & seqs_identities$run_name == run, "peak_identity"])/length(seqs_identities[seqs_identities$peak_cluster_ID == clust & seqs_identities$run_name == run, "peak_identity"])
    identity_table$min_identity[index] <- min(seqs_identities[seqs_identities$peak_cluster_ID == clust & seqs_identities$run_name == run, "peak_identity"])
    identity_table$max_identity[index] <- max(seqs_identities[seqs_identities$peak_cluster_ID == clust & seqs_identities$run_name == run, "peak_identity"])
    # increment the index counter
    index <- index + 1
  }
}
is.na(identity_table)<-sapply(identity_table, is.infinite)
identity_table[is.na(identity_table)]<-0

# set zeros to NA
identity_table_subset <- filter(identity_table, avg_identity > 0, min_identity > 0, max_identity > 0)

# sort cluster data
identity_table_subset <- identity_table_subset[order(identity_table_subset$family_ID, decreasing = FALSE),]  

# line plot with average identity per round for each of the families
family_identities_subset_plot <- ggplot(data=identity_table_subset, aes(x=as.character(run_num), y=100-avg_identity, group=family_ID, color=family_color))+
  geom_line(size = 1.25) +
  geom_point(size = 2.75) +
  geom_point() +
  theme_classic(base_size = 16) +
  scale_color_identity(name = "Family", labels = identity_table_subset$family_ID, breaks = identity_table_subset$family_color, guide = "legend") +
  #scale_y_continuous(labels = function(x) paste0(x, "%")) +
  guides(y = guide_axis(cap = "upper"), x = guide_axis(cap = "upper")) +
  ylab("Percent Unique") +
  xlab("Round")
# save the plot
exportFile <- paste(out_dir, "/persistent_family_avg_identities.png", sep = "")
png(exportFile, units="in", width=5, height=4, res=300)
print(family_identities_subset_plot)
dev.off()

# export plotting data
write.csv(cluster_data, file = paste(out_dir, "/family_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(peaks_data, file = paste(out_dir, "/peaks_abundances.csv", sep = ""), row.names = FALSE, quote = FALSE)
write.csv(peaks_table, file = paste(out_dir, "/peaks_table.csv", sep = ""), row.names = FALSE, quote = FALSE)
