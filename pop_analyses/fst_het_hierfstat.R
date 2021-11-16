# Looking at diversity index, using help from Matt Thorstensen's script

library(hierfstat)
library(vcfR)
library(tidyverse)
library(adegenet)

setwd("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/diversity")

# Reading in vcf
vcf <- read.vcfR("NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.vcf") 

# Reformatting the vcfR object into a genlight object
vcf <- vcfR2genlight(vcf)

# Reformatting this to a matrix for input into hierfstat
mat <- as.matrix(vcf)

# Prepare 0, 1, 2 matrix in hierfstat format
# Use locations to test for population differentiation
snps <- mat
snps[snps==0] <- 11
snps[snps==1] <- 12
snps[snps==2] <- 22

# Using partial data as test run
#subset <- snps[,1:1000]

df <- rownames_to_column(as.data.frame(snps), var="Ind")

# Bring in sample info
sample_info <- read.csv("C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/Hyperoodon_popinfo4_reseq_Oct2021.csv", header=T)
sample_info_subset <- subset(sample_info, remove!="Y") #removing individuals with high miss and kin pair

# Extract pop info
pops <- tibble(ID=sample_info_subset$Genome_ID, region=sample_info_subset$region)
colnames(pops)[1] <- "Ind"

# Make pops info numeric (automatically goes in alphabetical order)
# 1=Arctic, 2=Iceland, 3=Labrador, 4=Newfoundland, 5=Scotian_shelf
pops$region<- as.numeric(as.factor(pops$region))

# Left_join info to snps
df <- left_join(pops, df) %>%
  column_to_rownames(var="Ind")

save(df, file="NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.37.RData")

# Look at pairwise Fst with hierfstat
#wc <- pairwise.WCfst(df, diploid = TRUE)
#wc
#write.csv(wc,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.csv")

# Estimate confidence intervals
#wc_bootstrp_confidence <- boot.ppfst(dat = df, nboot = 1000)
#wc_bootstrp_confidence$ll
#wc_bootstrp_confidence$ul

# save outputs
#write.csv(wc_bootstrp_confidence$ll,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ll.csv")
#write.csv(wc_bootstrp_confidence$ul,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ul.csv")
#save(wc_bootstrp_confidence,file="NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.wc_bootstrp_CI.RData")

# Look at heterozygosity stats
# Re-sort the data 
df <- df[order(df$region),]
df[,1:3] # check

# Basic stats for all pops overall
basic_stats <- basic.stats(data = df)
basic_stats

# Save environment:
save.image(file = "overall_basic.stats.RData")

# Basic stats for population 1 
basic_1 <- basic.stats(df[df$region==1,], diploid = TRUE)
basic_1$overall

# Basic stats for population 2 
basic_2 <- basic.stats(df[df$region==2,], diploid = TRUE)
basic_2$overall

# Basic stats for population 3 
basic_3 <- basic.stats(df[df$region==3,], diploid = TRUE)
basic_3$overall

# Basic stats for population 4 
basic_4 <- basic.stats(df[df$region==4,], diploid = TRUE)
basic_4$overall

# Basic stats for population 5 
basic_5 <- basic.stats(df[df$region==5,], diploid = TRUE)
basic_5$overall

# Save environment:
save.image(file = "population_level_basic.stats.RData")
