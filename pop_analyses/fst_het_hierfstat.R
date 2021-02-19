# Looking at FST & heterozygosity between populations. Help from Matt Thorstensen to for code to run Hierfstat.
# 2) Making fst heatmap
# 3) Examining isolation-by-distance pattern & running mantel test


# 1) estimate FST 
library(hierfstat)
library(vcfR)
library(tidyverse)
library(adegenet)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/fst")

# Read in vcf
vcf <- read.vcfR("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.vcf") 

# Reformat vcfR object into a genlight object
vcf <- vcfR2genlight(vcf)

# Reformat to a matrix for input into hierfstat
mat <- as.matrix(vcf)

# Prepare 0, 1, 2 matrix in hierfstat format
snps <- mat
snps[snps==0] <- 11
snps[snps==1] <- 12
snps[snps==2] <- 22

# using partial data as test run
#subset <- snps[,1:1000]

df <- rownames_to_column(as.data.frame(snps), var="Ind")

# bring in sample info
sample_info <- read.csv("C:/Users/Evelien de Greef/Dropbox/NBW-me/sample_info.csv", header=T)
sample_info_subset <- subset(sample_info, remove_indiv!="Y") #removing individuals with high miss and kin pair

# extract pop info
pops <- tibble(ID=sample_info_subset$sample, region=sample_info_subset$region)
colnames(pops)[1] <- "Ind"

# make pops info numeric (automatically goes in alphabetical order)
# 1=Arctic, 2=Iceland, 3=Labrador, 4=Newfoundland, 5=Scotian_shelf
pops$region<- as.numeric(as.factor(pops$region))

# Left_join info to snps
df <- left_join(pops, df) %>%
  column_to_rownames(var="Ind")

save(df, file="NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.df.leftjoined.RData")

# Look at pairwise Fst with hierfstat
wc <- pairwise.WCfst(df, diploid = TRUE)
wc
write.csv(wc,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.csv")

# Estimate confidence intervals
wc_bootstrp_confidence <- boot.ppfst(dat = df, nboot = 1000)
wc_bootstrp_confidence$ll
wc_bootstrp_confidence$ul

# save outputs
write.csv(wc_bootstrp_confidence$ll,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ll.csv")
write.csv(wc_bootstrp_confidence$ul,"NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.FST.ul.csv")
save(wc_bootstrp_confidence,file="NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.wc_bootstrp_CI.RData")

# look at heterozygosity stats
# re-sort the data 
df <- df[order(df$region),]
df[,1:3] # check

# basic stats for all pops overall
basic_stats <- basic.stats(data = df)
basic_stats

# save environment:
save.image(file = "overall_basic.stats.RData")

# basic stats for population 1 
basic_1 <- basic.stats(df[df$region==1,], diploid = TRUE)
basic_1$overall

# basic stats for population 2 
basic_2 <- basic.stats(df[df$region==2,], diploid = TRUE)
basic_2$overall

# basic stats for population 3 
basic_3 <- basic.stats(df[df$region==3,], diploid = TRUE)
basic_3$overall

# basic stats for population 4 
basic_4 <- basic.stats(df[df$region==4,], diploid = TRUE)
basic_4$overall

# basic stats for population 5 
basic_5 <- basic.stats(df[df$region==5,], diploid = TRUE)
basic_5$overall

# save environment:
save.image(file = "population_level_basic.stats.RData")
