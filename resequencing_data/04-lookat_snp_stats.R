# Examining some stats on unfiltered snps/vcf to determine filtering thresholds

library(tidyverse)
library(readr)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps")

# Import unfiltered vcf information (.info file from bcftools)
vcfInfo <- read_table2("NBW_platypus_SNPs.info")
colnames(vcfInfo) <- c("CHROM", "POS", "REF", "ALT", "QUAL", "MQ", "QD", "TC", "FR")

# allele freq, missingness, hwe from vcftools
var_frq <- read_delim("NBW_platypus_SNPs.frq", delim="\t",col_names=c("CHR", "POS", "N_ALLELES", "N_CHR", "A1", "A2"), skip = 1)
var_miss <- read_delim("NBW_platypus_SNPs.lmiss", delim="\t")
ind_miss <- read_delim("NBW_platypus_SNPs.imiss", delim="\t")
hwe <- read_delim("NBW_platypus_SNPs.hwe", delim="\t", col_names=c("CHR", "POS", "OBS.HOM1.HET.HOM2", "E.HOM1.HET.HOM2", "ChiSq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS"), skip=1)

## View distributions. Red line setting filtering threshold
# QUAL
min(vcfInfo$QUAL)
max(vcfInfo$QUAL)
hist(vcfInfo$QUAL, breaks=seq(1,2000000,50), xlim=c(0,3500), main="QUAL")
abline(v=20, col="red")

# MQ
min(vcfInfo$MQ, na.rm=TRUE)
max(vcfInfo$MQ, na.rm=TRUE)
hist(vcfInfo$MQ, main="MQ")
abline(v=30, col="red")

# QD
min(vcfInfo$QD, na.rm=TRUE)
max(vcfInfo$QD, na.rm=TRUE)
summary(vcfInfo$QD)
hist(vcfInfo$QD, breaks=100, main="QD")
abline(v=2, col="red")
hist(vcfInfo$QD, breaks=1000, xlim=c(0,100))
abline(v=2, col="red")

# Minor allele frequency
var_frq$MAF <- var_frq %>% select(A1, A2) %>% apply(1, function(z) min(z))
summary(var_frq$MAF)
hist(var_frq$MAF, breaks=100, main="minor allele frequency")
abline(v=0.041, col="red") #2/49, 49 being the sample size

# Variant missingness
summary(var_miss$F_MISS)
hist(var_miss$F_MISS, breaks=100, main="variant missingness")
abline(v=0.4, col="red") #would be 0.6 for vcftools (number is reversed; i.e. 1=no missing)

# HWE for heterozygosity
summary(-log10(hwe$P_HWE))
hist(-log10(hwe$P_HWE), breaks=100, ylim=c(0,1000000), main="hwe p-value")

# Individual missingness
summary(ind_miss$N_MISS)
bar <- barplot(ind_miss$F_MISS, main="individual missingness", names.arg=ind_miss$INDV, las=2,cex.names=0.5)
abline(h=0.5, col="red")

# loking at individual missingness after filters (and LD pruned) -- just curious 
ind_miss2 <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.vcf.imiss", delim="\t")
summary(ind_miss2$N_MISS)
bar <- barplot(ind_miss2$F_MISS, main="individual missingness postfilter(autosomes, LDpruned)", names.arg=ind_miss2$INDV, las=2,cex.names=0.5)
abline(h=0.5, col="red")

