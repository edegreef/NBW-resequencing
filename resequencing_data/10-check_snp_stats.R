# Examining some stats on unfiltered snps/vcf to determine filtering thresholds

library(tidyverse)
library(readr)

setwd("C:/Users/Evelien de Greef/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/stats")

# Import unfiltered vcf information (.info file from bcftools)
vcfInfo <- read_table2("NBW2_SNPS_2M.info")
colnames(vcfInfo) <- c("CHROM", "POS", "REF", "ALT", "QUAL", "MQ", "QD", "TC", "FR")

# allele freq, missingness, hwe from vcftools
var_frq <- read_delim("NBW2_SNPS_2M.frq", delim="\t",col_names=c("CHR", "POS", "N_ALLELES", "N_CHR", "A1", "A2"), skip = 1)
var_miss <- read_delim("NBW2_SNPS_2M.lmiss", delim="\t")
#var_miss_updatedset <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.MISSINGNESS.lmiss", delim="\t")
ind_miss <- read_delim("NBW2_SNPS_2M.imiss", delim="\t")
#ind_miss_updatedset <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.MISSINGNESS.imiss", delim="\t")
hwe <- read_delim("NBW2_SNPS_2M.hwe", delim="\t", col_names=c("CHR", "POS", "OBS.HOM1.HET.HOM2", "E.HOM1.HET.HOM2", "ChiSq_HWE", "P_HWE", "P_HET_DEFICIT", "P_HET_EXCESS"), skip=1)

## View distributions. Red line setting filtering threshold
# QUAL
min(vcfInfo$QUAL)
max(vcfInfo$QUAL)
library(ggplot2)
qual <- ggplot(vcfInfo, aes(x=QUAL))+
  geom_histogram(fill="gray60", bins=50)+
  theme_classic()+
  ggtitle("Quality Score")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  geom_vline(xintercept=20, col="red")+
  ylab("Count")
  
qual


QUAL <- hist(vcfInfo$QUAL, breaks=seq(1,2000000,50), xlim=c(0,3100), main="QUAL", xlab="") 
abline(v=20, col="red")

# MQ
min(vcfInfo$MQ, na.rm=TRUE)
max(vcfInfo$MQ, na.rm=TRUE)
hist(vcfInfo$MQ, main="MQ")
abline(v=30, col="red")

mq <- ggplot(vcfInfo, aes(x=MQ))+
  geom_histogram(fill="gray60", bins=40)+
  theme_classic()+
  ggtitle("Mapping Quality Score")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  geom_vline(xintercept=30, col="red")+
  ylab("Count")

mq


# QD
min(vcfInfo$QD, na.rm=TRUE)
max(vcfInfo$QD, na.rm=TRUE)
summary(vcfInfo$QD)
hist(vcfInfo$QD, breaks=100, main="QD")
abline(v=2, col="red")
hist(vcfInfo$QD, breaks=1000, xlim=c(0,100))
abline(v=2, col="red")

qd <- ggplot(vcfInfo, aes(x=QD))+
  geom_histogram(fill="gray60", bins=100)+
  theme_classic()+
  ggtitle("Quality By Depth")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  geom_vline(xintercept=2, col="red")+
  ylab("Count")

qd

# Minor allele frequency
var_frq$MAF <- var_frq %>% select(A1, A2) %>% apply(1, function(z) min(z))
summary(var_frq$MAF)
hist(var_frq$MAF, breaks=100, main="minor allele frequency")
abline(v=0.041, col="red") #2/49, 49 being the sample size

maf <- ggplot(var_frq, aes(x=MAF))+
  geom_histogram(fill="gray60", bins=80)+
  theme_classic()+
  ggtitle("Minor Allele Frequency")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  geom_vline(xintercept=0.041, col="red")+
  ylab("Count")
  
maf

# Variant missingness
summary(var_miss$F_MISS)
hist(var_miss$F_MISS, breaks=100, main="variant missingness")
abline(v=0.4, col="red") #would be 0.6 for vcftools (number is reversed; i.e. 1=no missing)

variant_miss <- ggplot(var_miss, aes(x=F_MISS))+
  geom_histogram(fill="gray60", bins=50)+
  theme_classic()+
  ggtitle("Variant Missingness")+
  theme(plot.title=element_text(hjust=0.5, face="bold"))+
  geom_vline(xintercept=0.4, col="red")+
  xlab("Missingness")+
  ylab("Count")
variant_miss

# HWE for heterozygosity
summary(-log10(hwe$P_HWE))
hist(-log10(hwe$P_HWE), breaks=100, ylim=c(0,1000000), main="hwe p-value")

# Individual missingness
summary(ind_miss$N_MISS)
bar <- barplot(ind_miss$F_MISS, main="Individual missingness", names.arg=ind_miss$INDV, las=2,cex.names=0.5)
abline(h=0.5, col="red")

indiv_miss <- ggplot(data=ind_miss, aes(x=factor(INDV), y=F_MISS))+
  geom_bar(stat="identity", fill="gray60")+
  theme_classic()+
  xlab("Individual")+
  ylab("Missingness")+
  ggtitle("Individual Sample Missingness")+
  theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_blank())+
  geom_hline(yintercept=0.4, col="red")
#scale_fill_manual( values = c( "no"="gray", "yes"="orange" ), guide = FALSE )
indiv_miss

library(patchwork)
(qual | mq ) / (qd | maf ) / (variant_miss | indiv_miss)
ggsave("NBW2_SNPS_2M_prefilter.png", width = 9, height = 9, dpi = 300)

# loking at individual missingness after filters (and LD pruned) -- just curious 
#ind_miss2 <- read_delim("NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.vcf.imiss", delim="\t")
#summary(ind_miss2$N_MISS)
#bar <- barplot(ind_miss2$F_MISS, main="individual missingness postfilter(autosomes, LDpruned)", names.arg=ind_miss2$INDV, las=2,cex.names=0.5)
#abline(h=0.5, col="red")

#
#variant_miss_updated <- ggplot(var_miss_updatedset, aes(x=F_MISS))+
#  geom_histogram(fill="gray60", bins=38)+
#  theme_classic()+
#  ggtitle("Variant Missingness (post-filter)")+
#  theme(plot.title=element_text(hjust=0.5, face="bold"))+
#  ylab("Count")+
#  xlab("Missingness")+
#  xlim(0,1)
#variant_miss_updated

#indiv_miss_updated <- ggplot(data=ind_miss_updatedset, aes(x=factor(INDV), y=F_MISS))+
#  geom_bar(stat="identity", fill="gray60")+
#  theme_classic()+
#  xlab("individual")+
#  ylab("Missingness")+
#  ggtitle("Individual Sample Missingness (post-filter)")+
#  theme(plot.title=element_text(hjust=0.5, face="bold"),axis.text.x=element_blank())+
#  ylim(0,1)
#indiv_miss_updated

#variant_miss_updated + indiv_miss_updated
#ggsave("platypus_snps_postfilter_plots.png", width = 9, height = 3, dpi = 300)

#sum(var_miss_updatedset$N_MISS) / sum(var_miss_updatedset$N_DATA)

