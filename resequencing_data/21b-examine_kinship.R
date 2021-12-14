# Looking at kinship estimates with plink output

library(plinkQC)

directory="C:/Users/eveli/Dropbox/NBW-me/NBW_oct2021_updated/snps_2M/kinship"
files="NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDpruned.42" ##for .genome and .imiss
  
# Quick plot

png("plinkQC_NBW2_SNPS_2M_42.png", w=1000, h=500, res=100)
evaluate_check_relatedness(
  qcdir=directory,
  name=files,
  highIBDTh = 0.1875,
  imissTh = 0.03,
  interactive = FALSE,
  verbose = FALSE
)
dev.off()
 
# Extract values for the related individuals from the .genome file (filtering by pi_hat values)
IBD <- read.table(paste(directory,"/", files,".genome", sep=""), header=T)
IBD_order_PH <- IBD[order(IBD[,10], decreasing=TRUE),]
write.csv(IBD_order_PH, "IBD_NBW2_SNPS_2M_42_PIHAT.csv")

# Make nicer histogram, also remove duplicate so we are just looking at real kins:
# Remove duplicate: Hyam-2016-06-pooled & NBW-2016-06 --remove row123. 
# IBD_2 <- IBD[-123,]

# Order by pi_hat and then add nrow (for pair#s)
IBD_2 <- IBD[order(IBD[,10], decreasing=TRUE),]
IBD_2$pair <- 1:nrow(IBD_2)

library(ggplot2)
ggplot(data=IBD_2, aes(x=PI_HAT))+
  geom_histogram(bins=125, color="steelblue", fill="steelblue")+
  theme_bw()+
  ylab("number of pairs")+
  xlab("pi hat")
ggsave("IBD_NBW2_SNPS_2M_42_geomhistogram.png", width=6,height=4,dpi=300)



# PI_HAT value reference:
# 1st degree relative: 0.5
# 2nd degree relative: 0.25
# 3rd degree relative: 0.125

# Z score references:
# Z0=P(IBD=0)
# Z1=P(IBD=1)
# Z2=P(IBD=2)
# unrelated indiv, Z0=1
# unrelated indiv, Z2=0
# unrelated indiv, ratio <2 ish


# Can also use SNPRelate to look at relatedness (code below from Matt Thorstensen)
# Code to install SNPRelate
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("SNPRelate")
BiocManager::install("gdsfmt")
BiocManager::install("SeqArray")
library(gdsfmt)
library(SNPRelate)

vcf.fn <- "C:/Users/Evelien de Greef/Dropbox/NBW-me/reseq_newsnps/kinship/NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.vcf"
snpgdsVCF2GDS(vcf.fn, "nbw.gds", method="biallelic.only")

snpgdsSummary("nbw.gds")

genofile <- SNPRelate::snpgdsOpen("nbw.gds")

# Getting a list of the individuals in the SNP relate dataset
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Using the LD-pruned SNPs for kinship analysis using the method of moments (MoM)
# Based on the guide provided in the vignette here: https://bioconductor.org/packages/release/bioc/vignettes/SNPRelate/inst/doc/SNPRelateTutorial.html#estimating-ibd-using-plink-method-of-moments-mom
ibd <- snpgdsIBDMoM(genofile, sample.id=NULL, snp.id=NULL,
                    maf=0.05, missing.rate=0.05, num.thread=2, autosome.only = FALSE)
# Make a data.frame. This contains pairwise relatedness measures by the method of moments.
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)

snpgdsClose("nbw.gds")
