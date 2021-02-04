#!/bin/bash

# 1) First round of filtering, for QUAL, MQ, QD

# Gatk seemed to run better when going through each filtering criteria one at a time and on command line (crashed if I combined them, not sure why)

# filtering out QUAL < 20
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -R /home/degreefe/NBW/reference/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW_platypus_SNPs.vcf.gz -O NBW_platypus_SNPs.QUAL.vcf.gz --filter-name "QUALlt20" --filter-expression "QUAL < 20" 

# switching to SelectVariants function for next filter steps because the VariantFiltration does not exclude NAs
# fitering out MQ < 30
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R /home/degreefe/NBW/reference/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW_platypus_SNPs.QUAL.vcf.gz -O NBW_platypus_SNPs.QUAL.MQ.vcf.gz -select "MQ >= 30.0"

# filtering out QD < 2
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R /home/degreefe/NBW/reference/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW_platypus_SNPs.QUAL.MQ.vcf.gz -O NBW_platypus_SNPs.QUAL.MQ.QD.vcf.gz -select "QD >= 2.0"

# clean up files
mv NBW_platypus_SNPs.QUAL.MQ.QD.vcf.gz NBW_platypus_SNPs.filter1.vcf.gz
mv NBW_platypus_SNPs.QUAL.MQ.QD.vcf.gz.tbi NBW_platypus_SNPs.filter1.vcf.gz.tbi
rm NBW_platypus_SNPs.QUAL*


# 2) Second round of filtering, for minor allele count, max missingness, and non-biallelic sites

# minor allele count
vcftools --gzvcf NBW_platypus_SNPs.filter1.vcf.gz --mac 2 --recode --recode-INFO-all --out NBW_platypus_SNPs.filter1.mac

# missingness (in vcftools 1=no missing, 0=all missing)
vcftools --vcf NBW_platypus_SNPs.filter1.mac.recode.vcf --max-missing 0.6 --recode --recode-INFO-all --out NBW_platypus_SNPs.filter1.mac.miss

# remove non-biallelic sites
vcftools --vcf NBW_platypus_SNPs.filter1.mac.miss.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out NBW_platypus_SNPs.filter1.mac.miss.biallel

mv NBW_platypus_SNPs.filter1.mac.miss.biallel.recode.vcf NBW_platypus_SNPs.filter1.filter2.vcf
rm NBW_platypus_SNPs.filter1.mac*
