#!/bin/bash

# 1) First round of filtering, for QUAL, MQ, QD

# Gatk seemed to run better when going through each filtering criteria one at a time and on command line (crashed if I combined them, not sure why)
# need reference genome .dict file in same directory as reference fasta

# Filter out QUAL < 20
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantFiltration -R /home/degreefe/NBW/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW2_SNPS_2M.vcf.gz -O NBW2_SNPS_2M.QUAL.vcf.gz --filter-name "QUALlt20" --filter-expression "QUAL < 20" 

# Switching to SelectVariants function for next filter steps because the VariantFiltration does not exclude NAs
# Filter out MQ < 30
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R /home/degreefe/NBW/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW2_SNPS_2M.QUAL.vcf.gz -O NBW2_SNPS_2M.QUAL.MQ.vcf.gz -select "MQ >= 30.0"

# Filter out QD < 2
java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar SelectVariants -R /home/degreefe/NBW/ref_genome/Northern_bottlenose_whale_051018_shortLabel.fasta -V NBW2_SNPS_2M.QUAL.MQ.vcf.gz -O NBW2_SNPS_2M.QUAL.MQ.QD.vcf.gz -select "QD >= 2.0"

# Clean up files
mv NBW2_SNPS_2M.QUAL.MQ.QD.vcf.gz NBW2_SNPS_2M.filter1.vcf.gz
mv NBW2_SNPS_2M.QUAL.MQ.QD.vcf.gz.tbi NBW2_SNPS_2M.filter1.vcf.gz.tbi
rm NBW2_SNPS_2M.QUAL*


# 2) Second round of filtering, for max missingness, non-biallelic sites

# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing)
vcftools --gzvcf NBW2_SNPS_2M.filter1.vcf.gz --max-missing 0.6 --recode --recode-INFO-all --out NBW2_SNPS_2M.filter1.miss

# Remove non-biallelic sites
vcftools --vcf NBW2_SNPS_2M.filter1.miss.recode.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out NBW2_SNPS_2M.filter1.miss.biallel

mv NBW2_SNPS_2M.filter1.miss.biallel.recode.vcf NBW2_SNPS_2M.filter1.miss.biallel.vcf


# Going to filter Minor allele count later, will need to remove autosomes and structural variants next (need to keep MAC for smcpp)
# Filter out minor allele count 2
#vcftools --vcf $snps.vcf --mac 2 --recode --recode-INFO-all --out $snps.mac

