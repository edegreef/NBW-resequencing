#!/bin/bash

# looking at snp statistics

all_variants=/home/degreefe/NBW/snps_withmerged/NBW2_allvariantcalls_2mergedbam
snps=/home/degreefe/NBW/snps_withmerged/NBW2_SNPS_2M

# filter out indels 
vcftools --vcf $all_variants.vcf --remove-indels --recode --recode-INFO-all --out $snps
mv $snps.recode.vcf $snps.vcf

# prep vcf to look at quality of snps
bgzip $snps.vcf
tabix -p vcf $snps.vcf.gz

# getting snp info using bcftools 
/home/degreefe/programs/bcftools-1.9/bcftools query -H -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/MQ\t%INFO/QD\t%INFO/TC\t%INFO/FR\n' $snps.vcf.gz -o $snps.info

# allele frequency (outputs .frq file)
vcftools --gzvcf $snps.vcf.gz --out $snps --freq2 --max-alleles 2

# proportion of missing data per individual (outputs .imiss file)
vcftools --gzvcf $snps.vcf.gz --out $snps --missing-indv 

# propotion of missing data per site (outputs .lmiss file)
vcftools --gzvcf $snps.vcf.gz --out $snps --missing-site

# assess sites for Hardy-Weinberg Equilibrium (outputs .hwe file)
vcftools --gzvcf $snps.vcf.gz --out $snps --hardy
