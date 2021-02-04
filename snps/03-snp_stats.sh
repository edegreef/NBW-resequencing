#!/bin/bash

# looking at snp stats

all_variants=/home/degreefe/NBW/reseq_newsnps/NBW_variantcalls_platypus
snps=/home/degreefe/NBW/reseq_newsnps/NBW_platypus_SNPs

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


## some notes below for GATK snp info, but ended up using outputs from bcftools and vcftools
#reference=/home/degreefe/NBW/reference/Northern_bottlenose_whale_051018_shortLabel.fasta

# make .dict and .fai file for reference genome
#java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar CreateSequenceDictionary -R $reference
#samtools faidx $reference

# make index (.tbi) file for vcf
#java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar IndexFeatureFile -I $snps.vcf.gz

# create variants table of snp info prior to filtering
#java -jar /home/degreefe/programs/gatk-4.1.9.0/gatk-package-4.1.9.0-local.jar VariantsToTable -V $snps.vcf.gz -F CHROM -F POS -F QUAL -F MQ -F QD -F TC -F FR --show-filtered -O $snps.tab