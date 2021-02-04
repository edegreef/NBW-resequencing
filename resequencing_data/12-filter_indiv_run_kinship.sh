#!/bin/bash

# remove selected individuals from vcf, convert to bfiles, then run plink pairwise estimates

snps=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned
snps_out=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.rm8highmiss
list=remove_8highmiss

# remove individuals from list
vcftools --vcf $snps.vcf --remove $list --recode --recode-INFO-all --out $snps_out
mv $snps_out.recode.vcf $snps_out.vcf

# convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps_out.vcf --set-missing-var-ids @:#\$1,\$2 --out $snps_out

# run plink pairwise estimates, will use the .genome and .imiss files for plinkQC in R
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps_out --genome --out $snps_out
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps_out --missing --out $snps_out
