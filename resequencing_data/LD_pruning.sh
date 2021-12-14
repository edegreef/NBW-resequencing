#!/bin/bash

snps=$1

# LD pruning SNPs
/home/degreefe/programs/plink --vcf $snps.vcf --allow-extra-chr --indep-pairwise 50 5 0.5 --set-missing-var-ids @:#\$1,\$2 --out pruned.50window.r0.5

# preparing the prune-in file to filter vcf
sed 's/:/\t/g' pruned.50window.r0.5.prune.in > temp.in
sed 's/...$//' temp.in > temp2.in
mv temp2.in list.pruned.50window.r0.5.prune.in
rm temp.in

# prune the snps
vcftools --vcf $snps.vcf --positions list.pruned.50window.r0.5.prune.in --recode --recode-INFO-all --out $snps.LDprunedr05

# renaming to remove the "recode" thing
mv $snps.LDprunedr05.recode.vcf $snps.LDprunedr05.vcf
