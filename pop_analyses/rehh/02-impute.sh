#!/bin/bash

cd /home/degreefe/NBW/snps_2M/rehh/impute

# Filter out the single snp scaffolds
vcftools --vcf NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.vcf --exclude-positions single_snp_windows_min50kb_CHROM_POS --recode --recode-INFO-all --out NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles

mv NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.recode.vcf NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.vcf

# Zip vcf
bgzip NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.vcf
tabix -p vcf NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.vcf.gz

# Impute! (and phase)
java -Xmx70g -jar beagle.28Jun21.220.jar gt=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.vcf.gz iterations=20 out=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.imputed

