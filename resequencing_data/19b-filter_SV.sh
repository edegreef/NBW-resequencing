snps=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes
SV_list=SV_withpart2_CHROMPOS

# Filter out structural variants
vcftools --gzvcf $snps.vcf.gz --exclude-positions $SV_list --recode --recode-INFO-all --out $snps.SV
mv $snps.SV.recode.vcf $snps.SV.vcf
