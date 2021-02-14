#!/bin/bash

# Identify private and shared SNPs between populations

# first need to prep vcfs for each pop as well as 'all pops minus target pop'

cd /home/degreefe/NBW/reseq_newsnps/samplesizes/isec

snps=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36

# making vcfs by pop
vcftools --vcf $snps.vcf --recode --recode-INFO-all --keep pop_arctic --out $snps.arctic
vcftools --vcf $snps.vcf --recode --recode-INFO-all --keep pop_iceland --out $snps.iceland
vcftools --vcf $snps.vcf --recode --recode-INFO-all --keep pop_labrador --out $snps.labrador
vcftools --vcf $snps.vcf --recode --recode-INFO-all --keep pop_newfoundland --out $snps.newfoundland
vcftools --vcf $snps.vcf --recode --recode-INFO-all --keep pop_scotianshelf --out $snps.scotianshelf

# making vcfs by not-pop
vcftools --vcf $snps.vcf --recode --recode-INFO-all --remove pop_arctic --out $snps.no_arctic
vcftools --vcf $snps.vcf --recode --recode-INFO-all --remove pop_iceland --out $snps.no_iceland
vcftools --vcf $snps.vcf --recode --recode-INFO-all --remove pop_labrador --out $snps.no_labrador
vcftools --vcf $snps.vcf --recode --recode-INFO-all --remove pop_newfoundland --out $snps.no_newfoundland
vcftools --vcf $snps.vcf --recode --recode-INFO-all --remove pop_scotianshelf --out $snps.no_scotianshelf

# remove 'recode' in each vcf title
mv $snps.arctic.recode.vcf $snps.arctic.vcf
mv $snps.iceland.recode.vcf $snps.iceland.vcf
mv $snps.labrador.recode.vcf $snps.labrador.vcf
mv $snps.newfoundland.recode.vcf $snps.newfoundland.vcf
mv $snps.scotianshelf.recode.vcf $snps.scotianshelf.vcf
mv $snps.no_arctic.recode.vcf $snps.no_arctic.vcf
mv $snps.no_iceland.recode.vcf $snps.no_iceland.vcf
mv $snps.no_labrador.recode.vcf $snps.no_labrador.vcf
mv $snps.no_newfoundland.recode.vcf $snps.no_newfoundland.vcf
mv $snps.no_scotianshelf.recode.vcf $snps.no_scotianshelf.vcf

# zip files
bgzip $snps.arctic.vcf
bgzip $snps.iceland.vcf
bgzip $snps.labrador.vcf
bgzip $snps.newfoundland.vcf
bgzip $snps.scotianshelf.vcf
bgzip $snps.no_arctic.vcf
bgzip $snps.no_iceland.vcf
bgzip $snps.no_labrador.vcf
bgzip $snps.no_newfoundland.vcf
bgzip $snps.no_scotianshelf.vcf

# index files
tabix -p vcf $snps.arctic.vcf.gz 
tabix -p vcf $snps.iceland.vcf.gz
tabix -p vcf $snps.labrador.vcf.gz
tabix -p vcf $snps.newfoundland.vcf.gz
tabix -p vcf $snps.scotianshelf.vcf.gz
tabix -p vcf $snps.no_arctic.vcf.gz
tabix -p vcf $snps.no_iceland.vcf.gz
tabix -p vcf $snps.no_labrador.vcf.gz
tabix -p vcf $snps.no_newfoundland.vcf.gz
tabix -p vcf $snps.no_scotianshelf.vcf.gz

# next use bcftools isec to look at overlaps between vcfs to identify private and shared snps. -p flag creates directory with listed name, and in directory contains 4 vcfs (private to vcf1, private to vcf2, shared to vcf1, and shared to vcf2) and readme file.
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_arctic $snps.arctic.vcf.gz $snps.no_arctic.vcf.gz
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_iceland $snps.iceland.vcf.gz $snps.no_iceland.vcf.gz
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_labrador $snps.labrador.vcf.gz $snps.no_labrador.vcf.gz
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_newfoundland $snps.newfoundland.vcf.gz $snps.no_newfoundland.vcf.gz
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_scotianshelf $snps.scotianshelf.vcf.gz $snps.no_scotianshelf.vcf.gz