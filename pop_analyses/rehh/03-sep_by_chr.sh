#!/bin/bash

# Preparing vcf for rehh program, separating vcf for each chromosome and also dividing by northern (NOR=iceland, arctic, labrador newfoundland) and scotian shelf (SS), also making a sep file for Iceland (IC)

# Using files (chr1.txt, chr2.txt...) which has lists of scaffolds in chromosome

# Renamed 'NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.miss01.mac.min50kb.rmsingles.imputed.vcf.gz' to 'NBW2_SNPS_2M.abbr.37.min50kb.imputed.vcf.gz'


cd /home/degreefe/NBW/snps_2M/rehh

snps=NBW2_SNPS_2M.abbr.37.min50kb.imputed
list=chr_scafs_list_of_list 
pop_SS=pop_scotianshelf
pop_IC=pop_iceland

# Loop for each chromosome

while read chr
do

# Add word 'contig' because this was added to names in vcf (original contig names were numbers)
sed 's/^/contig/' $chr.txt > $chr.contigname

# Convert list to a one-liner to use in bcftools
awk '{print $1}' $chr.contigname | paste -s -d, - > temp.listline

contig_keep=`cat temp.listline`

# Filter vcf for these scaffolds
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $contig_keep $snps.vcf.gz > $snps.$chr.vcf

# Make NOR, SS, and IC vcfs
vcftools --vcf $snps.$chr.vcf --remove $pop_SS --recode --recode-INFO-all --out $snps.$chr.NOR
mv $snps.$chr.NOR.recode.vcf $snps.$chr.NOR.vcf

vcftools --vcf $snps.$chr.vcf --keep $pop_SS --recode --recode-INFO-all --out $snps.$chr.SS
mv $snps.$chr.SS.recode.vcf $snps.$chr.SS.vcf

vcftools --vcf $snps.$chr.vcf --keep $pop_IC --recode --recode-INFO-all --out $snps.$chr.IC
mv $snps.$chr.IC.recode.vcf $snps.$chr.IC.vcf

done < $list

# Clean up a little
rm *contigname
rm temp.listline
