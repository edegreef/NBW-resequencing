#!/bin/bash

# Script to make a vcf containing scaffolds with a min length cut-off, here using min 100kb
# For smc++ starting with vcf that has base quality filtering done and bialleliic/autosomes only
vcf=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37
new_vcf=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb

# Zip and index vcf if not already done
bgzip $vcf.vcf
tabix -p vcf $vcf.vcf.gz

# Have to add word 'contig' because this was added to names in vcf (original contig names were numbers)
sed 's/^/contig/' scafs_min100kb_autosomes.txt > scafs_min100kb_autosomes.contigname

# Split if list if super long, otherwise not necessary, and if using 'split' command it will spit out chunks "xaa", "xab", "xac", "xad"...

# Convert list to a one-liner to use in bcftools
awk '{print $1}' scafs_min100kb_autosomes.contigname | paste -s -d, - > scafs_min100kb_autosomes.contigname.listline

list1=`cat scafs_min100kb_autosomes.contigname.listline`

# Filter vcf for these scaffolds
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list1 $vcf.vcf.gz > $new_vcf.vcf

# Zip and index file
bgzip $new_vcf.vcf
tabix -p vcf $new_vcf.vcf.gz
