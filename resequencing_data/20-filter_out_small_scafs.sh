#!/bin/bash

# Script to make a vcf containing scaffolds with a min length cut-off, here using min 50kb

vcf=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV
scaf_list=scafs_min50kb

# index vcf if not indexed already
bgzip $vcf.vcf
tabix -p vcf $vcf.vcf.gz

# Add word 'contig' because this was added to names in vcf (original contig names were numbers)
sed 's/^/contig/' $scaf_list.txt > $scaf_list.contigname

# Split contig list it's super long. The 'split' command will spit out chunks "xaa", "xab", "xac", "xad"..., The -l 5000 parameter means splitting every 5000 lines
split -l 5000 $scaf_list.contigname

# Convert lists to one-liners to use in bcftools
awk '{print $1}' xaa | paste -s -d, - > xaa.listline
awk '{print $1}' xab | paste -s -d, - > xab.listline
awk '{print $1}' xac | paste -s -d, - > xac.listline

# Set up lists
list1=`cat xaa.listline`
list2=`cat xab.listline`
list3=`cat xac.listline`


# Filter vcf for these scaffolds
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list1 $vcf.vcf.gz > xaa.vcf
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list2 $vcf.vcf.gz > xab.vcf
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list3 $vcf.vcf.gz > xac.vcf

# Create list of CHROM and POS for the snps, also removing the ##header stuff with sed
awk '{print $1, $2}' xaa.vcf > xaa_CHROM_POS
sed '/^#/d' xaa_CHROM_POS > xaa_CHROM_POS_list
rm xaa_CHROM_POS

awk '{print $1, $2}' xab.vcf > xab_CHROM_POS
sed '/^#/d' xab_CHROM_POS > xab_CHROM_POS_list
rm xab_CHROM_POS

awk '{print $1, $2}' xac.vcf > xac_CHROM_POS
sed '/^#/d' xac_CHROM_POS > xac_CHROM_POS_list
rm xac_CHROM_POS

# Merge the CHROM_POS_lists to make one file for filtering in vcftools
cat xaa_CHROM_POS_list xab_CHROM_POS_list xac_CHROM_POS_list  > xaabc_CHROM_POS_list

# Filter to keep xaa-c snps
vcftools --gzvcf $vcf.vcf.gz --positions xaabc_CHROM_POS_list --recode --recode-INFO-all --out $vcf.min50kb
mv $vcf.min50kb.recode.vcf $vcf.min50kb.vcf

# Can make a list of contigs to double check
grep -v "^#" $vcf.min50kb.vcf | cut -f1 | sort | uniq > filtered_contig_list_check.txt
