#!/bin/bash

# script for editing labels in the vcf, shortening sample ID and adding "contig" to each contig number

snps=/home/degreefe/NBW/reseq_newsnps/NBW_platypus_SNPs.filter1.filter2

#1) change sample ID in vcf to just contain the individual ID (removing the path that is currently in sample name).
# make index file
bgzip -c $snps.vcf > $snps.vcf.gz
tabix -p vcf $snps.vcf.gz

# make list of current sample ID names
vcf-query -l $snps.vcf.gz > vcf_sampleID

# shorten sample IDs by removing text up to the last "/"
awk '{print $NF}' FS=/ vcf_sampleID > vcf_sampleID_shortened

# modify header in vcf file to change sample ID to the shortened ID list
/home/degreefe/programs/bcftools-1.9/bcftools reheader -s vcf_sampleID_shortened $snps.vcf.gz -o $snps.sample.vcf.gz
tabix -p vcf $snps.sample.vcf.gz

# 2) add the word "contig" to each contig number (i.e. "3" will turn into "contig3"). This will help for some downstream programs (such as plink).
# make list of current contig numbers
/home/degreefe/programs/bcftools-1.9/bcftools index --stats $snps.sample.vcf.gz > $snps.contigstat
awk '{print $1}' $snps.contigstat > contig_list
sed -e 's/^/contig/' contig_list > contig_list_newname
paste contig_list contig_list_newname > contig_rename.txt
rm contig_list
rm contig_list_newname

# annotate vcf to update contig names
/home/degreefe/programs/bcftools-1.9/bcftools annotate --rename-chrs contig_rename.txt $snps.sample.vcf.gz -o $snps.IDtemp.vcf
rm $snps.sample*

# also remove the ##contig lines since contig names are updated now and they take up a lot of unecessary space in header info 
egrep -v "^##contig" $snps.IDtemp.vcf > $snps.ID.vcf
rm $snps.IDtemp.vcf