#!/bin/bash

# script to make a vcf containing scaffolds with a min length cut-off, here using min 100kb
# for smc++ starting with vcf that has base quality filtering done and bialleliic/autosomes only; don't use MAF filter here 
vcf=NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss04.n36.min5kb
new_vcf=NBW_platypus_SNPs.filter1.ID.biallel.autosomes.miss04.n36.min100kb

# index vcf if not indexed already
bgzip $vcf.vcf
tabix -p vcf $vcf.vcf.gz

# have to add word 'contig' because this was added to names in vcf (original contig names were numbers)
sed 's/^/contig/' scafs_min100kb_autosomes.txt > scafs_min100kb_autosomes.contigname

# split if list if super long, otherwise not necessary, and if using 'split' command it will spit out chunks "xaa", "xab", "xac", "xad"...

# convert list to a one-liner to use in bcftools
awk '{print $1}' scafs_min100kb_autosomes.contigname | paste -s -d, - > scafs_min100kb_autosomes.contigname.listline

list1=`cat scafs_min100kb_autosomes.contigname.listline`

# filter vcf for these scaffolds
/home/degreefe/programs/bcftools-1.9/bcftools filter --regions $list1 $vcf.vcf.gz > $new_vcf.vcf

# zip and index file
bgzip $new_vcf.vcf
tabix -p vcf $new_vcf.vcf.gz