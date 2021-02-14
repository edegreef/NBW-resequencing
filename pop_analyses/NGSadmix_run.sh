#!/bin/bash

# 1) filter down genotype likelihoods (beagle file) to same loci as snps for LEA and ADMIXTURE, and 2) run admixture analyses with NGSadmix

# extract CHROM and POS info from vcf
awk '{print $1, $2}' NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.miss01.vcf > snp_list_auto_LDpruned_miss01_temp

# remove lines with #'s
sed '/^#/d' snp_list_auto_LDpruned_miss01_temp > snp_list_auto_LDpruned_miss01

# delete temp file
rm snp_list_auto_LDpruned_miss01_temp

# remove "contig" in scaffold names (to match beagle file)
cut -c7- snp_list_auto_LDpruned_miss01 > snp_list_auto_LDpruned_miss01_cut

# replace space with "_" (also to match beagle file)
sed 's/[[:blank:]]/_/g' snp_list_auto_LDpruned_miss01_cut > snp_list_auto_LDpruned_miss01_for_GLfilter

# use nano to add "marker" as row1 (this is so filtering beagle file keeps header)
nano snp_list_auto_LDpruned_miss01_for_GLfilter
#(then press enter, up, and add "marker")

# filter beagle file based on snp list
awk 'NR==FNR {a[$1]++; next} $1 in a' snp_list_auto_LDpruned_miss01_for_GLfilter NBW_angsd.beagle > NBW_angsd.beagle.filtered.miss01

# lower sample size to 36 indivs to compare with snp set (in beagle file, each individual has 3 columns, filtering them out by removing corresponding columns)
cut --complement -f31,32,33,34,35,36,40,41,42,46,47,48,55,56,57,64,65,66,76,77,78,100,101,102,106,107,108,124,125,126,133,134,135,139,140,141,145,146,147 NBW_angsd.beagle.filtered.miss01 > NBW_angsd.beagle.filtered.miss01.n36

# zip it up again for NGSadmix input
gzip NBW_angsd.beagle.filtered.miss01.n36


# 2) run NGS admix on filtered beagle file

#install program
#wget https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp
#g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

cd /home/degreefe/NBW/NGSadmix

# run NGSadmix. using minMaf 0 because loci are already filtered in previous step with maf. plot qopt output file in R
/home/degreefe/programs/NGSadmix -likes NBW_angsd.beagle.filtered.miss01.n36.gz -K 2 -minMaf 0 -o NBW_GL_filtered_miss01_n36_K2 -P 6
