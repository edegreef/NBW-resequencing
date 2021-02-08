#!/bin/bash

# running NGS admix on GL file

#install program
#wget https://raw.githubusercontent.com/ANGSD/angsd/master/misc/ngsadmix32.cpp
#g++ ngsadmix32.cpp -O3 -lpthread -lz -o NGSadmix

cd /home/degreefe/NBW/NGSadmix

# running NGSadmix on full beagle file
/home/degreefe/programs/NGSadmix -likes NBW_angsd.beagle.gz -K 2 -o NBW_GL_K3 -P 6

# lowering sample size to 36 indivs to compare with snp set (in beagle file, each individual has 3 columns, filtering them out by removing corresponding columns)
gunzip NBW_angsd.beagle.gz
cut --complement -f31,32,33,34,35,36,40,41,42,46,47,48,55,56,57,64,65,66,76,77,78,100,101,102,106,107,108,124,125,126,133,134,135,139,140,141,145,146,147 NBW_angsd.beagle > NBW_n36_beagle

# zip it up again for NGSadmix input
gzip NBW_n36_beagle

# choosing 0.0408163 for minMaf filter to represent minor allele count of 2 (keeping consistent with snp set). plot the .qopt output file in R to see admixture results
/home/degreefe/programs/NGSadmix -likes NBW_n36_beagle.gz -K 2 -minMaf 0.0408163 -o NBW_GL_K2_MAC2_n36 -P 6
