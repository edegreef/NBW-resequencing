# Notes on installing SMC++ and preparing input files, code guidance from Matt Thorstensen

cd /scratch/edegreef/SMC

# download smc++
wget https://github.com/popgenmethods/smcpp/releases/download/v1.15.2/smcpp-1.15.2-Linux-x86_64.sh

chmod +x smcpp-1.15.2-Linux-x86_64.sh

# install
./smcpp-1.15.2-Linux-x86_64.sh

# go through reading and acceptng license - answer 'yes' at the end, then program will be installed under /home/edegreef/smcpp

# after installation finishes, then run this line to use it
source /home/edegreef/smcpp/bin/activate

# check it works
smc++ {options}

# prep input files
# use vcf that is filtered for quality but not MAF yet. Also using autosomes only.
# planning to use NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf.gz file

# make vcf of snps that were filted from the raw file for creating a mask file
# will need to run the ID.sh script to match contig and sample IDs
# going back to biology-01 server
./04-edit_label_ID.sh #NBW_platypus_allvariantcalls.vcf

# zip # oops can just skip this since will zip later..
bgzip NBW2_allvariantcalls_2mergedbam.ID.vcf
tabix -p vcf NBW2_allvariantcalls_2mergedbam.ID.vcf.gz

# filter individuals from raw file to keep consistent
vcftools --gzvcf NBW2_allvariantcalls_2mergedbam.ID.vcf.gz --remove 7highmiss_1dup_4kin --recode --recode-INFO-all --out NBW2_allvariantcalls_2mergedbam.ID.37

mv NBW2_allvariantcalls_2mergedbam.ID.37.recode.vcf NBW2_allvariantcalls_2mergedbam.ID.37.vcf

# zip 
bgzip NBW2_allvariantcalls_2mergedbam.ID.37.vcf
tabix -p vcf NBW2_allvariantcalls_2mergedbam.ID.37.vcf.gz

# now make vcf of filtered out sites
/home/degreefe/programs/bcftools-1.9/bcftools isec -p isec_filteredout NBW2_allvariantcalls_2mergedbam.ID.37.vcf.gz NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf.gz

# isec_filteredout/0000.vcf is the one we want (snps that were filtered out)
cp isec_filteredout/0000.vcf 0000.vcf

# rename for smc input
mv 0000.vcf smcpp_removed_sites.vcf

# make list of contigs from vcf file
#gunzip -c NBW_platypus_allvariantcalls.ID.n36.vcf.gz 
gunzip -c NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf.gz > NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf
#grep -v "^#" NBW_platypus_allvariantcalls.ID.n36.vcf | cut -f1 | sort | uniq > unfiltered_contig_list.txt​
grep -v "^#" NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf | cut -f1 | sort | uniq > filtered_contig_list.txt​

# making bed file for masking
# install bedops
cd /home/degreefe/programs
mkdir bedops
cd bedops
wget https://github.com/bedops/bedops/releases/download/v2.4.39/bedops_linux_x86_64-v2.4.39.tar.bz2
tar jxvf bedops_linux_x86_64-v2.4.39.tar.bz2

cd /home/degreefe/NBW/reseq_newsnps/for_smc
export PATH=$PATH:/home/degreefe/programs/bedops/bin

# convert to vcf to bed
switch-BEDOPS-binary-type --megarow

vcf2bed --do-not-sort --deletions < smcpp_removed_sites.vcf > nbw_deletions.bed
vcf2bed --snvs < smcpp_removed_sites.vcf > nbw_snps.bed

# merge bed files
bedops --everything nbw_{deletions,snps}.bed > smcpp_removed_sites.bed

# sort bed
sort-bed smcpp_removed_sites.bed > smcpp_removed_sites_sorted.bed

# zip & index
bgzip smcpp_removed_sites_sorted.bed
tabix -p bed smcpp_removed_sites_sorted.bed.gz

# use this final bed file in --mask for SMC++ vcf2smc. 
# transfering input files to cedar cluster now to run SMC++ things on arrays

# need to add contig lengths to vcf header (I guess the bcftools I installed on biology-01 doesn't have the -f option so just going to do this on graham)
module load StdEnv/2020 gcc/9.3.0 bcftools/1.13

# -f file is the reference genome index file. "contig" because I just added the word 'contig' to each contig to match vcf names
bcftools reheader -f contig_fai.fai NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.37.min100kb.vcf.gz -o smcpp_snps_input.vcf.gz
tabix -p vcf smcpp_snps_input.vcf.gz

# next run vcf2smc script
# max jobs at a time on cedar is 1000, also finding out that max array range, cannot go past 9999.. So will just use a loop instead. 
