# Notes on filtering out SVs and running kinship

snps=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes
SV_list=SV_withpart2_CHROMPOS
high_miss=7highmiss_1dup
n=41

# Filter out structural variants
vcftools --gzvcf $snps.vcf.gz --exclude-positions $SV_list --recode --recode-INFO-all --out $snps.SV
mv $snps.SV.recode.vcf $snps.SV.vcf

# SNPs for pop structure need to filter MAC2 and LD pruned
vcftools --vcf $snps.SV.vcf --mac 2 --recode --recode-INFO-all --out $snps.SV.mac
mv $snps.SV.mac.recode.vcf $snps.SV.mac.vcf

# Make a file with LD pruned snps using LD_pruning.sh
./LD_pruning.sh $snps.SV.mac

# Also need to remove individuals with high-missingness before running kinship, using list of samples in "indiv_highmiss" file
vcftools --vcf $snps.SV.mac.LDpruned.vcf --remove $high_miss --recode --recode-INFO-all --out $snps.SV.mac.LDpruned.$n
mv $snps.SV.mac.LDpruned.$n.recode.vcf $snps.SV.mac.LDpruned.$n.vcf

# Convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps.SV.mac.LDpruned.$n.vcf --set-missing-var-ids @:#\$1,\$2 --out $snps.SV.mac.LDpruned.$n

# Run plink pairwise estimates, will use the .genome and .imiss files for plinkQC in R
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps.SV.mac.LDpruned.$n --genome --out $snps.SV.mac.LDpruned.$n
/home/degreefe/programs/plink --allow-extra-chr --bfile $snps.SV.mac.LDpruned.$n --missing --out $snps.SV.mac.LDpruned.$n

