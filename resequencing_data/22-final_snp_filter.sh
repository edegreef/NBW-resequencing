# Prep snp files for analyses
# After identifying kinship pairs and deciding which individuals of each pair to remove

snps=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV
snps2=NBW2_SNPS_2M.filter1.miss.biallel.ID.autosomes.SV.mac.LDprunedr05
list=7highmiss_1dup_4kin
n=37

# -----------------
# Making a file for smcp++ that still retains MAC and LD
vcftools --vcf $snps.vcf --remove $list --recode --recode-INFO-all --out $snps.$n
mv $snps.$n.recode.vcf $snps.$n.vcf

# -----------------
# Prep for pop structure analyses
vcftools --vcf $snps2.vcf --remove $list --recode --recode-INFO-all --out $snps2.$n
mv $snps2.$n.recode.vcf $snps2.$n.vcf

# -----------------
# Create separate vcf files by group in prep for Reich's fst
# Arctic
vcftools --vcf $snps2.$n.vcf --keep pop_arctic --recode --recode-INFO-all --out $snps2.$n.arctic
mv $snps2.$n.arctic.recode.vcf $snps2.$n.arctic.vcf

# Iceland
vcftools --vcf $snps2.$n.vcf --keep pop_iceland --recode --recode-INFO-all --out $snps2.$n.iceland
mv $snps2.$n.iceland.recode.vcf $snps2.$n.iceland.vcf

# Labrador
vcftools --vcf $snps2.$n.vcf --keep pop_labrador --recode --recode-INFO-all --out $snps2.$n.labrador
mv $snps2.$n.labrador.recode.vcf $snps2.$n.labrador.vcf

# Newfoundland
vcftools --vcf $snps2.$n.vcf --keep pop_newfoundland --recode --recode-INFO-all --out $snps2.$n.newfoundland
mv $snps2.$n.newfoundland.recode.vcf $snps2.$n.newfoundland.vcf

# Scotian Shelf
vcftools --vcf $snps2.$n.vcf --keep pop_scotianshelf --recode --recode-INFO-all --out $snps2.$n.scotianshelf
mv $snps2.$n.scotianshelf.recode.vcf $snps2.$nscotianshelf.vcf

# Zip
bgzip $snps2.$n.arctic.vcf
bgzip $snps2.$n.iceland.vcf
bgzip $snps2.$n.labrador.vcf
bgzip $snps2.$n.newfoundland.vcf
bgzip $snps2.$n.scotianshelf.vcf

# -----------------
# Prep files for admixture
# Filter out snps with high missingness (in vcftools 1=no missing, 0=all missing)
vcftools --vcf $snps2.$n.vcf --max-missing 0.9 --recode --recode-INFO-all --out $snps2.$n.miss01

mv $snps2.$n.miss01.recode.vcf $snps2.$n.miss01.vcf

# Convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps2.$n.miss01.vcf --set-missing-var-ids @:#\$1,\$2 --out $snps2.$n.miss01

# Prep for SNMF
# Convert to ped/map 
/home/degreefe/programs/bcftools-1.9/bcftools view -H $snps2.$n.miss01.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
vcftools --vcf $snps2.$n.miss01.vcf --plink --chrom-map chrom-map.txt --out $snps2.$n.miss01

