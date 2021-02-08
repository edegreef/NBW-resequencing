#!/bin/bash

# 1) Remove individuals in kin pairs (one in each pair) and prep files for analyses
snps=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned
snps_out=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36
list=remove_8highmiss_5kinpair

# remove individuals from list (high missingness and kin pair)
vcftools --vcf $snps.vcf --remove $list --recode --recode-INFO-all --out $snps_out
mv $snps_out.recode.vcf $snps_out.vcf

# convert to bim/bed/fam, using these files for most pop analyses
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $snps_out.vcf --set-missing-var-ids @:#\$1,\$2 --out $snps_out

# increase max missingness filter to 0.1 (0.9 for vcftools) to use in LEA and Rtsne
vcftools --vcf $snps_out.vcf --max-missing 0.9 --recode --recode-INFO-all --out $snps_out.miss01
mv $snps_out.miss01.recode.vcf $snps_out.miss01.vcf

# convert to ped/map (for LEA program)
/home/degreefe/programs/bcftools-1.9/bcftools view -H $snps_out.miss01.vcf | cut -f 1 | uniq | awk '{print $0"\t"$0}' > chrom-map.txt
vcftools --vcf $snps_out.miss01.vcf --plink --chrom-map chrom-map.txt --out $snps_out.miss01

# imputing remaining snps (for Rtsne program). LD-pruned file probably ok for tsne, but if want to impute snps for other things may need to do so with not-pruned file.
java -jar $EBROOTBEAGLE/beagle.jar gt=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.miss01.vcf iterations=20 out=NBW_platypus_SNPs.filter1.filter2.ID.autosomes.LDpruned.n36.miss01.imputed


# 2) sex-linked snps
Xsnps=NBW_platypus_SNPs.filter1.filter2.ID.Xchr
Ysnps=NBW_platypus_SNPs.filter1.filter2.ID.Ychr
list_females=females
list_males=males

# prune snps
./LD_pruning.sh $Xsnps.vcf
./LD_pruning.sh $Ysnps.vcf

# make datasets to separate females and males for X and Y
vcftools --vcf $Xsnps.LDpruned.vcf --keep $list_females --recode --recode-INFO-all --out $Xsnps.LDpruned.femalesX
mv $Xsnps.LDpruned.femalesX.recode.vcf $Xsnps.LDpruned.femalesX.vcf

vcftools --vcf $Ysnps.LDpruned.vcf --keep $list_males --recode --recode-INFO-all --out $Ysnps.LDpruned.malesY
mv $Ysnps.LDpruned.malesY.recode.vcf $Ysnps.LDpruned.malesY.vcf

# convert to bim/bed/fam
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $Xsnps.LDpruned.femalesX.vcf --set-missing-var-ids @:#\$1,\$2 --out $Xsnps.LDpruned.femalesX
/home/degreefe/programs/plink --allow-extra-chr --make-bed --vcf $Ysnps.LDpruned.malesY.vcf --set-missing-var-ids @:#\$1,\$2 --out $Ysnps.LDpruned.malesY
