#!/bin/bash

#SBATCH --time=48:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=%x-%j.out
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=vcf2smc_run3

# array jobs was a disaster with so many contigs, so this time I'm just making 1 job with a loop

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/smcpp
mkdir masked_data_run3

list=scafs_min100kb_autosomes.contigname

# running SMC++ with data separated as three groups (LB+NF, AR+IC, and SS). SS already done in the 5group run. Running the LB+NF and AR+IC twice with each subregion as distinguished lineage

while read contig
do
# Labrador & Newfoundland. LB first
smc++ vcf2smc -d HamLB03-1 HamLB03-1 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run3/LBNF.LB.HamLB03-1.$contig.smc $contig LBNF:HamLB03-1,HamLB03-2,HamLB03-3,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505

# Labrador & Newfoundland. NF first
smc++ vcf2smc -d Hyam-2016-06-pooled Hyam-2016-06-pooled --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run3/LBNF.NF.Hyam-2016-06-pooled.$contig.smc $contig LBNF:Hyam-2016-06-pooled,HamLB03-1,HamLB03-2,HamLB03-3,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505

# Arctic & Iceland. AR first
smc++ vcf2smc -d NBW-12-710-509 NBW-12-710-509 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run3/ARIC.AR.NBW-12-710-509.$contig.smc $contig ARIC:NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,Ha14-01-711-503,Ha14-03-717-504,Ha14-04-712-502,Ha14-08-710-503,NBW-2019-B10-718-504,NBW-2019-B11-718-505,Ha14-13-716-505,Ha14-18-717-505,NBW-2019-B13-711-505,NBW-2019-B14-712-505,NBW-2019-B15-712-503,NBW-2019-B16-717-503,NBW-2019-B17-716-503,Ha14-06-716-502-2,Ha14-10-711-502-2

# Arctic & Iceland. IC first
smc++ vcf2smc -d Ha14-01-711-503 Ha14-01-711-503 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run3/ARIC.IC.Ha14-01-711-503.$contig.smc $contig ARIC:Ha14-01-711-503,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,Ha14-03-717-504,Ha14-04-712-502,Ha14-08-710-503,NBW-2019-B10-718-504,NBW-2019-B11-718-505,Ha14-13-716-505,Ha14-18-717-505,NBW-2019-B13-711-505,NBW-2019-B14-712-505,NBW-2019-B15-712-503,NBW-2019-B16-717-503,NBW-2019-B17-716-503,Ha14-06-716-502-2,Ha14-10-711-502-2

# Scotian Shelf
smc++ vcf2smc -d NBW-2016-09 NBW-2016-09 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run3/SS.NBW-2016-09.$contig.smc $contig SS:NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,NBW1904-718-503,NBW1906-710-504,NBW19-03-712-504,NBW1908-711-504

done < $list
