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
#SBATCH --job-name=vcf2smc_ICSS

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/smcpp
mkdir masked_data_run_ICSS

list=scafs_min100kb_autosomes.contigname

while read contig
do

smc++ vcf2smc -d Ha14-01-711-503 Ha14-01-711-503 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run_ICSS/IC_SS.$contig.smc $contig IC:Ha14-01-711-503,Ha14-03-717-504,Ha14-04-712-502,Ha14-08-710-503,Ha14-13-716-505,Ha14-18-717-505,Ha14-06-716-502-2,Ha14-10-711-502-2 SS:NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,NBW1904-718-503,NBW1906-710-504,NBW19-03-712-504,NBW1908-711-504

smc++ vcf2smc -d NBW-2016-09 NBW-2016-09 --mask /scratch/edegreef/smcpp/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/smcpp/smcpp_snps_input.vcf.gz /scratch/edegreef/smcpp/masked_data_run_ICSS/SS_IC.$contig.smc $contig SS:NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,NBW1904-718-503,NBW1906-710-504,NBW19-03-712-504,NBW1908-711-504 IC:Ha14-01-711-503,Ha14-03-717-504,Ha14-04-712-502,Ha14-08-710-503,Ha14-13-716-505,Ha14-18-717-505,Ha14-06-716-502-2,Ha14-10-711-502-2

done < $list