#!/bin/bash
#SBATCH --mail-user=edegreef@myucdavis.edu
#SBATCH --mail-type=ALL
#SBATCH --time=48:00:00
#SBATCH --mem=10GB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --job-name=vcf2smc_run2

source /home/edegreef/smcpp/bin/activate

cd /scratch/edegreef/SMC
mkdir masked_data_run2

list=scafs_min100kb_autosomes_contigname.txt

# running SMC with all samples as one population. Running it 5 times to switch which region is listed at the distinguished lineage

while read contig
do
# Labrador as main
smc++ vcf2smc -d HamLB03-1 HamLB03-1 --mask /scratch/edegreef/SMC/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/SMC/smcpp_input_min100kb_contigheader.vcf.gz /scratch/edegreef/SMC/masked_data_run2/ALL.LB-HamLB03-1.$contig.smc $contig ALL:HamLB03-1,HamLB03-2,HamLB03-3,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,O-BKW1305-710-504,O-BKW1312-711-504,O-BKW1313-718-505,O-BKW1315-716-504,O-BKW1317-710-503,O-BKW1318-711-503,O-BKW1319-712-503,O-BKW1328-711-505,O-BKW1330-716-505,O-BKW1332-716-502,O-BKW1333-717-502,O-BKW1335-710-505,O-BKW1340-717-505,O-BKW1345-718-504,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505,NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,O-BKW0957-712-505,O-BKW0958-716-503,O-BKW0960-717-503,O-BKW0962-718-503

# Arctic as main
smc++ vcf2smc -d NBW-12-710-509 NBW-12-710-509 --mask /scratch/edegreef/SMC/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/SMC/smcpp_input_min100kb_contigheader.vcf.gz /scratch/edegreef/SMC/masked_data_run2/ALL.AR-NBW-12-710-509.$contig.smc $contig ALL:HamLB03-1,HamLB03-2,HamLB03-3,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,O-BKW1305-710-504,O-BKW1312-711-504,O-BKW1313-718-505,O-BKW1315-716-504,O-BKW1317-710-503,O-BKW1318-711-503,O-BKW1319-712-503,O-BKW1328-711-505,O-BKW1330-716-505,O-BKW1332-716-502,O-BKW1333-717-502,O-BKW1335-710-505,O-BKW1340-717-505,O-BKW1345-718-504,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505,NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,O-BKW0957-712-505,O-BKW0958-716-503,O-BKW0960-717-503,O-BKW0962-718-503

# Iceland as main
smc++ vcf2smc -d O-BKW1328-711-505 O-BKW1328-711-505 --mask /scratch/edegreef/SMC/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/SMC/smcpp_input_min100kb_contigheader.vcf.gz /scratch/edegreef/SMC/masked_data_run2/ALL.IC-O-BKW1328-711-505.$contig.smc $contig ALL:HamLB03-1,HamLB03-2,HamLB03-3,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,O-BKW1305-710-504,O-BKW1312-711-504,O-BKW1313-718-505,O-BKW1315-716-504,O-BKW1317-710-503,O-BKW1318-711-503,O-BKW1319-712-503,O-BKW1328-711-505,O-BKW1330-716-505,O-BKW1332-716-502,O-BKW1333-717-502,O-BKW1335-710-505,O-BKW1340-717-505,O-BKW1345-718-504,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505,NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,O-BKW0957-712-505,O-BKW0958-716-503,O-BKW0960-717-503,O-BKW0962-718-503

# Newfoundland as main
smc++ vcf2smc -d Hyam-2016-06-pooled Hyam-2016-06-pooled --mask /scratch/edegreef/SMC/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/SMC/smcpp_input_min100kb_contigheader.vcf.gz /scratch/edegreef/SMC/masked_data_run2/ALL.NF-Hyam-2016-06-pooled.$contig.smc $contig ALL:HamLB03-1,HamLB03-2,HamLB03-3,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,O-BKW1305-710-504,O-BKW1312-711-504,O-BKW1313-718-505,O-BKW1315-716-504,O-BKW1317-710-503,O-BKW1318-711-503,O-BKW1319-712-503,O-BKW1328-711-505,O-BKW1330-716-505,O-BKW1332-716-502,O-BKW1333-717-502,O-BKW1335-710-505,O-BKW1340-717-505,O-BKW1345-718-504,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505,NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,O-BKW0957-712-505,O-BKW0958-716-503,O-BKW0960-717-503,O-BKW0962-718-503

# Scotian Shelf as main
smc++ vcf2smc -d NBW-2016-09 NBW-2016-09 --mask /scratch/edegreef/SMC/smcpp_removed_sites_sorted.bed.gz /scratch/edegreef/SMC/smcpp_input_min100kb_contigheader.vcf.gz /scratch/edegreef/SMC/masked_data_run2/ALL.SS-NBW-2016-09.$contig.smc $contig ALL:HamLB03-1,HamLB03-2,HamLB03-3,NBW-12-710-509,NBW-13-702-507,NBW-15-708-507,NBW-16-705-509,O-BKW1305-710-504,O-BKW1312-711-504,O-BKW1313-718-505,O-BKW1315-716-504,O-BKW1317-710-503,O-BKW1318-711-503,O-BKW1319-712-503,O-BKW1328-711-505,O-BKW1330-716-505,O-BKW1332-716-502,O-BKW1333-717-502,O-BKW1335-710-505,O-BKW1340-717-505,O-BKW1345-718-504,Hyam-2016-06-pooled,NBW-2016-02-706-509,NBW-2017-01,NBW-2017-03,NBW-2017-2-705-505,NBW-2017-4-704-505,NBW-2016-09,NBW-2016-11,NBW-2016-13-707-505,NBW-2016-14,NBW-2019-01-709-505,O-BKW0957-712-505,O-BKW0958-716-503,O-BKW0960-717-503,O-BKW0962-718-503

done < $list