#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=run_busco
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --cpus-per-task=48
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

# Running BUSCO 4.1.4, and using Blast 2.10.1+
# Need to edit busco/config/config.ini file with appropriate paths 
# Using mammalia_odb10 for lineage dataset (https://busco-data.ezlab.org/v4/data/lineages/mammalia_odb10.2020-09-10.tar.gz)

cd /scratch/edegreef
module load gcc/5.4.0 hmmer/3.1b2 emboss/6.6.0 r/3.5.0 augustus/3.2.3/ python/3.5.4

export AUGUSTUS_CONFIG_PATH=/home/edegreef/augustus_config
export BUSCO_CONFIG_FILE=/home/edegreef/busco/config/config.ini

/home/edegreef/busco_env/bin/busco --offline --in /scratch/edegreef/Northern_bottlenose_whale_051018_shortLabel.fasta --mode genome --lineage_dataset /scratch/edegreef/mammalia_odb10 -o NBW_busco_run -c 48 --augustus_species human
