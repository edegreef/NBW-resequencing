#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=blastp
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /scratch/edegreef/maker

# download the UniProtKB protein data from uniprot.org (as fastas)
# get for blue whale, sperm whale, and cow
# combine the 3 fastas into 1

#cat uniprot-taxonomy_9755.SPWH.fasta uniprot-taxonomy_9771.BLWH.fasta uniprot-taxonomy_9913.COW.fasta > uniprot_SPWH_BLWH_COW.fasta

module load StdEnv/2020 gcc/9.3.0 blast+/2.11.0

makeblastdb -in uniprot_SPWH_BLWH_COW.fasta -dbtype prot

blastp -query NBW_ref_051018_round4.proteins.fasta -db uniprot_SPWH_BLWH_COW.fasta -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt 6 -out NBW_uniprot.blastp