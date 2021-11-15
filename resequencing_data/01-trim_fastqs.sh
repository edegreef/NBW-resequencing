#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=trim_fastq
#SBATCH --output=%x-%j.out

# Help with code from Tony Einfeldt

# Load modules
module load nixpkgs/16.09 trimmomatic/0.36

# Use TruSeq3-PE-2.fa to account for read-through, obtained at https://github.com/timflutre/trimmomatic/blob/master/adapters/TruSeq3-PE-2.fa


# Added to TruSeq3-PE-2.fa to account for Nextera transposase sequence:
	# >transposase1
	# GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAG
	# >transposase1_rc
	# CTGTCTCTTATACACATCTCTGATGGCGCGAGGGAGGC
	# >transposase2
	# GCCTTGCCAGCCCGCTCAGAGATGTGTATAAGAGACAG
	# >transposase2_rc
	# CTGTCTCTTATACACATCTCTGAGCGGGCTGGCAAGGC	
	
cd /scratch/edegreef/NBW_oct2021

# Run trimmomatic	

ls ./*_R1_001.fastq.gz | sed 's/_R1_001.fastq.gz$//' | parallel 'java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.36.jar PE {}_R1_001.fastq.gz {}_R2_001.fastq.gz /scratch/edegreef/NBW/{}_R1_paired.fastq.gz /scratch/edegreef/NBW/{}_R1_unpaired.fastq.gz /scratch/edegreef/NBW/{}_R2_paired.fastq.gz /scratch/edegreef/NBW/{}_R2_unpaired.fastq.gz ILLUMINACLIP:/scratch/edegreef/illumina_adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36'
