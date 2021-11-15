#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=dedup
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 picard/2.20.6

cd /scratch/edegreef/NBW/aligned_bams

# Remove duplicate reads
ls *.sorted.bam | sed 's/.sorted.bam$//' | parallel --jobs 4 'java -Xmx28000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.sorted.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true'