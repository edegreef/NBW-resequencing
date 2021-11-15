#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=170G
#SBATCH --cpus-per-task=32
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=dedup_large
#SBATCH --output=%x-%j.out

# Load modules
module load nixpkgs/16.09 picard/2.20.6

#cd /scratch/edegreef/NBW/aligned_bams
cd /scratch/edegreef/NBW/aligned_bams/part2

# Ran second chunk of files separately because program would keep crashing at my >20GB bam files
# Increased memory and also added these two parameters: "USE_JDK_DEFLATER=true USE_JDK_INFLATER=true"

# Remove duplicate reads
ls *.sorted.bam | sed 's/.sorted.bam$//' | parallel --jobs 1 'java -Xmx50000m -jar $EBROOTPICARD/picard.jar MarkDuplicates I={}.sorted.bam O={}.deDup.bam M={}_deDupMetrics.txt REMOVE_DUPLICATES=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'