#!/bin/bash

#SBATCH --time=10:00:00
#SBATCH --account=def-coling
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu
#SBATCH --nodes=1
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --kill-on-invalid=yes
#SBATCH --job-name=bam_coverage
#SBATCH --output=%x-%j.out

cd /scratch/edegreef/NBW/aligned_bams/dedupRG

# Load modules
module load StdEnv/2020 samtools/1.12

# Check modal coverage for all deDupRG.bam files in folder
for i in *deDupRG.bam
do
echo "Running bam_coverage: "$i        
./bam_coverage.sh $i
done

# Can also do one file at a time (example)
#./bam_coverage.sh Ha14-13-716-505.deDupRG.bam
