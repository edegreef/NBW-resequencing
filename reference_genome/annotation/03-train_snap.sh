#!/bin/bash

# merge maker outputs and train snap

maker_out=/scratch/edegreef/ref_genome/annotation/nbw/maker_round1/NBW_ref_051018_round1.gff 
snap_round=1

# load modules
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 maker/2.31.10

# set up snap directory
mkdir snap/round$snap_round
cd snap/round$snap_round

# convert .gff file to SNAP format, outputs genome.ann and genome.dna
maker2zff -n $maker_out

# check info, can look at these outputs to see how many genes, etc
fathom -gene-stats genome.ann genome.dna > stats_output
fathom -validate genome.ann genome.dna > validate_output

# check model with errors
cat validate_output | grep "error" > model_errors

# make a file with model IDs that have errors and remove them, making an updated genome.ann file
awk '{print $2}' model_errors > model_err_num
grep -vwE -f model_err_num genome.ann > genome.ann2

# check again that all models with errors removed
fathom -gene-stats genome.ann2 genome.dna

# train snap
fathom -categorize 1000 genome.ann2 genome.dna
fathom -export 1000 -plus uni.ann uni.dna
forge export.ann export.dna
hmm-assembler.pl train_snap$snap_round . > train_snap$snap_round.hmm
