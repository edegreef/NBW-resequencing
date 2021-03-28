#!/bin/bash
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --job-name=interproscan
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=8G
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
#SBATCH --mail-type=ALL
#SBATCH --mail-user=edegreef@ucdavis.edu 

cd /home/edegreef/maker/allappl

module load nixpkgs/16.09 python/2.7.14 gcc/5.4.0 interproscan/5.23-62.0

# test
#interproscan.sh -i test_all_appl.fasta -f tsv -dp -o test_ipr.out

# if using just pfam database, it will take much less time to run (here took about 2 hours)
#interproscan.sh -appl Pfam -dp -f TSV -goterms -iprlookup -pa -t p -i NBW_ref_051018_round4.proteins.fasta -o NBW_iprscan_Pfam.output

# adding more applications here to include infro from multiple databases (~2 days run time). I think this way the functional info is more complete.
interproscan.sh -appl TIGRFAM,PIRSF,SMART,ProSitePatterns,Hamap,Pfam,PRINTS,SUPERFAMILY,Coils -dp -f TSV -goterms -iprlookup -pa -t p -i NBW_ref_051018_round4.proteins.fasta -o NBW_iprscan.allappl.output