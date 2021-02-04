#!/bin/bash

#Running stats on reference genome assembly
cd home/degreefe/NBW/reference

#1) Assemblathon (need FAlite.pm in working directory, downloaded from KorfLab github page along with assemblathon_stats.pl script)
assemblathon_stats.pl Northern_bottlenose_whale_051018_shortLabel.fasta > NBW_051018.shortlabel.assemblathon.txt

#2) QUAST 
/home/degreefe/prog/quast-5.0.2/quast.py -t 8 Northern_bottlenose_whale_051018_shortLabel.fasta

#3) Scaffold lengths
cat Northern_bottlenose_whale_051018_shortLabel.fasta | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > scaffoldlengths.csv
