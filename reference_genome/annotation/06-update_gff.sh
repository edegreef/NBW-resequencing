#!/bin/bash

# Sciprt for adding GMOD IDs to annotations, and updatin gff and protein files with blastp and interproscan outputs. Also a couple final filtering steps at the end.

# load modules
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 maker/2.31.10

cd /home/edegreef/maker

# make map file with more standardized gene IDs (start with the full gff file, then can re-filter AED)
maker_map_ids --prefix GMOD_ --justify 8 NBW_ref_051018_round4.gff > GMOD.map

# copy and rename files before changing IDs since the 'map_...' functions will overwrite
cp NBW_ref_051018_round4.gff NBW_ref_051018_round4.gmod.gff
cp NBW_ref_051018_round4.proteins.fasta NBW_ref_051018_round4.proteins.gmod.fasta
cp NBW_ref_051018_round4.transcripts.fasta NBW_ref_051018_round4.transcripts.gmod.fasta
cp NBW_iprscan.allappl.output NBW_iprscan.allappl.gmod.output
cp NBW_uniprot.blastp NBW_uniprot.gmod.blastp

# use map file to rename IDs in the *gmod files
map_gff_ids GMOD.map NBW_ref_051018_round4.gmod.gff
map_fasta_ids GMOD.map NBW_ref_051018_round4.proteins.gmod.fasta
map_fasta_ids GMOD.map NBW_ref_051018_round4.transcripts.gmod.fasta
map_data_ids GMOD.map NBW_iprscan.allappl.gmod.output
map_data_ids GMOD.map NBW_uniprot.gmod.blastp

# update the maker files with blast and interproscan outputs

# add blast info to gff file (make sure the whole uniprot_SPWH_BLWH_COW.fasta database is in same directory)
maker_functional_gff uniprot_SPWH_BLWH_COW.fasta NBW_uniprot.gmod.blastp NBW_ref_051018_round4.gmod.gff > NBW_ref_051018_round4.gmod.blastp.gff

# add blast info to protiens fasta file
maker_functional_fasta uniprot_SPWH_BLWH_COW.fasta NBW_uniprot.gmod.blastp NBW_ref_051018_round4.proteins.gmod.fasta > NBW_ref_051018_round4.proteins.gmod.blastp.fasta

# add blast info to transcripts fasta file
maker_functional_fasta uniprot_SPWH_BLWH_COW.fasta NBW_uniprot.gmod.blastp NBW_ref_051018_round4.transcripts.gmod.fasta > NBW_ref_051018_round4.transcripts.gmod.blastp.fasta

# add interproscan info to updated gff file
ipr_update_gff NBW_ref_051018_round4.gmod.blastp.gff NBW_iprscan.allappl.gmod.output > NBW_ref_051018_round4.gmod.blastp.iperscan.gff

# not entirely sure what this one adds but let's run it too
iprscan2gff3 NBW_iprscan.allappl.gmod.output NBW_ref_051018_round4.gmod.gff > NBW_ref_051018_round4.gmod.visible_iprscan_domains.gff


##### next finalizing the gff file
# anything in the format "match/match_part" is evidence alignment or rejected model and is there for reference purposes (http://gmod.827538.n3.nabble.com/Re-Maker-accessory-scripts-td4033291.html)
# could keep it in, but uneccessary. The annotations we want are "gene/mRNA/exon/CDS"

# first double check what's in the list
cat NBW_ref_051018_round4.gmod.blastp.iperscan.gff | awk '{a[$3]++}END{for(k in a){print k,a[k]}}'

# can also print 3rd column to see
awk '{print $3}' NBW_ref_051018_round4.gmod.blastp.iperscan.gff  > see_3rdcol.txt

# remove "protein_match/match_part/match" and keep the rest
awk '{ if ($3!="protein_match" && $3!="match_part" && $3!="match") { print } }' NBW_ref_051018_round4.gmod.blastp.iperscan.gff > NBW_ref_051018_round4.gmod.blastp.iperscan.annot.gff

# b/c many contigs in genome and no rna, the annotation won't be complete so I am using a strict AED filter to remove false-positives
perl quality_filter.pl -a 0.25 NBW_ref_051018_round4.gmod.blastp.iperscan.annot.gff > NBW_ref_051018_round4.gmod.blastp.iperscan.annot.filter.gff

# look at gene/stats count again
cat NBW_ref_051018_round4.gmod.blastp.iperscan.annot.filter.gff | awk '{a[$3]++}END{for(k in a){print k,a[k]}}'
