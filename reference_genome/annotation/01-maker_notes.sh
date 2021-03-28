# Notes for prepping runs on MAKER on reference genome
# (running this on graham cluster on Compute Canada)

# modules (will need to email CC to get initial access to maker)
module load nixpkgs/16.09 gcc/7.3.0 openmpi/3.1.2 maker/2.31.10

# make copies of the control files
maker -CTL

# had some issues with repbase, even with it installed for maker version3; ended up using older maker (2.31.10) version and added this for repeatmasker in exe.ctl file
# change the repeatmasker thing to this in the exe.ctl files
/cvmfs/soft.computecanada.ca/easybuild/software/2017/avx2/Compiler/gcc5.4/repeatmasker/4-0-7/RepeatMasker

# download protein files from ensembl

# Blue Whale: https://uswest.ensembl.org/Balaenoptera_musculus/Info/Index
wget http://ftp.ensembl.org/pub/release-103/fasta/balaenoptera_musculus/pep/Balaenoptera_musculus.mBalMus1.v2.pep.all.fa.gz

# Sperm whale: https://useast.ensembl.org/Physeter_catodon/Info/Index
wget http://ftp.ensembl.org/pub/release-103/fasta/physeter_catodon/pep/Physeter_catodon.ASM283717v2.pep.all.fa.gz

# Cattle: https://uswest.ensembl.org/Bos_taurus/Info/Index
wget http://ftp.ensembl.org/pub/release-103/fasta/bos_taurus/pep/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz

# unzip files
gunzip Balaenoptera_musculus.mBalMus1.v2.pep.all.fa.gz
gunzip Physeter_catodon.ASM283717v2.pep.all.fa.gz
gunzip Bos_taurus.ARS-UCD1.2.pep.all.fa.gz

# prepare reference genome chunks (splitting into 20 parts to run simultaneously to save time)
mv Northern_bottlenose_whale_051018_shortLabel.fasta* reference_chunks/NBW_ref_051018.fasta
fasta_tool -chunks 20 NBW_ref_051018.fasta

# need to index each part
samtools faidx NBW_ref_051018_00.fasta
samtools faidx NBW_ref_051018_01.fasta
samtools faidx NBW_ref_051018_02.fasta
samtools faidx NBW_ref_051018_03.fasta
samtools faidx NBW_ref_051018_04.fasta
samtools faidx NBW_ref_051018_05.fasta
samtools faidx NBW_ref_051018_06.fasta
samtools faidx NBW_ref_051018_07.fasta
samtools faidx NBW_ref_051018_08.fasta
samtools faidx NBW_ref_051018_09.fasta
samtools faidx NBW_ref_051018_10.fasta
samtools faidx NBW_ref_051018_11.fasta
samtools faidx NBW_ref_051018_12.fasta
samtools faidx NBW_ref_051018_13.fasta
samtools faidx NBW_ref_051018_14.fasta
samtools faidx NBW_ref_051018_15.fasta
samtools faidx NBW_ref_051018_16.fasta
samtools faidx NBW_ref_051018_17.fasta
samtools faidx NBW_ref_051018_18.fasta
samtools faidx NBW_ref_051018_19.fasta
samtools faidx NBW_ref_051018_20.fasta

# Edit opts.ctl file for each maker run
###### Maker round 1:
protein=/scratch/edegreef/ref_genome/annotation/nbw/model_proteins/Balaenoptera_musculus.mBalMus1.v2.pep.all.fa,/scratch/edegreef/ref_genome/annotation/nbw/model_proteins/Physeter_catodon.ASM283717v2.pep.all.fa,/scratch/edegreef/ref_genome/annotation/nbw/model_proteins/Bos_taurus.ARS-UCD1.2.pep.all.fa
model_org=all
protein2genome=1
cpus=20
mincontig=10000

# then run the maker array script (maker_run_array.sh)

# merge gff chunks
gff3_merge *.all.gff -o NBW_ref_051018_round1.gff  

# move all the maker files into a new folder (maker_round1)

# train snap in between maker rounds(train_snap.sh) and add the output to snaphmm in opts.ctl file


###### Maker round 2:
maker_gff=/scratch/edegreef/ref_genome/annotation/nbw/maker_round1/NBW_ref_051018_round1.gff
protein_pass=1
rm_pass=1
snaphmm=/scratch/edegreef/ref_genome/annotation/nbw/snap/round1/train_snap1.hmm
protein=#remove
model_org=#remove
repeat_protein=#remove
protein2genome=0
pred_stats=1

# maker_run_array.sh
# merge
gff3_merge *.all.gff -o NBW_ref_051018_round2.gff 

# move all the maker files into a new folder (maker_round2)

# train_snap.sh with gff file


###### Maker round 3:
maker_gff=/scratch/edegreef/ref_genome/annotation/nbw/maker_round2/NBW_ref_051018_round2.gff
snaphmm=/scratch/edegreef/ref_genome/annotation/nbw/snap/round2/train_snap2.hmm
augustus_species=human

# maker_run_array.sh
# merge
gff3_merge *.all.gff -o NBW_ref_051018_round3.gff 

# move all the maker files into a new folder (maker_round3)

# train_snap.sh with gff file

###### Maker round 4:
maker_gff=/scratch/edegreef/ref_genome/annotation/nbw/maker_round3/NBW_ref_051018_round3.gff
snaphmm=/scratch/edegreef/ref_genome/annotation/nbw/snap/round3/train_snap3.hmm
AED_threshold=0.5

# merge the output files
gff3_merge *.all.gff -o NBW_ref_051018_round4.gff  
cat *all.maker.proteins.fasta > NBW_ref_051018_round4.proteins.fasta
cat *all.maker.transcripts.fasta > NBW_ref_051018_round4.transcripts.fasta
cat *all.maker.augustus_masked.transcripts.fasta > NBW_ref_051018_round4.augustus_masked.transcripts.fasta
cat *all.maker.non_overlapping_ab_initio.proteins.fasta > NBW_ref_051018_round4.non_overlapping_ab_initio.proteins.fasta 
cat *all.maker.non_overlapping_ab_initio.transcripts.fasta > NBW_ref_051018_round4.non_overlapping_ab_initio.transcripts.fasta
cat *all.maker.snap_masked.proteins.fasta > NBW_ref_051018_round4.snap_masked.proteins.fasta
cat *all.maker.snap_masked.transcripts.fasta > NBW_ref_051018_round4.snap_masked.transcripts.fasta

# move all the maker files into a new folder (maker_round4)

# can use this for a quick gene count:
cat NBW_ref_051018_round4.gff | awk '{a[$3]++}END{for(k in a){print k,a[k]}}'

# keep in mind this genome has a lot of tiny contigs so the gene count will be inflated b/c there will be genes split across contigs.
# for whale we should get around ~20,000 genes

# AED = Annotation Edit Distance (0-1) 0=best, 1=bad
# look at AED distribution
perl AED_cdf_generator.pl -b 0.025 NBW_ref_051018_round4.gff > AED_distribution.txt

# Before doing final AED filter, going to do interprocan and blast+ steps