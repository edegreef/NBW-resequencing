#!/bin/bash

# Running Satsuma Synteny on NBW genome using the blue whale genome as model (from NCBI)

# download
#wget ftp://ftp.broadinstitute.org/distribution/software/spines/satsuma-3.0.tar.gz
#tar xvzf satsuma-3.0.tar.gz

# Because genome is larger than 1.5GB, going to do chromosome by chromosome
# t=target genome fasta; q=query genome fasta; n=number of cpus; o=output directory

cd /home/degreefe/NBW/reference/satsuma

model_ref=/home/degreefe/NBW/reference/bluewhale/bluewhale_chr1.fasta 
NBW_genome=/home/degreefe/NBW/reference/Northern_bottlenose_whale_051018_shortLabel.fasta
out_dir=satsuma_out_chr1
block_display_out=BlockDisplaySatsuma_chr1

# run SatsumaSynteny
/home/degreefe/programs/satsuma-code-0/SatsumaSynteny -t $model_ref -q $NBW_genome -n 12 -o $out_dir

# convert satsuma summary output (satsuma_summary.chained.out) to MizBee format
/home/degreefe/programs/satsuma-code-0/BlockDisplaySatsuma -i $out_dir/satsuma_summary.chained.out -t $model_ref -q $NBW_genome > $block_display_out


# I ended up using the satsuma_summary.chained outputs into R, but below are some old notes for editing the blockdisplay outputs to run in MizBee:

# may need to adjust chr name from target/model ref so it says "chr1", or "chr2" etc, instead of "NC_......."
# running this part manually on command line for each chr b/c each one has different names
# can also use sed to extract full chr name on line 3: 
#sed -n '3p' BlockDisplaySatsuma_chr1 | awk '{print $1}'

# Then replace the name with "chr#"
#sed -i 's/NC_045785.1_Balaenoptera_musculus_isolate_JJ_BM4_2016_0621_chromosome_1,_mBalMus1.pri.v3,_whole_genome_shotgun_sequence/chr1/g' BlockDisplaySatsuma_chr1