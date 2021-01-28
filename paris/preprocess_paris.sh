#/bin/bash

# Preprocess PARIS data from Lu et al., 2016
# A. M. Chakrabarti
# 28th January 2021

conda activate comp-hiclip-dev

cd /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/paris

DATADIR=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/data/paris

#        543210987654321
#        123456789012345
# /5phos/DDDNNXXXXNNNNNNTACCCTTCGCTTCACACACAAG/iSp18/GGATCC/iSp18/TACTGAACCGCNNNNNN

# where the DDD indicates non-cytosine bases
# XXXX indicates the barcode
# N indicates any of the 4 bases
# Sp18 is spacer

# Barcode is in positions 7-10, so adapter is reverse complemented
# Therefore UMI format will be NNNNNNXXXXNNNNN

#                        543210987654321
#                        123456789012345
# iCLIP_ddRT_BC1: /5phos/DDDNNAACCNNNNAGATCGGAAGAGCGTCGTGAT/iSp18/GG ATCC/iSp18/TACTGAACCGC

for i in SRR2814763 SRR2814764 SRR2814765; do

    sbatch -t 12:00:00 -c 8 -o $i.log --wrap="umi_tools extract -p NNNNNNXXXXNNNNN -I $DATADIR/$i.fastq.gz -S $i.umi.fastq.gz && cutadapt -j 8 -m 12 -u 4 -a AGATCGGAAGAGC -o $i.preprocessed.fastq.gz $i.umi.fastq.gz"

done