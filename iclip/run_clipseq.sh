#!/bin/bash

# Process iCLIP data
# 14th December 2022

# ==========
# nf-core/clipseq pipeline run
# ==========

# LOAD REQUIRED MODULES
ml purge
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref
GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip
cd $GITHUBDIR/iclip

## UPDATE PIPELINE
nextflow pull nf-core/clipseq -r 1.0.0

## RUN PIPELINE FOR HuR
nextflow run nf-core/clipseq -r 1.0.0 \
-resume \
-profile crick \
--input hur_samplesheet.csv \
--outdir /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_hur \
--umi_separator : \
--smrna_org human \
--fasta $REFDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz

## RUN PIPELINE FOR TDP-43
nextflow run nf-core/clipseq -r 1.0.0 \
-resume \
-profile crick \
--input tdp43_samplesheet.csv \
--outdir /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_tdp43 \
--umi_separator _ \
--smrna_org human \
--fasta $REFDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz

# ==========
# iCount peak calling
# ==========

conda activate icount-mini

mkdir -p /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/tmp
export ICOUNT_TMP_ROOT=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/tmp

# ==========
# HuR
# ==========

cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_hur

# Group crosslinks from replicates
zcat xlinks/HuR_rep1.xl.bed.gz xlinks/HuR_rep2.xl.bed.gz xlinks/HuR_rep3.xl.bed.gz | \
sort -k1,1 -k2,2n -k3,3n -k6,6 | \
bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | \
awk '{OFS="\t"}{print $1, $2, $3, ".", $5, $4}' | \
pigz > xlinks/HuR_merged.xl.bed.gz

# Run peak calling for merged crosslinks
iCount-Mini sigxls \
$REFDIR/icount_mini_utr3/gencode.v33.annotation.utr3seg.gtf.gz \
xlinks/HuR_merged.xl.bed.gz \
HuR_merged.3nt.sigxls.bed.gz \
--half_window 3 \
--fdr 0.05

iCount-Mini peaks \
--dist 3 \
--slop 1 \
xlinks/HuR_merged.xl.bed.gz \
HuR_merged.3nt.sigxls.bed.gz \
HuR_merged.3nt_3nt.peaks.bed.gz

rm iCount.log

# ==========
# TDP-43
# ==========

cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_tdp43

# Group crosslinks from replicates
zcat xlinks/TDP43_rep1.xl.bed.gz xlinks/TDP43_rep2.xl.bed.gz | \
sort -k1,1 -k2,2n -k3,3n -k6,6 | \
bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | \
awk '{OFS="\t"}{print $1, $2, $3, ".", $5, $4}' | \
pigz > xlinks/TDP43_merged.xl.bed.gz

# Run peak calling for merged crosslinks
iCount-Mini sigxls \
$REFDIR/icount_mini_utr3/gencode.v33.annotation.utr3seg.gtf.gz \
xlinks/TDP43_merged.xl.bed.gz \
TDP43_merged.10nt.sigxls.bed.gz \
--half_window 10 \
--fdr 0.05

iCount-Mini peaks \
--dist 10 \
--slop 1 \
xlinks/TDP43_merged.xl.bed.gz \
TDP43_merged.10nt.sigxls.bed.gz \
TDP43_merged.10nt_10nt.peaks.bed.gz

rm iCount.log