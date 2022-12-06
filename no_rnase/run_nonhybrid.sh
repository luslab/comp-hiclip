#!/bin/sh

# Process nonhybrid data from Sugimoto et al., 2015
# 21 May 2021
# Updated 6th December 2022

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
cd $GITHUBDIR/no_rnase

## UPDATE PIPELINE
nextflow pull nf-core/clipseq -r 1.0.0

## RUN PIPELINE
nextflow run nf-core/clipseq -r 1.0.0 \
-resume \
-profile crick \
--input nornase.csv \
--outdir /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase \
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

cd /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase

# Group crosslinks from high and low RNase conditions
zcat xlinks/stau1_high.xl.bed.gz xlinks/stau1_low.xl.bed.gz | \
sort -k1,1 -k2,2n -k3,3n -k6,6 | \
bedtools groupby -i stdin -g 1,2,3,6 -c 5 -o sum | \
awk '{OFS="\t"}{print $1, $2, $3, ".", $5, $4}' | \
pigz > xlinks/stau1.xl.bed.gz

# Run peak calling for merged crosslinks
iCount-Mini sigxls \
$REFDIR/icount_mini_utr3/gencode.v33.annotation.utr3seg.gtf.gz \
xlinks/stau1.xl.bed.gz \
stau1.10nt.sigxls.bed.gz \
--half_window 10 \
--fdr 0.05

iCount-Mini peaks \
--dist 10 \
xlinks/stau1.xl.bed.gz \
stau1.10nt.sigxls.bed.gz \
stau1.10nt_10nt.peaks.bed.gz

rm iCount.log