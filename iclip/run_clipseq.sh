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