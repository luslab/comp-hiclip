#!/bin/sh

# Process hiCLIP not nolinker hybrids
# A. M. Chakrabarti
# 28th January 2021

#SBATCH --job-name="comp_not_nolinker"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=comp_not_nolinker-%A.out

ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

REFDIR=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/ref

nextflow pull nf-core/clipseq -r dev

# nextflow run nf-core/clipseq -r dev \
# -resume \
# -profile crick \
# --input not_nolinker.csv \
# --outdir /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/nornase \
# --smrna_org human \
# --fasta $REFDIR/GRCh38.primary_assembly.genome.fa.gz \
# --gtf $REFDIR/gencode.v33.annotation.gtf.gz \
# --umi_separator _ \
# --save_index true

nextflow run nf-core/clipseq -r dev \
-resume \
-profile crick \
--input not_nolinker.csv \
--outdir /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/nornase \
--smrna_org human \
--fasta $REFDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz \
--umi_separator _ \
--star_index /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/nornase/STAR_index