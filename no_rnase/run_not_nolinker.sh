#!/bin/sh

# Process hiCLIP not nolinker hybrids
# A. M. Chakrabarti
# 28th January 2021

#SBATCH --job-name="comp_not_nolinker"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=comp_not_nolinker-%A.out

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

REFDIR=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/ref

nextflow pull nf-core/clipseq -r dev

nextflow run nf-core/clipseq -r dev \
-resume \
-profile crick \
--input not_nolinker.csv \
--outdir /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/nornase \
--smrna_org human \
--fasta $REFDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz \
--save_index true