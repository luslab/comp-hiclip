#!/bin/sh

# Process hiCLIP nolinker hybrids
# A. M. Chakrabarti
# 28th January 2021
# 2nd December 2022

#SBATCH --job-name="comp_hiclip"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=comp_hiclip-%A.out

ml purge
ml Nextflow/21.10.3
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref
# Regenerate index as comp-hiclip-dev is 2.7.7a and tosca-1.0.0 is 2.7.4a
# cd $REFDIR && STAR --runMode genomeGenerate --runThreadN 8 --genomeDir STAR_GRCh38_GencodeV33 --genomeFastaFiles GRCh38.gencode_v33.fa 
GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip
cd $GITHUBDIR/no_linker

nextflow pull amchakra/tosca -r v1.0.0

nextflow run amchakra/tosca -r v1.0.0 \
-resume \
-profile crick \
--input nolinker.csv \
--outdir /camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nolinker \
--genome_fai $REFDIR/GRCh38.primary_assembly.genome.fa.fai \
--star_genome $REFDIR/STAR_GRCh38_GencodeV33 \
--transcript_fa $REFDIR/GRCh38.gencode_v33.fa \
--transcript_fai $REFDIR/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $REFDIR/GRCh38.gencode_v33.tx.gtf.gz \
--regions_gtf $REFDIR/icount_mini_utr3/regions.gtf.gz \
--percent_overlap 0.5 \
--analyse_structures true \
--clusters_only true \
--shuffled_energies true