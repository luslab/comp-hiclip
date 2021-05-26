#!/bin/sh

# Process hiCLIP nolinker hybrids
# A. M. Chakrabarti
# 28th January 2021

#SBATCH --job-name="comp_hiclip"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=comp_hiclip-%A.out

ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

PROJ=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip

mkdir -p $PROJ/results_nolinker
mkdir -p /camp/lab/luscomben/scratch/chakraa2/projects/comp_hiclip
export NXF_WORK=/camp/lab/luscomben/scratch/chakraa2/projects/comp_hiclip

nextflow pull amchakra/tosca

nextflow run amchakra/tosca \
-resume \
-profile crick,conda \
--org comp_hiclip \
--input nolinker.csv \
--outdir $PROJ/results_nolinker \
--umi_separator _ \
--percent_overlap 0.5 \
--analyse_structure true \
--clusters_only true \
--shuffled_mfe true \
--transcript_fa $PROJ/ref/GRCh38.gencode_v33.fa \
--transcript_fai $PROJ/ref/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $PROJ/ref/GRCh38.gencode_v33.tx.gtf.gz


nextflow run amchakra/tosca \
-resume \
-profile conda \
--org comp_hiclip \
--input silico.csv \
--outdir $PROJ/results_silico \
--umi_separator _ \
--percent_overlap 0.5 \
--atlas false \
--transcript_fa $PROJ/ref/GRCh38.gencode_v33.fa \
--transcript_fai $PROJ/ref/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $PROJ/ref/GRCh38.gencode_v33.tx.gtf.gz


mkdir -p $PROJ/results_nolinker_qc

nextflow pull amchakra/tosca -r multiqc

nextflow run amchakra/tosca -r multiqc \
-resume \
-profile conda \
--org comp_hiclip \
--input nolinker.csv \
--outdir $PROJ/results_nolinker_qc \
--umi_separator _ \
--percent_overlap 0.5 \
--atlas false \
--transcript_fa $PROJ/ref/GRCh38.gencode_v33.fa \
--transcript_fai $PROJ/ref/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $PROJ/ref/GRCh38.gencode_v33.tx.gtf.gz


ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

PROJ=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip

nextflow run amchakra/tosca -r multiqc \
-resume \
-profile conda \
--org comp_hiclip \
--input silico.csv \
--umi_separator _ \
--percent_overlap 0.5 \
--atlas false \
--transcript_fa $PROJ/ref/GRCh38.gencode_v33.fa \
--transcript_fai $PROJ/ref/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $PROJ/ref/GRCh38.gencode_v33.tx.gtf.gz


ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

cd /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/silico/
PROJ=/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip

nextflow run /camp/home/chakraa2/.nextflow/assets/amchakra/tosca/main.nf \
-resume \
-profile conda \
--org comp_hiclip \
--input silico.csv \
--umi_separator _ \
--percent_overlap 0.5 \
--atlas false \
--transcript_fa $PROJ/ref/GRCh38.gencode_v33.fa \
--transcript_fai $PROJ/ref/GRCh38.gencode_v33.fa.fai \
--transcript_gtf $PROJ/ref/GRCh38.gencode_v33.tx.gtf.gz