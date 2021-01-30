#!/bin/sh

# Process PARIS hybrids
# A. M. Chakrabarti
# 28th January 2021

#SBATCH --job-name="comp_paris"
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --output=comp_paris-%A.out

ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b

nextflow pull amchakra/tosca

nextflow run amchakra/tosca \
-resume \
-profile crick,conda \
--org comp_hiclip \
--input paris.csv \
--outdir /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/paris/results \
--split_size 1000000 \
--umi_separator _ \
--dedup_method none


