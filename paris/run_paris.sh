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

nextflow pull amchakra/tosca -r star_options

nextflow run amchakra/tosca -r star_options \
-resume \
-profile crick,conda \
-N anob.chakrabarti@crick.ac.uk \
--org comp_hiclip \
--input paris.csv \
--outdir /camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/paris \
--umi_separator _ \
--star_args '--limitOutSJcollapsed 5000000' \
--split_size 1000000
