#!/bin/sh

# Merge linker and nolinker nonhybrid data from Sugimoto et al., 2015 and run nf-core clipseq
# 21 May 2021


DATADIR_LINKER=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_linker/linker
DATADIR_NOLINKER=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nolinker/nonhybrids
RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

cd $RESULTSDIR

# Output files after merging nonhybrid reads from the linker and nolinker data
HIGH_MERGED_FASTQ=$RESULTSDIR/stau1_high.nonhybrid.merged.fastq.gz
LOW_MERGED_FASTQ=$RESULTSDIR/stau1_low.nonhybrid.merged.fastq.gz


# Merge linker and nolinker nonhybrids 

if [ -f "$HIGH_MERGED_FASTQ" ]; then
    echo "$HIGH_MERGED_FASTQ exists."
else 
    echo "$HIGH_MERGED_FASTQ does not exist. Merging.."
    cat $DATADIR_LINKER/LigPlusHigh.linker.nonhybrids.fastq.gz $DATADIR_NOLINKER/stau1_high.nonhybrid.fastq.gz > $HIGH_MERGED_FASTQ
fi


if [ -s "$LOW_MERGED_FASTQ" ]; then
    echo "$LOW_MERGED_FASTQ exists."
else 
    echo "$LOW_MERGED_FASTQ does not exist. Merging.."
    cat $DATADIR_LINKER/LigPlusLow.linker.nonhybrids.fastq.gz $DATADIR_NOLINKER/stau1_low.nonhybrid.fastq.gz > $LOW_MERGED_FASTQ
fi


# Download primary assembly fasta (not in REFDIR)

GENOME_FASTA=$RESULTSDIR/GRCh38.primary_assembly.genome.fa.gz

if [ -s "$GENOME_FASTA" ]; then
    echo "$GENOME_FASTA exists."
else 
    echo "$GENOME_FASTA does not exist. Downloading from GENCODE v33.."
    wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz
fi


# Run clipseq pipeline

# LOAD REQUIRED MODULES
ml purge
ml Nextflow/20.10.0
ml Singularity/3.6.4
ml Graphviz/2.38.0-foss-2016b


## iCount folder assignment
export ICOUNT_TMP_ROOT=/camp/lab/luscomben/home/users/iosubi/tmp
export NXF_SINGULARITY_CACHEDIR=/camp/lab/luscomben/home/shared/singularity

## UPDATE PIPELINE
nextflow pull nf-core/clipseq -r dev 

## RUN PIPELINE
nextflow run nf-core/clipseq -r dev \
-resume \
-profile crick \
--input nonhybrid_metadata.csv \
--outdir $RESULTSDIR \
--umi_separator _ \
--smrna_org human \
--fasta $RESULTSDIR/GRCh38.primary_assembly.genome.fa.gz \
--gtf $REFDIR/gencode.v33.annotation.gtf.gz \
--peakcaller icount --half_window 10 --merge_window 10 \
-N ira.iosub@crick.ac.uk

