#!/bin/sh

REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid

RNASEQ_DIR=$RESULTSDIR/rnaseq/results_demux/salmon

cd $RESULTSDIR

# Get highest TPM transcript per gene
#./get_max_tpm_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --tpm=$RNASEQ_DIR/salmon.merged.transcript_tpm.tsv --out=$RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz

# Get longest protein-coding transcripts
#./get_longest_pcoding_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --out=$RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# Annotate peaks with transcript id
./annotate_peaks.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite

# RNAplfold (3'UTRs of STAU1-bound transcripts)
./run_rnaplfold.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/stau1 --shuffle

# Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
./get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.bed,$RESULTSDIR/stau1_threeutrs.rnaplfold.shuffled.bed --prefix=$RESULTSDIR/stau1

