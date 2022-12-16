#!/bin/sh

GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip/no_rnase
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref
RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase
RNASEQDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/results_demux/salmon

HURDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_hur
TARDPDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/iclip/results_tdp43

cd $RESULTSDIR

# # Get highest TPM transcript per gene
# $GITHUBDIR/get_max_tpm_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --tpm=$RNASEQDIR/salmon.merged.transcript_tpm.tsv --out=$RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz

# # Get longest protein-coding transcripts
# $GITHUBDIR/get_longest_pcoding_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --out=$RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# # ==========
# # STAU1 analyses
# # ==========

# # Annotate peaks with transcript id
# $GITHUBDIR/annotate_peaks.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --max_tpm_bed $RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz --longest_tx $RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# # RNAplfold (3'UTRs of STAU1-bound transcripts)
# $GITHUBDIR/run_rnaplfold.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/stau1 --shuffle

# # Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
# $GITHUBDIR/get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/stau1

# # Profile for shuffled control
# $GITHUBDIR/get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/stau1_shuffled

# ==========
# HuR analyses
# ==========

# Annotate peaks with transcript id
$GITHUBDIR/annotate_peaks.R --bed=$HURDIR/HuR_merged.3nt_3nt.peaks.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --max_tpm_bed $RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz --longest_tx $RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# RNAplfold (3'UTRs of HuR-bound transcripts)
$GITHUBDIR/run_rnaplfold.R --bed=$HURDIR/HuR_merged.3nt_3nt.peaks.annot.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/hur --shuffle

# Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
$GITHUBDIR/get_structure_metaprofile.R --bed=$HURDIR/HuR_merged.3nt_3nt.peaks.annot.bed.gz --prob=$RESULTSDIR/HuR_merged_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/hur

# Profile for shuffled control
$GITHUBDIR/get_structure_metaprofile.R --bed=$HURDIR/HuR_merged.3nt_3nt.peaks.annot.bed.gz --prob=$RESULTSDIR/HuR_merged_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/hur_shuffled

# ==========
# TDP43 analyses
# ==========

# Annotate peaks with transcript id
$GITHUBDIR/annotate_peaks.R --bed=$TARDPDIR/TDP43_merged.10nt_10nt.peaks.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --max_tpm_bed $RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz --longest_tx $RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# RNAplfold (3'UTRs of TDP43-bound transcripts)
$GITHUBDIR/run_rnaplfold.R --bed=$TARDPDIR/TDP43_merged.10nt_10nt.peaks.annot.bed.gz --txdb=$REFDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/tdp43 --shuffle

# Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
$GITHUBDIR/get_structure_metaprofile.R --bed=$TARDPDIR/TDP43_merged.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/TDP43_merged_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/tdp43

# Profile for shuffled control
$GITHUBDIR/get_structure_metaprofile.R --bed=$TARDPDIR/TDP43_merged.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/TDP43_merged_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/tdp43_shuffled
