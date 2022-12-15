#!/bin/sh

GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip/no_rnase
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref
RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase
RNASEQ_DIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/rnaseq/results_demux/salmon

HURDIR=$RESULTSDIR/hur/flora

cd $RESULTSDIR

# Get highest TPM transcript per gene
$GITHUBDIR/get_max_tpm_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --tpm=$RNASEQ_DIR/salmon.merged.transcript_tpm.tsv --out=$RESULTSDIR/gencode.v33_max_tpm_transcripts.bed.gz

# Get longest protein-coding transcripts
$GITHUBDIR/get_longest_pcoding_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --out=$RESULTSDIR/gencode.v33.longest_pcoding_transcripts.tsv.gz

# ## STAU1 analyses

# # Annotate peaks with transcript id
# ./annotate_peaks.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite

# # RNAplfold (3'UTRs of STAU1-bound transcripts)
# ./run_rnaplfold.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/stau1 --shuffle

# # Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
# #./get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.bed.gz,$RESULTSDIR/stau1_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/stau1
# ./get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/stau1

# # Profile for shuffled control
# ./get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/stau1_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/stau1_shuffled


# ## HUR analyses

# # Annotate peaks with transcript id
# ./annotate_peaks.R --bed=$HURDIR/ELAVL1_merged.3nt_3nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite

# # RNAplfold (3'UTRs of STAU1-bound transcripts)
# ./run_rnaplfold.R --bed=$HURDIR/ELAVL1_merged.3nt_3nt.peaks.annot.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/ELAVL1 --shuffle

# # Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
# ./get_structure_metaprofile.R --bed=$HURDIR/ELAVL1_merged.3nt_3nt.peaks.annot.bed.gz --prob=$RESULTSDIR/ELAVL1_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/ELAVL1

# # Profile for shuffled control
# ./get_structure_metaprofile.R --bed=$HURDIR/ELAVL1_merged.3nt_3nt.peaks.annot.bed.gz --prob=$RESULTSDIR/ELAVL1_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/ELAVL1_shuffled


# ## TDP43 analyses

# # Annotate peaks with transcript id
# ./annotate_peaks.R --bed=$TARDP_DIR/tardp.10nt_10nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite

# # RNAplfold (3'UTRs of STAU1-bound transcripts)
# ./run_rnaplfold.R --bed=$TARDP_DIR/tardp.10nt_10nt.peaks.annot.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite --prefix=$RESULTSDIR/TARDP --shuffle

# # Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
# ./get_structure_metaprofile.R --bed=$TARDP_DIR/tardp.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/TARDP_threeutrs.rnaplfold.bed.gz --prefix=$RESULTSDIR/TARDP

# # Profile for shuffled control
# ./get_structure_metaprofile.R --bed=$TARDP_DIR/tardp.10nt_10nt.peaks.annot.bed.gz --prob=$RESULTSDIR/TARDP_threeutrs.rnaplfold.shuffled.bed.gz --prefix=$RESULTSDIR/TARDP_shuffled
