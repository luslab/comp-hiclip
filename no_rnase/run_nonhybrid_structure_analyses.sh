#!/bin/sh

REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

RESULTSDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid
STRUCTURE_DIR=$RESULTSDIR/forgi
RNASEQ_DIR=$RESULTSDIR/rnaseq/results/salmon

cd $STRUCTURE_DIR

# Get highest TPM transcript per gene
./get_max_tpm_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --tpm=$RNASEQ_DIR/salmon.merged.transcript_tpm.tsv --out=$STRUCTURE_DIR/gencode.v33_max_tpm_transcripts.bed.gz

# Get longest protein-coding transcripts
./get_longest_pcoding_transcripts.R --gtf=$REFDIR/gencode.v33.annotation.gtf.gz --out=$STRUCTURE_DIR/gencode.v33.longest_pcoding_transcripts.tsv.gz


# Annotate peaks with transcript id
./annotate_peaks.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite

# RNAplfold (3'UTRs of STAU1-bound tranacripts)
./run_rnaplfold.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --txdb=$RESULTSDIR/gencode.v33.txdb.sqlite -prefix=stau1

# Generate RNAplfold unpaired probability metaprofile (-100nt and +100nt window)
./get_structure_metaprofile.R --bed=$RESULTSDIR/stau1.10nt_10nt.peaks.bed.gz --prob=stau1_threeutrs.rnaplfold.bed,stau1_threeutrs.rnaplfold.shuffled.bed --prefix=stau1

# Cluster RNAplfold probability profiles
./cluster_probability_profiles.R -p stau1_threeutrs.rnaplfold_prob.df.txt -s stau1_threeutrs.rnaplfold.shuffled_prob.df.txt -c 5

# Run RNAfold and annotate structures with forgi for all STAU1 peaks
./run_rnafold.R -b stau1_high.10nt_10nt.peaks.bed.gz.annot.bed.gz,stau1_low.10nt_10nt.peaks.bed.gz.annot.bed.gz -p stau1

# Analyse structure annotions for the 3'UTR peaks with predicted stems downstream according to RNAplfold (in the"stau1_threeutrs.rnaplfold_prob_clusters.df.txt" output from steps above)
./analyse_nonhybrid_duplexes.R -f stau1.forgi.tsv.gz -r stau1.rnafold.tsv.gz -d=15