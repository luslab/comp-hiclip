#/bin/bash

# Obtain MFE and annotate structures for STAU1, PARIS and RIC-seq interactions

GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip
STAU1DIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/atlas
PARISDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/paris/results_paris/atlas_clusters 
RICSEQDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/ricseq/results_ricseq/atlas_clusters
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

cd $GITHUBDIR/merged_clustered

# ==========
# STAU1
# ==========

# Structures were-precomputed, so annotate them using forgi

./analyse_structure.R --input=$STAU1DIR/merged.atlas.clusters.tsv.gz --output=$STAU1DIR/merged.atlas.clusters.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --structure_annotation

# ==========
# PARIS
# ==========

# ./analyse_structure.R --input=$PARISDIR/all.atlas_clusters.gc.annotated.tsv.gz --output=$PARISDIR/paris.atlas.clusters.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --shuffled_mfe
./analyse_structure.R --input=$PARISDIR/all.atlas_clusters.gc.annotated.tsv.gz --output=$PARISDIR/paris.atlas.clusters.utr3.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --shuffled_mfe --threeutr

# ==========
# RIC-seq
# ==========

# ./analyse_structure.R --input=$RICSEQDIR/all.atlas_clusters.gc.annotated.tsv.gz --output=$RICSEQDIR/ricseq.atlas.clusters.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --shuffled_mfe
./analyse_structure.R --input=$RICSEQDIR/all.atlas_clusters.gc.annotated.tsv.gz --output=$RICSEQDIR/ricseq.atlas.clusters.utr3.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --shuffled_mfe --threeutr