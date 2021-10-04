#!/bin/sh

DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/merged_clusters

cd $DATADIR

 # Obtain hybridisation enery and forgi annotations for all predicted short-range structures
./analyse_clusters_structure.R --clusters=$DATADIR/short_range_duplexes.tsv.gz --fasta=$DATADIR/../ref/GRCh38.gencode_v33.fa --output=$DATADIR/short_range_duplexes.mfe.tsv.gz --shuffled_mfe

# Get atlas, and append short range structures that were not clustered
./get_merged_atlas.R

# Obtain hybridisation enery and forgi annotations for collapsed duplexes
./analyse_clusters_structure.R --clusters=$DATADIR/stau1_atlas/merged_atlas.clusters.collapsed_plus_nonhybrids.tsv.gz --fasta=$DATADIR/../ref/GRCh38.gencode_v33.fa --output=$DATADIR/merged_atlas.clusters.collapsed.mfe.tsv.gz --clusters_only --shuffled_mfe


