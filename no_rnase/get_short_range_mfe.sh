#!/bin/bash

DATADIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase
GITHUBDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip/merged_clustered
REFDIR=/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref

$GITHUBDIR/analyse_structure.R --input=$DATADIR/short_range_duplexes.tsv.gz --output=$DATADIR/short_range_duplexes.mfe.tsv.gz --fasta=$REFDIR/GRCh38.gencode_v33.fa --shuffled_mfe 

