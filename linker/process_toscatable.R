# Script to process deduplicated linker hybrids
# A. M. Chakrabarti
# 5th December 2022

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(rtracklayer))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 4) stop("Run command: process_toscatable.R <input hybrids table> <transcript gtf> <regions gtf> <output prefix>")
if(!file.exists(args[1])) stop("Please provide an input hybrids table.")

genes.gr <- rtracklayer::import.gff2(args[2])
regions.gr <- rtracklayer::import.gff2(args[3])

hybrids.dt <- fread(args[1])
hybrids.dt <- reorient_hybrids(hybrids.dt)
hybrids.gc.dt <- convert_coordinates(hybrids.dt, genes.gr)
hybrids.gc.annotated.dt <- annotate_hybrids(hybrids.gc.dt, regions.gr)

fwrite(hybrids.gc.annotated.dt, paste0(args[4], ".hybrids.gc.annotated.tsv.gz"), sep = "\t")