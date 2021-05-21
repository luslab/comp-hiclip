#!/usr/bin/env Rscript 

library(rtracklayer)
library(GenomicFeatures)
library(stringr)


# ==========
# Define functions
# ==========

annotate_intervals <- function(bed.filename, max.tpm.transcripts.gr, longest.pc.transcripts.gr) {
  
  xl.gr <- import.bed(paste0(bed.filename))
  xl.gr <- keepStandardChromosomes(xl.gr, pruning.mode = "coarse")
  
  # Annotate based on most highly expressed transcript per gene
  tpm_overlap <- findOverlaps(xl.gr, max.tpm.transcripts.gr)
  xl.gr$name <- as.numeric(NA)
  xl.gr[queryHits(tpm_overlap)]$name <- max.tpm.transcripts.gr[subjectHits(tpm_overlap)]$name
  
  xl_mapped.gr <- xl.gr[!is.na(xl.gr$name)]
  
  # Annotate the remaining unmapped intervals with the longest protein-coding transcript per gene
  xl_unmapped.gr <- xl.gr[is.na(xl.gr$name)]
  
  if (!(length(xl_unmapped.gr) > 0)) {
    
    all.gr <- xl_mapped.gr
  } else {
    
    longest_pc_overlap <- findOverlaps(xl_unmapped.gr, longest.pc.transcripts.gr)
    xl_unmapped.gr[queryHits(longest_pc_overlap)]$name <- longest.pc.transcripts.gr[subjectHits(longest_pc_overlap)]$tx_name
    xl_unmapped.gr <- xl_unmapped.gr[!is.na(xl_unmapped.gr$name)]
    
    all.gr <- append(xl_mapped.gr, xl_unmapped.gr)
  }
  
  prefix <- str_remove_all(bed.filename, ".bed.gz")
  export.bed(all.gr, paste0(prefix, ".annot.bed.gz"))
  return(all.gr)
  
}


# ==========
# Files and parameters
# ==========

genomes.dir <- "/camp/lab/luscomben/home/users/iosubi/genomes"
data.dir <- "/camp/lab/luscomben/home/users/iosubi/projects/comp_hiclip/nonhybrids/rnafold_forgi_peaks_10nt_10nt"

files.list <- list.files(path = data.dir, pattern = "peaks.bed.gz", full.names = TRUE)
longest_pcoding.df <- read.csv(paste0(data.dir, "/longest_pcoding_transcripts.tsv"), sep = "\t")
max_tpm_transcripts.gr <- import.bed(paste0(data.dir, "/gencode_v33_max_tpm_transcripts.bed.gz"))

# Load TxDb object
TxDb  <- loadDb(paste0(genomes.dir,"/gencode_V33_txdb.sqlite"))
TxDb <- keepStandardChromosomes(TxDb, pruning.mode="coarse")

# ==========
# Annotate peaks
# ==========

transcripts.gr <- transcripts(TxDb)
transcripts.gr <- transcripts.gr[transcripts.gr$tx_name %in% unique(longest_pcoding.df$transcript_id)]


lapply(files.list, annotate_intervals, max.tpm.transcripts.gr = max_tpm_transcripts.gr, longest.pc.transcripts.gr = transcripts.gr)


