#!/usr/bin/env Rscript 

library(rtracklayer)
library(GenomicFeatures)
library(stringr)
library(optparse)



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


# =========
# Options and paths
# =========

option_list <- list(make_option(c("", "--bed"), action = "store", type = "character", default=NA, help = "Peaks bed file"),
                    make_option(c("", "--txdb"), action = "store", type = "character", default=NA, help = "TxDb object"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/forgi"



# ==========
# Files and parameters
# ==========


#files.list <- list.files(path = opt$bed, pattern = "peaks.bed.gz", full.names = TRUE)
files.list <- opt$bed

longest_pcoding.df <- read.csv(paste0(data.dir, "/gencode.v33.longest_pcoding_transcripts.tsv.gz"), sep = "\t")
max_tpm_transcripts.gr <- import.bed(paste0(data.dir, "/gencode.v33_max_tpm_transcripts.bed.gz"))

# Load TxDb object
TxDb  <- loadDb(opt$txdb)
TxDb <- keepStandardChromosomes(TxDb, pruning.mode="coarse")

# ==========
# Annotate peaks
# ==========

transcripts.gr <- transcripts(TxDb)
transcripts.gr <- transcripts.gr[transcripts.gr$tx_name %in% unique(longest_pcoding.df$transcript_id)]


lapply(files.list, annotate_intervals, max.tpm.transcripts.gr = max_tpm_transcripts.gr, longest.pc.transcripts.gr = transcripts.gr)


