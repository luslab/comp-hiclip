#!/usr/bin/env Rscript

library(rtracklayer)
library(tidyverse)
library(robustbase)
library(reshape)
library(stringr)
library(optparse)

# ==========
# Define functions
# ==========

resize_peaks <- function(bedfiles.list, left = 100, right = 100) {
  
  w <- left + right + 1  # width of interval: xl site + flanks
  
  grl <- GRangesList(lapply(bedfiles.list, import.bed))
  gr <- unlist(grl)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- dropSeqlevels(gr, c("chrM", "chrY"), pruning.mode = "coarse")

  gr <- resize(gr, width = 1, fix = "start") # resize the peaks, start of peak = 1
  gr <- unique(gr)  # keep unique xl positions
  gr <- resize(resize(gr, width = right + 1, fix = "start"), width = w, fix = "end") # add +/- flanks
  gr$id <- paste0("ID", seq(1, length(gr)))

  return(gr)
  
}

get_overlaps <- function(gr, element.gr, left = 100, right = 100) {
  
  w <- left + right + 1  # width of interval: xl site + flanks
  
  gr.nt <- unlist(tile(gr, width = 1))
  overlap <- findOverlaps(gr.nt, element.gr) 
  gr.nt$structure_prob <- as.numeric(NA)
  gr.nt[queryHits(overlap)]$structure_prob <- element.gr[subjectHits(overlap)]$score
  
  overlap.df <- as.data.frame(gr.nt)
  overlap.df$id <- 1 + seq(0, nrow(overlap.df) - 1) %/% w  # add peak ids #one ID every w nt
  overlap.df$id <- paste0("ID",overlap.df$id)
  
  stopifnot(unique(unique(overlap.df$id) == unique(gr$id)) == TRUE)
  
  overlap.df <- rowid_to_column(overlap.df, "nt_id") # record nt order in nt_id column
  
  plus <- overlap.df %>% dplyr::filter(overlap.df$strand == "+") # separate by strands to assign nt position
  plus$pos <- seq(1:w)
  minus <- overlap.df %>% dplyr::filter(overlap.df$strand == "-")
  minus$pos <- rev(seq(1:w))
  overlap.df <- rbind(plus, minus) %>%
    arrange(nt_id) %>%
    dplyr::select(-nt_id) # order by nt_id and remove the nt_id col
  
  # write.table(overlap.df, paste0(prefix,"_overlap.df.txt"), quote = FALSE, sep = "\t")
  
  pos.df <- unstack(overlap.df, structure_prob ~ pos) # reshape df and keep only the nt positions and scores
  pos.df <- rowid_to_column(pos.df, var = "id")
  pos.df$id <- paste0("ID", pos.df$id)
  rownames(pos.df) <- pos.df$id # make the id column the index
  pos.df <- select(pos.df, -id)
  colnames(pos.df) <- seq(-left, right)
  pos.df <- pos.df[rowSums(is.na(pos.df)) != ncol(pos.df), ] # remove peaks with all NAs (i.e. no overlaps found)
  
  return(pos.df)
  
}

get_metaprofile <- function(gr, element.gr) {
  
  structure.gr <- import.bed(element.gr, which = gr)
  transcript.list <- unique(unlist(gr$name))
  
  pref <- str_split(element.gr, pattern = ".bed")[[1]][1]
  structure.gr <- structure.gr[structure.gr$name %in% transcript.list]
  #stopifnot(length(element.gr) == length(unique(element.gr)))
  
  pos.df <- get_overlaps(gr, structure.gr, left = opt$left, right = opt$right)
  data.table::fwrite(pos.df, paste0(pref,"_prob.tsv.gz"), sep = "\t", row.names = TRUE)
  
}

# ==========
# Define options and params
# ==========

option_list <- list(make_option(c("-b", "--bed"), action = "store", type = "character", default=NA, help = "Comma separated transcript (ENST)-annotated peaks/xlink bed files"),
                    make_option(c("", "--prob"), action = "store", type = "character", default=NA, help = "Comma separated structure probability bed files"),
                    make_option(c("", "--prefix"), action = "store", type = "character", default=NA, help = "Prefix for output files"),
                    make_option(c("-l", "--left"), action = "store", type = "integer", default = 100, help = "Number of nt upstream [default: %default]"),
                    make_option(c("-r", "--right"), action = "store", type = "integer", default = 100, help = "Number of nt downstream [default: %default]"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

prefix <- opt$prefix

# ==========
# Obtain structure metaprofiles
# ==========

# Load the peak/xl bed files; they need to be annotated with ENST in a 'name' metadata column
files.list <- unlist(strsplit(opt$bed, ","))
peaks.gr <- resize_peaks(files.list, left = opt$left, right = opt$right)
message("Merging and resizing peaks...")

# Load the structure probability bed files
structure.files.list <- unlist(strsplit(opt$prob, ","))
structure.files.list <- structure.files.list[str_detect(structure.files.list, pattern = ".bed")]

message("Calculating profiles...")

lapply(structure.files.list, get_metaprofile, gr = peaks.gr)
message("Metaprofiles generated.")



