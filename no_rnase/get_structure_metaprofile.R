#!/usr/bin/env Rscript

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(tidyverse)
library(robustbase)
library(reshape)
library(stringr)
library(tictoc)

resizePeaks <- function(bedfiles.list, left = 100, right = 100) {
  
  w <- left + right +1  # width of interval: xl site + flanks
  
  grl <- GRangesList(lapply(bedfiles.list, import.bed))
  gr <- unlist(grl)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- dropSeqlevels(gr, c("chrM", "chrY"), pruning.mode = "coarse")
  gr <- resize(gr, width = 1, fix = "start") # resize the peaks, start of peak = 1
  gr <- unique(gr)  # keep unique xl positions
  gr <- resize(resize(gr, width = right+1, fix = "start"), width = w, fix = "end") # add +/- flanks
  return(gr)
  
}

getOverlaps <- function(gr, element.gr, left = 100, right = 100) {
  
  w <- left + right +1  # width of interval: xl site + flanks
  
  gr.nt <- unlist(tile(gr, width = 1))
  overlap <- findOverlaps(gr.nt, element.gr) 
  gr.nt$structure_prob <- as.numeric(NA)
  gr.nt[queryHits(overlap)]$structure_prob <- element.gr[subjectHits(overlap)]$score
  
  overlap.df <- as.data.frame(gr.nt)
  overlap.df$id <- 1 + seq(0, nrow(overlap.df) - 1) %/% w  # add peak ids #one ID every w
  overlap.df$id <- paste0("ID",overlap.df$id)
  overlap.df <- rowid_to_column(overlap.df, "nt_id") # record nt order in nt_id column
  
  plus <- overlap.df %>% dplyr::filter(overlap.df$strand == "+") # separate by strands to assign nt position
  plus$pos <- seq(1:w)
  minus <- overlap.df %>% dplyr::filter(overlap.df$strand == "-")
  minus$pos <- rev(seq(1:w))
  overlap.df <- rbind(plus, minus) %>% arrange(nt_id) %>% dplyr::select(-nt_id) # order by nt_id and remove the nt_id col
  #write.table(overlap.df, paste0(args[3],"_overlap.df.txt"), quote = FALSE, sep = "\t")
  
  pos.df <- unstack(overlap.df, structure_prob ~ pos) #reshape df and keep only the nt positions and scores
  pos.df <- rowid_to_column(pos.df, var = "id")
  pos.df$id <- paste0("ID",pos.df$id)
  rownames(pos.df) <-pos.df$id #make the id column the index
  pos.df <- select(pos.df, -id)
  colnames(pos.df) <- seq(-left, right)
  pos.df <- pos.df[rowSums(is.na(pos.df)) != ncol(pos.df), ] #remove peaks with all NAs (i.e. no overlaps found)
  return(pos.df)
  
}

getMetaprofile <- function(gr, element.gr, prefix) {
  
  structure.gr <- import.bed(element.gr, which = gr)
  transcript.list <- unique(unlist(gr$name))
  
  prefix <- paste0(prefix, "_",str_split(element.gr, pattern = ".bed")[[1]][1])
  structure.gr <- structure.gr[structure.gr$name %in% transcript.list]
  #stopifnot(length(element.gr) == length(unique(element.gr)))
  
  pos.df <- getOverlaps(gr, structure.gr)
  print(paste0(nrow(pos.df), " peaks analysed"))
  write.table(pos.df, paste0(prefix,"_prob.df.txt"), quote = FALSE, sep = "\t", row.names = TRUE) #this df = input for plotHeatmaps.R
  scores.mean <- as.data.frame(colMeans(data.matrix(pos.df), na.rm = TRUE))
  colnames(scores.mean) <- paste0(prefix,"_mean_prob")
  write.table(scores.mean, paste0(prefix,"_mean_prob.txt"), quote = FALSE, sep = "\t") #this df = input for plotting average profiles w plotProfiles.R
  
}


###########
# Obtain structure  metaprofiles
###########

tic()
# load the peak/xl bed files; they need to be annotated with ENST in a name metadata column
data.dir <- getwd()

files.list <- list.files(path = data.dir, pattern = "annot.bed.gz", full.names = FALSE)
files.list

peaks.gr <- resizePeaks(files.list)
message("Merging and resizing peaks...")

# load the structure probability bed files
structure.files.list <- list.files(path = data.dir, pattern = "threeutrs.rnaplfold", full.names = FALSE)
structure.files.list <- structure.files.list[str_detect(structure.files.list, pattern = ".bed")]
structure.files.list
message("Calculating profiles...")


lapply(structure.files.list, getMetaprofile, gr = peaks.gr, prefix = "stau1")
message("Metaprofiles generated.")
toc()


