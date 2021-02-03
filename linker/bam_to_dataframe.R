library(GenomicAlignments)
library(data.table)
library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) stop("Run command: bam_to_dataframe.R <input bam> <output prefix>")
if(!file.exists(args[1])) stop("Please provide an input bam file.")


bam_to_dataframe <- function(bam.file) {
  bam <- readGAlignments(bam.file,
                         use.names = TRUE)
  bam <- bam[njunc(bam) == 0] # filter out sj
  bam.df <- as.data.frame(bam)
  bam.df$seqnames <- as.character(bam.df$seqnames)
  bam.df$qname <- row.names(bam.df)
  bam.df <- bam.df %>% 
    dplyr::arrange(seqnames) %>%
    dplyr::filter((seqnames != "chrM") & !str_detect(seqnames,"KI"))
  # unique(bam.df$seqnames)
  
  bam.df$read <- gsub("^.*?\\.", "", bam.df$qname)
  bam.df$arm <- gsub("\\..*$", "", bam.df$qname)
  left.df <- bam.df %>%
    dplyr::filter(arm == "L")
  right.df <- bam.df %>%
    dplyr::filter(arm == "R")
  hybrids.df <- inner_join(left.df, right.df, by = "read")
  hybrids.df <- hybrids.df %>% arrange(read)
  print(paste0("There are ", nrow(hybrids.df), " valid hybrid reads (with both arms mapped)."))
  
  #create L and R columns
  hybrids.df <- hybrids.df %>%
    mutate(L_start = start.x, R_start = start.y,
           L_seqnames = seqnames.x, R_seqnames = seqnames.y,
           L_strand = as.character(strand.x), R_strand = as.character(strand.y),
           L_width = width.x, R_width = width.y,
           L_qname = qname.x, R_qname = qname.y) %>% 
    select(read, L_seqnames, L_start, L_strand, L_width, L_qname, R_seqnames, R_start, R_strand, R_width, R_qname)
  return(data.frame(hybrids.df))
}

reorient_hybrids <- function(hybrids.dt) {
  
  # First do starts
  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]
  
  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)
  
  setnames(incorrect.dt, renamed)
  
  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  
  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))
  
  # Then do subject (to make sure intergenics in same order)
  correct.dt <- reoriented.dt[L_seqnames <= R_seqnames]
  incorrect.dt <- reoriented.dt[L_seqnames > R_seqnames]
  
  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)
  
  setnames(incorrect.dt, renamed)
  
  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_subject <= reoriented.dt$R_subject))
  stopifnot(nrow(reoriented.dt) == nrow(hybrids.dt))
  
  return(reoriented.dt)
  
}


hybrids.df <- bam_to_dataframe(args[1])
hybrids.dt <- data.table(hybrids.df)
hybrids_reoriented.dt <- reorient_hybrids(hybrids.dt)

write.csv(hybrids_reoriented.dt, sep = "\t", paste0(args[2],"_hybrids.txt"))


