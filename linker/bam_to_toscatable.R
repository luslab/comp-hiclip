suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(toscatools))

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) stop("Run command: bam_to_toscatable.R <input bam> <transcript gtf> <regions gtf> <output prefix>")
if(!file.exists(args[1])) stop("Please provide an input bam file.")

bam_to_toscatable <- function(bam.file) {

  bam <- readGAlignments(bam.file,
                         use.names = TRUE)
  bam <- bam[njunc(bam) == 0] # filter out sj
  bam.df <- as.data.frame(bam)
  bam.df$seqnames <- as.character(bam.df$seqnames)
  bam.df$qname <- row.names(bam.df)
  bam.df <- bam.df %>% 
    dplyr::arrange(seqnames) %>%
    dplyr::filter((seqnames != "chrM") & !str_detect(seqnames, "KI"))
  # unique(bam.df$seqnames)
  
  bam.df$name <- gsub("^.*?\\.", "", bam.df$qname)
  bam.df$arm <- gsub("\\..*$", "", bam.df$qname)
  left.df <- bam.df %>%
    dplyr::filter(arm == "L")
  right.df <- bam.df %>%
    dplyr::filter(arm == "R")
  hybrids.df <- inner_join(left.df, right.df, by = "name")
  hybrids.df <- hybrids.df %>% arrange(name)
  print(paste0("There are ", nrow(hybrids.df), " valid hybrid reads (with both arms mapped)."))
  
  #create L and R columns
  hybrids.df <- hybrids.df %>%
    mutate(L_start = start.x, R_start = start.y,
           L_seqnames = seqnames.x, R_seqnames = seqnames.y,
           L_strand = as.character(strand.x), R_strand = as.character(strand.y),
           L_width = width.x, R_width = width.y,
           L_qname = qname.x, R_qname = qname.y) %>% 
    select(name, L_seqnames, L_start, L_strand, L_width, L_qname, R_seqnames, R_start, R_strand, R_width, R_qname)
  
  # Calculations and annotations for downstream
  hybrids.dt <- data.table(hybrids.df)
  hybrids.dt[, `:=` (L_end = L_start + L_width - 1,
                     R_end = R_start + R_width - 1)]
  hybrids.dt[, type := ifelse(L_seqnames == R_seqnames, "intragenic", "intergenic")]
  hybrids.dt[, orientation := "linker"]
  hybrids.dt[, hybrid_selection := "single"]

  return(hybrids.dt)

}

hybrids.dt <- bam_to_toscatable(args[1])
# hybrids_reoriented.dt <- reorient_hybrids(hybrids.dt) # Do not reorient before deduplication

fwrite(hybrids.dt, sep = "\t", file = paste0(args[2], ".hybrids.tsv.gz"))