#!/usr/bin/env Rscript

library(rtracklayer)
library(tidyverse)
library(robustbase)
library(reshape)
library(stringr)
library(optparse)
library(rslurm)
library(data.table)

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

get_overlaps <- function(transcript_id, gr, rnaplfold.gr, left = 100, right = 100) {
  
  w <- left + right + 1  # width of interval: xl site + flanks
  
  print(transcript_id)
  transcript_id <- unlist(transcript_id)
  #tx.gr <- gr[gr$name %in% transcript_id]
  tx.gr <- gr[(elementMetadata(gr)[, "name"] %in% transcript_id)]
  #rnaplfold.tx.gr <- rnaplfold.gr[rnaplfold.gr$name %in% transcript_id]
  rnaplfold.tx.gr <- rnaplfold.gr[(elementMetadata(rnaplfold.gr)[, "name"] %in% transcript_id)]
  gr.nt <- unlist(tile(tx.gr, width = 1))
  overlap <- findOverlaps(gr.nt, rnaplfold.tx.gr)
  gr.nt$structure_prob <- rep(as.numeric(NA),  times = length(gr.nt))
  gr.nt[queryHits(overlap)]$structure_prob <- rnaplfold.tx.gr[subjectHits(overlap)]$score
  gr.nt$id <- rep(tx.gr$id, each = w)
  return(gr.nt)
  
}


# ==========
# Define options and params
# ==========

option_list <- list(make_option(c("-b", "--bed"), action = "store", type = "character", default=NA, help = "(ENST)-annotated peaks/xlink bed files"),
                    make_option(c("", "--prob"), action = "store", type = "character", default=NA, help = "Structure probability bed file"),
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

message("Resizing peaks...")
peaks.gr <- resize_peaks(opt$bed, left = opt$left, right = opt$right)
w <- opt$left + opt$right + 1

# Load the structure probability bed files

rnaplfold.file <- opt$prob
structure.gr <- import.bed(rnaplfold.file, which = peaks.gr)
transcript.list <- intersect(peaks.gr$name, structure.gr$name)

structure.gr <- structure.gr[(elementMetadata(structure.gr)[, "name"] %in% transcript.list)]
peaks.gr <- peaks.gr[(elementMetadata(peaks.gr)[, "name"] %in% transcript.list)]

transcript.ls <- as.list(c(transcript.list))

sjob <- rslurm::slurm_map(transcript.ls,
                          get_overlaps, left = opt$left, right = opt$right,
                          gr=peaks.gr, rnaplfold.gr=structure.gr,
                          nodes = 100, cpus_per_node = 1,
                          jobname = "RNAplfold_metaprofiles",
                          slurm_options = list(time = "24:00:00"),
                          submit = TRUE)

Sys.sleep(60)
message("Profiles generating..")

# check job is finished
status <- FALSE
while(status == FALSE) {
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left, = all jobs are finished
  Sys.sleep(60)
}

output <- get_slurm_out(sjob, outtype = "raw")
saveRDS(output, paste0(opt$prefix,"_threeutrs.rnaplfold.profiles.rds"))

cleanup_files(sjob)

# ==========
# Convert tiled GRanges to profile dataframe
# ==========

rnaplfold.ls <- readRDS("threeutrs.rnaplfold.profiles.rds") # load the GRanges List

rnaplfold.grl <- GRangesList(rnaplfold.ls)
rnaplfold.gr <- unlist(rnaplfold.grl)

overlap.df <- as.data.frame(rnaplfold.gr)
overlap.df <- rowid_to_column(overlap.df, "nt_id") # record nt order in nt_id column

plus <- overlap.df %>% dplyr::filter(overlap.df$strand == "+") # separate by strands to assign nt position

plus$pos <- seq(1:w)
minus <- overlap.df %>% dplyr::filter(overlap.df$strand == "-")
minus$pos <- rev(seq(1:w))

overlap.df <- rbind(plus, minus) %>%
  arrange(nt_id) %>%
  dplyr::select(-nt_id) # order by nt_id and remove the nt_id col

id.list <- unique(overlap.df$id)
pos.df <- unstack(overlap.df, structure_prob ~ pos) # reshape df and keep only the nt positions and scores
pos.df$id <- id.list

rownames(pos.df) <- pos.df$id # make the id column the index
pos.df <- dplyr::select(pos.df, -id)
colnames(pos.df) <- seq(-opt$left, opt$right)

pos.df <- pos.df[rowSums(is.na(pos.df)) != ncol(pos.df), ] # remove peaks with all NAs (i.e. no overlaps found)
fwrite(pos.df, paste0(opt$prefix,"_prob_profiles.tsv.gz"), sep = "\t", row.names = TRUE)

if (nrow(pos.df) > 0) {message("Profiles generated.")}



