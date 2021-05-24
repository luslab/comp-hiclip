#!/usr/bin/env Rscript

library(rtracklayer)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(stringr)
library(rslurm)
library(tictoc)
library(tidyverse)
#library(primavera)
library(optparse)


# ==========
# Define functions
# ==========

run_rnaplfold <- function(sequence, args =  c("-u 1", "-W 100", paste0("< ",sequence))) {
  
  system2(command = "RNAplfold", args = args, stdout = TRUE)
}


get_rnaplfold_probability <- function(id, sequence) {
  
  temp.fasta <- paste0(id,".fasta")
  fasta.ds <- DNAStringSet(sequence)
  names(fasta.ds) <- id
  writeXStringSet(fasta.ds, filepath = paste0(id,".fasta"))
  
  run_rnaplfold(temp.fasta)
  rnaplfold.df <- read.csv(paste0(id,"_lunp"), comment.char = "#",sep="\t", header = FALSE, row.names = NULL, skip = 2, col.names = c("pos","prob"))
  
  nt <- unlist(tile(grl[[id]], width = 1))
  stopifnot(nrow(rnaplfold.df) == length(nt)) # check that size of the transcript from RNAplfold output is identical to GRanges annotation
  if(unique(strand(nt)) == "+") {
    nt$score <- as.numeric(rnaplfold.df$prob)
    nt$name <- noquote(id)
  } else {
    nt$score <- rev(as.numeric(rnaplfold.df$prob)) 
    nt$name <- noquote(id)
  }
  
  invisible(file.remove(temp.fasta))
  return(nt)
}

shuffle_sequence <- function(sequence, number = 1, klet = 2, seed = 42) {
  system(paste0("ushuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)
}


get_rnaplfold_shuffled <- function(id, sequence) {
  
  fasta.df <- data.frame(id = id, sequence = sequence)
  shuffled <- sapply(seq_along(id), function(i) shuffle_sequence(fasta.df[i,]$sequence, number = 100, klet = 2))
  
  # Write a multi fasta per ID containg the shuffled sequences
  shuff_fasta <- sapply(1:ncol(shuffled), function(i) DNAStringSet(shuffled[,i]))
  shuff_fasta <- unlist(DNAStringSetList(shuff_fasta))
  
  names(shuff_fasta) <- sapply(1:length(shuff_fasta), function(i) paste0(id, "_",i))
  
  writeXStringSet(shuff_fasta, filepath = paste0(id,".fa"))
  run_rnaplfold(paste0(id,".fa"))
  
  # Read all output files that start with id
  rnaplfold.out.ls <- vector(mode = "list", length = 1) # create empty list/ re-initialise list
  rnaplfold.out.ls <- list.files(pattern = paste0("^",id))
  rnaplfold.out.ls <- rnaplfold.out.ls[str_detect(rnaplfold.out.ls, pattern = "_lunp")]

  rnaplfold.df.ls <- lapply(rnaplfold.out.ls, read.csv, comment.char = "#",sep="\t", header = FALSE, row.names = NULL, skip = 2, colClasses = as.numeric(), col.names = c("pos","prob"))
  rnaplfold.df <- bind_cols(rnaplfold.df.ls)
  rnaplfold.df <- rnaplfold.df[!duplicated(as.list(rnaplfold.df))]
  rnaplfold.df <- rnaplfold.df %>% remove_rownames %>% column_to_rownames(var="pos...1")
  
  rnaplfold.df$prob <- rowMeans(rnaplfold.df, na.rm=TRUE)
  rnaplfold.df$prob_sd <- apply(rnaplfold.df, 1, function(x) { sd(x, na.rm = TRUE) })
  rnaplfold.df <- rnaplfold.df %>%
    select(prob, prob_sd)
  rnaplfold.out.ls <- list()
  
  nt <- unlist(tile(threeutr.grl[[id]], width = 1)) # make GRanges for each nucleotide position
  stopifnot(nrow(rnaplfold.df) == length(nt)) # check that size of the transcript from RNAplfold output is identical to GRanges annotation
  
  if(unique(strand(nt)) == "+") {
    nt$score <- as.numeric(rnaplfold.df$prob)
    nt$sd <- as.numeric(rnaplfold.df$prob_sd)
    nt$name <- noquote(id)
  } else {
    
    nt$score <- rev(as.numeric(rnaplfold.df$prob)) 
    nt$sd <- rev(as.numeric(rnaplfold.df$prob_sd))
    nt$name <- noquote(id)
  }
  
  invisible(file.remove(paste0(id,".fa")))
  #invisible(file.remove(rnaplfold.out.ls))
  return(nt)
}


# ==========
# Define options and params
# ==========

option_list <- list(make_option(c("-b", "--bed"), action = "store", type = "character", default=NA, help = "Comma separated trancript(ENST) annotated bed files list"),
                    make_option(c("-p", "--prefix"), action = "store", type = "character", default=NA, help = "Prefix for output files"),
                    make_option(c("-t", "--txdb"), action = "store", type = "character", default=NA, help = "TxDb (optional, but required if --threeutrs not provided)"),
                    make_option(c("", "--threeutrs"), action = "store", type = "character", default="threeutrs.grl.rds", help = "RDS of 3'UTRs GRangesList (optional)"),
                    make_option(c("-s", "--shuffle"), action = "store_true",  default=FALSE, type = "character", help = "Generate shuffled control"),
                    make_option(c("-n", "--nodes"), action = "store", type = "integer", default = 100, help = "Number of nodes to allocate [default: %default]"))

opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# genome and annotations
# genomes.dir <- "/camp/lab/luscomben/home/users/iosubi/genomes/"
# txdb <- paste0(genomes.dir,"gencode_V33_txdb.sqlite")
txdb <- opt$txdb
Hsapiens <- BSgenome.Hsapiens.UCSC.hg38 # load hg38 genome

# 3UTR plfold database (previous RNAplfold run on all 3'UTRs annotated in Gencode V33)
threeutr.plfold.db <- "threeutrs.rnaplfold.bed"

#inputs and outputs
#files.list <- unlist(strsplit(opt$bed, ","))
files.list <- opt$bed
prefix <- opt$prefix # RBP name of interest, will be added as prefix to output filenames

# ==========
# Obtain/load 3UTR coordinates
# ==========

if (!file.exists(opt$threeutrs)) {
  
  message("Getting 3UTR coordinates...")
  
  if(!file.exists(txdb)) {
    stop("Error: cannot extract coordinates, please provide a TxDb file")
    
  } else {
    
    gencode_V33.txdb <- loadDb(txdb) # load TxDb object
    gencode_V33.txdb <- keepStandardChromosomes(gencode_V33.txdb, pruning.mode="coarse")
    threeutr.grl <- threeUTRsByTranscript(gencode_V33.txdb, use.names=TRUE)
  }
  
} else {
  
  threeutr.grl <- readRDS(opt$threeutrs) # load the GRanges List of 3UTRs
}

# ==========
# Extract 3UTR fasta sequences
# ==========

# Load annotated bed files, join them, get transcript names

bed.grl <- GRangesList(lapply(files.list, import.bed))
bed.gr <- unlist(bed.grl)
bed.gr <- unique(bed.gr)
bed.gr <- keepStandardChromosomes(bed.gr, pruning.mode = "coarse")
transcript.ls <- unique(bed.gr$name)



# select only the transcripts from the peaks annotation, and get fasta
transcript.ls <- intersect(names(threeutr.grl), transcript.ls)
grl <- threeutr.grl[transcript.ls]

message("Extracting 3UTR sequences...")
fasta <- extractTranscriptSeqs(Hsapiens, grl)
fasta.df <- data.frame(id = names(fasta), sequence = as.character(fasta))

# ==========
# Run RNAplfold
# ==========

# only run rnaplfold if a 3UTR transcriptome-wide rnaplfold probability database isn't already available
# if unavailable, run rnaplfold only for the 3UTRs bound by the RBP of interest

tic()

if (!file.exists(threeutr.plfold.db)) {
  
  
  sjob <- slurm_apply(get_rnaplfold_probability, fasta.df, jobname = "RNAplfold", add_objects = c("grl","run_rnaplfold"),
                      nodes = opt$nodes, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), 
                      submit = TRUE)
  Sys.sleep(60)
  message("RNAplfold is running..")
  
  # check job is finished
  status <- FALSE
  while(status == FALSE) {
    squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left, = all jobs are finished
    Sys.sleep(60)
  }

  output <- get_slurm_out(sjob, outtype = "raw")
  saveRDS(output, "threeutrs.rnaplfold.grl.rds")
  
  #cleanup_files(sjob)

  rnaplfold.ls <- readRDS("threeutrs.rnaplfold.grl.rds") # load the GRanges List
  
  print(length(rnaplfold.ls))
  
  rnaplfold.grl <- GRangesList(rnaplfold.ls)
  export.bed(unlist(rnaplfold.grl), paste0(prefix, "_threeutrs.rnaplfold.bed"), format = "BED")
  message("Probability scores (BED) exported.")
  
} else {
  
  # filter 3UTR db for transcripts of interest
  message("Filtering 3UTR db for transcripts of interest...")
  rnaplfold.gr <- import.bed(threeutr.plfold.db)
  rnaplfold.gr <- rnaplfold.gr[rnaplfold.gr$name %in% transcript.ls]
  export.bed(rnaplfold.gr, paste0(prefix, "_threeutrs.rnaplfold.bed"), format = "BED")
  
}

toc()



# ==========
# Run RNAplfold for shuffled sequences
# ==========

tic()

if (opt$shuffle) {
  sjob <- slurm_apply(get_rnaplfold_shuffled, fasta.df, jobname = "RNAplfold_shuffled", add_objects = c("threeutr.grl","run_rnaplfold","shuffle_sequence"),
                      nodes = opt$nodes, cpus_per_node = 1, slurm_options = list(time = "24:00:00"),
                      submit = TRUE)
  
  Sys.sleep(60)
  message("RNAplfold is running for shuffled control...")
  
  # check job is finished
  status <- FALSE
  while(status == FALSE) {
    squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left, = all jobs are finished
    Sys.sleep(60)
  }
  
  output <- get_slurm_out(sjob, outtype = "raw")
  saveRDS(output, "threeutrs.rnaplfold.shuffled.grl.rds")
  
  #cleanup_files(sjob)
  
  rnaplfold.shuffled.ls <- readRDS("threeutrs.rnaplfold.shuffled.grl.rds") #load the GRanges List
  rnaplfold.shuffled.grl <- GRangesList(rnaplfold.shuffled.ls)
  export.bed(unlist(rnaplfold.shuffled.grl), paste0(prefix, "_threeutrs.rnaplfold.shuffled.bed"), format = "BED")
  message("Probability scores for shuffled control (BED) exported.")
} else {
  message("A shuffled control was not generated.")
}

toc()
