#!/usr/bin/env Rscript

library(rtracklayer)
library(BSgenome.Hsapiens.UCSC.hg38)
library(tictoc)
library(data.table)
library(pbapply)
library(rslurm)
library(optparse)
library(dplyr,warn.conflicts = FALSE)

run_rnafold <- function(sequence, args =  c("-p", "--noLP", "--noPS", "--MEA", paste0("< ",sequence))) {
  
  rnafold.out <- system2(command = "RNAfold", args = args, stdout = T)
  rnafold.dt <- data.table(sequence = rnafold.out[1],
                           mfe_structure = strsplit(rnafold.out[2], " ")[[1]][1],
                           partition_structure = strsplit(rnafold.out[3], " ")[[1]][1],
                           centroid_structure = strsplit(rnafold.out[4], " ")[[1]][1],
                           mea_structure = strsplit(rnafold.out[5], " ")[[1]][1],
                           mfe = as.numeric(gsub("\\(|\\)", "", strsplit(rnafold.out[2], " ")[[1]][2])),
                           partition_fe = as.numeric(gsub("\\[|\\]", "", strsplit(rnafold.out[3], " ")[[1]][2])),
                           centroid_fe = as.numeric(gsub("\\{|\\})", "", strsplit(rnafold.out[4], " ")[[1]][2])),
                           mea_fe = as.numeric(gsub("\\{|\\})", "", strsplit(rnafold.out[5], " ")[[1]][2])),
                           ensemble_diversity = as.numeric(gsub(".*diversity", "", rnafold.out[6])))
  
  return(rnafold.dt)
  
}


get_rnafold_out <- function(id, sequence) {
  
  temp.fasta <- paste0(id,".fasta")
  fwrite(list(sequence), file = paste0(id,".fasta"))
  run_rnafold(temp.fasta)
  
}


get_forgi <- function(db, script_path = ".") {
  
  # stopifnot(all(strsplit(db, "")[[1]] %in% c("(", ")", ".")))
  forgi.out <- system(command = paste0(script_path, "/get_elementstring.py"), input = db, intern = TRUE)
  
  forgi.elements <- grep("define", forgi.out, value = TRUE)
  # If only one entry then fread tries to find a file
  if(length(forgi.elements) == 1) {
    
    forgi.string <- strsplit(forgi.elements, " ")[[1]]
    forgi.dt <- data.table(V2 = forgi.string[2], V3 = forgi.string[3], V4 = forgi.string[4], V5 = forgi.string[5], V5 = forgi.string[6])
    
  } else {
    
    forgi.dt <- fread(text = forgi.elements, fill = TRUE)[, -1]
    
  }
  setnames(forgi.dt, c("element", "L_start", "L_end", "R_start", "R_end"))
  forgi.dt$L_start <- as.integer(forgi.dt$L_start)
  forgi.dt$L_end <- as.integer(forgi.dt$L_end)
  forgi.dt$R_start <- as.integer(forgi.dt$R_start)
  forgi.dt$R_end <- as.integer(forgi.dt$R_end)
  
  forgi.dt$L_width <- with(forgi.dt, L_end - L_start + 1)
  forgi.dt$R_width <- with(forgi.dt, R_end - R_start + 1)
  forgi.dt$element_type <- with(forgi.dt, substr(element, 1, 1))
  forgi.dt$element_number <- with(forgi.dt, substr(element, 2, nchar(element)))
  
  return(forgi.dt)
  
}

# =========
# Options and paths
# =========

option_list <- list(make_option(c("-b", "--bed"), action = "store", type = "character", default=NA, help = "Space separated bed file list"),
                    make_option(c("-p", "--prefix"), action = "store", type = "character", default=NA, help = "Prefix for output files"),
                    make_option(c("-n", "--nodes"), action = "store", type = "integer", default = 100, help = "Number of nodes to allocate [default: %default]"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

get_elementstring.py_path <- getwd()

# =========
# Load crosslinks bed file and get sequence
# =========

bed_list <-c("~/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1_low.10nt_10nt.peaks.bed.gz",
                "~/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1_high.10nt_10nt.peaks.bed.gz")

bed_list <- strsplit(opt$filelist, " ")
bed.grl <- GRangesList(lapply(bed_list, import.bed))
bed <- unlist(bed.grl)
bed <- keepStandardChromosomes(bed, pruning.mode = "coarse")

# mcols(bed) = NULL

# get start of peaks  +100
bed <- resize(bed, width = 1, fix = "start") # resize the peaks, start of peak = 1
# keep unique positions
bed <- unique(bed)
bed$id <- paste0("ID", 1:length(bed))
bed <- resize(bed, width = 100+1, fix = "start") # add + flank

fasta <- getSeq(Hsapiens, bed)

names(fasta) <- bed$id
fasta.df <- data.frame(id = names(fasta), sequence = as.character(fasta))


# ==========
# Run RNAfold
# ==========

tic()
sjob <- slurm_apply(get_rnafold_out, fasta.df, jobname = "RNAfold", add_objects = c("run_rnafold"),
                    nodes = opt$nodes, cpus_per_node = 1, slurm_options = list(time = "24:00:00"),
                    submit = TRUE)
Sys.sleep(60)
message("RNAfold is running..")

#check job is finished
status <- FALSE
while(status == FALSE) {
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left, = all jobs are finished
  Sys.sleep(60)
}

rnafold <- get_slurm_out(sjob, outtype = 'raw')
rnafold.dt <- rbindlist(rnafold)
cleanup_files(sjob)
toc()

fwrite(rnafold.dt, file = paste0(opt$prefix, ".rnafold.tsv.gz"), sep = "\t")

# ==========
# Run forgi
# ==========

forgi_input.df <- data.frame(db = as.character(rnafold.dt$mea_structure), script_path = get_elementstring.py_path)
forgi_input.df <- mutate(forgi_input.df, db = as.character(forgi_input.df$db))


tic()
sjob <- slurm_apply(get_forgi, forgi_input.df, jobname = paste0("forgi_", paste0(sample(c(LETTERS, letters, 0:9), 10), collapse = "")),
                    nodes = opt$nodes, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)

status <- FALSE
while(status == FALSE) {
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
  Sys.sleep(60)
}

forgi <- get_slurm_out(sjob, outtype = 'raw')
forgi.dt <- rbindlist(forgi)
cleanup_files(sjob)
toc()

forgi.dt$id <- rep(bed$name, elementNROWS(forgi))
fwrite(forgi.dt, file = paste0(opt$prefix, ".forgi.tsv.gz"), sep = "\t")
