#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(rslurm))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(stringr))


shuffle_sequence <- function(sequence, number = 1, klet = 2, seed = 42) {
  system(paste0("ushuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)
}

get_mfe <- function(name, L_sequence, R_sequence) {
  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)
  
  # Get MFE
  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5)
  if (mfe == "(") mfe <- sapply(strsplit(rnaduplex, "_"), "[", 6) # Positives < 10 have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))
  
  return(data.table(name = name, mfe = mfe))
}

get_shuffled_mfe <- function(name, L_sequence, R_sequence) {
  L <- shuffle_sequence(L_sequence, number = 100, klet = 2)
  R <- shuffle_sequence(R_sequence, number = 100, klet = 2)
  
  shuffled_mfe.dt <- rbindlist(lapply(seq_along(L), function(i) get_mfe(name, L[i], R[i])))
  shuffled_mfe.dt[, `:=`(
    mean_shuffled_mfe = mean(mfe, na.rm = TRUE),
    sd_shuffled_mfe = sd(mfe, na.rm = TRUE)
  ),
  by = name
  ]
  
  return(unique(shuffled_mfe.dt[, .(name, mean_shuffled_mfe, sd_shuffled_mfe)]))
}

# ==========
# Define functions
# ==========

get_forgi <- function(id, db, script_path = ".") {
  
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
  
  forgi.dt$id <- id
  
  return(forgi.dt)
  
}

# ==========
# Files and parameters
# ==========

option_list <- list(make_option(c("", "--input"), action = "store", type = "character", help = "Hybrids or clusters file"),
                    make_option(c("", "--fasta"), action = "store", type = "character", help = "Transcript fasta"),
                    make_option(c("", "--output"), action = "store", type = "character", help = "Output file"),
                    make_option(c("", "--nodes"), action = "store", type = "integer", default = 100, help = "Number of nodes to allocate [default: %default]"),
                    make_option(c("", "--shuffled_mfe"), action = "store_true", type = "logical", help = "Calculate shuffled binding energy (100 iterations)", default = FALSE),
                    make_option(c("", "--clusters_only"), action = "store_true", type = "logical", help = "Calculate only for duplexes in clusters", default = FALSE),
                    make_option(c("", "--intragenic"), action = "store_true", type = "logical", help = "Consider only intragenic clusters", default = FALSE),
                    make_option(c("", "--threeutr"), action = "store_true", type = "logical", help = "Consider only 3'UTR-3'UTR clusters (for forgi annotation only)", default = FALSE))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


get_elementstring.py_path <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/comp-hiclip/no_rnase"

fa.dss <- readDNAStringSet(opt$fasta)
clusters.dt <- fread(opt$input)

# ==========
# Get genomic sequence
# ==========

# Filter intragenic clusters if option

if(opt$intragenic) {
  clusters.dt <- clusters.dt %>%
    dplyr::filter(L_gene_id == R_gene_id)
} else {
  clusters.dt <- clusters.dt
}

# Extract genomic sequence

fa.dt <- data.table(gene_id = names(fa.dss), sequence = as.character(fa.dss))

# focus on reads that have been clustered, if analysing the collapsed clusters

if (opt$clusters_only) {
  
  clusters.dt <- data.table(dplyr::filter(clusters.dt, str_detect(cluster, "C")))
  
} else {
  
  clusters.dt <- clusters.dt
}


if (c("mfe", "structure") %in% colnames(clusters.dt)) {
  
  clusters.dt <- dplyr::select(clusters.dt, -mfe, -structure)
} else {
  clusters.dt <- clusters.dt
}

message(paste0("Analysing ", nrow(clusters.dt), " clusters"))

clusters.dt <- primavera::get_sequence(clusters.dt, fa.dt)
stopifnot(!any(is.na(c(clusters.dt$L_sequence, clusters.dt$L_sequence))))

# ==========
# Get MFE
# ==========

# Cluster jobs
tic()
sjob <- slurm_apply(analyse_structure, clusters.dt[, .(name, L_sequence, R_sequence)],
                    jobname = sapply(strsplit(basename(opt$input), "\\."), "[[", 1),
                    nodes = opt$nodes,
                    cpus_per_node = 1,
                    slurm_options = list(time = "24:00:00"),
                    submit = TRUE)

Sys.sleep(60)

message("RNAduplex is running..")

status <- FALSE
while(status == FALSE) {
  
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
  Sys.sleep(60)
  
}

structure.list <- get_slurm_out(sjob)
# saveRDS(structure.list, file = "mfe.rds")
toc()

cleanup_files(sjob) # Remove temporary files


# structure.list <- readRDS("mfe.rds")
structure.dt <- rbindlist(structure.list)
clusters.dt <- merge(clusters.dt, structure.dt, by = "name")
clusters.dt <- data.table(distinct(clusters.dt))


# ==========
# Shuffled 
# ==========

if(opt$shuffled_mfe) {
  sjob <- slurm_apply(get_shuffled_mfe, clusters.dt[, .(name, L_sequence, R_sequence)], 
                      jobname = sapply(strsplit(basename(opt$input), "\\."), "[[", 1), 
                      nodes = opt$nodes, 
                      add_objects = c("shuffle_sequence", "get_mfe"),
                      cpus_per_node = 1, 
                      slurm_options = list(time = "24:00:00"), 
                      submit = TRUE)
  
  Sys.sleep(60) # To give it enough time to submit before the first check
  status <- FALSE
  while(status == FALSE) {
    
    squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
    if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
    Sys.sleep(60)
    
  }
  
  mfe.list <- get_slurm_out(sjob)
  #saveRDS(mfe.list, file = "shuffled_mfe.rds")
  cleanup_files(sjob) # Remove temporary files
  
  mfe.dt <- rbindlist(mfe.list)
  clusters.dt <- merge(clusters.dt, mfe.dt, by = "name")
  
}


clusters.dt <- data.table(distinct(clusters.dt))
fwrite(clusters.dt, opt$output, sep = "\t")


# Filter intragenic 3"UTR hybrids

if(opt$threeutr) {
  clusters.dt <- clusters.dt %>%
    dplyr::filter(L_gene_id == R_gene_id & L_region == "UTR3" & R_region == "UTR3")
} else {
  clusters.dt <- clusters.dt
}

# ==========
# Annotate structures
# ==========

# Separate db structure into L_db and R_db
clusters.dt <- clusters.dt %>%
  separate(structure, into = c("L_db", "R_db"), sep = "&", remove = FALSE)

forgi_db.df <- clusters.dt %>%
  unite("forgi_db", c(L_db, R_db), sep="...", remove = FALSE) %>%  # replace & with "..."
  dplyr::select(name, forgi_db)

forgi_input.df <- data.frame(id = forgi_db.df$name, db = forgi_db.df$forgi_db, script_path = get_elementstring.py_path)
forgi_input.df <- mutate(forgi_input.df, db = as.character(forgi_input.df$db))

tic()
sjob <- slurm_apply(get_forgi, forgi_input.df, jobname = paste0("forgi_", paste0(sample(c(LETTERS, letters, 0:9), 10), collapse = "")),
                    nodes = opt$nodes, cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)

message("forgi is running..")

status <- FALSE
while(status == FALSE) {
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
  Sys.sleep(60)
}

forgi <- get_slurm_out(sjob, outtype = 'raw')
# saveRDS(forgi, file = "forgi.rds")
forgi.dt <- rbindlist(forgi)
cleanup_files(sjob)
toc()

# remove the hairpins
forgi.dt <- forgi.dt %>%
  dplyr::filter(element_type != "h")

forgi.output.filename <- str_replace(opt$output, "mfe", "forgi")
fwrite(forgi.dt, file = forgi.output.filename, sep = "\t")
