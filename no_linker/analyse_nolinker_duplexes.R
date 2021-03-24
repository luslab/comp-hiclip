library(data.table)
library(rtracklayer)
library(primavera)
library(tidyverse)
library(stringr)
library(Biostrings)
library(tictoc)
library(rslurm)

# ==========
# Define functions
# ==========

GetMFE <- function(L_sequence, R_sequence) {
  
  input <- paste0(L_sequence, "\n", R_sequence)
  rnaduplex <- system("RNAduplex --noLP", input = input, intern = TRUE)
  
  rnaduplex <- gsub("\\s+", "_", rnaduplex)
  mfe <- sapply(strsplit(rnaduplex, "_"), "[", 5) # Need to fix for positives < 10 ? as they have an extra space
  mfe <- as.numeric(gsub("\\(|\\)", "", mfe))
  
  # Get dot.bracket
  db <- sapply(strsplit(rnaduplex, "_"), "[", 1)
  
  # Get from/to positions
  l <- sapply(strsplit(rnaduplex, "_"), "[", 2)
  r <- sapply(strsplit(rnaduplex, "_"), "[", 4)
  l_from <- as.numeric(sapply(strsplit(l, ","), "[", 1))
  l_to <- as.numeric(sapply(strsplit(l, ","), "[", 2))
  r_from <- as.numeric(sapply(strsplit(r, ","), "[", 1))
  r_to <- as.numeric(sapply(strsplit(r, ","), "[", 2))
  
  l_width <- nchar(L_sequence)
  r_width <- nchar(R_sequence)
  
  # Pad db with . according to from/to
  l_db <- sapply(strsplit(db, "&"), "[", 1)
  r_db <- sapply(strsplit(db, "&"), "[", 2)
  l_db <- unlist(strsplit(l_db, "")) # convert string to vector of characters
  r_db <- unlist(strsplit(r_db, ""))
  
  l_db <- c(rep(".", l_from - 1L), l_db, rep(".", l_width - l_to))
  
  r_db <- c(rep(".", r_from - 1L), r_db, rep(".", r_width - r_to))
  
  stopifnot(all(length(l_db == l_width), length(r_db) == r_width))
  db <- paste0(c(l_db, "&", r_db), collapse = "")
  
  #if(mfe > 0) mfe <- 0
  
  return(list(mfe = mfe, db = db))
  
}


ShuffleSequence <- function(sequence, number = 1, klet = 2, seed = 42) {
  
  system(paste0("ushuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)
}


slurm_GetMFE <- function(id, L_sequence, R_sequence) {
  
  binding <- GetMFE(L_sequence, R_sequence)
  mfe <- binding$mfe
  db <- binding$db
  names(mfe) <- id
  names(db) <- id
  
  #return(list(mfe = mfe, structure = structure))
  
  mfe.dt <- as.data.table(list(id = id, mfe = mfe, db = db))
  return(mfe.dt)
  
}



get_shuffled_mfe <- function(id, L_sequence, R_sequence) {
  
  #t <- data[data$id == id,]
  
  L <- ShuffleSequence(L_sequence, number = 100, klet = 2)
  R <- ShuffleSequence(R_sequence, number = 100, klet = 2)
  
  # make a df with the shuffled sequences for each id
  
  shuffled.dt <- as.data.table(list(L_sequence = L, R_sequence = R))
  
  # Get MFE mean and sd
  shuffled.df <- shuffled.dt %>%
    rowwise() %>%
    mutate(binding = list(GetMFE(L_sequence, R_sequence))) %>%
    mutate(mfe = binding$mfe) %>%
    ungroup() %>%
    dplyr::select(-binding) %>%
    summarise(mean_mfe = mean(mfe), mean_sd = sd(mfe))
  
  results.ls <- list(id = id, shuffled_mfe = shuffled.df$mean_mfe,
                     shuffled_mfe_sd = shuffled.df$mean_sd)
  shuff.dt <- as.data.table(results.ls)
  return(shuff.dt)
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

# Processing structure functions:

# find position of the first "(" and ")" in L and R respectively
# find position of the last "(" and ")" in L and R respectively
# trim external dots..to get all paired (find pos of first and last ())
# trim same positions from seq

get_paired_regions <- function(id, L_db, L_sequence, R_db, R_sequence) {
  
  l_brackets <- str_locate_all(L_db, pattern = "\\(")[[1]][,1]
  l_start <- l_brackets[1]
  l_end <- l_brackets[length(l_brackets)]
  l_trimmed_db <- str_sub(L_db, l_start, l_end)
  l_trimmed_seq <- str_sub(L_sequence, l_start, l_end)
  
  r_brackets <- str_locate_all(R_db, pattern = "\\)")[[1]][,1]
  r_start <- r_brackets[1]
  r_end <- r_brackets[length(r_brackets)]
  r_trimmed_db <- str_sub(R_db, r_start, r_end)
  r_trimmed_seq <- str_sub(R_sequence, r_start, r_end)
  
  stopifnot(length(l_brackets) == length(r_brackets)) # check L and R arms of the duplex have an equal number of bp
  
  paired_residues_count <- length(l_brackets)
  l_paired_residues <- paste(str_sub(L_sequence, start = l_brackets, end = l_brackets), collapse = '')
  r_paired_residues <- paste(str_sub(R_sequence, start = r_brackets, end = r_brackets), collapse = '')
  
  results.ls <- list(id = id, l_start =l_start, l_end = l_end, l_trimmed_db=l_trimmed_db, l_trimmed_seq=l_trimmed_seq,
                     r_start=r_start, r_end=r_end, r_trimmed_db=r_trimmed_db, r_trimmed_seq=r_trimmed_seq,
                     paired_residues_count = paired_residues_count,
                     l_paired_residues=l_paired_residues, r_paired_residues=r_paired_residues)
  
  return(results.ls)
  
}


# ==========
# Files and parameters
# ==========

# Annotation files

prefix <- "stau1_3utr_minmfe_cluster_hybrids."

get_elementstring.py_path <- getwd()

genes.gr <- import.gff2("/camp/lab/luscomben/home/users/iosubi/genomes/comp_hiclip_annotations/human_GencodeV33.gtf.gz")
regions.gr <- import.gff2("/camp/lab/luscomben/home/users/iosubi/genomes/comp_hiclip_annotations/regions.gtf.gz")

# hybrid table containing cluster information

hybrids.dt <- fread("all.clusters.tsv.gz")
hybrids_mfe.dt <- fread("all.mfe.clusters.tsv.gz")

# ==========
# Get genomic coordinates and annotate hybrids
# ==========

hybrids.dt <- hybrids.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
hybrids.dt <- hybrids.dt[grep("Mt", L_seqnames, invert = TRUE)]

# Reformat names of hybrids seqnames so they match the annotation gtf
hybrids.dt[, c("L_seqnames") := tstrsplit(L_seqnames, "::", fixed=TRUE)[1]]
hybrids.dt[, c("R_seqnames") := tstrsplit(L_seqnames, "::", fixed=TRUE)[1]]

# Convert to genomic coordinates
hybrids.dt <- convert_coordinates(hybrids.dt, genes.gr)
#fwrite(hybrids.dt, "clusters.gc.tsv.gz", sep = "\t")

# Annotate hybrids
hybrids.dt <- annotate_hybrids(hybrids.dt, regions.gr)
#fwrite(hybrids.dt, "annotated_clusters.gc.tsv.gz", sep = "\t")


hybrids_mfe.dt <- hybrids_mfe.dt[grep("^rRNA", L_seqnames, invert = TRUE)] # Remove rRNA
hybrids_mfe.dt <- hybrids_mfe.dt[grep("Mt", L_seqnames, invert = TRUE)]

hybrids_mfe.dt <- hybrids_mfe.dt %>%
  select(-umi) %>%
  mutate(Experiment = "STAU1 - No linker")

hybrids_mfe.dt <- hybrids_mfe.dt %>% dplyr::filter(str_detect(cluster, "C")) # keep hybrids within clusters

# ==========
# Get 3'UTR intragenic hybrids that belong to clusters and have the min MFE of the cluster
# ==========

hybrids_mfe.dt <- hybrids_mfe.dt %>%
  dplyr::filter(str_detect(cluster, "C")) %>%
  group_by(cluster) %>%  #had to group also by gene id because redundant cluster names in PARIS data
  dplyr::slice(which.min(mfe)) %>%
  ungroup()


# get list of 3'UTR cluster names from the anntated clusters table
threeutr_clusters.dt <- hybrids.dt %>%
  dplyr::filter(L_region == "UTR3" & R_region == "UTR3") 

threeutr_clusters.ls <- unique(threeutr_clusters.dt$name)

# filter hybrids table by 3'UTR cluster names
threeutr.dt <- hybrids_mfe.dt[hybrids_mfe.dt$cluster %in% threeutr_clusters.ls,]

# ==========
# Get MFE and shuffled control MFE
# ==========

seq.df <- threeutr.dt %>%
  dplyr::select(id, L_sequence, R_sequence)

sjob <- slurm_apply(slurm_GetMFE, seq.df, jobname = "mfe", nodes = 100, add_objects = c("GetMFE", "ShuffleSequence"),
                    cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)

Sys.sleep(60) # To give it enough time to submit before the first check

status <- FALSE
while(status == FALSE) {
  
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
  Sys.sleep(60)
  
}

mfe.ls <- get_slurm_out(sjob, outtype = 'raw')
mfe.dt <- rbindlist(mfe.ls)

threeutr.dt <- left_join(threeutr.dt, mfe.dt, by = c("id","mfe"))

saveRDS(mfe.ls, file = "mfe.rds")

# Remove temporary files
cleanup_files(sjob) 

# ==========
# Get shuffled control MFE
# ==========

sjob <- slurm_apply(get_shuffled_mfe, seq.df, jobname = "shuffled_mfe", nodes = 100, add_objects = c("GetMFE", "ShuffleSequence"),
                    cpus_per_node = 1, slurm_options = list(time = "24:00:00"), submit = TRUE)

Sys.sleep(60) # To give it enough time to submit before the first check

status <- FALSE
while(status == FALSE) {
  
  squeue.out <- system(paste("squeue -n", sjob$jobname), intern = TRUE) # Get contents of squeue for this job
  if(length(squeue.out) == 1) status <- TRUE # i.e. only the header left
  Sys.sleep(60)
  
}

shuffled.ls <- get_slurm_out(sjob, outtype = 'raw')
shuffled.dt <- rbindlist(shuffled.ls)
threeutr.dt <- left_join(threeutr.dt, shuffled.dt, by = "id")

saveRDS(shuffled.ls, file = "shuffled.rds")
cleanup_files(sjob)

# ==========
# Get paired positions, db, sequences, total sum of paired residues within a duplex
# ==========

# separate L_db and R_db
threeutr.dt <- threeutr.dt %>%
  separate(db, into = c("L_db", "R_db"), sep = "&", remove = FALSE)

threeutr.dt <- threeutr.dt %>%
  rowwise() %>%
  mutate(paired_regions = list(get_paired_regions(id, L_db, L_sequence, R_db, R_sequence))) %>%
  mutate(L_duplex_start= paired_regions$l_start, 
         L_duplex_end = paired_regions$l_end,
         L_duplex_db =  paired_regions$l_trimmed_db,
         L_duplex_sequence =  paired_regions$l_trimmed_seq,
         R_duplex_start= paired_regions$r_start, 
         R_duplex_end = paired_regions$r_end,
         R_duplex_db =  paired_regions$r_trimmed_db,
         R_duplex_seq =  paired_regions$r_trimmed_seq,
         total_paired = paired_regions$paired_residues_count,
         l_paired_residues = paired_regions$l_paired_residues,
         r_paired_residues = paired_regions$r_paired_residues) %>%
  ungroup() %>%
  dplyr::select(-paired_regions)

fwrite(threeutr.dt, paste0(prefix, "mfe_analyses.tsv.gz"), sep = "\t")

# ==========
# Annotate structures
# ==========

forgi_db <- threeutr.dt %>%
  unite("forgi_db", c(L_duplex_db, R_duplex_db), sep="...", remove = FALSE) %>%  # replace & with "..."
  dplyr::select(id, forgi_db)

forgi_db.ls <- forgi_db$forgi_db

forgi_input.df <- data.frame(db = forgi_db.ls, script_path = get_elementstring.py_path)

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
forgi.dt <- rbindlist(forgi)
cleanup_files(sjob)
toc()

forgi.dt$id <- rep(forgi_db$id, elementNROWS(forgi))
# remove the hairpins

forgi.dt <- forgi.dt %>%
  dplyr::filter(element_type != "h")

fwrite(forgi.dt, file = paste0(prefix, "forgi.tsv.gz"), sep = "\t")



