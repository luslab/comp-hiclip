#!/usr/bin/env Rscript

# Script to cluster linker hybrid atlas and get structures
# A. M. Chakrabarti
# 16th December 2022

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rslurm))

setDTthreads(8)
set.seed(42)

# ==========
# Data
# ==========

linker.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_linker"
nolinker.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nolinker/hybrids"
nornase.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase"

ref.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref"
results.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/atlas"


# Linker
linker.dt <- rbindlist(list(fread(paste0(linker.dir, "/stau1_high.hybrids.gc.annotated.tsv.gz"))[, sample := "stau1_linker"],
                             fread(paste0(linker.dir,"/stau1_low.hybrids.gc.annotated.tsv.gz"))[, sample := "stau1_linker"]),
                       use.names = TRUE)


# No linker
nolinker.dt <- rbindlist(list(fread(paste0(nolinker.dir, "/stau1_high.hybrids.gc.annotated.tsv.gz"))[, sample := "stau1_nolinker"],
                             fread(paste0(nolinker.dir,"/stau1_low.hybrids.gc.annotated.tsv.gz"))[, sample := "stau1_nolinker"]),
                       use.names = TRUE)

# Non-hybrids, file produced in Figure_3.Rmd
nonhybrid.dt <- fread(paste0(nornase.dir, "/short_range_duplexes_min8bp.tsv.gz"))
nonhybrid.dt$sample <- "stau1_derived"

# ==========
# Annotation files
# ==========

genes.gr <- rtracklayer::import.gff2(paste0(ref.dir, "/GRCh38.gencode_v33.tx.gtf.gz"))
regions.gr <- rtracklayer::import.gff2("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/icount_mini_utr3/regions.gtf.gz")

# ==========
# Check data
# ==========

# Make sure all tables have the same cols
colnames.ls <- intersect(colnames(nonhybrid.dt), colnames(linker.dt))
nolinker.dt[, c(setdiff(colnames(nolinker.dt), colnames.ls)):=NULL]
linker.dt[, c(setdiff(colnames(linker.dt), colnames.ls)):=NULL]
nonhybrid.dt[, c(setdiff(colnames(nonhybrid.dt), colnames.ls)):=NULL]

all.hybrids.dt <- rbind(linker.dt, nolinker.dt)
nrow(all.hybrids.dt) # linker and no linker
all.hybrids.dt <- rbind(all.hybrids.dt, nonhybrid.dt)
nrow(all.hybrids.dt)

atlas.hybrids.dt <- all.hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]


# ==========
# Remove rRNA intramolecular hybrids
# ==========

atlas.hybrids.dt <- atlas.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
nrow(atlas.hybrids.dt )

# ==========
# Cluster hybrids
# ==========

atlas.hybrids.dt <- atlas.hybrids.dt[total_count > 1]

# Keep ones not clustered to add back in later
unclustered.hybrids.dt <- all.hybrids.dt[!name %in% atlas.hybrids.dt$name]
stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(all.hybrids.dt))

# Split into list to parallelise
atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))

# Cluster
# atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = 0.5, fraction = FALSE, mc.cores = 8)
atlas.clusters.list <- lapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = 0.5)
clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)

# Merge back and sanity check
clusters.dt <- rbindlist(list(clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)
setorder(clusters.dt, name)

stopifnot(nrow(clusters.dt) == nrow(all.hybrids.dt))
stopifnot(all(clusters.dt$name %in% all.hybrids.dt$name))

collapsed.dt <- collapse_clusters(clusters.dt, mode = "median")
collapsed.dt <- convert_coordinates(collapsed.dt, genes.gr)

# ==========
# Add the rest of the non-hybrid reads
# ==========

# Non-hybrids, file produced in Figure_3.Rmd
nonhybrid.dt <- fread(paste0(nornase.dir, "/short_range_duplexes_min8bp.tsv.gz"))
nonhybrid.dt$sample <- "stau1_derived"

nonhybrid.dt$cluster <- nonhybrid.dt$name

# sanity check
test  <-  clusters.dt %>%
  dplyr::filter(sample == "stau1_derived")
stopifnot(nrow(test) == nrow(nonhybrid.dt))


short_range_clustered.df <- clusters.dt %>%
  dplyr::filter(str_detect(cluster, "C") & sample == "stau1_derived")

# Get non-clustered derived 
nonhybrid.dt <- nonhybrid.dt %>%
  dplyr::filter(!(name %in% short_range_clustered.df$name)) %>%
  # mutate(total_count = as.integer(NA)) %>% # add in these columns so can merge later with the clustered
  dplyr::select(colnames(collapsed.dt))
nonhybrid.dt <- as.data.table(nonhybrid.dt)

nrow(collapsed.dt)
collapsed.dt <- rbind(collapsed.dt, nonhybrid.dt)
nrow(collapsed.dt)

# ==========
# Annotate collapsed clusters
# ==========

collapsed.dt <- annotate_hybrids(collapsed.dt, regions.gr)

# ==========
# Analyse structure of collapsed clusters
# ==========

genome.fa <- Biostrings::readDNAStringSet(paste0(ref.dir,"/GRCh38.gencode_v33.fa"))
genome.dt <- data.table(gene_id = names(genome.fa),
                        sequence = as.character(genome.fa))


collapsed.dt  <- get_sequence(hybrids.dt = collapsed.dt , genome.dt = genome.dt)
stopifnot(!any(is.na(c(collapsed.dt$L_sequence, collapsed.dt$R_sequence))))

sjob <- slurm_apply(analyse_structure, collapsed.dt[, .(name, L_sequence, R_sequence)], 
                  jobname = "structure", 
                  nodes = 100, 
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
structure.list <- get_slurm_out(sjob)
structure.dt <- rbindlist(structure.list, use.names = TRUE)
cleanup_files(sjob)

# ==========
# Shuffled
# ==========
sjob <- slurm_apply(get_shuffled_mfe, collapsed.dt[, .(name, L_sequence, R_sequence)], 
                  jobname = "shuffled_mfe", 
                  nodes = 100, 
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
shuffled.list <- get_slurm_out(sjob)
shuffled.dt <- rbindlist(shuffled.list, use.names = TRUE)
cleanup_files(sjob)

# ==========
# Export tables
# ==========

structures.collapsed.dt <- merge(collapsed.dt, structure.dt, by = "name", all.x = TRUE)
shuffled.collapsed.dt <- merge(structures.collapsed.dt, shuffled.dt, by = "name", all.x = TRUE)

fwrite(shuffled.collapsed.dt, paste0(results.dir,"/merged", ".", "atlas", ".clusters.tsv.gz"), sep = "\t")

fwrite(clusters.dt, paste0(results.dir,"/merged.atlas.clustered.tsv.gz"), sep = "\t")
