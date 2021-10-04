#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))

# ==========
# Data
# ==========

results.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/merged_clusters/stau1_atlas"

# Linker
linker.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_linker/linker.clusters.mfe.tsv.gz")
linker.dt$sample <- "stau1_linker"

# No linker
nolinker.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nolinker/atlas/all.hybrids.tsv.gz")
#nolinker.dt[, sample := tstrsplit(sample, "\\.")[[1]]]
nolinker.dt$sample <- "stau1_nolinker"

# Non-hybrids, files produced in Figure_3.Rmd
nonhybrid.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/short_range_duplexes_min8bp.tsv.gz")
nonhybrid.dt$sample <- "stau1_nonhybrid"


# ==========
# Annotation files
# ==========

ref.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/"

genes.gr <- rtracklayer::import.gff2(paste0(ref.dir, "/GRCh38.gencode_v33.tx.gtf.gz"))
regions.gr <- rtracklayer::import.gff2(paste0(ref.dir,"/regions.gtf.gz"))


# ==========
# Check data
# ==========

colnames.ls <- intersect(colnames(nonhybrid.dt), colnames(linker.dt))
nolinker.dt[, c(setdiff(colnames(nolinker.dt), colnames.ls)):=NULL]
linker.dt[, c(setdiff(colnames(linker.dt), colnames.ls)):=NULL]
nonhybrid.dt[, c(setdiff(colnames(nonhybrid.dt), colnames.ls)):=NULL]

nrow(nolinker.dt)
nrow(linker.dt)
nrow(nonhybrid.dt)

all.hybrids.dt <- rbind(linker.dt, nolinker.dt)
nrow(all.hybrids.dt) # linker and no linker

all.hybrids.dt <- rbind(all.hybrids.dt, nonhybrid.dt)
nrow(all.hybrids.dt)
#all.hybrids.dt[, c("cluster", "cluster_hybrid_count"):=NULL]

# ==========
# Remove rRNA and tRNA hybrids
# ==========

#all.hybrids.dt[, total_count := NULL]

all.hybrids.dt <- all.hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
all.hybrids.dt <- all.hybrids.dt[L_seqnames == R_seqnames]
all.hybrids.dt <- all.hybrids.dt[!(L_seqnames == "rRNA_45S" & R_seqnames == "rRNA_45S")]
all.hybrids.dt <- all.hybrids.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
all.hybrids.dt <- all.hybrids.dt[!(L_seqnames == "rRNA_5S" & R_seqnames == "rRNA_5S")]
all.hybrids.dt <- all.hybrids.dt[!(L_seqnames == "rRNA5S" & R_seqnames == "rRNA5S")]
all.hybrids.dt <- all.hybrids.dt[!grepl("tRNA", L_seqnames)]
all.hybrids.dt <- all.hybrids.dt[!grepl("tRNA", R_seqnames)]

nrow(all.hybrids.dt)

# ==========
# Cluster hybrids
# ==========

all.list <- split(all.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
all.clusters.list <- parallel::mclapply(all.list, cluster_hybrids, percent_overlap = 0.5, mc.cores = 4)
all.clusters.dt <- rbindlist(all.clusters.list, use.names = TRUE, fill = TRUE)
all.clusters.dt <- convert_coordinates(all.clusters.dt, genes.gr)
all.clusters.dt <- annotate_hybrids(all.clusters.dt, regions.gr)
#all.clusters.dt[L_seqnames == R_seqnames][!L_seqnames %in% c("tRNA", "rDNA", "rRNA_5S")][grep("C", cluster), .N, by = .(cluster, L_seqnames, R_seqnames)]

all.collapsed  <- collapse_clusters(all.clusters.dt, mode = "median")
all.collapsed.dt <- convert_coordinates(all.collapsed, genes.gr)

# ==========
# Add the rest of the non-hybrid reads
# ==========

nonhybrid.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/short_range_duplexes_min8bp.tsv.gz")
nonhybrid.dt$sample <- "stau1_nonhybrid"
nonhybrid.dt$cluster <- nonhybrid.dt$name

test  <-  all.clusters.dt %>%
  dplyr::filter(sample == "stau1_nonhybrid")
stopifnot(nrow(test) == nrow(nonhybrid.dt))
short_range_clustered.df <- all.clusters.dt %>%
  dplyr::filter(str_detect(cluster, "C") & sample == "stau1_nonhybrid")

nonhybrid.dt <- nonhybrid.dt %>%
  dplyr::filter(!(name %in% short_range_clustered.df$name)) %>%
  dplyr::select(colnames(all.collapsed.dt))

nonhybrid.dt <- as.data.table(nonhybrid.dt)

nrow(all.collapsed.dt)
all.collapsed.dt <- rbind(all.collapsed.dt, nonhybrid.dt)
nrow(all.collapsed.dt)

# ==========
# Annotate hybrids and clusters
# ==========

all.collapsed.dt <- annotate_hybrids(all.collapsed.dt, regions.gr)

# ==========
# Export tables
# ==========

fwrite(all.clusters.dt, paste0(results.dir,"/merged_atlas.clusters.tsv.gz"), sep = "\t")
fwrite(all.collapsed.dt, paste0(results.dir,"/merged_atlas.clusters.collapsed_plus_nonhybrids.tsv.gz"), sep = "\t")
