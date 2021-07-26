#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(primavera))
suppressPackageStartupMessages(library(tictoc))
suppressPackageStartupMessages(library(parallel))

# ==========
# Data
# ==========

results.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/merged_clusters/stau1_atlas"

# Linker
linker.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_linker/linker.clusters.mfe.tsv.gz")

# No linker
nolinker.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nolinker/atlas/all.hybrids.tsv.gz")
nolinker.dt[, sample := tstrsplit(sample, "\\.")[[1]]]

# Non-hybrids, files produced in Figure_2.Rmd
nonhybrid.dt <- fread("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/stau1_nonhybrid.gc.txt")
nonhybrid.dt$sample <- "stau1_nonhybrid"


# ==========
# Annotation files
# ==========

genes.gr <- rtracklayer::import.gff2("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/GRCh38.gencode_v33.tx.gtf.gz")
regions.gr <- rtracklayer::import.gff2("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/regions.gtf.gz")


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
# Remove rRNA and tRNA hybrids
# ==========

all.list <- split(all.hybrids.dt, by = c("L_seqnames", "R_seqnames"))
all.clusters.list <- parallel::mclapply(all.list, cluster_hybrids, percent_overlap = 0.5, mc.cores = 4)
all.clusters.dt <- rbindlist(all.clusters.list, use.names = TRUE, fill = TRUE)

#all.clusters.dt[L_seqnames == R_seqnames][!L_seqnames %in% c("tRNA", "rDNA", "rRNA_5S")][grep("C", cluster), .N, by = .(cluster, L_seqnames, R_seqnames)]

all.collapsed  <- collapse_clusters(all.clusters.dt, mode = "median")
nrow(all.collapsed)

# ==========
# Annotate clusters
# ==========

all.collapsed.dt <- convert_coordinates(all.collapsed, genes.gr)
all.collapsed.dt <- annotate_hybrids(all.collapsed.dt, regions.gr)

all.clusters.dt <- convert_coordinates(all.clusters.dt, genes.gr)
all.clusters.dt <- annotate_hybrids(all.clusters.dt, regions.gr)

# ==========
# Export tables
# ==========

fwrite(all.clusters.dt, paste0(results.dir,"/merged_atlas.clusters.tsv.gz"), sep = "\t")
fwrite(all.collapsed.dt, paste0(results.dir,"/merged_atlas.clusters.collapsed.tsv.gz"), sep = "\t")





