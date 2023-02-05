#!/usr/bin/env Rscript

# Script to make static visualisation for Fig. 4 and Suppl. Fig. 4
# A. M. Chakrabarti
# Last updated: 24th June 2022

library(rtracklayer)
library(GenomicFeatures)
library(ggplot2)
library(scales)
library(cowplot)
library(patchwork)
library(data.table)
library(ggthemes)
library(toscatools)

results.path <- "/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions"
plot.path <- "~/Dropbox (The Francis Crick)/comp_hiclip/revisions/plots/figure_4"
if(!dir.exists(plot.path)) dir.create(plot.path)

# Based off clipplotr code

# ARF1
region.gr <- GRanges(seqnames = "chr1",
                     ranges = IRanges(228097850, 228099296),
                     strand = "+")

# ==========
# Arcs
# ==========

# stau1.dt <- fread("~/Dropbox (The Francis Crick)/comp_hiclip/stau1_atlas/merged_atlas.clusters.collapsed_plus_nonhybrids.tsv.gz")
stau1.dt <- fread(file.path(results.path, "atlas/merged.atlas.clusters.tsv.gz"))

sel.gene.dt <- stau1.dt[L_gene_name == "ARF1"]
sel.gene.dt[, `:=` (L_start = L_genomic_start,
                    L_end = L_genomic_end,
                    R_start = R_genomic_start,
                    R_end = R_genomic_end)]

# Adjust so L and R arms are same width
sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]
sel.gene.dt[, adjust := R_width - L_width]

sel.gene.dt[adjust > 0, `:=` (R_start = R_start + floor(adjust/2),
                              R_end = R_end - ceiling(adjust/2))]

sel.gene.dt[adjust < 0, `:=` (L_start = L_start + floor(abs(adjust)/2),
                              L_end = L_end - ceiling(abs(adjust)/2))]

sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]

stopifnot(all(sel.gene.dt$L_width == (sel.gene.dt$R_width)))

# Expand for plotting
sel.expanded.gene.dt <- sel.gene.dt[rep(seq(1, nrow(sel.gene.dt)), sel.gene.dt$L_width)]
sel.expanded.gene.dt[, `:=` (L_adjust = 1:.N,
                             R_adjust = .N:1), by = cluster]
sel.expanded.gene.dt[, `:=` (L = L_start + L_adjust - 1,
                             R = R_start + R_adjust - 1)]

setorder(sel.expanded.gene.dt, count)

p1 <- ggplot(sel.expanded.gene.dt[grepl("^ID", cluster)]) +
  geom_curve(aes(x = L, xend = R, y = 0, yend = 0, colour = count), lineend = "round", curvature = 1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(-10, 0), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 36)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

p2 <- ggplot(sel.expanded.gene.dt[grepl("^C", cluster)]) +
  geom_curve(aes(x = L, xend = R, y = 0, yend = 0, colour = count), lineend = "round", curvature = -1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(0, 10), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 36)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top")

# ==========
# Peaks
# ==========

# peaks.gr <- import.bed("~/Dropbox (The Francis Crick)/comp_hiclip/short_range/sept_2021/stau1.10nt_10nt.peaks.bed.gz")
# peaks.gr <- subsetByOverlaps(peaks.gr, region.gr, ignore.strand = FALSE)
# auxiliary.dt <- as.data.table(peaks.gr)
# auxiliary.dt$centre <- with(auxiliary.dt, start + width/2)
# 
# p.auxiliary <- ggplot(auxiliary.dt, aes(x = centre, width = width, y = factor(1))) +
#   geom_tile() +
#   scale_y_discrete(labels = "") +
#   xlim(start(region.gr), end(region.gr)) +
#   labs(y = "",
#        x = "") +
#   theme_minimal_grid() + theme(legend.position = "none")

# ==========
# xlink tracks
# ==========

SubsetBedgraph <- function(gr, selected.region.gr) {
  
  xlinks.gr <- unlist(tile(selected.region.gr, width = 1))
  
  ol <- findOverlaps(xlinks.gr, gr)
  xlinks.gr$score <- 0
  xlinks.gr[queryHits(ol)]$score <- gr[subjectHits(ol)]$score
  # xlinks.gr$score[is.na(xlinks.gr$score)] <- 0
  
  return(xlinks.gr)
  
}

xlinks.gr <- import.bed(file.path(results.path, "results_nornase/xlinks/stau1.xl.bed.gz"))
seqlevelsStyle(xlinks.gr) <- "UCSC"
xlinks.gr <- SubsetBedgraph(gr = xlinks.gr, selected.region.gr = region.gr)
xlinks.dt <- as.data.table(xlinks.gr)[, sample := "STAU1"]
xlinks.dt[, norm := score]
xlinks.dt[, smoothed := zoo::rollmean(norm, 10, fill = 0)]
y.label <- "Crosslink signal"

p.iclip <- ggplot(xlinks.dt) +
  geom_line(aes(x = start, y = smoothed, group = sample, color = sample)) +
  labs(x = "",
       y = y.label,
       colour = "") +
  scale_colour_tableau(palette = "Tableau 10") +
  theme_minimal_grid() + theme(legend.position = "top") +
  xlim(start(region.gr),end(region.gr))

# ==========
# Annotation
# ==========

# TxDb <- loadDb("~/Dropbox (The Francis Crick)/comp_hiclip/ref/gencode.v33.txdb.sqlite")
TxDb <- loadDb("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/gencode.v33.txdb.sqlite")
seqlevelsStyle(TxDb) <- "UCSC"

# gtf <- import.gff2("~/Dropbox (The Francis Crick)/comp_hiclip/ref/gencode.v33.annotation.gtf.gz")
gtf <- import.gff2("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/gencode.v33.annotation.gtf.gz")
seqlevelsStyle(gtf) <- "UCSC"
genes.gr <- gtf[gtf$type == "gene"]

rosetta.dt <- as.data.table(mcols(genes.gr))[, .(gene_id, gene_name)]
setkey(rosetta.dt, gene_id)

# Get transcripts that overlap region and order for plotting
sel.tx_genes <- transcriptsByOverlaps(TxDb, region.gr, columns = c("gene_id", "tx_name"))
sel.tx_genes <- sel.tx_genes[order(width(sel.tx_genes), decreasing = TRUE)]
sel.tx_genes <- sel.tx_genes[sel.tx_genes$tx_name == "ENST00000272102.10"]
tx.order.dt <- data.table(transcript_id = sel.tx_genes$tx_name,
                          gene_id = unlist(sel.tx_genes$gene_id),
                          centre = start(sel.tx_genes) + width(sel.tx_genes)/2)[, group := 1:.N]
setkey(tx.order.dt, gene_id)
tx.order.dt <- rosetta.dt[tx.order.dt]
tx.order.dt[, gene := paste0(gene_name, " | ", gene_id)]

# Region (for arrows)
region <- region.gr
region$gene_id <- NULL
region.tiled <- rep(region, length(sel.tx_genes))
region.tiled$tx_name <- sel.tx_genes$tx_name
region.tiled <- tile(region.tiled, width = round(width(region)/15, -1))

sel.region.tiled <- GRangesList(lapply(seq_along(sel.tx_genes), function(i) {
  gr <- subsetByOverlaps(region.tiled[[i]], sel.tx_genes[i], type = "any")
  start(gr[1]) <- start(sel.tx_genes[i])
  end(gr[length(gr)]) <- end(sel.tx_genes[i])
  return(gr)
}))

names(sel.region.tiled) <- sel.tx_genes$tx_name
sel.region.tiled.dt <- as.data.table(sel.region.tiled)
sel.region.tiled.dt[, group := NULL]
setnames(sel.region.tiled.dt, "group_name", "transcript_id")
if(as.character(strand(region.gr)) == "-") setnames(sel.region.tiled.dt, c("start", "end"), c("end", "start"))

sel.region.tiled.dt <- merge(sel.region.tiled.dt, tx.order.dt, by = "transcript_id")

# CDS
cds_tx <- cdsBy(TxDb, by = "tx", use.names = TRUE)
sel.cds_tx <- cds_tx[names(cds_tx) %in% sel.tx_genes$tx_name]
sel.cds_tx.dt <- as.data.table(sel.cds_tx)
sel.cds_tx.dt[, group := NULL]
setnames(sel.cds_tx.dt, "group_name", "transcript_id")
sel.cds_tx.dt <- merge(sel.cds_tx.dt, tx.order.dt, by = "transcript_id")

# UTR5
utr5_tx <- fiveUTRsByTranscript(TxDb, use.names = TRUE)
sel.utr5_tx <- utr5_tx[names(utr5_tx) %in% sel.tx_genes$tx_name]
sel.utr5_tx.dt <- as.data.table(sel.utr5_tx)
sel.utr5_tx.dt[, group := NULL]
setnames(sel.utr5_tx.dt, "group_name", "transcript_id") 
sel.utr5_tx.dt <- merge(sel.utr5_tx.dt, tx.order.dt, by = "transcript_id")

# UTR3
utr3_tx <- threeUTRsByTranscript(TxDb, use.names = TRUE)
sel.utr3_tx <- utr3_tx[names(utr3_tx) %in% sel.tx_genes$tx_name]
sel.utr3_tx.dt <- as.data.table(sel.utr3_tx)
sel.utr3_tx.dt[, group := NULL]
setnames(sel.utr3_tx.dt, "group_name", "transcript_id")
sel.utr3_tx.dt <- merge(sel.utr3_tx.dt, tx.order.dt, by = "transcript_id")

# Exons
exons_tx <- exonsBy(TxDb, by = "tx", use.names = TRUE)
sel.exons_tx <- exons_tx[names(exons_tx) %in% sel.tx_genes$tx_name[!sel.tx_genes$tx_name %in% names(sel.cds_tx)]] # Don't want genes with CDS, just e.g. ncRNA
sel.exons_tx.dt <- as.data.table(sel.exons_tx)
sel.exons_tx.dt[, group := NULL]
setnames(sel.exons_tx.dt, "group_name", "transcript_id")
sel.exons_tx.dt <- merge(sel.exons_tx.dt, tx.order.dt, by = "transcript_id")

# Plot
p.annot <- ggplot() +
  geom_segment(data = sel.region.tiled.dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
  geom_rect(data = sel.exons_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
  geom_rect(data = sel.cds_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
  geom_rect(data = sel.utr5_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
  geom_rect(data = sel.utr3_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
  scale_fill_tableau() +
  theme_minimal_vgrid() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
  labs(x = "",
       y = "",
       fill = "") +
  coord_cartesian(xlim = c(start(region), end(region)))

# ==========
# Assemble
# ==========

ggsave(p2 / p1 / p.iclip / p.annot + plot_layout(heights = c(10, 3, 2, 1)), 
       filename = file.path(plot.path, "arf1_arc.pdf"), 
       width = 297, height = 210, units = "mm")



# ==============================

# ==========
# Contact map
# ==========

# hybrids.dt <- sfread(file.path(results.path, "atlas/merged.atlas.clusters.tsv.gz"))
# sel.gene.dt <- stau1.dt[L_gene_name == "ARF1"]
# sel.gene.dt <- sel.gene.dt[rep(seq(1, nrow(sel.gene.dt)), sel.gene.dt$count)]

fai.dt <- fread("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/GRCh38.gencode_v33.fa.fai", select = 1:2, col.names = c("gene", "length"))
genome.size <- as.integer(fai.dt[grep("^ARF1", gene)]$length)
# mat <- get_contact_map(hybrid.dt = sel.gene.dt, genome.size = genome.size)
# 
# binned.mat <- bin_matrix(mat, bin.size = 25)
# rm(mat)
# 
# binned.dt <- data.table(reshape2::melt(binned.mat))
# binned.dt <- binned.dt[value != 0]
# 
# p1 <- ggplot(binned.dt, aes(x = Var1, y = Var2, fill = log10(value))) +
#   geom_tile() +
#   geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
#   scale_fill_viridis_c(na.value = "white") +
#   theme_minimal_grid() + theme(legend.position = "top") +
#   # coord_equal(xlim = c(0, utr3.width), ylim = c(0, utr3.width)) +
#   labs(title = "ARF1", x = "STAU2", y = "STAU1", fill = expression(log[10]~counts)) +
#   coord_equal()
# 
# p1

# ==========
# Split approach
# ==========

# linker.dt <- rbindlist(list(
#   fread(file.path(results.path, "/results_linker/lph.hybrids.dedup.tsv.gz")[, sample := "stau1_high"],
#   fread("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_linker/lpl.hybrids.dedup.tsv.gz")[, sample := "stau1_low"]),
#   use.names = TRUE)
# 
# # Add calculations
# linker.dt[, `:=` (L_end = L_start + L_width - 1,
#                   R_end = R_start + R_width - 1)]
# linker.dt[, type := ifelse(L_seqnames == R_seqnames, "intragenic", "intergenic")]
# linker.dt[, orientation := "linker"]
# 
# genes.gr <- import.gff2("~/Dropbox (The Francis Crick)/comp_hiclip/ref/GRCh38.gencode_v33.tx.gtf.gz")
regions.gr <- import.gff2("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/icount_mini_utr3/regions.gtf.gz")
# 
# linker.dt <- convert_coordinates(hybrids.dt = linker.dt, genes.gr = genes.gr)
# linker.dt <- annotate_hybrids(hybrids.dt = linker.dt, regions.gr = regions.gr)
# 
# # No linker
# nolinker.dt <- fread("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nolinker/atlas/all.hybrids.tsv.gz")
# nolinker.dt[, sample := tstrsplit(sample, "\\.")[[1]]]
# 
# hybrids.dt <- rbindlist(list(nolinker.dt, linker.dt),
#                         use.names = TRUE,
#                         fill = TRUE)
# hybrids.dt[, exp := ifelse(orientation == "linker", "Linker", "Direct")]
# hybrids.dt[sample == "stau1_high", sample := "High RNase"]
# hybrids.dt[sample == "stau1_low", sample := "Low RNase"]
# 
# hybrids.dt[, type := ifelse(type == "intragenic", "intra-transcript", "inter-transcript")]


linker.dt <- fread(file.path(results.path, "results_linker/all.atlas.clustered.tsv.gz"))
nolinker.dt <- fread(file.path(results.path, "results_nolinker/atlas/all.atlas.clustered.tsv.gz"))
nolinker.dt[, `:=` (cluster_hybrid_count = NULL,
                    cluster_hybrid_count.x = NULL,
                    sample = tstrsplit(sample, "\\.")[[1]])]
setnames(nolinker.dt, "cluster_hybrid_count.y", "cluster_hybrid_count")

hybrids.dt <- rbindlist(list(linker.dt, nolinker.dt),
                        use.names = TRUE, 
                        fill = TRUE)

# Adjust parameter names for plots
hybrids.dt[, exp := ifelse(orientation == "linker", "Linker", "Direct")]
hybrids.dt[, type := ifelse(type == "intragenic", "intra-transcript", "inter-transcript")]
hybrids.dt[sample == "stau1_high", sample := "High RNase"]
hybrids.dt[sample == "stau1_low", sample := "Low RNase"]

# ARF1 Map

arf1.dt <- hybrids.dt[type == "intra-transcript"][L_gene_name == "ARF1"]
arf1.utr3.size <- width(regions.gr[regions.gr$gene_name == "ARF1" & regions.gr$type == "UTR3"])

# fai.dt <- fread("~/Dropbox (The Francis Crick)/comp_hiclip/ref/GRCh38.gencode_v33.fa.fai", select = 1:2, col.names = c("gene", "length"))
# genome.size <- as.integer(fai.dt[grep("^ARF1", gene)]$length)

# hiCLIP
mat <- get_contact_map(hybrid.dt = arf1.dt, genome.size = genome.size)
binned.dt <- data.table(reshape2::melt(mat))
rm(mat)
hiclip.dt <- binned.dt[value != 0]

# Derived Short range
hybrids.dt <- fread(file.path(results.path, "atlas/merged.atlas.clusters.tsv.gz"))
sel.gene.dt <- stau1.dt[L_gene_name == "ARF1"][grepl("^ID", cluster)]
sel.gene.dt <- sel.gene.dt[rep(seq(1, nrow(sel.gene.dt)), sel.gene.dt$count)]

mat <- get_contact_map(hybrid.dt = sel.gene.dt, genome.size = genome.size)
binned.dt <- data.table(reshape2::melt(mat[15000:genome.size, 15000:genome.size]))
rm(mat)
derived.dt <- binned.dt[value != 0]
setnames(derived.dt, c("Var2", "Var1", "value"))
derived.dt[, `:=` (Var2 = Var2 + 15000 - 1,
                   Var1 = Var1 + 15000 - 1)]         

map.dt <- rbindlist(list(hiclip.dt, derived.dt), use.names = TRUE)
map.dt[, `:=` (Var2 = Var2 + min(start(regions.gr[regions.gr$gene_name == "ARF1"])) - 1,
               Var1 = Var1 + min(start(regions.gr[regions.gr$gene_name == "ARF1"])) - 1)]   

# Plot

p1 <- ggplot(map.dt, aes(x = Var1, y = Var2, fill = log10(value))) +
  geom_tile() +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  scale_fill_viridis_c(option = "G", direction = -1, breaks = c(0, 1, 2)) +
  scale_x_continuous(label = comma, n.breaks = 4) +
  scale_y_continuous(label = comma, n.breaks = 4) +
  coord_equal(xlim = c(max(end(regions.gr[regions.gr$gene_name == "ARF1"])) - arf1.utr3.size, 
                       max(end(regions.gr[regions.gr$gene_name == "ARF1"]))), 
              ylim = c(max(end(regions.gr[regions.gr$gene_name == "ARF1"])) - arf1.utr3.size, 
                       max(end(regions.gr[regions.gr$gene_name == "ARF1"])))) +
  labs(title = "ARF1", 
       x = "Derived duplexes", 
       y = "hiCLIP duplexes", 
       fill = expression(log[10]~counts)) +
  theme_minimal_grid() + theme(legend.position = "top")

ggsave(file.path(plot.path, "arf1_map.pdf"),
       p1,
       width = 150, 
       height = 150,
       units = "mm")

# ==========
# SRSF1
# ==========

region.gr <- GRanges(seqnames = "chr17",
                     ranges = IRanges(58000600, 58006400),
                     strand = "-")

# ==========
# Arcs
# ==========

stau1.dt <- fread(file.path(results.path, "atlas/merged.atlas.clusters.tsv.gz"))

sel.gene.dt <- stau1.dt[L_gene_name == "SRSF1"]
sel.gene.dt[, `:=` (L_start = L_genomic_start,
                    L_end = L_genomic_end,
                    R_start = R_genomic_start,
                    R_end = R_genomic_end)]

# Adjust so L and R arms are same width
sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]
sel.gene.dt[, adjust := R_width - L_width]

sel.gene.dt[adjust > 0, `:=` (R_start = R_start + floor(adjust/2),
                              R_end = R_end - ceiling(adjust/2))]

sel.gene.dt[adjust < 0, `:=` (L_start = L_start + floor(abs(adjust)/2),
                              L_end = L_end - ceiling(abs(adjust)/2))]

sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]

stopifnot(all(sel.gene.dt$L_width == (sel.gene.dt$R_width)))

# Expand for plotting
sel.expanded.gene.dt <- sel.gene.dt[rep(seq(1, nrow(sel.gene.dt)), sel.gene.dt$L_width)]
sel.expanded.gene.dt[, `:=` (L_adjust = 1:.N,
                             R_adjust = .N:1), by = cluster]
sel.expanded.gene.dt[, `:=` (L = L_start + L_adjust - 1,
                             R = R_start + R_adjust - 1)]

setorder(sel.expanded.gene.dt, count)

# p1 <- ggplot(sel.expanded.gene.dt[grepl("^ID", cluster)]) +
#   geom_curve(aes(x = L, xend = R, y = 0, yend = 0, colour = count), lineend = "round", curvature = 1, size = 1, ncp = 50) +
#   coord_cartesian(ylim = c(-10, 0), xlim = c(start(region.gr), end(region.gr))) +
#   scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 36)) +
#   theme_minimal_vgrid() +
#   labs(x = "",
#        colour = "Counts") +
#   theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
#         legend.position = "none")

# Flip as on negative strand

p1 <- ggplot(sel.expanded.gene.dt[grepl("^C", cluster)]) +
  geom_curve(aes(x = R, xend = L, y = 0, yend = 0, colour = count), lineend = "round", curvature = -1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(0, 10), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 28)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top")

p4 <- ggplot(sel.expanded.gene.dt[grepl("C003", cluster)]) +
  geom_curve(aes(x = R, xend = L, y = 0, yend = 0, colour = count), lineend = "round", curvature = -1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(0, 10), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 28)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "top")

paris.dt <- fread(file.path(results.path, "paris/results_paris/atlas_clusters/all.atlas_clusters.gc.annotated.tsv.gz"))
# paris.dt <- fread(file.path(results.path, "paris/results_paris/atlas/all.atlas.clustered.tsv.gz"))
# paris.dt <- fread("~/Dropbox (The Francis Crick)/comp_hiclip/paris_atlas/paris.all.atlas_clusters.gc.annotated.mfe.tsv.gz")

sel.gene.dt <- paris.dt[L_gene_name == R_gene_name][L_gene_name == "SRSF1"]
sel.gene.dt[, `:=` (L_start = L_genomic_start,
                    L_end = L_genomic_end,
                    R_start = R_genomic_start,
                    R_end = R_genomic_end)]

# Adjust so L and R arms are same width
sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]
sel.gene.dt[, adjust := R_width - L_width]

sel.gene.dt[adjust > 0, `:=` (R_start = R_start + floor(adjust/2),
                              R_end = R_end - ceiling(adjust/2))]

sel.gene.dt[adjust < 0, `:=` (L_start = L_start + floor(abs(adjust)/2),
                              L_end = L_end - ceiling(abs(adjust)/2))]

sel.gene.dt[, `:=` (L_width = L_end - L_start + 1,
                    R_width = R_end - R_start + 1)]

stopifnot(all(sel.gene.dt$L_width == (sel.gene.dt$R_width)))

# Expand for plotting
sel.expanded.gene.dt <- sel.gene.dt[rep(seq(1, nrow(sel.gene.dt)), sel.gene.dt$L_width)]
sel.expanded.gene.dt[, `:=` (L_adjust = 1:.N,
                             R_adjust = .N:1), by = cluster]
sel.expanded.gene.dt[, `:=` (L = L_start + L_adjust - 1,
                             R = R_start + R_adjust - 1)]

setorder(sel.expanded.gene.dt, count)

p2 <- ggplot(sel.expanded.gene.dt) +
  geom_curve(aes(x = R, xend = L, y = 0, yend = 0, colour = count), lineend = "round", curvature = 1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(-10, 0), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 28)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")

p3 <- ggplot(sel.expanded.gene.dt[cluster == "C019"]) +
  geom_curve(aes(x = R, xend = L, y = 0, yend = 0, colour = count), lineend = "round", curvature = 1, size = 1, ncp = 50) +
  coord_cartesian(ylim = c(-10, 0), xlim = c(start(region.gr), end(region.gr))) +
  scale_colour_viridis_c(option = "G", direction = -1, limits = c(1, 28)) +
  theme_minimal_vgrid() +
  labs(x = "",
       colour = "Counts") +
  theme(axis.line.y = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none")


# ==========
# xlink tracks
# ==========

SubsetBedgraph <- function(gr, selected.region.gr) {
  
  xlinks.gr <- unlist(tile(selected.region.gr, width = 1))
  
  ol <- findOverlaps(xlinks.gr, gr)
  xlinks.gr$score <- 0
  xlinks.gr[queryHits(ol)]$score <- gr[subjectHits(ol)]$score
  # xlinks.gr$score[is.na(xlinks.gr$score)] <- 0
  
  return(xlinks.gr)
  
}

xlinks.gr <- import.bed(file.path(results.path, "results_nornase/xlinks/stau1.xl.bed.gz"))
seqlevelsStyle(xlinks.gr) <- "UCSC"
xlinks.gr <- SubsetBedgraph(gr = xlinks.gr, selected.region.gr = region.gr)
xlinks.dt <- as.data.table(xlinks.gr)[, sample := "STAU1"]
xlinks.dt[, norm := score]
xlinks.dt[, smoothed := zoo::rollmean(norm, 25, fill = 0)]
y.label <- "Crosslink signal"

p.iclip <- ggplot(xlinks.dt) +
  geom_line(aes(x = start, y = smoothed, group = sample, color = sample)) +
  labs(x = "",
       y = y.label,
       colour = "") +
  scale_colour_tableau(palette = "Tableau 10") +
  theme_minimal_grid() + theme(legend.position = "top") +
  xlim(start(region.gr),end(region.gr))

# ==========
# Annotation
# ==========

# TxDb <- loadDb("~/Dropbox (The Francis Crick)/comp_hiclip/ref/gencode.v33.txdb.sqlite")
# seqlevelsStyle(TxDb) <- "UCSC"
# 
# gtf <- import.gff2("~/Dropbox (The Francis Crick)/comp_hiclip/ref/gencode.v33.annotation.gtf.gz")
# seqlevelsStyle(gtf) <- "UCSC"
# genes.gr <- gtf[gtf$type == "gene"]

# rosetta.dt <- as.data.table(mcols(genes.gr))[, .(gene_id, gene_name)]
# setkey(rosetta.dt, gene_id)

# Get transcripts that overlap region and order for plotting
sel.tx_genes <- transcriptsByOverlaps(TxDb, region.gr, columns = c("gene_id", "tx_name"))
sel.tx_genes <- sel.tx_genes[order(width(sel.tx_genes), decreasing = TRUE)]
sel.tx_genes <- sel.tx_genes[sel.tx_genes$tx_name == "ENST00000258962.5"]
tx.order.dt <- data.table(transcript_id = sel.tx_genes$tx_name,
                          gene_id = unlist(sel.tx_genes$gene_id),
                          centre = start(sel.tx_genes) + width(sel.tx_genes)/2)[, group := 1:.N]
setkey(tx.order.dt, gene_id)
tx.order.dt <- rosetta.dt[tx.order.dt]
tx.order.dt[, gene := paste0(gene_name, " | ", gene_id)]

# Region (for arrows)
region <- region.gr
region$gene_id <- NULL
region.tiled <- rep(region, length(sel.tx_genes))
region.tiled$tx_name <- sel.tx_genes$tx_name
region.tiled <- tile(region.tiled, width = round(width(region)/15, -1))

sel.region.tiled <- GRangesList(lapply(seq_along(sel.tx_genes), function(i) {
  gr <- subsetByOverlaps(region.tiled[[i]], sel.tx_genes[i], type = "any")
  start(gr[1]) <- start(sel.tx_genes[i])
  end(gr[length(gr)]) <- end(sel.tx_genes[i])
  return(gr)
}))

names(sel.region.tiled) <- sel.tx_genes$tx_name
sel.region.tiled.dt <- as.data.table(sel.region.tiled)
sel.region.tiled.dt[, group := NULL]
setnames(sel.region.tiled.dt, "group_name", "transcript_id")
if(as.character(strand(region.gr)) == "-") setnames(sel.region.tiled.dt, c("start", "end"), c("end", "start"))

sel.region.tiled.dt <- merge(sel.region.tiled.dt, tx.order.dt, by = "transcript_id")

# CDS
cds_tx <- cdsBy(TxDb, by = "tx", use.names = TRUE)
sel.cds_tx <- cds_tx[names(cds_tx) %in% sel.tx_genes$tx_name]
sel.cds_tx.dt <- as.data.table(sel.cds_tx)
sel.cds_tx.dt[, group := NULL]
setnames(sel.cds_tx.dt, "group_name", "transcript_id")
sel.cds_tx.dt <- merge(sel.cds_tx.dt, tx.order.dt, by = "transcript_id")

# UTR5
utr5_tx <- fiveUTRsByTranscript(TxDb, use.names = TRUE)
sel.utr5_tx <- utr5_tx[names(utr5_tx) %in% sel.tx_genes$tx_name]
sel.utr5_tx.dt <- as.data.table(sel.utr5_tx)
sel.utr5_tx.dt[, group := NULL]
setnames(sel.utr5_tx.dt, "group_name", "transcript_id") 
sel.utr5_tx.dt <- merge(sel.utr5_tx.dt, tx.order.dt, by = "transcript_id")

# UTR3
utr3_tx <- threeUTRsByTranscript(TxDb, use.names = TRUE)
sel.utr3_tx <- utr3_tx[names(utr3_tx) %in% sel.tx_genes$tx_name]
sel.utr3_tx.dt <- as.data.table(sel.utr3_tx)
sel.utr3_tx.dt[, group := NULL]
setnames(sel.utr3_tx.dt, "group_name", "transcript_id")
sel.utr3_tx.dt <- merge(sel.utr3_tx.dt, tx.order.dt, by = "transcript_id")

# Exons
exons_tx <- exonsBy(TxDb, by = "tx", use.names = TRUE)
sel.exons_tx <- exons_tx[names(exons_tx) %in% sel.tx_genes$tx_name[!sel.tx_genes$tx_name %in% names(sel.cds_tx)]] # Don't want genes with CDS, just e.g. ncRNA
sel.exons_tx.dt <- as.data.table(sel.exons_tx)
sel.exons_tx.dt[, group := NULL]
setnames(sel.exons_tx.dt, "group_name", "transcript_id")
sel.exons_tx.dt <- merge(sel.exons_tx.dt, tx.order.dt, by = "transcript_id")

# Plot
p.annot <- ggplot() +
  geom_segment(data = sel.region.tiled.dt, mapping = aes(x = start, xend = end, y = group, yend = group), arrow = arrow(length = unit(0.1, "cm")), colour = "grey50") +
  geom_rect(data = sel.exons_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
  geom_rect(data = sel.cds_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.25, ymax = group + 0.25, fill = gene)) +
  geom_rect(data = sel.utr5_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
  geom_rect(data = sel.utr3_tx.dt, mapping = aes(xmin = start, xmax = end, ymin = group - 0.15, ymax = group + 0.15, fill = gene)) +
  scale_fill_tableau() +
  theme_minimal_vgrid() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "bottom") +
  labs(x = "",
       y = "",
       fill = "") +
  coord_cartesian(xlim = c(start(region), end(region)))

# ==========
# Assemble
# ==========

ggsave(p1 / p2 / p.iclip / p.annot + plot_layout(heights = c(4, 4, 3, 0.5)), 
       filename = file.path(plot.path, "srsf1_arc.pdf"), 
       width = 297, height = 210, units = "mm")
