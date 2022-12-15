suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(toscatools))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(rslurm))

setDTthreads(8)
set.seed(42)

# ============================== #
# Linker
# ============================== #

setwd("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nolinker/atlas")

hybrids.dt <- fread("all.hybrids.tsv.gz")
hybrids.dt[, `:=` (mfe = NULL,
                   structure = NULL,
                   mean_shuffled_mfe = NULL,
                   sd_shuffled_mfe = NULL,
                   cluster = NULL,
                   cluster_hybrid_count = NULL,
                   sample = tstrsplit(sample, "\\.")[[1]])]
atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
atlas.hybrids.dt <- atlas.hybrids.dt[total_count > 1]

# Keep ones not clustered to add back in later
unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt$name]
stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

# Split into list to parallelise
atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))

# Cluster
atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = 0.5, fraction = TRUE, mc.cores = 8)
clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)

# Merge back and sanity check
clusters.dt <- rbindlist(list(clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)
setorder(clusters.dt, name)

stopifnot(nrow(clusters.dt) == nrow(hybrids.dt))
stopifnot(all(clusters.dt$name %in% hybrids.dt$name))

fwrite(clusters.dt, paste0("all_inc", ".", "atlas", ".clustered.tsv.gz"), sep = "\t")

# ==========
# Get structures
# ==========

sel.clusters.dt <- clusters.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
sel.clusters.dt <- sel.clusters.dt[!is.na(cluster)][grep("^C", cluster)]

genome.fa <- Biostrings::readDNAStringSet("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/GRCh38.gencode_v33.fa")
genome.dt <- data.table(gene_id = names(genome.fa),
                        sequence = as.character(genome.fa))
sel.clusters.dt  <- get_sequence(hybrids.dt = sel.clusters.dt , genome.dt = genome.dt)
stopifnot(!any(is.na(c(sel.clusters.dt$L_sequence, sel.clusters.dt$R_sequence))))

structure.list <- parallel::mclapply(seq_len(nrow(sel.clusters.dt)), function(i) {

        analyse_structure(name = sel.clusters.dt$name[i], L_sequence = sel.clusters.dt$L_sequence[i], R_sequence = sel.clusters.dt$R_sequence[i])
    
}, mc.cores = 8)
structure.dt <- rbindlist(structure.list, use.names = TRUE)

sjob <- slurm_apply(get_shuffled_mfe, sel.clusters.dt[, .(name, L_sequence, R_sequence)], 
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

# shuffled.list <- parallel::mclapply(seq_len(nrow(sel.clusters.dt)), function(i) {

#       get_shuffled_mfe(name = sel.clusters.dt$name[i], L_sequence = sel.clusters.dt$L_sequence[i], R_sequence = sel.clusters.dt$R_sequence[i])

# }, mc.cores = 8)
# shuffled.dt <- rbindlist(shuffled.list, use.names = TRUE)

structures.hybrids.dt <- merge(clusters.dt, structure.dt, by = "name", all.x = TRUE)
shuffled.hybrids.dt <- merge(structures.hybrids.dt, shuffled.dt, by = "name", all.x = TRUE)
fwrite(shuffled.hybrids.dt, paste0("all_inc", ".", "atlas", ".gc.annotated.mfe.shuffled.tsv.gz"), sep = "\t")

# ============================== #
# Linker
# ============================== #

setwd("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_linker")

hybrids.dt <- rbindlist(list(fread("stau1_high.hybrids.gc.annotated.tsv.gz")[, sample := "stau1_high"],
                             fread("stau1_low.hybrids.gc.annotated.tsv.gz")[, sample := "stau1_low"]),
                       use.names = TRUE)
atlas.hybrids.dt <- hybrids.dt[, total_count := .N, by = .(L_seqnames, R_seqnames)]
atlas.hybrids.dt <- atlas.hybrids.dt[total_count > 1]

# Keep ones not clustered to add back in later
unclustered.hybrids.dt <- hybrids.dt[!name %in% atlas.hybrids.dt$name]
stopifnot(nrow(unclustered.hybrids.dt) + nrow(atlas.hybrids.dt) == nrow(hybrids.dt))

# Split into list to parallelise
atlas.hybrids.list <- split(atlas.hybrids.dt, by = c("L_seqnames", "R_seqnames"))

# Cluster
atlas.clusters.list <- parallel::mclapply(atlas.hybrids.list, cluster_hybrids, percent_overlap = 0.5, fraction = TRUE, mc.cores = 8)
clusters.dt <- rbindlist(atlas.clusters.list, use.names = TRUE, fill = TRUE)

# Merge back and sanity check
clusters.dt <- rbindlist(list(clusters.dt, unclustered.hybrids.dt), use.names = TRUE, fill = TRUE)
setorder(clusters.dt, name)

stopifnot(nrow(clusters.dt) == nrow(hybrids.dt))
stopifnot(all(clusters.dt$name %in% hybrids.dt$name))

fwrite(clusters.dt, paste0("all_inc", ".", "atlas", ".clustered.tsv.gz"), sep = "\t")

# ==========
# Get structures
# ==========

sel.clusters.dt <- clusters.dt[!(L_seqnames == "rDNA" & R_seqnames == "rDNA")]
sel.clusters.dt <- sel.clusters.dt[!is.na(cluster)][grep("^C", cluster)]

genome.fa <- Biostrings::readDNAStringSet("/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/GRCh38.gencode_v33.fa")
genome.dt <- data.table(gene_id = names(genome.fa),
                        sequence = as.character(genome.fa))
sel.clusters.dt  <- get_sequence(hybrids.dt = sel.clusters.dt , genome.dt = genome.dt)
stopifnot(!any(is.na(c(sel.clusters.dt$L_sequence, sel.clusters.dt$R_sequence))))

structure.list <- parallel::mclapply(seq_len(nrow(sel.clusters.dt)), function(i) {

        analyse_structure(name = sel.clusters.dt$name[i], L_sequence = sel.clusters.dt$L_sequence[i], R_sequence = sel.clusters.dt$R_sequence[i])
    
}, mc.cores = 8)
structure.dt <- rbindlist(structure.list, use.names = TRUE)

shuffled.list <- parallel::mclapply(seq_len(nrow(sel.clusters.dt)), function(i) {

      get_shuffled_mfe(name = sel.clusters.dt$name[i], L_sequence = sel.clusters.dt$L_sequence[i], R_sequence = sel.clusters.dt$R_sequence[i])

}, mc.cores = 8)
shuffled.dt <- rbindlist(shuffled.list, use.names = TRUE)

structures.hybrids.dt <- merge(clusters.dt, structure.dt, by = "name", all.x = TRUE)
shuffled.hybrids.dt <- merge(structures.hybrids.dt, shuffled.dt, by = "name", all.x = TRUE)
fwrite(shuffled.hybrids.dt, paste0("all_inc", ".", "atlas", ".gc.annotated.mfe.shuffled.tsv.gz"), sep = "\t")

