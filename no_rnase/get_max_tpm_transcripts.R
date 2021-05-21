#!/usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


# ==========
# Load annotations and RNAseq data
# ==========

genomes.dir <- "/camp/lab/luscomben/home/users/iosubi/genomes"
data.dir <- "/camp/lab/luscomben/home/users/chakraa2/projects/comp_hiclip/rnaseq/results/salmon"

if (file.exists(paste0(genomes.dir, "/gencode_V33_txdb.sqlite"))) {
  TxDb <- loadDb(paste0(genomes.dir, "/gencode_V33_txdb.sqlite"))
  TxDb <- keepStandardChromosomes(TxDb , pruning.mode="coarse")
} else {
  TxDb <- makeTxDbFromGFF(paste0(genomes.dir, "/gencode.v33.primary_assembly.annotation.gtf.gz"), format="gtf",
                          organism = "Homo sapiens", chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene) )
  saveDb(TxDb, file=paste0(genomes.dir, "/gencode_V33_txdb.sqlite"))
  TxDb <- loadDb(paste0(genomes.dir, "/gencode_V33_txdb.sqlite"))
}


rnaseq.dt <- fread(paste0(data.dir, "/salmon.merged.transcript_tpm.tsv.gz"))


# ==========
# Get all transcripts
# ==========

transcripts.df <- as.data.frame(transcriptsBy(TxDb))
gene_tx_map.df <- transcripts.df %>%
  dplyr::select(group_name, tx_name)

# ==========
# Expression profile
# ==========

rnaseq.dt <- rnaseq.dt %>%
  dplyr::filter(ERX575536_R1 != 0 & ERX575537_R1 != 0 & ERX575538_R1 != 0)
rnaseq.dt <- left_join(rnaseq.dt, gene_tx_map.df, by = c("transcript_id" = "tx_name"))
rnaseq.dt <- dplyr::filter(rnaseq.dt, !is.na(group_name))
rnaseq.dt <- dplyr::rename(rnaseq.dt, gene_id = group_name)
rnaseq.dt <- dplyr::select(rnaseq.dt, gene_id, transcript_id, ERX575536_R1, ERX575537_R1, ERX575538_R1)

rnaseq.dt$tpm_mean <- rowMeans(subset(rnaseq.dt, select = c(ERX575536_R1, ERX575537_R1, ERX575538_R1)), na.rm = TRUE)

# ==========
# For each gene, get the transcript with max TPM
# ==========

max_tpm.df <- rnaseq.dt %>%
  group_by(gene_id) %>%
  dplyr::slice(which.max(tpm_mean)) #in case of ties, first transcript is taken; with this data there were no ties

max_tpm_transcripts.df <- semi_join(transcripts.df, max_tpm.df, by = c("tx_name" = "transcript_id"))
max_tpm_transcripts.df <- max_tpm_transcripts.df %>%
  dplyr::select(seqnames, start, end, width, strand, tx_name, group_name) %>%
  dplyr::rename(gene_id = group_name)

max_tpm_transcripts.gr <- GRanges(max_tpm_transcripts.df)
max_tpm_transcripts.gr$name <- max_tpm_transcripts.gr$tx_name

export.bed(max_tpm_transcripts.gr, con="gencode_v33_max_tpm_transcripts.bed.gz")
