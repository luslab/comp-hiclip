#!/usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)
library(optparse)


# =========
# Options and paths
# =========

option_list <- list(make_option(c("", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"),
                    make_option(c("", "--tpm"), action = "store", type = "character", default=NA, help = "Transcript quantification table"),
                    make_option(c("", "--out"), action = "store", type = "character", default=NA, help = "Output name for bed file with tranascript coordinates"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# ==========
# Load annotations and RNAseq data
# ==========

ref.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref"
# data.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions/results_nornase"


# Generate TxDb object
if (file.exists(paste0(ref.dir, "/gencode.v33.txdb.sqlite"))) {
  TxDb <- loadDb(paste0(ref.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- keepStandardChromosomes(TxDb , pruning.mode="coarse")
} else {
  TxDb <- makeTxDbFromGFF(opt$gtf, format="gtf",
                          organism = "Homo sapiens", chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene) )
  saveDb(TxDb, file=paste0(ref.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- loadDb(paste0(ref.dir, "/gencode.v33.txdb.sqlite"))
}

# TPM data from nf-core rnaseq pipeline run
rnaseq.dt <- fread(opt$tpm)

# ==========
# Expression profile
# ==========

samples.ls <- colnames(rnaseq.dt)[str_detect(colnames(rnaseq.dt), "untreated")]
rnaseq.dt <- dplyr::select(rnaseq.dt, c(tx, gene_id, all_of(samples.ls)))

rnaseq.dt <- rnaseq.dt %>%
  dplyr::filter(if_any(matches("untreated"), ~ . != 0)) # filter out rows with all samples with TPM = 0
  

rnaseq.dt$tpm_mean <- rowMeans(subset(rnaseq.dt, select = samples.ls), na.rm = TRUE)

# ==========
# For each gene, get the transcript with max TPM
# ==========

#  Get all transcripts from TxDb 
transcripts.df <- as.data.frame(transcriptsBy(TxDb))

max_tpm.df <- rnaseq.dt %>%
  group_by(gene_id) %>%
  dplyr::slice(which.max(tpm_mean)) # in case of ties, first transcript is taken; with this dataset there were no ties

# Get coordinates of max TPM transcripts
max_tpm_transcripts.df <- semi_join(transcripts.df, max_tpm.df, by = c("tx_name" = "tx"))
max_tpm_transcripts.df <- max_tpm_transcripts.df %>%
  dplyr::select(seqnames, start, end, width, strand, tx_name, group_name) %>%
  dplyr::rename(gene_id = group_name)

max_tpm_transcripts.gr <- GRanges(max_tpm_transcripts.df)
max_tpm_transcripts.gr$name <- max_tpm_transcripts.gr$tx_name

export.bed(max_tpm_transcripts.gr, con=opt$out)
