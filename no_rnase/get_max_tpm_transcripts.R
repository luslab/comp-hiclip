#!/usr/bin/env Rscript

library(rtracklayer)
library(dplyr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(stringr)

# =========
# Options and paths
# =========

option_list <- list(make_option(c("", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"),
                    make_option(c("", "--tpm"), action = "store", type = "character", default=NA, help = "Transcript quantification table"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)


# ==========
# Load annotations and RNAseq data
# ==========

#ref.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref"
data.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/"

if (file.exists(paste0(data.dir, "/gencode.v33.txdb.sqlite"))) {
  TxDb <- loadDb(paste0(data.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- keepStandardChromosomes(TxDb , pruning.mode="coarse")
} else {
  TxDb <- makeTxDbFromGFF(opt$gtf, format="gtf",
                          organism = "Homo sapiens", chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene) )
  saveDb(TxDb, file=paste0(data.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- loadDb(paste0(data.dir, "/gencode.v33.txdb.sqlite"))
}

# TPM data 
rnaseq.dt <- fread(opt$tpm)

# ==========
# Get all transcripts
# ==========

transcripts.df <- as.data.frame(transcriptsBy(TxDb))

# ==========
# Expression profile
# ==========

samples.ls <- colnames(rnaseq.dt)[str_detect(colnames(rnaseq.dt), "ERR")]

rnaseq.dt <- rnaseq.dt %>%
  dplyr::filter(if_any(matches("ERR"), ~ . != 0)) # filter out rows with all samples with TPM = 0
  
rnaseq.dt <- dplyr::select(rnaseq.dt, tx, gene_id, samples.ls)
rnaseq.dt$tpm_mean <- rowMeans(subset(rnaseq.dt, select = samples.ls), na.rm = TRUE)

# ==========
# For each gene, get the transcript with max TPM
# ==========

max_tpm.df <- rnaseq.dt %>%
  group_by(gene_id) %>%
  dplyr::slice(which.max(tpm_mean)) #in case of ties, first transcript is taken; with this data there were no ties

max_tpm_transcripts.df <- semi_join(transcripts.df, max_tpm.df, by = c("tx" = "transcript_id"))
max_tpm_transcripts.df <- max_tpm_transcripts.df %>%
  dplyr::select(seqnames, start, end, width, strand, tx_name, group_name) %>%
  dplyr::rename(gene_id = group_name)

max_tpm_transcripts.gr <- GRanges(max_tpm_transcripts.df)
max_tpm_transcripts.gr$name <- max_tpm_transcripts.gr$tx_name

export.bed(max_tpm_transcripts.gr, con=paste0(data.dir,"/gencode.v33_max_tpm_transcripts.bed.gz"))
