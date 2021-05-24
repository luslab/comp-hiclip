#!/usr/bin/env Rscript

library(rtracklayer)
library(GenomicFeatures)
library(dplyr)
library(data.table)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(optparse)



# =========
# Options and paths
# =========

option_list <- list(make_option(c("", "--gtf"), action = "store", type = "character", default=NA, help = "GTF annotation file"),
                    make_option(c("", "--out"), action = "store", type = "character", default=NA, help = "Output name"))
opt_parser = OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

data.dir <- "/camp/lab/luscomben/home/shared/projects/ira-nobby/comp_hiclip/results_nonhybrid/"


# ==========
# Load annotations
# ==========

gtf <- import.gff2(opt$gtf)


if (file.exists(paste0(data.dir, "/gencode.v33.txdb.sqlite"))) {
  TxDb <- loadDb(paste0(data.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- keepStandardChromosomes(TxDb , pruning.mode="coarse")
} else {
  TxDb <- makeTxDbFromGFF(opt$gtf, format="gtf",
                          organism = "Homo sapiens", chrominfo = seqinfo(TxDb.Hsapiens.UCSC.hg38.knownGene) )
  saveDb(TxDb, file=paste0(data.dir, "/gencode.v33.txdb.sqlite"))
  TxDb <- loadDb(paste0(data.dir, "/gencode.v33.txdb.sqlite"))
}


# ==========
# Obtain longest protein coding transcripts
# ==========


txlengths <- transcriptLengths(TxDb, with.cds_len = TRUE,
                                  with.utr5_len = TRUE,
                                  with.utr3_len = TRUE)

txlengths.dt <- data.table(txlengths, key = c("tx_name", "gene_id"))

pc = c("protein_coding", "IG_V_gene", "TR_V_gene", "IG_C_gene", "IG_J_gene", "TR_J_gene", "TR_C_gene", "IG_D_gene", "TR_D_gene")


gtf.df <- as.data.frame(gtf)
gtf.dt <- data.table(gtf.df, key = c("transcript_id", "gene_id"))
gtf.dt <- gtf.dt[txlengths.dt]

longest.pc.dt <- gtf.dt[gene_type %in% pc & transcript_type %in% pc, longest := max(tx_len), by = gene_id] # select out where both are protein coding as sometimes a processed transcript is the longest

longest.pc.dt <- longest.pc.dt[gene_type %in% pc & transcript_type %in% pc & tx_len == longest] # selects longest
longest.pc.df <- longest.pc.df %>% arrange(desc(longest), desc(nexon), desc(utr3_len), desc(cds_len), desc(utr5_len))

unique.longest.pc.df <- longest.pc.df[ !duplicated(longest.pc.df$gene_id), ] 

fwrite(unique.longest.pc.df, opt$out, sep = "\t")
