#!/usr/bin/env Rscript

# Script to get human custom fasta
# A. M. Chakrabarti
# 29th October 2020
# Updated 28th December 2020 to merge overlapping genes
# Updated 28th January 2021 to add tRNA

suppressPackageStartupMessages(library(rentrez))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(rtracklayer))

# setwd("/camp/lab/luscomben/home/users/chakraa2/projects/flora/ref/human")

# =========
# Get rRNA
# =========

get_ncbi_sequence <- function(accession) {
  
  fa <- entrez_fetch(db="nucleotide", id=accession, rettype="fasta")
  lines <- str_split(fa, "\n")[[1]]
  id <- gsub(">", "", lines[1])
  seq <- paste0(lines[-1], collapse = "")
  
  # Convert to DNAStringSet
  dss <- DNAStringSet(seq)
  names(dss) <- id
  
  return(dss)
  
}

rDNA.ds <- get_ncbi_sequence("U13369.1") # NCBI rDNA
rRNA_5S.ds <- get_ncbi_sequence("NR_023363.1") # NCBI rRNA 5S

rrna.ds <- c(rDNA.ds, rRNA_5S.ds)
names(rrna.ds) <- c("rDNA", "rRNA5S")

# =========
# Get tRNA
# =========

trna.ds <- DNAStringSet(readRNAStringSet("hg38-mature-tRNAs.fa"))
names(trna.ds) <- gsub("^Homo_sapiens_", "", sapply(strsplit(names(trna.ds), " "), "[[", 1))

# =========
# Get genes
# =========

ref.gff <- "gencode.v33.annotation.gff3.gz"

if(!file.exists(ref.gff)) {

    download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz", 
                destfile = ref.gff,
                method = "wget",
                quiet = TRUE)

}

# Load annotation
genes.gr <- import.gff3(ref.gff)
genes.gr <- keepStandardChromosomes(genes.gr, pruning.mode = "coarse")
genes.gr <- dropSeqlevels(genes.gr, "chrY", pruning.mode = "coarse")
genes.gr <- genes.gr[genes.gr$type == "gene"]

# Get protein coding genes
biotypes <- unique(genes.gr$gene_type)
pc <- grep("protein_coding|IG_[A-Z]_gene|TR_[A-Z]_gene", biotypes, value = TRUE)
sel.biotypes <- c(pc, "lncRNA", "vault_RNA")
sel.genes.gr <- genes.gr[genes.gr$gene_type %in% sel.biotypes]
sel.genes.gr$ol <- countOverlaps(sel.genes.gr, drop.self = TRUE)

# Collapse overlapping
unique.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol == 0]
mcols(unique.sel.genes.gr) <- mcols(unique.sel.genes.gr)[, c("gene_name", "gene_id", "gene_type")]

multi.sel.genes.gr <- sel.genes.gr[sel.genes.gr$ol > 0]
multi.reduced.sel.genes.gr <- reduce(multi.sel.genes.gr, with.revmap = TRUE, min.gapwidth = 0) # gapwidth so doesn't merge immediately juxtaposed ranges

# Stitch together names
revmap <- mcols(multi.reduced.sel.genes.gr)$revmap
multi.sel.genes.gr.grl <- relist(multi.sel.genes.gr[unlist(revmap)], revmap)

multi.reduced.sel.genes.gr$gene_name <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$gene_name), collapse = "_"))
multi.reduced.sel.genes.gr$gene_id <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(x$gene_id), collapse = "_"))
multi.reduced.sel.genes.gr$gene_type <- sapply(multi.sel.genes.gr.grl, function(x) paste0(sort(unique(x$gene_type)), collapse = "-"))
multi.reduced.sel.genes.gr$revmap <- NULL

reduced.sel.genes.gr <- sort(c(unique.sel.genes.gr, multi.reduced.sel.genes.gr))
reduced.sel.genes.gr$name <- paste0(reduced.sel.genes.gr$gene_name, ":", reduced.sel.genes.gr$gene_id)

invisible(file.remove(ref.gff))

# =========
# Mask snRNAs, rRNA and tRNA
# =========

# snRNA and rRNA
snrna.gr <- genes.gr[genes.gr$gene_type %in% c("rRNA", "rRNA_pseudogene", "snRNA")]
mcols(snrna.gr) <- NULL

# tRNA
# Get tRNA bed from gtRNAdb
if(!file.exists("hg38-tRNAs.bed")) {

  download.file(url = "http://gtrnadb.ucsc.edu/genomes/eukaryota/Hsapi38/hg38-tRNAs.tar.gz",
                destfile = "hg38-tRNAs.tar.gz",
                method = "wget",
                quiet = TRUE)
  
  untar(tarfile = "hg38-tRNAs.tar.gz",
        files = "hg38-tRNAs.bed")
  
  invisible(file.remove("hg38-tRNAs.tar.gz"))

}

trna.gr <- import.bed("hg38-tRNAs.bed")
trna.gr <- keepStandardChromosomes(trna.gr, pruning.mode = "coarse")
mcols(trna.gr) <- NULL

mask.gr <- c(snrna.gr, trna.gr)
mask.bed <- tempfile(tmpdir = ".", fileext = ".bed")
export.bed(mask.gr, mask.bed)

# =========
# Create custom reference fasta
# =========

ref.fasta <- "GRCh38.primary_assembly.genome.fa"
if(!file.exists(ref.fasta)) {

    download.file(url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/GRCh38.primary_assembly.genome.fa.gz", 
                destfile = "GRCh38.primary_assembly.genome.fa.gz",
                method = "wget",
                quiet = TRUE)

    system("pigz -d GRCh38.primary_assembly.genome.fa.gz")

}


masked.fasta <- tempfile(tmpdir = ".", fileext = ".fa")
cmd <- paste("bedtools maskfasta -fi", ref.fasta, "-bed", mask.bed, "-fo", masked.fasta)
# message(cmd)
system(cmd)

invisible(file.remove(mask.bed))

reduced.sel.genes.bed <- tempfile(tmpdir = ".", fileext = ".bed")
export.bed(reduced.sel.genes.gr, reduced.sel.genes.bed)

cmd <- paste("samtools faidx", masked.fasta)
# message(cmd)
system(cmd)

cmd <- paste("bedtools getfasta -s -name -fi", masked.fasta, "-bed", reduced.sel.genes.bed, "| pigz > human.fa.gz")
# message(cmd)
system(cmd)

invisible(file.remove(reduced.sel.genes.bed))
invisible(file.remove(ref.fasta))
invisible(file.remove(masked.fasta))
invisible(file.remove(paste0(masked.fasta, ".fai")))

# Add rRNA and tRNA
writeXStringSet(rrna.ds, "human.fa.gz", append = TRUE, compress = TRUE)
writeXStringSet(trna.ds, "human.fa.gz", append = TRUE, compress = TRUE)

message("Completed")