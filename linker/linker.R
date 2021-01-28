#!/usr/bin/env Rscript

# Script to extract linker hybrids
# A. M. Chakrabarti
# 28th January 2021 (modified from hiclipr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 2) stop("Run command: linker.R <input fastq> <output fastq>")
if(file.exists(args[2])) stop("Output file already exists.")

suppressPackageStartupMessages(library(ShortRead))
suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(GenomicRanges))

# setwd("~/projects/comp_hiclip")

split_read <- function(fq, adapter) {
  
  # Select reads with the adapter
  B <- vcountPattern(adapter, subject = sread(fq), max.mismatch = 2, with.indels = FALSE)
  fq.nolinker <- fq[B == 0]
  fq <- fq[B > 0]
  
  # Extract components of ShortRead
  reads.ds <- sread(fq)
  quality.bs <- quality(fq)@quality # extract quality as a BString
  names(reads.ds) <- id(fq)
  names(quality.bs) <- id(fq)
  
  # Find locations of adapter B and reduce instances of BB/BBB
  hybrid.ir <- unlist(vmatchPattern(adapter, subject = reads.ds, max.mismatch = 2, with.indels = FALSE)) # can unlist and still know which read it refers to as has name
  hybrid.irl <- split(hybrid.ir, names(hybrid.ir)) # split into list by read
  hybrid.ir <- unlist(reduce(hybrid.irl)) # reduce polyB and unlist back
  
  # Remove those that are ...B...B...
  duplicated <- names(hybrid.ir)[duplicated(names(hybrid.ir))] # identify those that are ...B...B... and remove from ongoing analysis (may remove a few hybrids that are ...B...B, but these are very few ~150)
  hybrid.ir <- hybrid.ir[!(names(hybrid.ir) %in% duplicated), ] # so these are the ranges of B for reads that are ...B...A or ...B...BA
  reads.ds <- reads.ds[names(reads.ds) %in% names(hybrid.ir), ] # and these are the reads for the above hybrids
  quality.bs <- quality.bs[names(quality.bs) %in% names(hybrid.ir)] # and there are the qualities
  
  # Order reads
  hybrid.ir <- hybrid.ir[order(names(hybrid.ir))] # order both by read name for trimming step
  reads.ds <- reads.ds[order(names(reads.ds))]
  quality.bs <- quality.bs[order(names(quality.bs))]
  
  #Remove reads that are too short for trimming i.e. adapter B at 3' end
  tooshort <- end(hybrid.ir) > width(reads.ds) # need to remove reads that are too short
  hybrid.ir <- hybrid.ir[!tooshort]
  reads.ds <- reads.ds[!(tooshort)]
  quality.bs <- quality.bs[!(tooshort)]
  
  # Remove reads that have a negative start (i.e. adapter B overhang)
  overhang <- start(hybrid.ir) < 1
  hybrid.ir <- hybrid.ir[!overhang]
  reads.ds <- reads.ds[!(overhang)]
  quality.bs <- quality.bs[!(overhang)]
  
  # Create R and L arms of hybrid reads
  R.reads.ds <- DNAStringSet(reads.ds, start = 1, end = (start(hybrid.ir)-1)) # R arm is from start of read to start of adapter B
  R.quality.bs <- BStringSet(quality.bs, start = 1, end = (start(hybrid.ir)-1))
  L.reads.ds <- DNAStringSet(reads.ds, start = (end(hybrid.ir)+1)) # L arm is from end of adapter B to end of read
  L.quality.bs <- BStringSet(quality.bs, start = (end(hybrid.ir)+1))
  
  # Only keep those that have both R and L arms > 16 nt
  R.reads.ds <- R.reads.ds[width(R.reads.ds) > 12]
  L.reads.ds <- L.reads.ds[width(L.reads.ds) > 12]
  R.reads.ds <- R.reads.ds[names(R.reads.ds) %in% names(L.reads.ds), ]
  L.reads.ds <- L.reads.ds[names(L.reads.ds) %in% names(R.reads.ds), ]
  
  # Get corresponding qualities
  R.quality.bs <- R.quality.bs[names(R.quality.bs) %in% names(R.reads.ds)]
  L.quality.bs <- L.quality.bs[names(L.quality.bs) %in% names(L.reads.ds)]
  
  # Create ShortReads
  if(length(R.reads.ds) == 0 | length(L.reads.ds) == 0) {
    
    message("No linker adapter hybrids remaining after filtering for read length")
    
  } else {
    
    R.hybrid <- ShortReadQ(sread = R.reads.ds[order(names(R.reads.ds))],
                           id = BStringSet(paste0("R.", names(R.reads.ds)[order(names(R.reads.ds))])),
                           quality = R.quality.bs[order(names(R.reads.ds))])
    L.hybrid <- ShortReadQ(sread = L.reads.ds[order(names(L.reads.ds))],
                           id = BStringSet(paste0("L.", names(L.reads.ds)[order(names(L.reads.ds))])),
                           quality = L.quality.bs[order(names(L.reads.ds))])
    hybrid <- append(R.hybrid, L.hybrid)
   
  }

  return(list(hybrid = hybrid,
              nolinker = fq.nolinker))
    
}

linker <- "CTGTAGGCACCATACAATG"

args <- commandArgs(trailingOnly = TRUE)

if(file.exists(args[2])) stop("Output file already exists.")

fq <- readFastq(args[1])
B.hybrid <- split_read(fq = fq, adapter = linker)
B_minus1.hybrid <- split_read(fq = B.hybrid$nolinker, adapter = substr(linker, 1, nchar(linker) - 1))
B_minus2.hybrid <- split_read(fq = B_minus1.hybrid$nolinker, adapter = substr(linker, 1, nchar(linker) - 2))

hybrid.fq <- append(B.hybrid$hybrid, 
                    B_minus1.hybrid$hybrid,
                    B_minus2.hybrid$hybrid)

writeFastq(hybrid.fq, compress = TRUE, file = args[2])