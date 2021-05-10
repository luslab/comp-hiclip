#!/usr/bin/env Rscript 

library(stringr)
library(rtracklayer)
library(GenomicFeatures)
library(dplyr)

# ==========
# Define functions
# ==========

resize_peaks <- function(bedfiles.list, left = 100, right = 100) {
  
  w <- left + right +1  # width of interval: xl site + flanks
  
  grl <- GRangesList(lapply(bedfiles.list, import.bed))
  gr <- unlist(grl)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- dropSeqlevels(gr, c("chrM", "chrY"), pruning.mode = "coarse")
  gr.df <- as.data.frame(gr)
  
  gr <- resize(gr, width = 1, fix = "start") # resize the peaks, start of peak = 1
  gr <- unique(gr)  # keep unique xl positions
  gr <- resize(resize(gr, width = right+1, fix = "start"), width = w, fix = "end") # add +/- flanks
  gr$id <- paste0("ID", seq(1,length(gr)))
  
  resized.gr.df <- as.data.frame(gr)

  return(gr)
  
}


merge_by_partial_string <- function(df1, df2) {
  
  df1$matchID = row.names(df1)
  df2$matchID = sapply(df2$gene_name, function(x) grep(x, df1$fasta_id)) # not very fast but ok
  
  df_merge = merge(df1, df2, by = "matchID")[-1]
  
  return(df_merge)

  }


# ==========
# Data and parameters
# ==========

data.dir <- "/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt"
first_stem_loop.df <- read.csv(paste0(data.dir, "/stau1.peaks.10nt_10nt_forgi_analyses.df.txt"), sep = "\t")
bed.files <- list.files(path = data.dir, pattern = "annot.bed.gz", full.names = TRUE)


# Params
xl <- 1 # crosslink site position 
prefix = "stau1.peaks.10nt_10nt"

# Annotation
txdb <- paste0(data.dir,"/gencode_V33_txdb.sqlite")
gencode.txdb <- AnnotationDbi::loadDb(txdb)
human.gtf <- paste0(data.dir, "/human_GencodeV33.gtf") # contains fasta_id used by Tosca

message(paste0("Analysing ", prefix, " data"))

# ==========
# Convert forgi stem-loops data to a format that is similar to the hybrids and cluster tables from Tosca
# ==========

# ==========
# Prepare gene IDs and coordinates
# ==========

# Load peaks data, resize exactly as for rnafold used to predict structures
peaks.gc.gr <- resize_peaks(bed.files, left = 0, right = 100)

# Add a new column containing gene ids by directly mapping the transcript ids
keys = peaks.gc.gr$name
mapped.ids <- AnnotationDbi::select(gencode.txdb, keys = keys, columns="GENEID", keytype="CDSNAME")
peaks.gc.gr$gene_name <- mapped.ids$GENEID

# Add fasta_id from annotation used for Tosca as a new column
human.gr <- import.gff2(human.gtf)
human.gr <- keepStandardChromosomes(human.gr, pruning.mode = "coarse")
human.gr <- dropSeqlevels(human.gr, c("chrM", "chrY"), pruning.mode = "coarse")
human.df <- as.data.frame(human.gr)

genes.ls <- unique(peaks.gc.gr$gene_name) # create list of unique gene names from peak file to filter the human annotation

human.df <- human.df %>%
  filter(str_detect(fasta_id, paste(genes.ls, collapse = "|"))) # filter for gene names in the peak gr

human.df$row <- row.names(human.df)

human.df <- dplyr::rename(human.df, fasta_id_start = start, fasta_id_end = end)

human.df <- as.data.frame(human.df[,c("row", "fasta_id", "fasta_id_start", "fasta_id_end")])

peaks.gc.df <- as.data.frame(peaks.gc.gr)

peaks.gc.df <- merge_by_partial_string(human.df, peaks.gc.df)
peaks.gc.df <- peaks.gc.df %>%
  dplyr::select(-row)

# Add transcriptomic coordinates
peaks.df <- peaks.gc.df %>%
  mutate(tx_start = case_when((strand == "+") ~ start - fasta_id_start + 1,
                                     (strand == "-") ~  fasta_id_end - end + 1),
         tx_end = case_when((strand == "+") ~ end - fasta_id_start + 1,
                                   (strand == "-") ~  fasta_id_end - start + 1))

# ==========
# Load predicted stem-loops arms data
# ==========

# Extract the start and end position of the stem-loop individually for each arm
arms.df <- first_stem_loop.df %>% 
  distinct(id, L_min, L_max, R_min, R_max, L_seq, R_seq, L_db, R_db)

# Merge df containing the peaks starts + 100 in genomic and transcriptomic coordinates to arms dataframe
arms.df <- as.data.frame(left_join(arms.df, peaks.df, by = "id"))

# ==========
# Obtain coordinates for each arm of the stem loop
# ==========

# Calculate L and R arms genomic coordinates
arms.df <- arms.df %>%
  mutate(L_genomic_start = case_when((strand == "+") ~ start + L_min - 1,
                                     (strand == "-") ~ end - L_max + 1),
         L_genomic_end = case_when((strand == "+") ~ start + L_max - 1,
                                   (strand == "-") ~ end - L_min + 1),
         R_genomic_start = case_when((strand == "+") ~ start + R_min - 1,
                                     (strand == "-") ~ end - R_max + 1),
         R_genomic_end = case_when((strand == "+") ~ start + R_max - 1,
                                   (strand == "-") ~ end - R_min + 1))

# Calculate L and R arms transcriptomic coordinates
arms.df <- arms.df %>%
    mutate(L_tx_start = tx_start + L_min - 1,
           L_tx_end = tx_start + L_max - 1,
           R_tx_start = tx_start + R_min - 1,
           R_tx_end = tx_start + R_max - 1)

# Check L and R arm coordinate widths are correct
stopifnot(abs(arms.df$L_genomic_end - arms.df$L_genomic_start) == 
            abs(arms.df$L_tx_end - arms.df$L_tx_start))
stopifnot(abs(arms.df$L_max - arms.df$L_min) == 
            abs(arms.df$L_genomic_end - arms.df$L_genomic_start))

stopifnot(abs(arms.df$R_genomic_end - arms.df$R_genomic_start) == 
            abs(arms.df$R_tx_end - arms.df$R_tx_start))
stopifnot(abs(arms.df$R_max - arms.df$R_min) == 
            abs(arms.df$R_genomic_end - arms.df$R_genomic_start))


# Write out tables

# Transcriptomic (similar to Tosca hybrids table):
forgi.reformated.tx.df <- data.frame(id=arms.df$id, name=arms.df$name, orientation = NA, type = "intragenic", hybrid_selection = NA, L_seqnames = arms.df$fasta_id,
                                     L_read_start = arms.df$L_min, L_read_end = arms.df$L_max, L_start = arms.df$L_tx_start, L_end = arms.df$L_tx_end, L_width = NA, L_strand = "+", 
                                     R_seqnames = arms.df$fasta_id, R_read_start = arms.df$R_min, R_read_end = arms.df$R_max, R_start = arms.df$R_tx_start, R_end = arms.df$R_tx_end, R_width = NA, R_strand = "+",
                                     umi= NA, L_sequence=arms.df$L_seq, R_sequence=arms.df$R_seq, mfe = NA)

forgi.reformated.tx.df <- forgi.reformated.tx.df %>%
  rowwise() %>%
  mutate(L_width = L_end - L_start + 1,
         R_width = R_end - R_start +1,
         L_sequence = str_replace_all(L_sequence, "U", "T"), # replace Us with Ts
         R_sequence = str_replace_all(R_sequence, "U", "T")) %>%
  ungroup()


# Genomic (similar to tosca clusters table):
forgi.reformated.gc.df <- data.frame(name = arms.df$id, count = NA, L_seqnames=arms.df$fasta_id,	
                                     L_start=arms.df$L_tx_start, L_end=arms.df$L_tx_end, L_strand="+", L_genomic_seqnames=arms.df$seqnames,
                                     L_genomic_start=arms.df$L_genomic_start, L_genomic_end= arms.df$L_genomic_end,
                                     L_genomic_strand=arms.df$strand, R_seqnames=arms.df$fasta_id, R_start=arms.df$R_tx_start, R_end=arms.df$R_tx_end,
                                     R_strand="+", R_genomic_seqnames=arms.df$seqnames,
                                     R_genomic_start=arms.df$R_genomic_start,  R_genomic_end= arms.df$R_genomic_end,
                                     R_genomic_strand =arms.df$strand)

stopifnot(abs(forgi.reformated.gc.df$L_genomic_end - forgi.reformated.gc.df$L_genomic_start) == 
            abs(forgi.reformated.gc.df$L_end - forgi.reformated.gc.df$L_start))
stopifnot(abs(forgi.reformated.gc.df$R_genomic_end - forgi.reformated.gc.df$R_genomic_start) == 
            abs(forgi.reformated.gc.df$R_end - forgi.reformated.gc.df$R_start))

# Export tables
write_tsv(forgi.reformated.gc.df, paste0(prefix, "_gc.txt"), quote = FALSE)
write_tsv(forgi.reformated.tx.df, paste0(prefix, "_tx.txt"), quote = FALSE)

