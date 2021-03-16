#!/usr/bin/env Rscript

library(ggplot2)
library(ggrepel)
library(forcats)
library(scales)
library(tidyverse, warn.conflicts = FALSE)
library(stringr)
library(stringi)
library(Biostrings)
library(tictoc)
library(ggpubr)
theme_set(theme_bw() +
            theme(legend.position = "top"))
library(optparse)

library(rtracklayer)
library(GenomicFeatures)


# ==========
# Define functions
# ==========

draw_pie_chart <- function(data.df, column) {
  data_counts.df <- data.df %>%
    select(id, {{column}}) %>%
    distinct() %>%
    group_by({{column}}) %>%
    summarise(counts = n()) %>%
    arrange(desc(counts)) %>%
    mutate(percentage = percent(counts / sum(counts))) 
  
  pie <- ggplot(data_counts.df, aes(x="", y=counts, fill={{column}})) +
    geom_bar(stat = "identity", width = 0.05, color="black") +
    coord_polar("y", start=0) +
    scale_fill_brewer(palette = "Set3", direction = -4) +
    geom_label_repel(aes(label = percentage), size=3.5, show.legend = F, position = position_stack(vjust = .5)) +
    guides(fill = guide_legend(title = "Type of stems near the xl sites")) +
    ggtitle(paste0("Stem types within", d," nt from peak starts\n(", sum(data_counts.df$counts)," peaks)")) +
    theme_void() # remove background, grid, numeric labels
  return(pie)
}


get_stem_length <- function(data.df) {
  data_sh.df <- data.df %>%
    group_by(id) %>%
    dplyr::filter(element == "s0" | element == "h0") %>%  #use these positions to estimate stem length
    mutate(stem_length_L = L_start - lag(L_start),
           stem_length_R = lag(R_end) - L_end) %>%
    dplyr::filter(element == "h0") %>%
    select(id, stem_length_L, stem_length_R)
  data.df <- left_join(data.df, data_sh.df, by = "id")
  return(data.df)
}


get_hairpin_length <- function(data.df) {
  data_h.df <- data.df %>% 
    group_by(id) %>%
    dplyr::filter(element == "h0") %>%
    mutate(hairpin_length = sum(L_width)) %>%
    select(id, hairpin_length)
  data.df <- left_join(data.df, data_h.df, by = "id")
  return(data.df)
}


get_paired_regions <- function(data.df) {
  duplex_max.df <- data.df %>%
    group_by(id) %>%
    dplyr::filter(element_type == "s") %>%
    summarise(max_duplex = max(L_width),
              sum_paired = sum(L_width)) %>%
    select(id, max_duplex, sum_paired)
  data.df <- left_join(data.df, duplex_max.df)
  
}


count_element <- function(data.df, structure_type) {
  counter <- data.df %>% 
    group_by(id) %>% 
    dplyr::filter(element_type == structure_type) %>%
    summarise(count=n())
  return(counter)
}


analyse_duplexes <- function(data.df) {
  data.df <- get_paired_regions(data.df)
  data.df <- get_stem_length(data.df)
  data.df <- get_hairpin_length(data.df)
  #data.df <- get_iloop_length_sum(data.df)
  
  iloop_counts.df <- count_element(data.df, "i") %>% 
    select(id, count) %>%
    dplyr::rename(iloop_counts = count)
  data.df <- left_join(data.df, iloop_counts.df, by = "id")
  return(data.df)
}


get_id_features <- function(data.df) {
  id_features.df <- data.df %>%
    select(id, stem_type, max_duplex, sum_paired, iloop_counts, hairpin_length) %>%
    distinct() %>%
    replace(is.na(.), 0)
  return(id_features.df)
}


find_residue <- function(residue, x) {
  pos <- str_locate_all(x, residue)
  return(unique(unlist(pos)))
}


extract_paired_sequences <- function(forgi.df, rnafold.df) {
  
  duplexes.df <- forgi.df %>%
    group_by(id) %>%
    dplyr::filter(element_type == "s") %>%
    select(1:10)
  
  duplexes.df <- left_join(duplexes.df, rnafold.df, by = "id") # add the sequence and structure to the forgi.df
  
  # Extract sequences corresponding to the left and right arms of the stem, respectively
  duplexes.df <- duplexes.df %>%
    mutate(stem_L_seq = substr(sequence, L_start, L_end),
           stem_R_seq = substr(sequence, R_start, R_end))
  
  stem_L_seqs <- RNAStringSet(duplexes.df$stem_L_seq)
  stem_R_seqs <- RNAStringSet(duplexes.df$stem_R_seq)
  stem_L_freq <- as.data.frame(oligonucleotideFrequency(stem_L_seqs, width = 1, as.prob=TRUE))
  stem_R_freq <- as.data.frame(oligonucleotideFrequency(stem_R_seqs, width = 1, as.prob=TRUE))
  
  duplexes.df <- dplyr::bind_cols(duplexes.df, stem_L_freq, stem_R_freq)
  colnames(duplexes.df) = gsub("...15|...16|...17|...18", "_L", colnames(duplexes.df)) # rename columns
  colnames(duplexes.df) = gsub("...19|...20|...21|...22", "_R", colnames(duplexes.df))
  duplexes.df <- duplexes.df %>%
    dplyr::select(-sequence, -mea_structure)
  return(duplexes.df)
  
}


get_nuc_frequencies <- function(forgi.df, rnafold.df) {
  
  duplexes.df <- extract_paired_sequences(forgi.df, rnafold.df)
  # Pur content - calculate max pur L vs. R arm 
  duplexes.df <- duplexes.df %>%
    rowwise() %>% 
    mutate(L_pur = sum(c(A_L, G_L)), # purine total per stem segment
           R_pur = sum(c(A_R, G_R)))
  
  duplexes.df <- duplexes.df %>%
    group_by(id) %>%
    mutate(L_pur_mean = mean(L_pur), # purine mean per arm 
           R_pur_mean = mean(R_pur))  %>%
    select(-c(L_pur, R_pur)) %>%
    rowwise() %>%
    mutate(max_pur = max(L_pur_mean, R_pur_mean)) %>%
    ungroup() 
  
  duplexes.df$arm <- ifelse(duplexes.df$max_pur == duplexes.df$L_pur_mean, "L", "R") # assign maximal purine arm
  
  # Indivdual nucleotide frequencies
  duplexes.df <- duplexes.df %>%
    group_by(id) %>%
    mutate(L_A_mean = mean(A_L),
           L_G_mean = mean(G_L),
           R_A_mean = mean(A_R),
           R_G_mean = mean(G_R),
           L_C_mean = mean(C_L),
           L_U_mean = mean(U_L),
           R_C_mean = mean(C_R),
           R_U_mean = mean(U_R))
  
  # reorient duplexes based on maximal purine content (column "arm")
  nuc_freq_max_pur_arm.df <- duplexes.df %>% select(id, tail(names(.), 10)) %>%
    distinct() %>% dplyr::filter(arm == "L") %>%
    mutate(A_freq_max_pur_arm = L_A_mean, A_freq_min_pur_arm = R_A_mean,
           G_freq_max_pur_arm = L_G_mean, G_freq_min_pur_arm = R_G_mean,
           C_freq_max_pur_arm = L_C_mean, C_freq_min_pur_arm = R_C_mean,
           U_freq_max_pur_arm = L_U_mean, U_freq_min_pur_arm = R_U_mean) %>%
    select(-L_A_mean, -R_A_mean, -L_G_mean, -R_G_mean,
           -L_C_mean, -R_C_mean, -L_U_mean, -R_U_mean)
  nuc_freq_min_pur_arm.df <- duplexes.df %>% select(id, tail(names(.), 10)) %>%
    distinct() %>% dplyr::filter(arm == "R") %>%
    mutate(A_freq_max_pur_arm = R_A_mean, A_freq_min_pur_arm = L_A_mean,
           G_freq_max_pur_arm = R_G_mean, G_freq_min_pur_arm = L_G_mean,
           C_freq_max_pur_arm = R_C_mean, C_freq_min_pur_arm = L_C_mean,
           U_freq_max_pur_arm = R_U_mean, U_freq_min_pur_arm = L_U_mean) %>%
    select(-L_A_mean, -R_A_mean, -L_G_mean, -R_G_mean,
           -L_C_mean, -R_C_mean, -L_U_mean, -R_U_mean)
  
  nuc_freq_reordered.df <- rbind(nuc_freq_max_pur_arm.df, nuc_freq_min_pur_arm.df)
  nuc_freq_reordered.df <- nuc_freq_reordered.df %>% arrange(desc(max_pur))
  duplexes.df <- duplexes.df %>%
    select(-tail(names(.), 12))
  duplexes.df <- left_join(duplexes.df, nuc_freq_reordered.df, by = "id")
  
  forgi.df <- left_join(forgi.df, 
                        select(duplexes.df, id, element, stem_L_seq, stem_R_seq, max_pur,arm, A_freq_max_pur_arm, G_freq_max_pur_arm, 
                               A_freq_min_pur_arm, G_freq_min_pur_arm,U_freq_max_pur_arm, C_freq_max_pur_arm,
                               U_freq_min_pur_arm, C_freq_min_pur_arm),
                        by = c("id", "element"))
  
  return(forgi.df)
  
}

get_gc_content <- function(data.df) {
  gc.df <- data.df %>%
    dplyr::filter(element_type == "s") %>%
    mutate_at("stem_R_seq", stri_reverse) %>%
    rowwise() %>%
    mutate(G_pos_L = list(find_residue(stem_L_seq, residue = "G")),
           C_pos_L = list(find_residue(stem_L_seq, residue = "C")),
           G_pos_R = list(find_residue(stem_R_seq, residue = "G")),
           C_pos_R = list(find_residue(stem_R_seq, residue = "C")))
  
  #obtain the positions of Gs and Cs in L and see if the same but complement in R
  gc.df <- gc.df %>%
    rowwise() %>%
    mutate(gc_pairs = list(intersect(G_pos_L, C_pos_R)),
           cg_pairs = list(intersect(C_pos_L, G_pos_R))) %>%
    mutate(GC_pos = list(union(gc_pairs, cg_pairs))) %>%
    mutate(GC_pos_count = length(GC_pos)) %>%
    group_by(id) %>%
    mutate(sum_width = sum(L_width),
           sum_gc = sum(GC_pos_count)) %>%
    mutate(GC_percentage = sum_gc*100/sum_width) %>%
    select(id, GC_percentage) %>%
    distinct()
  data.df <- left_join(data.df, gc.df, by = "id") 
  
  return(data.df)
  
}


plot_nuc_frequencies <- function(forgi.df) {
  # Prep for plotting
  high_pur.df <- forgi.df %>% 
    select(id, A_freq_max_pur_arm, G_freq_max_pur_arm, C_freq_max_pur_arm, U_freq_max_pur_arm) %>%
    ungroup() %>%
    distinct() %>%
    dplyr::arrange(G_freq_max_pur_arm)
  
  high_pur.df <-rowid_to_column(high_pur.df)
  long_high_pur.df <- high_pur.df %>% 
    gather(residue, frequency, A_freq_max_pur_arm:U_freq_max_pur_arm)
  long_high_pur.df$arm = "High purine arm"
  
  low_pur.df <- first_stem_loop.df %>% 
    select(id, A_freq_min_pur_arm, G_freq_min_pur_arm, C_freq_min_pur_arm, U_freq_min_pur_arm) %>%
    ungroup() %>%
    distinct() %>%
    dplyr::arrange(desc(U_freq_min_pur_arm))
  low_pur.df <- rowid_to_column(low_pur.df)
  long_low_pur.df <- low_pur.df %>% 
    gather(residue, frequency, A_freq_min_pur_arm:U_freq_min_pur_arm)
  long_low_pur.df$arm = "Low purine arm"
  
  nuc_freq.df <- rbind(long_high_pur.df, long_low_pur.df)
  nuc_freq.df$residue <- str_sub(nuc_freq.df$residue, end=-18) #trim string to get nucleotide letter
  
  # Plot individual nucleotide frequencies
  nuc_freq.gg <- ggplot(nuc_freq.df, aes(x=rowid, y=frequency, color = residue)) +
    geom_smooth(se=F) + coord_flip() +
    facet_grid(cols = vars(arm)) +
    ylab("Frequency") +
    xlab("ID") +
    scale_color_manual(values = c("darkgreen","dodgerblue4","goldenrod3","firebrick"))
  
  return(nuc_freq.gg)
  
}


run_rnaeval <- function(fasta.file, args =  c("-i", paste0(fasta.file,">temp.fa"))) {
  
  rnaeval.out <- system2(command = "RNAeval", args = args, stdout = TRUE)
  rnaeval.df <- as.data.frame(readBStringSet("temp.fa"))
  rnaeval.df <- rnaeval.df %>% rownames_to_column(var = "id")
  
  rnaeval.df <- rnaeval.df %>%
    rowwise() %>%
    mutate(eval_mfe = as.numeric(str_sub(x, -7,-2))) %>%
    ungroup() %>%
    select(id, eval_mfe)
  
  file.remove("temp.fa")
  return(rnaeval.df)
  
}


get_mfe <- function(forgi, rnafold, prefix) {
  
  data.df <- forgi %>%
    dplyr::filter(element_type == "s") %>%
    group_by(id) %>%
    mutate(L_limit = min(L_start),
           R_limit = max(R_end)) %>%
    select(id, L_limit, R_limit) %>%
    distinct() %>%
    ungroup()
  
  data.df <- left_join(data.df , rnafold, by = "id") # add the sequence and structure to the forgi
  data.df <- data.df %>%
    mutate(stem_loop_seq = substr(sequence, L_limit, R_limit),
           stem_loop_db = substr(mea_structure, L_limit, R_limit)) %>%
    select(-sequence, -mea_structure)
  # write multi-fasta with seq and db
  fasta <- character(nrow(data.df) * 3) # empty file
  fasta[c(TRUE, FALSE, FALSE)] <- paste0(">", data.df$id)
  fasta[c(FALSE, TRUE, FALSE)] <- data.df$stem_loop_seq
  fasta[c(FALSE, FALSE, TRUE)] <- data.df$stem_loop_db
  
  writeLines(fasta, paste0(prefix,"_rnaeval.fasta"))
  mfe.df <- run_rnaeval(paste0(prefix,"_rnaeval.fasta"))
  data.df <- left_join(data.df, mfe.df, by = "id")
  forgi <- left_join(forgi, data.df, by = "id")
  return(forgi)
  
}

shuffle_sequence <- function(sequence, number = 1, klet = 2, seed = 42) {
  
  system(paste0("ushuffle -seed ", seed, " -k ", klet, " -n ", number, " -s ", sequence), intern = TRUE)
}



resize_peaks <- function(bedfiles.list, left = 100, right = 100) {
  
  w <- left + right +1  # width of interval: xl site + flanks
  
  grl <- GRangesList(lapply(bedfiles.list, import.bed))
  gr <- unlist(grl)
  gr <- keepStandardChromosomes(gr, pruning.mode = "coarse")
  gr <- dropSeqlevels(gr, c("chrM", "chrY"), pruning.mode = "coarse")
  # gr$peak_id <- seq(1,length(gr))
  gr.df <- as.data.frame(gr)
  
  gr <- resize(gr, width = 1, fix = "start") # resize the peaks, start of peak = 1
  gr <- unique(gr)  # keep unique xl positions
  gr <- resize(resize(gr, width = right+1, fix = "start"), width = w, fix = "end") # add +/- flanks
  gr$id <- paste0("ID", seq(1,length(gr)))
  
  resized.gr.df <- as.data.frame(gr)
  peaks.df <- left_join(gr.df,resized.gr.df, by="peak_id", suffix = c("", ".unique_resized"))
  
  # write.table(peaks.df, paste0(prefix, "_gr.df.txt"), quote = FALSE, sep = "\t") # save this to be able to map processed data to the original peaks later on
  return(gr)
  
}


get_duplex_arms <- function(forgi, rnafold) {
  
  data.df <- forgi %>%
    dplyr::filter(element_type == "s") %>%
    group_by(id) %>%
    mutate(L_min = min(L_start), L_max = max(L_end),
           R_min = min(R_start), R_max = max(R_end)) %>%
    select(id, L_min, L_max, R_min, R_max) %>%
    distinct() %>%
    ungroup()
  
  data.df <- left_join(data.df , rnafold, by = "id") # add the sequence and structure
  
  data.df <- data.df %>%
    mutate(L_seq = substr(sequence, L_min, L_max), R_seq = substr(sequence, R_min, R_max),
           L_db = substr(mea_structure, L_min, L_max), R_db = substr(mea_structure, R_min, R_max)) %>%
    select(-sequence, -mea_structure)
  
  forgi <- left_join(forgi, data.df, by = "id")
  return(forgi)
}


# gr needs to have a "name" metadata column with transcript names
convert_to_transcriptomic_coordinates <- function(gr, txdb.sqlite) {
  txdb <- loadDb(txdb.sqlite)
  transcripts.gr <- transcripts(txdb)
  names(transcripts.gr) <- id2name(txdb, feature.type = "tx")
  transcripts.gr <- transcripts.gr[transcripts.gr$tx_name %in% gr$name]
  mapped.gr <- mapToTranscripts(x=gr, transcripts=transcripts.gr)
  return(mapped.gr)
}

merge_by_partial_string <- function(df1, df2) {
  df1$matchID = row.names(df1)
  df2$matchID = sapply(df2$gene_name, function(x) grep(x, df1$fasta_id)) # not very fast but ok
  
  df_merge = merge(df1, df2, by = "matchID")[-1]
  return(df_merge)
}


# ==========
# Define options
# ==========

option_list <- list(
  make_option(c("-f", "--forgi"), type="character", action="store",
              help="Provide an input forgi file"), 
  make_option(c("-r", "--rnafold"), type="character", action="store",
              help="Provide an input rnafold file"),
  make_option(c("-d", "--distance"), action="store", type="numeric", default=15,
              help="Distance from peak start position")
)

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)



# ==========
# Start the analysis
# ==========

if ( !is.na(opt$forgi )) {
  forgi.file <- opt$forgi
  forgi.df <- read.csv(opt$forgi, sep="\t")
} else {
  cat("Please provide the forgi input file\n", file = stderr())
}

if ( !is.na(opt$rnafold )) {
  rnafold.df <- read.csv(opt$rnafold, sep="\t")
  rnafold.df <- rowid_to_column(rnafold.df, "id")
  rnafold.df$id <- paste0("ID", rnafold.df$id, sep="")
  rnafold.df <- dplyr::select(rnafold.df, c(id, sequence, mea_structure)) # keep only the mea information
} else {
  cat("Please provide the rnafold input file\n", file = stderr())
}

if ( !is.na(opt$d )) {
  d <- opt$d
} else {
  d <- 15
  cat("No distance provided, the default of 15 nt will be used\n", file = stderr())
}

# ==========
# Load data
# ==========

forgi.df <- read.csv("/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1.forgi.tsv.gz", sep = "\t")
rnafold.df <- read.csv("/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1.rnafold.tsv.gz", sep = "\t")
rnafold.df <- rowid_to_column(rnafold.df, "id")
rnafold.df$id <- paste0("ID", rnafold.df$id, sep="")
rnafold.df <- dplyr::select(rnafold.df, c(id, sequence, mea_structure))
d = 15


txdb <- paste0(data.dir,"/gencode_V33_txdb.sqlite")

gencode.txdb <- loadDb(txdb)

human.gtf <- paste0(data.dir, "/human_GencodeV33.gtf")
data.dir <- "/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt"
rnaplfold.clusters.df <- read.csv(paste0(data.dir,"/stau1_threeutrs.rnaplfold_prob_clusters.df.txt"), sep = "\t")
xl <- 1

prefix <- str_c(str_split(forgi.file, "\\.")[[1]][1:3], collapse = ".")

prefix = "stau1.peaks.10nt_10nt_"

message(paste0("Analysing ", prefix, " data"))

tic()

# ==========
# Filter the forgi & rnafold files on IDs from clusters 1 and 3 from the cluster_probability_profiles.R output
# ==========

rnaplfold.clusters.df <- rnaplfold.clusters.df %>% dplyr::filter(.cluster == 1 | .cluster == 3 | .cluster == 5)

forgi.df <- semi_join(forgi.df, rnaplfold.clusters.df, by = "id")
rnafold.df <- semi_join(rnafold.df, rnaplfold.clusters.df, by = "id")

# ==========
# Process the forgi file
# ==========

# Classify structures based on their distance from xl site, whether they span multi-loops 
# and whether the first stem-loop is fused to another stem-loop

# Find the minimum location of the first non-null multiloop
first_non_null_m.df <- forgi.df %>%
  group_by(id) %>%
  dplyr::filter(element_type == "m") %>%
  dplyr::slice(which.min(L_start)) %>%
  select(id, L_start) %>%
  dplyr::rename(first_m = L_start)

forgi.df <- left_join(forgi.df, first_non_null_m.df, by="id")

# Classify structure types based on the position of first stem start relative to xl site and presence and position of multiloops
forgi.df <- forgi.df %>%
  group_by(id) %>%
  mutate(duplex_type = case_when(any(element == "f0" & L_end >= d) ~ "No stem",
                                 (!any(element == "f0" & L_end >= d)) & (any(element == "s0" & L_start == xl)) ~ "Stem and xl paired",
                                 (!any(element == "f0" & L_end >= d)) & (any(element == "s0" & L_start > xl & L_start <= d & is.na(first_m))) ~ "Stem and xl unpaired\n(no multiloops, first or all multiloops are null)",
                                 (!any(element == "f0" & L_end >= d)) & (any(element == "s0" & L_start > xl & L_start <= d & R_end > first_m)) ~ "Stem and xl unpaired\n(multiloops)",
                                 (!any(element == "f0" & L_end >= d)) & (any(element == "s0" & L_start > xl & L_start <= d & R_end < first_m)) &
                                   ((is.na(first_m) | any(element == "m0" & is.na(L_start)))) ~ "Stem and xl unpaired\n(no multiloops, first or all multiloops are null)",
                                 (!any(element == "f0" & L_end >= d)) & (any(element == "s0" & L_start > xl & L_start <= d & R_end < first_m)) &
                                   !((is.na(first_m) | any(element == "m0" & is.na(L_start)))) ~ "Stem and xl unpaired\n(no multiloops, single stems)"))

pie1.gg <- draw_pie_chart(forgi.df, duplex_type)
ggsave(paste0(prefix, "_pie1.pdf"), pie1.gg, height = 11, width = 7, dpi = 300)

# Analyse the single stems closest to the xl site

# select the first stem
first_stem_loop.df <- forgi.df %>%
  group_by(id) %>%
  dplyr::filter(duplex_type == "Stem and xl unpaired\n(no multiloops, single stems)") %>%
  dplyr::filter((R_end < first_m) | (element_type == "i" & is.na(R_end) & L_end < first_m)
                | (element_type == "h" & is.na(R_end) & L_end < first_m)) %>%
  select(-c(first_m, duplex_type))


# assign the types of first stems based on the presence of bulges and internal loops
first_stem_loop.df <- first_stem_loop.df %>% replace(is.na(.), 0) # replqce NA with zero to simplify filtering

first_stem_loop.df <- first_stem_loop.df %>%
  group_by(id) %>%
  mutate(stem_type = case_when(!any(element_type == "i") ~ "Uninterrupted stem",
                               any(element_type == "i") & (!any(element_type == "i" & (L_width > 1 | R_width > 1))) &
                                 (!any(element_type == "i" & L_width != R_width)) ~ "Stem with only symmetrical bulges",
                               (any(element_type == "i")) & (!any(element_type == "i" & (L_width > 1 | R_width > 1))) &
                                 (any(element_type == "i" & L_width != R_width)) ~ "Stem with non-symmetrical bulges",
                               (any(element_type == "i")) &
                                 (any(element_type == "i" & (L_width > 1) & (R_width == 0 | R_width >= 1))) |
                                 (any(element_type == "i" & (L_width = 1) & (R_width > 1))) ~ "Stem with internal loops"))


pie2.gg <- draw_pie_chart(first_stem_loop.df, stem_type)
ggsave(paste0(prefix, "_pie2.pdf"), pie2.gg, height = 11, width = 7, dpi = 300)

# First stem loop properties
first_stem_loop.df <- analyse_duplexes(first_stem_loop.df)

# Nucleotide frequency

first_stem_loop.df <- get_nuc_frequencies(first_stem_loop.df, rnafold.df)
# GC content
first_stem_loop.df <- get_gc_content(first_stem_loop.df)

# get mea of the first stem-loop
first_stem_loop.df <- get_mfe(first_stem_loop.df, rnafold.df, prefix)

# write_tsv(first_stem_loop.df, paste0(prefix, "_forgi_analyses.df.txt"), quote = F)

# ==========
# Plotting
# ==========

# Plot features

id_features.df <- get_id_features(first_stem_loop.df)

id_features.df$stem_type <- as.factor(id_features.df$stem_type)
id_features.df <- rowid_to_column(id_features.df)
id_features.df <- id_features.df %>% 
  gather("feature", "counts", max_duplex:hairpin_length)


id_features.gg <- ggplot(transform(id_features.df, feature = factor(feature, levels = c("sum_paired", "max_duplex","hairpin_length","iloop_counts"))),
                         aes(x = "", y=counts, fill=stem_type)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(. ~ feature, scales = "free") +
  theme(legend.position="right", legend.direction="vertical")

id_features_summed.gg <- ggplot(id_features.df, aes(x = fct_reorder(feature, counts, .desc=TRUE), y = counts, fill = feature)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_brewer(palette="Blues") +
  theme(legend.position="right", legend.direction="vertical")

ggsave(paste0(prefix, "_boxplot_each.pdf"), id_features.gg, height = 7, width = 7, dpi = 300)
ggsave(paste0(prefix, "_boxplot_all.pdf"), id_features_summed.gg, height = 7, width = 7, dpi = 300)

# Nuc frequency plots
nuc_freq.gg <- plot_nuc_frequencies(first_stem_loop.df)
ggsave(paste0(prefix, "_nuc_freq.pdf"), nuc_freq.gg, height = 11, width = 7, dpi = 300)

# MFE of the first stem loop

mfe.df <- first_stem_loop.df %>% dplyr::select(id, eval_mfe) %>% distinct()

# assign positive values to 0
mfe.df$Experiment = "Non-hybrids"
mfe.df$eval_mfe[mfe.df$eval_mfe > 0] <- 0

mfe_dens.gg <- ggplot(mfe.df, aes(x = eval_mfe, color = Experiment)) + 
  geom_density(alpha = 0.8) +
  scale_colour_brewer(palette="Spectral") +
  facet_grid(cols = vars(Experiment)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Density") +
  xlab("MFE (kcal/mol)") +
  xlim(-87.5, 0) +
  ggtitle("Non-hybrids")
mfe_dens.gg
ggsave(paste0(prefix, "_mfe_dens.pdf"), mfe_dens.gg, height = 4, width = 4, dpi = 300)

toc()
message("Analysis completed!")



# ==========
# Code below is in testing phase
# ==========

# ==========
# Filter for stem-loop strucures that have similar properties to the ones originating from hybrid reads
# ==========

### Many structures have high mfe; attempt filtering by the length of the stem + paired segment?

# quantile(unique(first_stem_loop.df$stem_length_L), 0.95)

# filter for structures at least x nt long and max y nt long, and also based on number of paired residues
test <- first_stem_loop.df %>%
  dplyr::filter(between(stem_length_L, 8, 44) | between(stem_length_R, 8, 44)) %>%
  dplyr::filter(sum_paired > 5 & max_duplex > 3 & max_duplex < 16)

mfe.df <- test %>% dplyr::select(id, eval_mfe) %>% distinct()

mfe.df$Experiment = "Non-hybrids"
mfe.df$eval_mfe[mfe.df$eval_mfe > 0] <- 0

mfe_dens.gg <- ggplot(mfe.df, aes(x = eval_mfe, color = Experiment)) + 
  geom_density(alpha = 0.8) +
  scale_colour_brewer(palette="Spectral") +
  facet_grid(cols = vars(Experiment)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Density") +
  xlab("MFE (kcal/mol)") +
  xlim(-87.5, 0) +
  ggtitle("Non-hybrids")


# ==========
# Convert forgi stem-loops data to a format that is similar to the hybrids table from Tosca
# ==========

# load peaks data, resize exactly as for rnafold annotationsu sed to predict structures

peaks.gc.gr <- resize_peaks(bed.files, left = 0, right = 100)

# add a new column containing gene ids by mapping the trancript ids
keys = peaks.gc.gr$name
mapped.ids <- select(gencode.txdb, keys = keys, columns="GENEID", keytype="CDSNAME")
peaks.gc.gr$gene_name <- mapped.ids$GENEID
genes.ls <- unique(peaks.gc.gr$gene_name)

# add fasta_id from annotation used for Tosca as a new column
human.gr <- import.gff2(human.gtf)
human.gr <- keepStandardChromosomes(human.gr, pruning.mode = "coarse")
human.gr <- dropSeqlevels(human.gr, c("chrM", "chrY"), pruning.mode = "coarse")
human.df <- as.data.frame(human.gr)

human.df <- human.df %>%
  filter(str_detect(fasta_id, paste(genes.ls, collapse = "|"))) #filter for gene names in the peak gr
human.df$row <- row.names(human.df)
human.df <- as.data.frame(human.df[,c("row", "fasta_id")])
peaks.gc.df <- as.data.frame(peaks.gc.gr)

peaks.gc.df <- merge_by_partial_string(human.df, peaks.gc.df)
peaks.gc.df <- peaks.gc.df %>%
  dplyr::select(-row)
  
# Get transcriptomic coordinates; Convert genomic coordinates to transcriptomic coordinates
peaks.tx.gr <- convert_to_transcriptomic_coordinates(peaks.gc.gr, txdb)
peaks.tx.gr$id <- paste0("ID",peaks.tx.gr$xHits)
peaks.tx.gr$name <- seqnames(peaks.tx.gr)

# make df containing the peaks starts + 100 in genomic and transcriptomic coordinates
peaks.df <- left_join(peaks.gc.df, as.data.frame(peaks.tx.gr), by = c("id", "name"), suffix=c("",".tx"))

# get the start and end position of the stem-loop individually for each arm
first_stem_loop.df <- get_duplex_arms(first_stem_loop.df, rnafold.df)

arms.df <- first_stem_loop.df %>% 
  distinct(id, L_min, L_max, R_min, R_max, L_seq, R_seq, L_db, R_db)

# merge GRanges info to forgi data frame
arms.df <- as.data.frame(left_join(arms.df, peaks.df, by = "id"))

# Calculate genomic coordinates, remember that xl = 1
arms.df <- arms.df %>%
  mutate(L_genomic_start = case_when((strand == "+") ~ start + L_min - xl,
                               (strand == "-") ~ end - L_max + xl),
         L_genomic_end = case_when((strand == "+") ~ start + L_max - xl,
                                     (strand == "-") ~ end - L_min + xl),
         R_genomic_start = case_when((strand == "+") ~ start + R_min - xl,
                               (strand == "-") ~ end - R_max + xl),
         R_genomic_end = case_when((strand == "+") ~ start + R_max - xl,
                                     (strand == "-") ~ end - R_min + xl))

# Calculate transcriptomic coordinates
arms.df <- arms.df %>%
  mutate(L_tx_start = start.tx + L_min - xl,
         L_tx_end = start.tx + L_max - xl,
         R_tx_start = start.tx + R_min - xl,
         R_tx_end = start.tx + R_max - xl)


# transcriptomic:
forgi.reformated.tx.df <- data.frame(id=arms.df$id, name=arms.df$name, orientation=NA, type = "intragenic", hybrid_selection = NA, L_seqnames = arms.df$fasta_id,
                L_read_start = arms.df$L_min, L_read_end = arms.df$L_max, L_start = arms.df$L_tx_start, L_end = arms.df$L_tx_end, L_width = NA, L_strand = "+", 
                R_seqnames = arms.df$fasta_id, R_read_start = arms.df$R_min, R_read_end = arms.df$R_max, R_start = arms.df$R_tx_start, R_end = arms.df$R_tx_end, R_width = NA, R_strand = "+",
                umi=NA, L_sequence=arms.df$L_seq, R_sequence=arms.df$R_seq, mfe = NA)

forgi.reformated.tx.df <- forgi.reformated.tx.df %>%
  rowwise() %>%
  mutate(L_width = L_end - L_start + 1,
         R_width = R_end - R_start +1,
         L_sequence = str_replace_all(L_sequence, "U", "T"), # replace Us with Ts
         R_sequence = str_replace_all(R_sequence, "U", "T")) %>%
  ungroup()


# genomic:
forgi.reformated.gc.df <- data.frame(name = arms.df$id, count = arms.df$score, L_seqnames=arms.df$fasta_id,	
                               L_start=arms.df$L_tx_start, L_end=arms.df$L_tx_end, L_strand="+", L_genomic_seqnames=arms.df$seqnames,
                               L_genomic_start=arms.df$L_genomic_start, L_genomic_end= arms.df$L_genomic_end,
                               L_genomic_strand=arms.df$strand, R_seqnames=arms.df$fasta_id, R_start=arms.df$R_tx_start, R_end=arms.df$L_tx_end,
                               R_strand="+", R_genomic_seqnames=arms.df$seqnames,
                               R_genomic_start=arms.df$R_genomic_start,  R_genomic_end= arms.df$R_genomic_end,
                               R_genomic_strand =arms.df$strand)

head(forgi.reformated.gc.df)
head(forgi.reformated.tx.df)

write.table(forgi.reformated.gc.df, paste0(prefix, "_gc.txt"), quote = FALSE, sep = "\t")
write.table(forgi.reformated.tx.df, paste0(prefix, "_tx.txt"), quote = FALSE, sep = "\t")

