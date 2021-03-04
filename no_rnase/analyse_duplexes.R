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


#' A function to count all occurrences of a structure annotation type residue within a structure/ID
#'
#' @description 
#' This function will count the structure element "structure_type" within a structure annotation dataframe "data.df"
#' 
#' @param data.df the dataframe you would like to search and count from
#' @param structure_type the string that denotes the structure element type
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
    ggtitle(paste0("Stem types within", d," nt from crosslink sites\n(", sum(data_counts.df$counts)," cross-link sites)")) +
    theme_void() # remove background, grid, numeric labels
  return(pie)
}

#' A function to calculate the length of a stem structure, including the bulges and internal loops, but not hairpins
#'
#' @description 
#' This function will calculate stem-loop lengths per ID on the left and right sides, respectively, from a structure annotation dataframe 
#' "data.df" and add them as new columns
#' 
#' @param data.df the dataframe you would like to analyse
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


# get_iloop_length_sum <- function(data.df) {
#   data_h.df <- data.df %>% 
#     group_by(id) %>%
#     dplyr::filter(element_type == "i") %>%
#     mutate(iloop_length_L = sum(L_width),
#            iloop_length_R = sum(R_width)) %>%
#     select(id, iloop_length_L, iloop_length_R)
#   data.df <- left_join(data.df, data_h.df, by = "id")
#   return(data.df)
# }


#' A function to calculate the longest duplex per ID
#'
#' @description 
#' This function will calculate the maximum continuous duplex stretch in each structure/ID
#' from a structure annotation dataframe "data.df" and add it as a new column
#' 
#' @param data.df the dataframe you would like to analyse
get_paired_regions <- function(data.df) {
  duplex_max.df <- data.df %>%
    group_by(id) %>%
    dplyr::filter(element_type == "s") %>%
    summarise(max_duplex = max(L_width),
              sum_paired = sum(L_width)) %>%
    select(id, max_duplex, sum_paired)
  data.df <- left_join(data.df, duplex_max.df)
  
}


#' A function to count all occurences of a structure annotation type residue within a structure/ID
#'
#' @description 
#' This function will count the structure element "structure_type" within a structure annotation dataframe "data.df"
#' 
#' @param data.df the dataframe you would like to search and count from
#' @param structure_type the string that denotes the structure element type
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

#' A function to find all locations of a specific residue within a sequence
#'
#' @description 
#' This function will find and locate the letter "residue" within a string sequence "x"
#' 
#' @param residue one letter you'd like to find
#' @param x the string you would like to find it in
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


get_cofold_mfe <- function(forgi, rnafold, prefix) {
  
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
  
  # fuse the sequences using "&" separator
  data.df <- data.df %>%
    unite("duplex_seq", L_seq:R_seq, remove = FALSE, sep = "&") %>%
    unite("duplex_db", L_db:R_db, remove = FALSE, sep = "&")
  
  # write multi-fasta with seq and db
  fasta <- character(nrow(data.df) * 3) # empty file
  fasta[c(TRUE, FALSE, FALSE)] <- paste0(">", data.df$id)
  fasta[c(FALSE, TRUE, FALSE)] <- data.df$duplex_seq
  fasta[c(FALSE, FALSE, TRUE)] <- data.df$duplex_db
  
  writeLines(fasta, paste0(prefix,"_duplex_rnaeval.fasta"))
  mfe.df <- run_rnaeval(paste0(prefix,"_duplex_rnaeval.fasta"))
  data.df <- left_join(data.df, mfe.df, by = "id")
  forgi <- left_join(forgi, data.df, by = "id")
  return(forgi)
  
}


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


forgi.df <- read.csv("/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1.10nt_10nt.peaks.forgi.tsv.gz", sep = "\t")
# rnafold.df <- read.csv("/Users/iosubi/Documents/projects/computational_hiCLIP/nonhybrids_10nt_10nt/stau1.10nt_10nt.peaks.rnafold.tsv.gz", sep = "\t")
# rnafold.df <- rowid_to_column(rnafold.df, "id")
# rnafold.df$id <- paste0("ID", rnafold.df$id, sep="")
# rnafold.df <- dplyr::select(rnafold.df, c(id, sequence, mea_structure))


xl <- 1
#d = 15
prefix <- str_c(str_split(forgi.file, "\\.")[[1]][1:3], collapse = ".")

#prefix = "test"

message(paste0("Analysing ", prefix, " data"))

tic()

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
#draw_stacked_barchart(first_stem_loop.df, stem_type)

# First stem loop properties
first_stem_loop.df <- analyse_duplexes(first_stem_loop.df)

# Nucleotide frequency

first_stem_loop.df <- get_nuc_frequencies(first_stem_loop.df, rnafold.df)
# GC content
first_stem_loop.df <- get_gc_content(first_stem_loop.df)

# get mea of the first stem-loop
#first_stem_loop.df <- get_mfe(first_stem_loop.df, rnafold.df, prefix)


nrow(first_stem_loop.df)
# get binding energy between the 2 arms of the first stem-loop
# MFE (arms of the first stem loop)
first_stem_loop.df <- get_cofold_mfe(first_stem_loop.df, rnafold.df, prefix)

write_tsv(first_stem_loop.df, "forgi_analyses.df.txt", quote = F)

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

# mfe.gg <- ggplot(mfe.df, aes(eval_mfe)) + stat_ecdf(geom = "step") +
#   labs(title=paste0(prefix, "MFE"),
#        y = "F(MFE)", x="MFE (kcal/mol")
# ggsave(paste0(prefix, "_mfe.pdf"), nuc_freq.gg, height = 11, width = 7, dpi = 300)

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


### Many structures have high mfe; attempt filtering by the length of the stem + paired segment
colnames(first_stem_loop.df)

median(unique(first_stem_loop.df$stem_length_L))

quantile(unique(first_stem_loop.df$stem_length_L), 0.95)


# filter for structures at least 6 nt long and max 36 nt long
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

mfe_dens.gg

median(unique(test$stem_length_L))
median(unique(test$eval_mfe))
mean(unique(test$eval_mfe))
median(first_stem_loop.df$stem_length_R)

