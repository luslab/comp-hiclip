library(GenomicAlignments)
library(data.table)
library(tidyverse)
library(stringr)

bam_to_dataframe <- function(bam.file) {
  bam <- readGAlignments(bam.file,
                         use.names = TRUE)
  bam <- bam[njunc(bam) == 0] # filter out sj
  bam.df <- as.data.frame(bam)
  bam.df$seqnames <- as.character(bam.df$seqnames)
  bam.df$qname <- row.names(bam.df)
  bam.df <- bam.df %>% 
    dplyr::arrange(seqnames) %>%
    dplyr::filter((seqnames != "chrM") & !str_detect(seqnames,"KI"))
  # unique(bam.df$seqnames)
  
  bam.df$read <- gsub("^.*?\\.", "", bam.df$qname)
  bam.df$arm <- gsub("\\..*$", "", bam.df$qname)
  left.df <- bam.df %>%
    dplyr::filter(arm == "L")
  right.df <- bam.df %>%
    dplyr::filter(arm == "R")
  hybrids.df <- inner_join(left.df, right.df, by = "read")
  hybrids.df <- hybrids.df %>% arrange(read)
  print(paste0("There are ", nrow(hybrids.df), " valid hybrid reads (with both arms mapped)."))
  
  # place L arm before the R arm
  hybrids.df <- hybrids.df %>%
    group_by(read) %>%
    mutate(L_start = case_when((start.x > start.y) ~ start.y, TRUE ~ start.x),
           R_start = case_when((start.x > start.y) ~ start.x, TRUE ~ start.y))
  
  # re-assign the positions of the other columns if swapped in previous step
  hybrids.df <- hybrids.df %>%
    mutate(L_seqnames = case_when((L_start == start.x) ~ seqnames.x, TRUE ~ seqnames.y),
           L_strand = case_when((L_start == start.x) ~ as.character(strand.x), TRUE ~ as.character(strand.y)),
           L_width = case_when((L_start == start.x) ~ width.x, TRUE ~ width.y),
           L_qname = case_when((L_start == start.x) ~ qname.x, TRUE ~ qname.y),
           R_seqnames = case_when((L_start == start.x) ~ seqnames.y, TRUE ~ seqnames.x),
           R_strand = case_when((L_start == start.x) ~ as.character(strand.y), TRUE ~ as.character(strand.x)),
           R_width = case_when((L_start == start.x) ~ width.y, TRUE ~ width.x),
           R_qname = case_when((L_start == start.x) ~ qname.y, TRUE ~ qname.x))
  hybrids.df <- hybrids.df %>% 
    select(read, L_seqnames, L_start, L_strand, L_width, L_qname, R_seqnames, R_start, R_strand, R_width, R_qname)  
  #write.csv(hybrids.df, paste0(bam.file, "_hybrids.df.txt"), sep = "\t" )
  return(data.frame(hybrids.df))
}

reorient_hybrids <- function(hybrids.dt) {
  
  # First do starts
  correct.dt <- hybrids.dt[L_start <= R_start]
  incorrect.dt <- hybrids.dt[L_start > R_start]
  
  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)
  
  setnames(incorrect.dt, renamed)
  
  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  
  stopifnot(all(reoriented.dt$L_start <= reoriented.dt$R_start))
  
  # Then do subject (to make sure intergenics in same order)
  correct.dt <- reoriented.dt[L_seqnames <= R_seqnames]
  incorrect.dt <- reoriented.dt[L_seqnames > R_seqnames]
  
  renamed <- gsub("^L_", "X_", names(incorrect.dt))
  renamed <- gsub("^R_", "L_", renamed)
  renamed <- gsub("^X_", "R_", renamed)
  
  setnames(incorrect.dt, renamed)
  
  reoriented.dt <- rbindlist(list(correct.dt, incorrect.dt), use.names = TRUE)
  stopifnot(all(reoriented.dt$L_subject <= reoriented.dt$R_subject))
  stopifnot(nrow(reoriented.dt) == nrow(hybrids.dt))
  
  return(reoriented.dt)
  
}

work.dir <- "/camp/lab/luscomben/home/users/iosubi/projects/comp_hiclip/linker"
# bam.ls <- list.files(path = work.dir, pattern = ".linker.Aligned.sortedByCoord.out.bam")
# bam.df.ls <- lapply(bam.ls, bam_to_dataframe)

high_hybrids.df <- bam_to_dataframe(paste0(work.dir,"/LigPlusHigh.linker.Aligned.sortedByCoord.out.bam"))
high_hybrids.dt <- data.table(high_hybrids.df)
high_hybrids_reoriented.dt <- reorient_hybrids(high_hybrids.dt)

low_hybrids.df <- bam_to_dataframe(paste0(work.dir,"/LigPlusLow.linker.Aligned.sortedByCoord.out.bam"))
low_hybrids.dt <- data.table(low_hybrids.df)
low_hybrids_reoriented.dt <- reorient_hybrids(low_hybrids.dt)


