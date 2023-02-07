

results.path <- "/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/revisions"

tic()


# Reference
TxDb <- loadDb("/Volumes/lab-luscomben/home/shared/projects/ira-nobby/comp_hiclip/ref/gencode.v33.txdb.sqlite")
txlengths.dt <- data.table(transcriptLengths(TxDb, with.cds_len = TRUE,
                                             with.utr5_len = TRUE,
                                             with.utr3_len = TRUE))
setorder(txlengths.dt, gene_id, -utr3_len)
utr3.grl <- threeUTRsByTranscript(TxDb, use.names = TRUE)

# Atlases
duplexes.dt <- fread(file.path(results.path, "atlas/merged.atlas.clusters.tsv.gz"))
duplexes.dt <- duplexes.dt[!grep("^ID", cluster)]

paris.dt <- fread(file.path(results.path, "paris/results_paris/atlas_clusters/all.atlas_clusters.gc.annotated.tsv.gz"))
ricseq.dt <- fread(file.path(results.path, "ricseq/results_ricseq/atlas_clusters/all.atlas_clusters.gc.annotated.tsv.gz"))

utr3.duplexes.dt <- ricseq.dt[L_seqnames == R_seqnames][L_region == "UTR3"][R_region == "UTR3"]
# utr3.duplexes.dt <- paris.dt[L_seqnames == R_seqnames][L_region == "UTR3"][R_region == "UTR3"]
# utr3.duplexes.dt <- duplexes.dt[L_seqnames == R_seqnames][L_region == "UTR3"][R_region == "UTR3"]
utr3.duplexes.dt[, ensg := sapply(strsplit(L_gene_id, "\\."), "[[", 1)]

L.gr <- convert_to_granges(utr3.duplexes.dt, arm = "L", genomic = TRUE)
R.gr <- convert_to_granges(utr3.duplexes.dt, arm = "R", genomic = TRUE)
stopifnot(all(countOverlaps(L.gr, utr3.grl) != 0))
stopifnot(all(countOverlaps(R.gr, utr3.grl) != 0))
duplexes.grl <- split(c(L.gr, R.gr), c(L.gr, R.gr)$name)

utr3.dt <- mclapply(1:length(duplexes.grl), function(i) {
  
  # utr3.dt <- lapply(1:length(duplexes.grl), function(i) {
  
  # message(i)
  x <- duplexes.grl[[i]]
  
  # First get 3' UTR for the transcript
  ol <- findOverlaps(x, utr3.grl)
  stopifnot(length(ol) > 0)
  
  if(length(unique(subjectHits(ol))) == 1) {
    
    utr3 <- utr3.grl[unique(subjectHits(ol))]
    
    utr3 <- unlist(utr3)
    utr3$tx_name <- names(utr3)
    names(utr3) <- NULL
    
  } else {
    
    sel.ol <- unique(subjectHits(ol))
    sel.ol <- sel.ol[sapply(sel.ol, function(y) all(queryHits(ol[subjectHits(ol) == y]) == c(1, 2)))] # Keep only ones where both arms overlapped by same 3' UTR overlapped
    # stopifnot(length(sel.ol) != 0)
    if(length(sel.ol) == 0) {
      return(NA) # these are intra-genic UTR3-UTR3 but across different transcripts
    }
    
    utr3 <- utr3.grl[sel.ol]
    
    utr3.tx.dt <- txlengths.dt[tx_name %in% names(utr3)] # select matching tx
    utr3.tx.dt[gene_id %in% x$L_gene_id] # keep only same genes
    utr3 <- utr3[names(utr3) %in% utr3.tx.dt[1]] # select longest 3' UTR 
    
    utr3 <- unlist(utr3)
    utr3$tx_name <- names(utr3)
    names(utr3) <- NULL
    
  }
  
  # Then work out introns
  utr3.dt <- data.table(name = unique(x$name),
                        ensg = unique(x$ensg),
                        utr3_start = min(start(utr3)),
                        utr3_end = max(end(utr3)),
                        utr3_width = sum(width(utr3)))
  
  g <- gaps(utr3, start = min(start(utr3)), end = max(end(utr3)))
  g <- g[seqnames(g) == unique(seqnames(utr3))]
  g <- g[strand(g) == unique(strand(utr3))]
  
  # Reversed order if strand is neg to match utr3
  if(unique(strand(utr3)) == "-") {
    g <- sort(g, decreasing = TRUE)
  } else {
    g <- sort(g)
  }
  
  utr3.dt[, utr3_intron := sum(width(g))]
  
  # Between arms   
  L.ol <- max(which(countOverlaps(utr3, x[1]) > 0)) # if intron in arm, select
  R.ol <- min(which(countOverlaps(utr3, x[2]) > 0))
  
  if(L.ol == R.ol) {
    
    utr3.dt[, intron := 0]
    
  } else {
    
    stopifnot(R.ol > L.ol) # Should be this for all
    k <- seq(L.ol, R.ol - 1, by = 1) # which introns to keep
    utr3.dt[, intron := sum(width(g[k]))]
    
  }
  
  return(utr3.dt)
  
}, mc.cores = 4)

stopifnot(all(elementNROWS(utr3.dt) == 1))
table(unlist(lapply(utr3.dt, length)))
table(unlist(lapply(utr3.dt, length)) == 1) # 1 is NA, i.e. across different 3' UTRs
utr3.dt <- utr3.dt[unlist(lapply(utr3.dt, length)) == 7]

utr3.dt <- rbindlist(utr3.dt)
utr3.duplexes.dt <- merge(utr3.duplexes.dt, utr3.dt, by = c("name", "ensg"), all.x = TRUE)

# Calculations
utr3.duplexes.dt[, `:=` (corrected_span = R_start - L_end - 1 - intron)]

# ==========

hiclip.utr3.dt <- copy(utr3.duplexes.dt)
paris.utr3.dt <- copy(utr3.duplexes.dt)
ricseq.utr3.dt <- copy(utr3.duplexes.dt)

hiclip.utr3.dt[!is.na(intron), span := R_start - L_end - 1]
paris.utr3.dt[!is.na(intron), span := R_start - L_end - 1]
ricseq.utr3.dt[!is.na(intron), span := R_start - L_end - 1]

fwrite(hiclip.utr3.dt, "~/Dropbox (The Francis Crick)/comp_hiclip/revisions/plots/figure_5/hiclip.utr3.dt.tsv.gz", sep = "\t")
fwrite(paris.utr3.dt, "~/Dropbox (The Francis Crick)/comp_hiclip/revisions/plots/figure_5/paris.utr3.dt.tsv.gz", sep = "\t")
fwrite(ricseq.utr3.dt, "~/Dropbox (The Francis Crick)/comp_hiclip/revisions/plots/figure_5/ricseq.utr3.dt.tsv.gz", sep = "\t")

table(hiclip.utr3.dt$span != hiclip.utr3.dt$corrected_span)
table(paris.utr3.dt$span != paris.utr3.dt$corrected_span)
table(ricseq.utr3.dt$span != ricseq.utr3.dt$corrected_span)

hiclip.utr3.dt[span <= 100 & span >= 25, grp := "medium range"]
hiclip.utr3.dt[span > 100, grp := "long range"]
hiclip.utr3.dt[span < 25, grp := "short range"]
hiclip.utr3.dt[corrected_span <= 100 & corrected_span >= 25, cor_grp := "medium range"]
hiclip.utr3.dt[corrected_span > 100, cor_grp := "long range"]
hiclip.utr3.dt[corrected_span < 25, cor_grp := "short range"]

paris.utr3.dt[span <= 100 & span >= 25, grp := "medium range"]
paris.utr3.dt[span > 100, grp := "long range"]
paris.utr3.dt[span < 25, grp := "short range"]
paris.utr3.dt[corrected_span <= 100 & corrected_span >= 25, cor_grp := "medium range"]
paris.utr3.dt[corrected_span > 100, cor_grp := "long range"]
paris.utr3.dt[corrected_span < 25, cor_grp := "short range"]

ricseq.utr3.dt[span <= 100 & span >= 25, grp := "medium range"]
ricseq.utr3.dt[span > 100, grp := "long range"]
ricseq.utr3.dt[span < 25, grp := "short range"]
ricseq.utr3.dt[corrected_span <= 100 & corrected_span >= 25, cor_grp := "medium range"]
ricseq.utr3.dt[corrected_span > 100, cor_grp := "long range"]
ricseq.utr3.dt[corrected_span < 25, cor_grp := "short range"]

table(hiclip.utr3.dt[!is.na(intron)]$grp != hiclip.utr3.dt[!is.na(intron)]$cor_grp)
table(paris.utr3.dt[!is.na(intron)]$grp != paris.utr3.dt[!is.na(intron)]$cor_grp)
table(ricseq.utr3.dt[!is.na(intron)]$grp != ricseq.utr3.dt[!is.na(intron)]$cor_grp)
