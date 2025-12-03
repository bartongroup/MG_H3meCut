library(DiffBind)

read_and_count_diffbind <- function(meta_db, min_overlap = 2) {
  
  ## Read all peaks & metadata 
  db <- DiffBind::dba(sampleSheet = meta_db)
  
  ## Build per-condition consensus peaksets
  peakset.cons <- DiffBind::dba.peakset(
    db,
    consensus = DBA_CONDITION,
    minOverlap = min_overlap
  )
  
  ## Restrict to consensus peaks and extract as GRanges
  peakset.cons <- DiffBind::dba(
    peakset.cons,
    mask = peakset.cons$masks$Consensus,
    minOverlap = 1
  )
  peaks.cons <- dba.peakset(peakset.cons, bRetrieve = TRUE)

  # Make sure multimappers are not silently dropped
  db$config$mapQCth <- 0

  db <- DiffBind::dba.count(
    db,
    peaks = peaks.cons,
    summits = FALSE,
    bRemoveDuplicates = FALSE,
    bUseSummarizeOverlaps = TRUE,
    bParallel = FALSE
  )
  db
}


meta <- read.csv("diffbind_metadata.csv")
db <- read_and_count_diffbind(meta, min_overlap = 2)
saveRDS(db, file = "diffbind_counted.rds")