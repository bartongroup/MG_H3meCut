make_diffbind_metadata <- function(meta, config, type = "broad", expm = "Me3") {
  if(type == "boad") {
    macs_path <- config$dir$macs_broad
    suffix <- "peaks.broadPeak"
  } else {
    macs_path <- config$dir$macs_narrow
    suffix <- "peaks.narrowPeak"
  }
  
  meta |> 
    filter(experiment == expm) |> 
    droplevels() |> 
    transmute(
      SampleID = sample,
      Tissue = time_point,
      Factor = "H3K9me3",
      Condition = group,
      Treatment = treatment,
      Replicate = replicate,
      bamReads = file.path(config$remote_path, config$dir$bam, str_glue("{raw_sample}.noMT.bam")),
      Spikein = file.path(config$remote_path, config$dir$bam_dm, str_glue("{raw_sample}.bam")),
      Peaks = file.path(config$remote_path, macs_path, str_glue("{raw_sample}_{suffix}")),
      PeakCaller = "bed",
      ScoreCol = 5
    )
}

make_diffbind_metadata_for_cluster <-  function(meta, config, type = "broad", expm = "Me3") {
  config$remote_path <- "."
  make_diffbind_metadata(meta, config, type, expm) |> 
    write_csv("diffbind_metadata.csv")
}

read_and_count_diffbind <- function(meta_db, min_overlap = 2) {
  
  ## Read all peaks & metadata 
  db <- DiffBind::dba(sampleSheet = meta_db)
  
  ## Build per-condition consensus peaksets
  peakset.cons <- DiffBind::dba.peakset(
    db,
    consensus = DiffBind::DBA_CONDITION,
    minOverlap = min_overlap
  )
  
  ## Restrict to consensus peaks and extract as GRanges
  peakset.cons <- DiffBind::dba(
    peakset.cons,
    mask = peakset.cons$masks$Consensus,
    minOverlap = 1
  )
  peaks.cons <- DiffBind::dba.peakset(peakset.cons, bRetrieve = TRUE)

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

