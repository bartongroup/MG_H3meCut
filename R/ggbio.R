plot_region_ggbio <- function(samples, config, meta, region, peaks_broad, peaks_narrow) {
  
  region_ucsc <- region
  seqlevels(region_ucsc) <- paste0("chr", seqlevels(region))
  genome(region_ucsc) <- "hg38"

  peaks_broad_reg <- peaks_broad |> 
    filter(chr == seqlevels(region) |> as.character(), start < end(region), end > start(region))
  peaks_narrow_reg <- peaks_narrow |> 
    filter(chr == seqlevels(region) |> as.character(), start < end(region), end > start(region))

  m <- map(samples, function(smpl) {
    rsmpl <- meta |> 
      filter(sample == smpl) |> 
      pull(raw_sample) |> 
      as.character()
    bam_file <- file.path(config$remote_path, config$dir$bam, str_glue("{rsmpl}.noMT.bam"))

    an <- tibble(text = smpl, x = -Inf, y = Inf, hj = -0.1, vj = 2)
    tr_bam <- autoplot(bam_file, which = region) +
      th() +
      geom_text(data = an, aes(label = text, x = x, y = y, hjust = hj, vjust = vj)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)))


    pks_broad <- peaks_broad_reg |> 
      filter(sample == smpl) |> 
      makeGRangesFromDataFrame()
    pks_narrow <- peaks_narrow_reg |> 
      filter(sample == smpl) |> 
      makeGRangesFromDataFrame()

    tr_pks_broad <- autoplot(pks_broad)
    tr_pks_narrow <- autoplot(pks_narrow)

    c(tr_bam, tr_pks_broad, tr_pks_narrow)
  })

  h <- rep(c(8, 1, 1), length(samples))
  tracks(unlist(m), heights = h)

  # Crashes
  #tr_genes <- autoplot(TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene, which = region_ucsc)

}
