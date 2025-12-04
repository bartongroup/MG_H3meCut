annotate_peaks_chipseeker <- function(bed) {
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

  pk <- annotatePeak(
    gr,
    tssRegion = c(0, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    level = "gene"
  ) |> 
    as_tibble() |> 
    clean_names() |> 
    rename(
      entrez_id = gene_id,
      gene_symbol = symbol,
      gene_id = ensembl,
      gene_name = genename
    )
  
  # Extract gene biotype from org.Hs.eg.db and add to the result
  gids <- pk$entrez_id
  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gids,
    columns = "GENETYPE",
    keytype = "ENTREZID"
  )

  pk |> 
    mutate(gene_biotype = gene_info$GENETYPE)
}

###############################################################

annotate_peaks_annotatr <- function(bed, annotations) {
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

  an <- annotate_regions(gr, annotations) |> 
    as_tibble() |> 
    clean_names()  |> 
    rename(
      chr = seqnames,
      gene_chr = annot_seqnames,
      entrez_id = annot_gene_id,
      transcript_id = annot_tx_id,
      gene_symbol = annot_symbol,
      annotation = annot_type
    ) |> 
    mutate(
      annotation = annotation |> str_remove("hg38_genes_") |> str_remove("s$"),
      peak_id = str_glue("{chr}:{start}-{end}{strand}")
    ) |> 
    relocate(peak_id, .before = 1)

  # Extract gene biotype from org.Hs.eg.db and add to the result
  gids <- an$entrez_id |> 
    unique()
  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db,
    keys = gids,
    columns = "GENETYPE",
    keytype = "ENTREZID"
  ) |> 
    as_tibble() |> 
    rename(
      entrez_id = ENTREZID,
      gene_biotype = GENETYPE
    )
 
  an |> 
    left_join(gene_info, by = c("entrez_id"))
}

######################################################################

annotate_all_chipseeker <- function(beds) {
  ids <- unique(beds$id)
  map(ids, function(i) {
    beds |> 
      filter(id == i) |> 
      annotate_peaks_bed()
  }) |> 
    list_rbind()
}

######################################################################

annotate_all_annotatr <- function(bed) {
  annotations <- build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
  seqlevels(annotations) <- seqlevels(annotations) |> str_remove("^chr")

  annotate_peaks_annotatr(bed, annotations)

}

######################################################################

# annotatr leaves a lot of transcripts without gene names
annotate_missing_genes <- function(annot, transcripts) {
  an <- annot |> 
    mutate(row_id = row_number()) |> 
    mutate(transcript_id = str_remove(transcript_id, "\\..*$")) |> 
    select(-gene_biotype)

  good_data <- an |> filter(!is.na(gene_symbol))
  bad_data <- an |> filter(is.na(gene_symbol))

  added <- bad_data |>   
    select(-c(gene_symbol)) |> 
    left_join(transcripts |> select(transcript_id, gene_symbol), by = "transcript_id")

  bind_rows(good_data, added) |> 
    left_join(transcripts |> select(transcript_id, gene_biotype), by = "transcript_id") |> 
    arrange(row_id) |> 
    select(-row_id)
}

######################################################################

combine_genes_annotatr <- function(annot) {
  annot |> 
    select(peak_id:cell_line, annotation, entrez_id, gene_symbol, gene_biotype) |> 
    distinct() |> 
    filter(annotation != "promoter")
}

######################################################################

combine_hits_rbp <- function(annot) {
  annot |> 
    group_split(rbp_symbol, cell_line) |> 
    map(function(an) {
      ns <- an |> 
        count(id) |> 
        arrange(desc(n))
      lst <- map(ns$id, function(selid) {
        an |> 
          filter(id == selid) |> 
          makeGRangesFromDataFrame(keep.extra.columns = TRUE)
      })
      Reduce(function(x, y) subsetByOverlaps(x, y), lst) |> 
        as_tibble()
    }) |> 
    list_rbind() |> 
    rename(chr = seqnames) |> 
    relocate(peak_id, .before = 1)
}