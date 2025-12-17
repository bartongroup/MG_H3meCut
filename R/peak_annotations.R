annotate_peaks_chipseeker <- function(bed) {
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  seqlevels(gr) <- paste0("chr",  seqlevels(gr))

  pk <- ChIPseeker::annotatePeak(
    gr,
    tssRegion = c(0, 3000),
    TxDb = TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene,
    annoDb = "org.Hs.eg.db",
    level = "gene"
  ) |> 
    as_tibble() |> 
    clean_names() |> 
    rename(
      chr = seqnames,
      entrez_id = gene_id,
      gene_symbol = symbol,
      gene_id = ensembl,
      gene_name = genename,
      annot_id = annotation
    ) |> 
    mutate(
      chr = str_remove(chr, "chr"),
      annotation = annot_id |> str_remove("(?<=^(Intron|Exon))\\s*\\([^)]*\\)"),
      annotation = fct_relevel(
        annotation, 
        c(
          "Promoter (2-3kb)",
          "Promoter (1-2kb)",
          "Promoter (<=1kb)",
          "5' UTR",
          "Exon",
          "Intron",
          "3' UTR",
          "Downstream (<1kb)",
          "Downstream (1-2kb)",
          "Downstream (2-3kb)",
          "Distal Intergenic"
        )
      )
    ) |> 
    relocate(peak_id, .before = 1)
  
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


#######################################################

annotate_peaks_annotatr <- function(bed, annotations) {
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)

  an <- annotatr::annotate_regions(gr, annotations) |> 
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
      annotation = fct_relevel(annotation, "1to5kb", "promoter", "5UTR", "exon", "intron", "3UTR"),
    ) |> 
    relocate(peak_id, .before = 1)

  # Extract gene biotype from org.Hs.eg.db and add to the result
  gids <- an$entrez_id |> 
    unique()
  gene_info <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
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


annotate_peaks_rmsk <- function(bed, rmsk) {
  if(!("peak_id" %in% colnames(bed)))
    bed <- bed |> mutate(peak_id = str_glue("{chr}:{start}-{end}"))
  
  gr <- makeGRangesFromDataFrame(bed, keep.extra.columns = TRUE)
  reps <- makeGRangesFromDataFrame(rmsk, keep.extra.columns = TRUE)

  an <- annotatr::annotate_regions(gr, reps) |> 
    as_tibble() |> 
    clean_names()  |> 
    rename(
      chr = seqnames,
      rep_id = annot_rep_id,
      rep_chr = annot_seqnames,
      rep_start = annot_start,
      rep_end = annot_end,
      rep_width = annot_width,
      rep_strand = annot_strand,
      rep_name = annot_rep_name,
      rep_class = annot_rep_class,
      rep_family = annot_rep_family,
      rep_fam = annot_rep_fam
    ) |> 
    relocate(peak_id, .before = 1) |> 
    relocate(rep_id, .before = rep_chr)
}

######################################################################

build_human_annotations <- function() {
  annotations <- annotatr::build_annotations(genome = 'hg38', annotations = 'hg38_basicgenes')
  seqlevels(annotations) <- seqlevels(annotations) |> str_remove("^chr")
  return(annotations)
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


######################################################################

get_repeat_masker <- function() {
  asm <- "hg38"

  file_path <- str_glue("https://hgdownload.soe.ucsc.edu/goldenPath/{asm}/database/rmsk.txt.gz")
  cols <- c(
    "bin","swScore","milliDiv","milliDel","milliIns",
    "genoName","genoStart","genoEnd","genoLeft",
    "strand","repName","repClass","repFamily",
    "repStart","repEnd","repLeft","id"
  )
  r <- read_tsv(file_path, col_names = cols, show_col_types = FALSE)

  r |> select(chr = genoName, start = genoStart, end = genoEnd, strand,
           rep_name = repName, rep_class = repClass, rep_family = repFamily) |> 
    mutate(
      rep_id = sprintf("R%07d", row_number()),
      chr = str_remove(chr, "chr"),
      across(c(start, end), as.integer)
    ) |> 
    unite(rep_fam, c(rep_class, rep_family), sep = ":", remove = FALSE)
}