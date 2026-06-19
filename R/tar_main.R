targets_main <- function() {

  get_annotations <- tar_plan(
    config_annot = list(ensembl = config$ensembl, species = config$species),
    mart = new_mart(config_annot),
    chromosomes = biomart_fetch_chromosomes(mart),
    #structures = biomart_fetch_structures(mart, gene_ids),
    gns = biomart_fetch_genes(mart),
    trans = biomart_fetch_transcripts(mart),
    terms = get_functional_terms(config_annot, gns),
    fterms = prepare_terms_fenr(terms),
    rmsk = get_repeat_masker()
  )

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    tar_target(sample_file, config$sample_file, format = "file"),
    metadata = read_metadata(sample_file),
    metadata3 = metadata |> dplyr::filter(experiment == "Me3") |> droplevels(),
    meta3 = metadata3 |> dplyr::mutate(dplyr::across(c(sample, group), ~stringr::str_remove(., "^Me3_")))
  )

  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),
    fig_ribo_prop = plot_ribo_proportion(fastq_len),

    fscrn = parse_fscreens(config, metadata),
    qcs = parse_qcs(config, metadata, paired = TRUE),
    idxstats = parse_idxstats(config, metadata),
    fig_chrom_reads = plot_chrom_proportion(idxstats, c(1:22, "X", "Y", "MT"), per_length = FALSE),
    fig_mitochondrial_prop = plot_mitochondrial_proportion(idxstats, metadata),
    
    fig_fscreen = plot_fscreen_map(fscrn),
    fig_read_qual = plot_qualities(qcs),
    fig_read_qual_clust = plot_cluster_qualities(qcs),

    bowtie_sum_human = read_bowtie_summaries(config, metadata),
    bowtie_sum_dm = read_bowtie_summaries(config, metadata, bam_dir = "bam_dm"),
    bowtie_sum_merged = merge_two_bowtie_summaries(bowtie_sum_human, bowtie_sum_dm),
    fig_bowtie_sum_human = plot_bowtie_summaries(bowtie_sum_human, metadata),
    fig_bowtie_sum_dm = plot_bowtie_summaries(bowtie_sum_dm, metadata),
    fig_spikein_ratios = plot_spikein_ratios(bowtie_sum_merged, metadata),

    picard_dup = read_picard_duplication(config, metadata),
    fig_picard_dup = plot_picard_duplication(picard_dup, metadata),

    spikein_norms = read_spikein_norms(config, metadata),
    fig_spikein_norms = plot_spikein_norms(spikein_norms |> dplyr::filter(stringr::str_detect(sample, "Me3")), metadata) 
  )

  read_data <- tar_plan(
    chrom_sizes = readr::read_tsv(file.path(config$data_path, config$dir$genome, "chrom_sizes.txt"), col_names = c("chr", "chrom_size"), show_col_types = FALSE),
    annotations = build_human_annotations()
  )

  map_peaks <- tar_map(
    values = tibble::tibble(TYPE = c("broad", "narrow")),
    names = TYPE,

    tar_plan(
      peaks = read_peak_files(config, metadata3, format = TYPE),
      fig_npeaks = plot_number_of_peaks(peaks, metadata3),
      fig_peak_width = plot_peak_width_distribution(peaks),

      # Annotate all peaks
      peaks_annotatr_raw = annotate_peaks_annotatr(peaks, annotations),
      peaks_annotatr = reduce_annotatr(peaks_annotatr_raw, trans),
      peaks_chipseeker = annotate_peaks_chipseeker(peaks),
      peaks_rmsk = annotate_peaks_rmsk(peaks, rmsk),

      # Random peaks
      prand = make_random_regions(peaks, chrom_sizes, n_regions = 100000),
      prand_annotatr = annotate_peaks_annotatr(prand, annotations) |> reduce_annotatr(trans),
      prand_chipseeker = annotate_peaks_chipseeker(prand),
      prand_rmsk = annotate_peaks_rmsk(prand, rmsk),

      fig_annotatr = plot_peak_annotatios(peaks_annotatr),
      fig_chipseeker = plot_peak_annotatios(peaks_chipseeker),

      # These are proper figures
      fig_groups_annotatr = plot_group_peak_annotations(peaks_annotatr, rand_an = prand_annotatr, palette = 2),
      fig_groups_chipseeker = plot_group_peak_annotations(peaks_chipseeker, rand_an = prand_chipseeker),
      fig_groups_rmsk = plot_group_peak_annotations(peaks_rmsk, annot_var = "rep_category", rand_an = prand_rmsk, palette = 3),

      png_groups_tss = plot_group_peak_tss(peaks_chipseeker, rand_an = prand_chipseeker, max_dist = 300) |> gs("groups_tss", 6, 4),

      # Prepare for DiffBind
      meta_db = make_diffbind_metadata(metadata3, config, type = TYPE),
      sav_meta_db = make_diffbind_metadata_for_cluster(metadata3, config, type = TYPE),
      # These steps are very slow, so we cache the result
      # db = read_and_count_diffbind(meta_db, min_overlap = 2),
      # db_spike = spikein_normalization(dbc),
    )
  )

  # Map over normalization and peak type
  map_diffbind <- tar_map(
    values = tidyr::expand_grid(
      NORM = c("libs", "spike"),
      TYPE = c("broad", "narrow")
    ) |> dplyr::mutate(
      NAME = stringr::str_glue("{NORM}_{TYPE}")
    ),
    names = NAME,

    tar_plan(
      db = read_cache(TYPE, NORM),

      fig_clust = plot_db_clustering(db, metadata3),
      fig_peaks_chr11 = plot_db_biding_region(db, metadata3, make_region("11", 40e6, 60e6)),
    ),

    # Map over two models
    tar_map(
      values = tibble::tibble(MODEL = c("full", "groups")),
      names = MODEL,

      tar_plan(
        # Perform DiffBind analysis
        dbc = de_model(db, MODEL),
        # Extract DE results in a useful format
        dec = extract_de_results(dbc),

        fig_volcano = plot_volcano(dec),
        fig_pdist = plot_pdist(dec),
        fig_up_down = plot_up_down(dec),

        # Annotate DE results
        dec_annotatr = annotate_peaks_annotatr(dec, annotations) |> reduce_annotatr(trans),
        dec_chipseeker = annotate_peaks_chipseeker(dec),
        dec_rmsk = annotate_peaks_rmsk(dec, rmsk),

        dec_annotatr_sig = dec_annotatr |> dplyr::filter(fdr < 0.01 & abs(log_fc) > 1),
        dec_chipseeker_sig = dec_chipseeker |> dplyr::filter(fdr < 0.01 & abs(log_fc) > 1),
        dec_rmsk_sig = dec_rmsk |> dplyr::filter(fdr < 0.01 & abs(log_fc) > 1),

        fig_dec_annotatr = plot_group_peak_annotations(dec_annotatr_sig, group_var = "contrast", palette = 2, n_top = NULL),
        fig_dec_chipseeker = plot_group_peak_annotations(dec_chipseeker_sig, group_var = "contrast", n_top = NULL),
        fig_dec_rmsk_da = plot_peak_rmsk(dec_rmsk_sig, fdr_limit = 0.01),
        fig_dec_rmsk_all = plot_peak_rmsk(dec_rmsk),

        fig_dec_rmsk = plot_group_peak_annotations(dec_rmsk_sig, group_var = "contrast", annot_var = "rep_category", n_top = NULL, palette = 3),

        sav_dec = write_dec(dec_chipseeker, paste0("diffbind_chipseeker_", MODEL, "_", NAME)),

        deg = prepare_de_genes(dec_chipseeker),
        fig_deg_biotype = plot_gene_biotype_distribution(deg, fdr_limit = 0.01),
        gse = fgsea_all_terms(deg, fterms, rank_expr = "-log_fc * log10(p_value)"),
        gse_tab = tabulate_gsea(gse, gns)
      )
    )

  )

  search_for_motifs <- tar_plan(
    tar_target(genome_file, file.path(config$data_path, config$dir$genome, "Homo_sapiens.canonical.fa"), format = "file"),
    db_ranges = dec_chipseeker_sig_groups_libs_broad |> distinct(peak_id, chr, start, end, strand),
    db_seqs = fetch_range_sequences(genome_file, db_ranges),
    tf_motifs = get_jaspar_matrix(),
    tf_res = search_sequences_for_motifs(db_seqs, tf_motifs),
    db_res = dec_chipseeker_sig_groups_libs_broad |> dplyr::left_join(tf_res |> dplyr::mutate(peak_id = as.character(peak_id)), by = "peak_id")
  )



 

  for_report <- tar_plan(
    fig_all_biotype = plot_all_genes_biotype(),

    consensus_peaks = db_libs_broad |> extract_consensus(),
    tbl_rmsk_peaks =  annotate_peaks_rmsk(consensus_peaks, rmsk),
    #tbl_rmsk_stats = rmsk_stats(rmsk, consensus_peaks, n_perm = 50, n_top = 20),
    tbl_rmsk_stats = readr::read_rds("cache/tbl_rmsk_stats.rds"),

    sav_rmsk_peaks = tbl_rmsk_peaks |> readr::write_csv("tab/rmask_peaks.csv"),

    db_none_broad = no_normalization(db_libs_broad),
    fig_norm_compare = plot_all_total_bindings(db_none_broad, db_libs_broad, db_spike_broad, metadata3),

    #fig_peak_example_1 = plot_region_ggbio("Me3_iPS_Ctrl_2", config, metadata3, make_region(17, 26615466, 26624227), peaks_broad, peaks_narrow),
    #fig_peak_example_2 = plot_region_ggbio("Me3_iPS_SET_2", config, metadata3, make_region(12, 34250000, 34415000), peaks_broad, peaks_narrow),

    fig_de_abu_example_1 = plot_peak_bindings(dbc_full_libs_broad, "5:70822898-70823402", metadata3),
    #fig_de_peak_example_1 = plot_region_ggbio_peak_id(dec_full_libs_broad |> dplyr::filter(contrast == "48h-iPS"), "5:70822898-70823402", config, metadata3, peaks_broad, peaks_narrow, normfac |> dplyr::mutate(norm = libs_norm)),

    normfac = tibble::tibble(
      sample = db_libs_broad$samples$SampleID,
      lib_size = db_libs_narrow$norm$DESeq2$lib.sizes,
      libs_norm = db_libs_narrow$norm$DESeq2$norm.facs,
      spike_norm = db_spike_narrow$norm$DESeq2$norm.facs
    ),
    fig_compare_spikein_norm = plot_compare_spikein_normalisation(normfac, spikein_norms),
    png_peaks_broad_ctrl = plot_peak_distribution_sample(peaks_broad, "Me3_iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_ctrl", width = 8, height = 8),
    png_peaks_narrow_ctrl = plot_peak_distribution_sample(peaks_narrow, "Me3_iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_ctrl", width = 8, height = 8),
    png_peaks_broad_set = plot_peak_distribution_sample(peaks_broad, "Me3_iPS_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_set", width = 8, height = 8),
    png_peaks_broad_48 = plot_peak_distribution_sample(peaks_broad, "Me3_48h_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_48", width = 8, height = 8),
    png_peaks_broad_11 = plot_peak_distribution_chrom(peaks_broad, "11", metadata) |> gs("peaks_broad_11", width = 8, height = 8)
  )

  for_manuscript <- tar_plan(
    mn_fdr_limit = 0.01,
    mn_logfc_limit = 0.76,
    mn_ctr = c("48h_Ctrl-iPS_Ctrl", "48h_SETDB1-48h_Ctrl", "iPS_SETDB1-iPS_Ctrl"),

    png_volcano = mn_plot_volcano(dec_chipseeker_groups_libs_broad |> filter(contrast %in% mn_ctr), fdr_limit = mn_fdr_limit, log_fc_limit = mn_logfc_limit) |> gp("volcano", width = 7, height = 3),

    tar_map(
      values = tibble::tibble(DIRECTION = c("up", "down")),
      names = DIRECTION,

      tar_plan(
        mn_dec = filter_dec(dec_chipseeker_groups_libs_broad, ctrs = mn_ctr, fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit, direction = DIRECTION),
        mn_dec_rmsk = filter_dec(dec_rmsk_groups_libs_broad, ctrs = mn_ctr, fdr_limit = mn_fdr_limit, logfc_limit = mn_logfc_limit, direction = DIRECTION),

        png_chipseeker_annotations = mn_plot_peak_annotatios(mn_dec) |> gp(paste0("chipseeker_annotations_", DIRECTION), width = 7, height = 3.5),
        png_rmsk_da = plot_peak_rmsk(mn_dec_rmsk , fdr_limit = mn_logfc_limit)|> gp(paste0("rmsk_da_", DIRECTION), width = 7, height = 4.5),
      )
    )
    
  )

  for_shiny <- tar_plan(
    sav_shiny = prepare_for_shiny(dbc_full_libs_broad, dbc_full_spike_broad, dbc_groups_libs_broad, dbc_groups_spike_broad,
                              metadata3, terms, gns)
  )


  get_environment <- tar_plan(
    tar_target(version_file, file.path(config$data_path, "environment.txt"), format = "file"),
    software_versions = get_software_versions(
      version_file,
      selection = c("fastqc", "fastq-screen", "fastp", "ribodetector", "multiqc", "bowtie2", "samtools", "bedtools", "snakemake", "drmaa", "picard", "macs3")
    )
  )


  c(
    get_environment,
    get_annotations,
    setup_experiment,
    qc,
    read_data,
    map_peaks,
    map_diffbind,
    search_for_motifs,
    for_report,
    for_manuscript,
    for_shiny
  )
  
}