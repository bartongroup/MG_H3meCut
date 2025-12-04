targets_main <- function() {

  get_annotations <- tar_plan(
    config_annot = list(ensembl = config$ensembl, species = config$species),
    mart = new_mart(config_annot),
    chromosomes = biomart_fetch_chromosomes(mart),
    #structures = biomart_fetch_structures(mart, gene_ids),
    gns = biomart_fetch_genes(mart)
    #terms = get_functional_terms(config_annot, gns),
    #fterms = prepare_terms_fenr(terms, counts$genes$id)
  )

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    tar_target(sample_file, config$sample_file, format = "file"),
    metadata = read_metadata(sample_file)
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
    peaks_broad = read_peak_files(config, metadata, format = "broad"),
    peaks_narrow = read_peak_files(config, metadata, format = "narrow")
  )

  map_diffbind <- tar_map(
    values = tibble::tibble(TYPE = c("broad", "narrow")),
    names = TYPE,

    tar_plan(
      meta_db = make_diffbind_metadata(metadata, config, type = TYPE, expm = "Me3"),
      sav_meta_db = make_diffbind_metadata_for_cluster(metadata, config, type = TYPE, expm = "Me3"),
      # These steps are very slow, so we cache the result
      # db = read_and_count_diffbind(meta_db, min_overlap = 2),
      # db_spike = spikein_normalization(dbc),
      db = read_correct_cache(TYPE, spike = FALSE),
      db_spike = read_correct_cache(TYPE, spike = TRUE),

      dbc_full = de_full_model(db_spike),
      dbc_groups = de_groups(db_spike),

      dec_full = extract_de_results(dbc_full),
      dec_groups = extract_de_results(dbc_groups)

    )


  )

 

  for_report <- tar_plan(
    png_peaks_broad_ctrl = plot_peak_distribution_sample(peaks_broad, "Me3_iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_ctrl", width = 8, height = 8),
    png_peaks_narrow_ctrl = plot_peak_distribution_sample(peaks_narrow, "Me3_iPS_Ctrl_1", chrom_sizes, bar_size = 1) |> gs("peaks_narrow_dist_ctrl", width = 8, height = 8),
    png_peaks_broad_set = plot_peak_distribution_sample(peaks_broad, "Me3_iPS_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_set", width = 8, height = 8),
    png_peaks_broad_48 = plot_peak_distribution_sample(peaks_broad, "Me3_48h_SET_1", chrom_sizes, bar_size = 1) |> gs("peaks_broad_dist_48", width = 8, height = 8),
    png_peaks_broad_11 = plot_peak_distribution_chrom(peaks_broad, "11", metadata) |> gs("peaks_broad_11", width = 8, height = 8)
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
    map_diffbind,
    for_report
  )
  
}