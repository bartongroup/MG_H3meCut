targets_main <- function() {

  get_annotations <- tar_plan(
    mart = new_mart(config),
    chromosomes = biomart_fetch_chromosomes(mart),
    #structures = biomart_fetch_structures(mart, gene_ids),
    gns = biomart_fetch_genes(mart)
    #terms = get_functional_terms(config, gns),
    #fterms = prepare_terms_fenr(terms, counts$genes$id)
  )

  setup_experiment <- tar_plan(
    tar_target(config_file, "config.json", format = "file"),
    config = jsonlite::read_json(config_file),
    metadata = read_metadata(config$sample_file)
  )

  qc <- tar_plan(
    fastq_len = read_fastq_lengths(config, metadata),
    fig_ribo_prop = plot_ribo_proportion(fastq_len),

    fscrn = parse_fscreens(config, metadata),
    qcs = parse_qcs(config, metadata, paired = TRUE),
    idxstats = parse_idxstats(config, metadata),
    
    fig_fscreen = plot_fscreen_map(fscrn),
    fig_read_qual = plot_qualities(qcs),
    fig_read_qual_clust = plot_cluster_qualities(qcs),

    bowtie_sum_human = read_bowtie_summaries(file.path(config$data_path, config$dir$bam_human), metadata),
    bowtie_sum_dm = read_bowtie_summaries(file.path(config$data_path, config$dir$bam_dm), metadata)
  )


  c(
    get_annotations,
    setup_experiment,
    qc
  )
  
}