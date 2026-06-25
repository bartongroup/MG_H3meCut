make_fastq_rename_script <- function(samples) {
  data_sets <- tribble(
    ~set, ~dir, ~infix,
    "Me3", "X204SC25110565-Z01-F001", "EKDL250033359-1A_237VTVLT4_L7",
    "Me2", "X204SC25114793-Z01-F001", "EKDL250033839-1A_237VV5LT4_L3"
  )

  map(1:nrow(data_sets), function(i) {
    r <- data_sets[i, ]

    sample_sel <- samples |> str_subset(r$set)
     
    expand_grid(
      sample = sample_sel,
      pair = c("1", "2")
    ) |> 
    mutate(
      old_fq = str_glue("./{r$dir}/01.RawData/{sample}/{sample}_{r$infix}_{pair}.fq.gz"),
      new_fq = str_glue("./fastq/{sample}_R{pair}.fastq.gz"),
      cmd = str_glue("cp {old_fq} {new_fq}")
    ) |> 
    select(cmd)
  }) |> 
  list_rbind() |> 
  add_row(cmd = "#!/usr/bin/env bash", .before = 1) |> 
  summarise(script = paste(cmd, collapse = "\n")) |> 
  pull(script) |> 
  write_lines("rename_fastq.sh")
}


read_metadata <- function(sample_file) {
  read_csv(sample_file, col_names = "raw_sample", show_col_types = FALSE) |> 
    separate_wider_delim(
      raw_sample,
      delim = "_", 
      names = c("experiment", "time_point", "treatment", "replicate"),
      cols_remove = FALSE
    ) |> 
    mutate(
      replicate = str_remove(replicate, "rep"),
      treatment = recode_values(
        treatment,
        "Nt" ~ "Ctrl",
        "set" ~ "SETDB1"
      ),
      time_point = recode_values(
        time_point,
        "48" ~ "48h",
        "iPS" ~ "iPS"
      )
    ) |> 
    unite(sample, c(experiment, time_point, treatment, replicate), remove = FALSE) |> 
    mutate(sample = str_remove(sample, "DB1")) |> 
    unite(group, c(experiment, time_point, treatment), remove = FALSE) |> 
    mutate(
      experiment = fct_relevel(experiment, "Me2"),
      time_point = fct_relevel(time_point, "iPS"),
      treatment = fct_relevel(treatment, "Ctrl"),
      replicate = fct_reorder(replicate, as.integer(replicate))
    ) |> 
    arrange(experiment, time_point, treatment, replicate) |> 
    mutate(sample = as_factor(sample), group = as_factor(group))
}


write_dec <- function(dec, name) {
  fname = file.path("tab", str_glue("{name}.csv"))
  dec |> 
    mutate(across(where(is.double), ~signif(.x, 4))) |> 
    write_csv(fname)
}



make_gr_regions <- function(df, chromosomes = c(1:22, "X", "Y")) {
  df |> 
    distinct(chr, start, end) |> 
    filter(chr %in% chromosomes) |> 
    mutate(chr = paste0("chr", chr)) |> 
    makeGRangesFromDataFrame()
}


random_overlap_test <- function(reg1, reg2, n_perm = 300) {
  pt <- regioneR::overlapPermTest(
    A = r1,
    B = r2,
    ntimes = n_perm,
    force.parallel = TRUE
  ) |> 
    suppressMessages() |> 
    suppressWarnings()
}

num_overlaps_once <- function(A, B) {
  regioneR::numOverlaps(A, B, count.once = TRUE)
}

rmsk_stats <- function(rmsk, pks, n_perm = 100, n_top = 20) {
  gr_pk <- make_gr_regions(pks)

  fams <- rmsk$rep_fam |> 
    unique()

  rc <- map(fams, function(fm) {
    gr_rm <- rmsk |> 
      filter(rep_fam == fm) |> 
      make_gr_regions()
    ov <- countOverlaps(gr_rm, gr_pk)
    n_over <- sum(ov > 0)
    tibble(
      rep_fam = fm,
      n_rmask = length(gr_rm),
      n_over = n_over,
      perc = 100 * n_over / length(gr_rm),
      median_width = median(width(gr_rm))
    )
  }) |> 
    list_rbind() |> 
    arrange(desc(perc)) 

  perm <- map(1:n_top, function(i) {
    r <- rc[i, ]
    gr_rm <- rmsk |> 
      filter(rep_fam == r$rep_fam) |> 
      make_gr_regions()

    print(paste("Permutation for", r$rep_fam))
    pt <- regioneR::permTest(
      A = gr_rm,
      B = gr_pk,
      randomize.function = regioneR::randomizeRegions,
      evaluate.function = num_overlaps_once,
      alternative = "greater",
      ntimes = n_perm,
      force.parallel = TRUE
    ) |> 
      suppressMessages() |> 
      suppressWarnings()

    tibble(
      rep_fam = r$rep_fam,
      p_value = pt$num_overlaps_once$pval,
      perm_over = pt$num_overlaps_once$permuted |> mean()
    )
  }) |> 
    list_rbind()

  rc |> 
    left_join(perm, by = join_by(rep_fam))
}




filter_dec <- function(dec, ctrs, fdr_limit = 0.01, logfc_limit = 0.76, direction = c("up", "down")) {
  direction <- match.arg(direction)
  dec |> 
    mutate(dir = case_when(
      log_fc > 0 ~ "up",
      log_fc <= 0 ~ "down"
    )) |>
    filter(contrast %in% ctrs & fdr < fdr_limit & abs(log_fc) > logfc_limit & dir == direction)
}



make_random_regions <- function(peaks, chrom_sizes, n_regions = 100000,
                                chromosomes = c(as.character(1:22), "X", "Y")) {
  cs <- chrom_sizes |>
    dplyr::filter(chr %in% chromosomes) |>
    dplyr::select(chr, chrom_size)

  widths <- peaks |>
    dplyr::filter(chr %in% chromosomes) |>
    dplyr::slice_sample(n = n_regions) |> 
    dplyr::mutate(width = end - start) |>
    dplyr::pull(width)
    
  # For each sampled width, choose a chromosome with probability proportional
  # to the number of valid start positions available on that chromosome.
  chr_id <- purrr::map_int(widths, \(w) {
    weights <- pmax(cs$chrom_size - w + 1, 0)
    sample.int(nrow(cs), size = 1, prob = weights)
  })

  max_start <- cs$chrom_size[chr_id] - widths

  tibble::tibble(
    peak_id = sprintf("peak_%d", 1:length(widths)),
    chr = cs$chr[chr_id],
    width = widths,
    start = purrr::map_int(max_start, \(m) sample.int(m + 1, size = 1) - 1L),
    end = start + width
  ) |>
    dplyr::select(peak_id, chr, start, end, width)
}

plot_peak_coverage <- function(peaks, chrom_sizes, meta, annot = NULL) {
  if(!is.null(annot)) {
    peaks <- peaks |> 
      filter(annotation == annot)
  }

  tot <- sum(chrom_sizes$chrom_size)
  
  d <- peaks |> 
    distinct(sample, chr, start, end) |> 
    filter(chr %in% chrom_sizes$chr) |> 
    mutate(width = end - start + 1) |> 
    group_by(sample) |> 
    summarise(peak_coverage = sum(width), .groups = "drop") |> 
    mutate(peak_perc = 100 * peak_coverage / tot) |> 
    left_join(meta, by = "sample")

  env <- new.env(parent = globalenv())
  env$d <- d

  with(env, {
    d |> 
      ggplot(aes(x = sample, y = peak_perc, fill = group)) +
      th() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none"
      ) +
      geom_col() +
      scale_fill_manual(values = okabe_ito_palette) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
      labs(x = "Sample", y = "Genome coverage by peaks (%)")
  })
}

get_motif_stats <- function(dec_res, dec_for_motifs) {
  d1 <- dec_for_motifs |> 
    distinct(contrast, peak_id) |> 
    count(contrast, name = "DB peaks")
  d2 <- dec_res |> 
    distinct(contrast, peak_id) |> 
    count(contrast, name = "Peaks with hits")
  d3 <- dec_res |> 
    distinct(contrast, motif_name, peak_id) |> 
    count(contrast, peak_id) |> 
    group_by(contrast) |> 
    summarise(`Motifs per peak` = mean(n), .groups = "drop")
  d4 <- dec_res |> 
    distinct(contrast, motif_name) |> 
    count(contrast, name = "Unique motifs")

  inner_join(d1, d2, by = "contrast") |>
    inner_join(d3, by = "contrast") |>
    inner_join(d4, by = "contrast")
}



tabulate_annotations <- function(an, group_var, annot_var, n_top = NULL, rand_an = NULL, order_annot = FALSE) {
  d <- an |> 
    mutate(grp = get(group_var), ann = get(annot_var))

  # Some data don't have signal
  if(!("signal" %in% names(d))) {
    d <- d |> mutate(signal = 0)
  }

  d <- d |> 
    distinct(grp, peak_id, ann, signal)

  if(!is.null(n_top)) {
    d <- d |> 
      group_by(grp) |> 
      slice_max(signal, n = n_top) |>
      ungroup() |> 
      select(-signal)
  }

  if(!is.null(rand_an)) {
    ra <- rand_an |> 
      mutate(ann = get(annot_var)) |> 
      select(peak_id, ann) |> 
      mutate(grp = "Random")
    d <- d |> 
      bind_rows(ra)
  }

  dc <- d |> 
    count(grp, ann, .drop = FALSE) |> 
    group_split(grp) |> 
    map(function(w) {
      w |> 
        mutate(tot = sum(n)) |> 
        select(grp, ann, n, tot) |> 
        nest(data = c(n, tot)) |> 
        mutate(
          tst = map(data, ~prop.test(.x$n, .x$tot)),
          tidied = map(tst, tidy)
        ) |> 
        unnest(tidied) |> 
        unnest(data) |> 
        select(-tst)
    }) |> 
    list_rbind() |> 
    rename(prop = estimate)

  if(order_annot)
    dc <- dc |> mutate(ann = fct_reorder(ann, -n))
  
  if(annot_var == "rep_category") {
    dc <- dc |>
      mutate(ann = factor(ann, levels = c("SINE", "LINE", "LTR", "Simple_repeat", "DNA", "Low_complexity", "Satellite", "Retroposon", 
"Other")))
  }

  dc
}

# Fisher's test of each annotation and group against all others
# dc created by tabulate_annotations, must contain columns 'ann', 'grp' and 'n'
annot_group_test <- function(dc) {
  annots <- unique(as.character(dc$ann))
  map(annots, function(ant) {
    # for each annotation, test against all others
    d1 <- dc |> 
      filter(ann == ant) |> 
      select(grp, n_sel = n)
    d2 <- dc |> 
      filter(ann != ant) |> 
      group_by(grp) |> 
      summarise(n_else = sum(n))
    d <- inner_join(d1, d2, by = "grp")
    # for each group test against others
    grps <- unique(as.character(d$grp))
    map(grps, function(gr) {
      x1 <- d |> 
        filter(grp == gr)
      x2 <- d |> 
        filter(grp != gr) |> 
        summarise(n_sel = sum(n_sel), n_else = sum(n_else)) |> 
        mutate(grp = "else")
      x <- bind_rows(x1, x2) |> 
        select(n_sel, n_else) |> 
        as.matrix()
      fisher.test(x, alternative = "two.sided")  |> 
        broom::tidy() |> 
        select(odds_ratio = estimate, ci_lo = conf.low, ci_up = conf.high, p_value = p.value) |> 
        add_column(ann = ant, grp = gr, .before = 1)
    }) |> 
      list_rbind()
  }) |> 
    list_rbind() |> 
    mutate(fdr = p.adjust(p_value, method = "BH"))
}