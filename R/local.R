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
  # This is used implicitly by regioneR::randomizeRegions()
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop("Package BSgenome.Hsapiens.UCSC.hg38 is required for rmsk_stats().")
  }

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
      randomize.function = function(A, ...) {
        regioneR::randomizeRegions(A = A, genome = "hg38")
      },
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


# Fisher's test of each annotation and group against pooled others
# dc created by tabulate_annotations, must contain columns 'ann', 'grp' and 'n'
annot_group_test <- function(dc) {
  dc <- dc |>
    mutate(
      ann = as.character(ann),
      grp = as.character(grp)
    )

  crossing(
    ann = unique(dc$ann),
    grp = unique(dc$grp)
  ) |>
    mutate(
      tidied = map2(ann, grp, function(ant, gr) {
        d <- dc |>
          summarise(
            a = sum(n[ann == ant & grp == gr]),
            b = sum(n[ann != ant & grp == gr]),
            c = sum(n[ann == ant & grp != gr]),
            d = sum(n[ann != ant & grp != gr])
          )
        x <- d |> 
          as.numeric() |> 
          matrix(nrow = 2, byrow = TRUE)

        depleted <- x[1, 1] / sum(x[1, ]) < x[2, 1] / sum(x[2, ])
        n_obs <- x[1, 1]

        fisher.test(x, alternative = "two.sided") |>
          broom::tidy() |>
          mutate(
            #estimate = if_else(depleted, -estimate, estimate),
            estimate = estimate,
            n_ann_grp = d$a,
            n_notann_grp = d$b,
            n_ann_notgrp = d$c,
            n_notann_notgrp = d$d
          )
      })
    ) |>
    unnest(tidied) |>
    select(
      ann, grp,
      starts_with("n_"),
      odds_ratio = estimate,
      ci_lo = conf.low,
      ci_up = conf.high,
      p_value = p.value
    ) |>
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      log_odds_ratio = if_else(odds_ratio == 0, 0, log2(odds_ratio)),
      p_in = n_ann_grp / (n_ann_grp + n_notann_grp),
      p_out = n_ann_notgrp / (n_ann_notgrp + n_notann_notgrp)
    )
}



annotation_fisher_example <- function(tst, annot = "Promoter", group = "48h_Ctrl-iPS_Ctrl_up") {
  d <- tst |> 
    filter(ann == annot & grp == group) |> 
    select(starts_with("n_"))

  tribble(
    ~., ~`With annotation`, ~`Without annotation`,
    "Significantly changing", d$n_ann_grp, d$n_notann_grp,
    "In background", d$n_ann_notgrp, d$n_notann_notgrp
  )
}




add_unannotated_peaks <- function(peaks, annotated, annot_var,
                                  none = "none",
                                  key = c("peak_id", "contrast")) {
  stopifnot(
    all(key %in% names(peaks)),
    all(key %in% names(annotated)),
    annot_var %in% names(annotated)
  )

  # The unannotated DE table should contain one row per peak and contrast.
  if(anyDuplicated(peaks[key])) {
    stop("`peaks` contains duplicated peak/contrast keys.")
  }

  annotated <- annotated |>
    mutate("{annot_var}" := coalesce(as.character(.data[[annot_var]]), none))

  annotated_keys <- annotated |>
    select(all_of(key)) |>
    distinct()

  unannotated <- peaks |>
    anti_join(annotated_keys, by = key) |>
    mutate("{annot_var}" := none)

  bind_rows(annotated, unannotated)
}

# dec contains all peaks, annotated and unannotated
annot_de_test <- function(dec, fdr_limit, logfc_limit, annot_var) {
  dc <- dec |> 
    mutate(
      ann = get(annot_var),
      sig = fdr < fdr_limit & abs(log_fc) > logfc_limit
    ) |> 
    select(peak_id, chr, start, end, width, strand, contrast, contrast_direction, sig, ann)

  ctrs <- unique(as.character(dc$contrast_direction))
  # "none" indicates a peak missing an annotation
  anns <- dc |>
    distinct(ann) |>
    filter(!is.na(ann), ann != "none") |>
    pull(ann) |>
    as.character()

  dct <- dc |> distinct(contrast, contrast_direction)

  map2(dct$contrast, dct$contrast_direction, function(ctr, ctrd) {
    d_sig <- dc |>
      filter(contrast_direction == ctrd & sig)
    d_nsig <- dc |>
      filter(contrast == ctr & !sig)

    if(n_distinct(d_sig$peak_id) < 10)
      return(NULL)

    # peaks have multiple annotations, now look at peaks only
    # we ignore GC content here, hence fake gc = 1 for make_matched_background()
    pks_sig <- d_sig |> distinct(peak_id, chr, start, end, width, strand) |> add_column(gc = 1)
    pks_nsig <- d_nsig |> distinct(peak_id, chr, start, end, width, strand) |> add_column(gc = 1)
    bkg <- make_matched_background(pks_sig, pks_nsig, gc_bins = 1)

    d_bkg <- d_nsig |> 
      filter(peak_id %in% bkg$background$peak_id)

    n_sig <- n_distinct(d_sig$peak_id)
    n_bkg <- n_distinct(d_bkg$peak_id)

    map(anns, function(an) {
      a <- d_sig |>
        filter(ann == an) |>
        summarise(n = n_distinct(peak_id)) |>
        pull(n)
      
      c <- d_bkg |>
        filter(ann == an) |>
        summarise(n = n_distinct(peak_id)) |>
        pull(n)
      
      d <- tibble(
        a = a,
        b = n_sig - a,
        c = c,
        d = n_bkg - c
      )

      x <- d |> unlist() |> matrix(nrow = 2, byrow = TRUE)

      n_obs <- x[1, 1]

      fisher.test(x, alternative = "two.sided") |>
        broom::tidy() |>
        mutate(
          estimate = estimate,
          n_ann_grp = d$a,
          n_notann_grp = d$b,
          n_ann_notgrp = d$c,
          n_notann_notgrp = d$d,
          ann = an,
          grp = ctrd
        )
    }) |>
    list_rbind() |> 
    select(
      ann,
      grp,
      starts_with("n_"),
      odds_ratio = estimate,
      ci_lo = conf.low,
      ci_up = conf.high,
      p_value = p.value
    ) |>
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      log_odds_ratio = log2(odds_ratio),
      p_in = n_ann_grp / (n_ann_grp + n_notann_grp),
      p_out = n_ann_notgrp / (n_ann_notgrp + n_notann_notgrp)
    )
  }) |> 
    list_rbind()
}



make_matched_background <- function(targets, universe, seed = 42,
                                    width_bins = 10, gc_bins = 10,
                                    n_per_target = 1) {
  set.seed(seed)

  stopifnot(all(c("peak_id", "width", "gc") %in% names(targets)))
  stopifnot(all(c("peak_id", "width", "gc") %in% names(universe)))

  target_ids <- unique(targets$peak_id)

  targets2 <- targets |>
    dplyr::ungroup() |>
    dplyr::select(peak_id, chr, start, end, strand, width, gc) |>
    dplyr::mutate(.set = "targets")

  universe2 <- universe |>
    dplyr::ungroup() |>
    dplyr::filter(!peak_id %in% target_ids) |>
    dplyr::select(peak_id, chr, start, end, strand, width, gc) |>
    dplyr::mutate(.set = "background")

  binned <- dplyr::bind_rows(targets2, universe2) |>
    dplyr::mutate(
      width_bin = dplyr::ntile(width, width_bins),
      gc_bin = dplyr::ntile(gc, gc_bins),
      bg_stratum = interaction(width_bin, gc_bin, drop = TRUE)
    )

  targets_strata <- binned |>
    dplyr::filter(.set == "targets") |>
    dplyr::count(bg_stratum, name = "n_targets")

  bg_pool <- binned |>
    dplyr::filter(.set == "background") |>
    dplyr::inner_join(targets_strata, by = "bg_stratum")

  bg_counts <- binned |>
    dplyr::filter(.set == "background") |>
    dplyr::count(bg_stratum, name = "n_background")

  diagnostics <- targets_strata |>
    dplyr::left_join(bg_counts, by = "bg_stratum") |>
    dplyr::mutate(
      n_background = dplyr::coalesce(n_background, 0L),
      n_requested = n_targets * n_per_target,
      n_sampled = pmin(n_background, n_requested),
      enough_background = n_background >= n_requested
    )

  background <- bg_pool |>
    dplyr::group_by(bg_stratum) |>
    dplyr::group_split() |>
    purrr::map_dfr(function(x) {
      n_slice <- min(nrow(x), x$n_targets[1] * n_per_target)
      dplyr::slice_sample(x, n = n_slice)
    }) |>
    dplyr::select(-.set, -n_targets)

  list(
    background = background,
    diagnostics = diagnostics
  )
}

