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
      treatment = case_match(
        treatment,
        "Nt" ~ "Ctrl",
        "set" ~ "SETDB1"
      ),
      time_point = case_match(
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


rmsk_stats <- function(rmsk, dec) {
  dc <- dec |> 
    distinct(peak_id, rep_id) |> 
    count(rep_id)

  fw <- rmsk |> 
    group_by(rep_fam) |> 
    mutate(width = end - start + 1) |> 
    summarise(median_width = median(width)) 

  rm <- rmsk |>
    left_join(dc, by = join_by(rep_id)) |> 
    mutate(n = replace_na(n, 0))

  rn <- rm |> 
    group_by(rep_fam) |> 
    summarise(mean_n_peaks = mean(n))

  rp <- rm |> 
    mutate(with_peak = if_else(n > 0, "Y", "N")) |> 
    group_by(rep_fam, with_peak) |> 
    summarise(n = n()) |> 
    pivot_wider(id_cols = rep_fam, names_from = with_peak, values_from = n) |> 
    drop_na() |> 
    mutate(
      n_rmask = Y + N,
      perc = 100 * Y / n_rmask
    ) |> 
    select(-c(Y, N)) |> 
    arrange(desc(perc)) |> 
    ungroup()

  rp |> 
    left_join(rn, by = join_by(rep_fam)) |> 
    left_join(fw, by = join_by(rep_fam))
}