make_fastq_rename_script <- function(samples) {
  expand_grid(
    sample = samples,
    pair = c("1", "2")
  ) |> 
  mutate(
    old_fq = str_glue("./X204SC25110565-Z01-F001/01.RawData/{sample}/{sample}_EKDL250033359-1A_237VTVLT4_L7_{pair}.fq.gz"),
    new_fq = str_glue("./fastq/{sample}_R{pair}.fastq.gz"),
    cmd = str_glue("cp {old_fq} {new_fq}")
  ) |> 
  select(cmd) |> 
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
    unite(group, c(time_point, treatment), remove = FALSE) |> 
    mutate(
      time_point = fct_relevel(time_point, "iPS"),
      treatment = fct_relevel(treatment, "Ctrl"),
      replicate = fct_reorder(replicate, as.integer(replicate))
    ) |> 
    arrange(experiment, time_point, treatment, replicate) |> 
    mutate(sample = as_factor(sample), group = as_factor(group))
}