mn_plot_peak_annotatios <- function(an, ncol = 4, group_var = "contrast", annot_var = "annotation") {

  d <- an |> 
    mutate(grp = get(group_var), ann = get(annot_var)) |> 
    distinct(grp, peak_id, ann) |> 
    count(grp, ann, .drop = FALSE) 
  
  d |> 
      ggplot(aes(x = ann, y = n, fill = ann)) +
      th() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      ) +
      geom_col() +
      facet_wrap(~grp, ncol = ncol, scales = "free_y") +
      scale_fill_brewer(palette = "Paired") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), labels = scientific_10) +
      labs(x = NULL, y = "Count")
}


mn_plot_volcano <- function(res, fdr_limit = 0.01, log_fc_limit = 0, point_size = 0.5,
                            point_alpha = 0.5, scales = "fixed") {
  
  r <- res |>
    mutate(
      x = log_fc,
      y = -log10(p_value),
      group = contrast,
      up = fdr < fdr_limit & log_fc > log_fc_limit,
      down = fdr < fdr_limit & log_fc < -log_fc_limit,
    ) |>
    select(x, y, up, down, group)

  r_up <- r |> filter(up)
  r_down <- r |> filter(down)
  r_nsel <- r |> filter(!up & !down)
  
  g <- res |>
    ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
              size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group, scale = scales) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))


    if (nrow(r_up) > 0) {
      g <- g + geom_point(data = r_up, aes(x = x, y = y), colour = "#0072B2",
                          size = point_size, alpha = point_alpha)
    }

    if (nrow(r_down) > 0) {
      g <- g + geom_point(data = r_down, aes(x = x, y = y), colour = "#D55E00",
                          size = point_size, alpha = point_alpha)
    }

  g
}
