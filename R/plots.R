okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#AAAAAA", "#000000")

th <- ggplot2::theme_bw() +
  ggplot2::theme(panel.grid = ggplot2::element_blank())

scientific_10 <- function(v, limit = 3.01) {
  sc <- map_chr(v, function(x) {
    if (is.na(x)) {
      as.character(x)
    } else if (x == 0 | abs(log10(x)) < limit) {
      format(x, scientific = FALSE)
    } else {
      format(x, scientific = TRUE) |> 
        str_replace("e", "%*%10^") |> 
        str_replace("\\+", "") |> 
        str_replace("\\^0+", "^")
    }
  })
  parse(text = sc)
}


gs <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("doc/images", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300, limitsize = FALSE)
}

gp <- function(gg, name, width, height) {
  ggplot2::ggsave(filename = file.path("fig", paste0(name, ".pdf")), plot = gg, device = "pdf",
                  width = width, height = height)
}



plot_volma <- function(res, point_size, point_alpha, scales = "fixed") {
  r <- res |>
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)

  rm(res, r)  # Minimise environment for serialisation

  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group, scale = scales)

  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

plot_ma <- function(res, expr_var = log_cpm, fc_var = log_fc, p_var = p_value, fdr_var = fdr, 
                    group_var = contrast, sel_genes = NULL, fdr_limit = 0.01, log_fc_limit = 0,
                    point_size = 0.5, point_alpha = 0.5) {
  
   if(is.null(sel_genes)) {
    res <- res |> mutate(sel = {{ fdr_var }} < fdr_limit & abs({{ fc_var }}) >= log_fc_limit)
  } else {
    res <- res |> mutate(sel = id %in% sel_genes)
  }
  
  res |>
    mutate(
      group = {{ group_var }},
      x = {{ expr_var }},
      y = {{ fc_var }},
    ) |>
    plot_volma(point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc_var = log_fc, p_var = p_value, fdr_var = fdr, group_var = contrast,
                         sel_genes = NULL, fdr_limit = 0.01, log_fc_limit = 0, point_size = 0.5,
                         point_alpha = 0.5, scales = "fixed") {
  if(is.null(sel_genes)) {
    res <- res |> mutate(sel = {{ fdr_var }} < fdr_limit & abs({{ fc_var }}) >= log_fc_limit)
  } else {
    res <- res |> mutate(sel = id %in% sel_genes)
  }
  
  res |>
    mutate(
      group = {{ group_var }},
      x = {{ fc_var }},
      y = -log10({{ p_var}})
    ) |> 
    plot_volma(point_size, point_alpha, scales) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))
}

plot_pdist <- function(res, p_var = p_value, group_var = contrast, n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |>
    mutate(p = {{ p_var }}, grp = {{ group_var }}) |>
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}

plot_vmp <- function(res, expr_var = log_cpm, fc_var = log_fc, p_var = p_value, fdr_var = fdr,
                     group_var = contrast, sel_genes = NULL, fdr_limit = 0.01, log_fc_limit = 0, point_size = 0.5, point_alpha = 0.5, n_bins = 50){
  list(
    volcano = plot_volcano(res, {{fc_var}}, {{p_var}}, {{fdr_var}}, {{group_var}}, sel_genes, fdr_limit, log_fc_limit, point_size, point_alpha, scales = "free_x"),
    ma = plot_ma(res, {{expr_var}}, {{fc_var}}, {{p_var}}, {{fdr_var}}, {{group_var}}, sel_genes, fdr_limit, log_fc_limit, point_size, point_alpha),
    pdist = plot_pdist(res, {{p_var}}, {{group_var}}, n_bins)
  )
}

