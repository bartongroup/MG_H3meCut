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
