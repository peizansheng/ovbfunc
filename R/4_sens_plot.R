# --------------------------- Helper Functions --------------------------------#

#' Set ggplot2 theme
#'
#' @return A ggplot2 theme object.
#' @noRd
base_theme <- function() {
  theme(
    legend.position = "bottom",
    text = element_text(size = 15, family = "CMU Serif"),
    plot.title = element_text(hjust = 0.5, size = 15),
    panel.grid.major = element_line(color = "grey95"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    axis.line = element_line(color = "black"),
    strip.background = element_rect(fill = "white", color = "white")
  )
}

#' Save ggplot object to PDF
#'
#' @param plot_obj A ggplot object.
#' @param filename Output filename.
#'
#' @return A saved PDF file.
#' @noRd
save_plot <- function(plot_obj, filename) {
  ggsave(
    filename, plot = plot_obj, device = CairoPDF, width = 20, height = 15, units = "cm"
  )
}
