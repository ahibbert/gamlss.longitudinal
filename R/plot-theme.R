#' Plot theme for gamlss.longitudinal
#'
#' A restrained, publication-ready theme used by all package plotting helpers.
#' Use it to make additional ggplot layers visually consistent with package
#' output.
#'
#' @param base_size Base font size in points.
#' @param base_family Base font family.
#'
#' @return A ggplot2 theme object.
#' @export
theme_gamlss_longitudinal <- function(base_size = 11, base_family = "") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", colour = "#1F2933", margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(colour = "#52606D", size = ggplot2::rel(0.9), margin = ggplot2::margin(b = 8)),
      plot.caption = ggplot2::element_text(colour = "#7B8794", hjust = 0),
      axis.title = ggplot2::element_text(colour = "#323F4B"),
      axis.text = ggplot2::element_text(colour = "#52606D"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "#E7ECF0", linewidth = 0.35),
      strip.background = ggplot2::element_rect(fill = "#F0F4F8", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", colour = "#323F4B"),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(face = "bold"),
      plot.margin = ggplot2::margin(8, 10, 8, 8)
    )
}

.gl_plot_colours <- function() {
  c(
    observed = "#3C4858",
    fitted = "#0072B2",
    interval = "#56B4E9",
    reference = "#D55E00",
    warning = "#E69F00",
    neutral = "#B8C2CC",
    light = "#E7ECF0"
  )
}

.gl_group_colours <- function(n) {
  n <- as.integer(n)
  if (!is.finite(n) || n < 1L) return(character())
  palette <- c(
    "#0072B2", "#E69F00", "#009E73", "#CC79A7",
    "#D55E00", "#56B4E9", "#F0E442", "#000000"
  )
  if (n <= length(palette)) palette[seq_len(n)] else grDevices::hcl.colors(n, "Dark 3")
}

.gl_new_plot_result <- function(x, primary = c("dashboard", "plot")) {
  primary <- match.arg(primary)
  structure(x, primary = primary, class = c("gamlss_longitudinal_plot", "list"))
}

#' @export
print.gamlss_longitudinal_plot <- function(x, ...) {
  primary <- attr(x, "primary") %||% "dashboard"
  plot_object <- x[[primary]]
  if (!is.null(plot_object)) print(plot_object)
  invisible(x)
}
