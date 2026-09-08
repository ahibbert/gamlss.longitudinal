jss_plot_colours <- function() {
  c(
    package = "#0072B2",
    fitted = "#0072B2",
    comparator = "#6B7280",
    observed = "#3C4858",
    interval = "#A7D8F0",
    neutral = "#C5CDD5",
    truth = "#111827",
    reference = "#52606D",
    warning = "#E69F00",
    negative = "#B35806",
    positive = "#2166AC",
    low = "#F7F9FB"
  )
}

jss_method_colours <- function(method = NULL) {
  colours <- jss_plot_colours()
  values <- c(
    "gamlss.longitudinal" = colours[["package"]],
    "ours_rs_joint" = colours[["package"]],
    "rs_joint" = colours[["package"]],
    "gamlss2" = colours[["comparator"]],
    "gamlss" = colours[["comparator"]],
    "rs_separate" = colours[["comparator"]]
  )
  if (is.null(method)) values else unname(values[method])
}

jss_method_fills <- function(method = NULL) {
  colours <- jss_plot_colours()
  values <- c(
    "gamlss.longitudinal" = colours[["interval"]],
    "ours_rs_joint" = colours[["interval"]],
    "rs_joint" = colours[["interval"]],
    "gamlss2" = colours[["neutral"]],
    "gamlss" = colours[["neutral"]],
    "rs_separate" = colours[["neutral"]]
  )
  if (is.null(method)) values else unname(values[method])
}

theme_jss_paper <- function(base_size = 9, legend_position = "top", show_title = FALSE) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.title = if (isTRUE(show_title)) {
        ggplot2::element_text(face = "bold", colour = "#1F2933", size = ggplot2::rel(1))
      } else {
        ggplot2::element_blank()
      },
      plot.subtitle = ggplot2::element_blank(),
      plot.caption = ggplot2::element_text(
        colour = "#52606D", size = ggplot2::rel(0.86), hjust = 0,
        margin = ggplot2::margin(t = 5)
      ),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "#E7ECF0", linewidth = 0.3),
      axis.title = ggplot2::element_text(colour = "#323F4B"),
      axis.text = ggplot2::element_text(colour = "#52606D"),
      strip.background = ggplot2::element_rect(fill = "#F0F4F8", colour = NA),
      strip.text = ggplot2::element_text(face = "bold", colour = "#323F4B"),
      strip.text.y = ggplot2::element_text(angle = 0),
      strip.text.y.right = ggplot2::element_text(angle = 0),
      legend.position = legend_position,
      legend.title = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(6, 8, 6, 6)
    )
}

jss_figure_dimensions <- function(size = c("standard", "shallow", "dashboard", "tall")) {
  size <- match.arg(size)
  dimensions <- list(
    standard = c(width = 7.1, height = 4.5),
    shallow = c(width = 7.1, height = 2.8),
    dashboard = c(width = 7.1, height = 5.5),
    tall = c(width = 7.1, height = 6.2)
  )
  dimensions[[size]]
}

jss_save_figure <- function(path, plot, size = "standard", dpi = 300, width = NULL, height = NULL) {
  dimensions <- jss_figure_dimensions(size)
  if (is.null(width)) width <- unname(dimensions[["width"]])
  if (is.null(height)) height <- unname(dimensions[["height"]])
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
  normalizePath(path, winslash = "/", mustWork = TRUE)
}
