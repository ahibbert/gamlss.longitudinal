#' @export
pithist.gamlss.longitudinal <- function(object, bins = 20, randomize = NULL, seed = 1L, plot = TRUE, by_time = FALSE, by = NULL, data = NULL, ...) {
  pit_out <- .gl_pit(object, randomize = randomize, seed = seed)

  pit <- pit_out$pit

  split_info <- .gl_diag_split_info(by_time, by, pit_out$diag, data = data, plot_name = "PIT histogram")

  pit_df <- data.frame(pit = pit, time = as.factor(pit_out$diag$time))

  if (split_info$split_by) {
    pit_df$split_group <- split_info$group
  }

  pit_df <- pit_df[is.finite(pit_df$pit), , drop = FALSE]

  if (!plot) {
    return(pit_df)
  }

  if (!split_info$split_by) {
    expected <- nrow(pit_df) / bins

    return(
      ggplot2::ggplot(pit_df, ggplot2::aes(x = pit)) +
        ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1L), closed = "right", fill = "#2c7fb8", color = "white") +
        ggplot2::geom_hline(yintercept = expected, linetype = "dashed", linewidth = 0.4, color = "#444444") +
        ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
        ggplot2::labs(x = "PIT", y = "Count", title = "PIT Histogram") +
        ggplot2::theme_minimal()
    )
  }

  expected_df <- aggregate(pit ~ split_group, data = pit_df, FUN = function(x) length(x) / bins)

  names(expected_df)[2] <- "expected"

  ggplot2::ggplot(pit_df, ggplot2::aes(x = pit)) +
    ggplot2::geom_histogram(breaks = seq(0, 1, length.out = bins + 1L), closed = "right", fill = "#2c7fb8", color = "white") +
    ggplot2::geom_hline(data = expected_df, ggplot2::aes(yintercept = expected), linetype = "dashed", linewidth = 0.4, color = "#444444") +
    ggplot2::facet_wrap(~split_group, scales = "free_y") +
    ggplot2::scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    ggplot2::labs(x = "PIT", y = "Count", title = paste0("PIT Histogram by ", split_info$by)) +
    ggplot2::theme_minimal()
}
