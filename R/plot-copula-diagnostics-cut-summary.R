.copula_v2_cut_summary <- function(df, plot2_cuts, split_name = NULL) {
  if (nrow(df) < 1) {
    return(data.frame())
  }

  # Use rank-based bins to avoid collapsed quantile cuts when many fitted tau values are tied.
  df <- df[is.finite(df$tau_fit), , drop = FALSE]
  if (nrow(df) < 1) {
    return(data.frame())
  }

  effective_cuts <- min(plot2_cuts, nrow(df))
  cut_labels <- paste0("C", seq_len(effective_cuts))
  tau_rank <- rank(df$tau_fit, ties.method = "first", na.last = "keep")
  df$cut_group <- cut(tau_rank, breaks = effective_cuts, include.lowest = TRUE, labels = cut_labels)

  out <- do.call(rbind, lapply(split(df, df$cut_group), function(x) {
    tau_emp <- suppressWarnings(stats::cor(x$u1, x$u2, method = "kendall", use = "complete.obs"))
    tau_fit <- mean(x$tau_fit, na.rm = TRUE)
    data.frame(
      cut_group = as.character(x$cut_group[1]),
      tau_emp = tau_emp,
      tau_fit = tau_fit,
      n_pairs = nrow(x),
      stringsAsFactors = FALSE
    )
  }))

  if (!is.null(split_name)) {
    out$split_group <- split_name
  }
  out
}

.copula_v2_cut_summary_by_group <- function(pair_data_plot, plot2_cuts, is_grouped) {
  if (is_grouped) {
    quartile_list <- lapply(split(pair_data_plot, pair_data_plot$split_group), function(x) {
      .copula_v2_cut_summary(x, plot2_cuts, split_name = as.character(x$split_group[1]))
    })
    return(do.call(rbind, quartile_list))
  }

  .copula_v2_cut_summary(pair_data_plot, plot2_cuts)
}
