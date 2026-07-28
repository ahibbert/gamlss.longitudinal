.gl_effect_counterfactual_values <- function(x) {
  if (is.factor(x)) {
    levels(x)
  } else if (is.character(x)) {
    sort(unique(x))
  } else {
    as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE))
  }
}

.gl_effect_add_factor_calibration_rows <- function(nd, template) {
  factor_cols <- names(nd)[vapply(nd, is.factor, logical(1))]
  extras <- list()
  for (fc in factor_cols) {
    levs <- levels(nd[[fc]])
    if (length(levs) < 2L) next
    present <- unique(as.character(nd[[fc]]))
    missing <- setdiff(levs, present)
    if (length(missing) == 0L) next
    for (lev in missing) {
      row <- template[1L, , drop = FALSE]
      row[[fc]] <- factor(lev, levels = levs, ordered = is.ordered(nd[[fc]]))
      for (other_fc in factor_cols) {
        if (other_fc == fc) next
        row[[other_fc]] <- factor(
          as.character(row[[other_fc]]),
          levels = levels(nd[[other_fc]]),
          ordered = is.ordered(nd[[other_fc]])
        )
      }
      extras[[length(extras) + 1L]] <- row
    }
  }

  if (length(extras) == 0L) {
    return(nd)
  }

  rbind(nd, do.call(rbind, extras))
}
