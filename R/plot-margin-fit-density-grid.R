#' @keywords internal
#' @noRd
.plot_margin_density_grid <- function(y, family, params, grid_n = 200, group = NULL) {
  y <- as.numeric(y)
  keep <- is.finite(y)
  y <- y[keep]
  family_args <- attr(params, "family_args", exact = TRUE) %||% list()
  params <- lapply(params, function(x) {
    x <- as.numeric(x)
    if (length(x) == length(keep)) {
      x[keep]
    } else {
      rep(x, length.out = length(y))
    }
  })
  family_args <- lapply(family_args, function(x) {
    if (length(x) == length(keep)) x[keep] else rep(x, length.out = length(y))
  })
  attr(params, "family_args") <- family_args
  if (!is.null(group)) {
    group <- as.character(group)[keep]
  }

  if (length(y) < 1L) {
    return(data.frame(response = numeric(), density = numeric(), split_group = character()))
  }

  family_name <- .plot_margin_family_name(family)
  is_discrete <- .plot_margin_is_discrete(family)
  support_bounds <- .plot_margin_support_bounds(family, params)
  response_range <- range(y, na.rm = TRUE)
  pad <- diff(response_range) * 0.04
  if (!is.finite(pad) || pad <= 0) {
    pad <- max(abs(response_range), 1) * 0.04
  }

  build_one <- function(idx, group_name = "All") {
    yy <- y[idx]
    if (is_discrete) {
      x_min <- floor(min(yy, na.rm = TRUE))
      x_max <- ceiling(max(yy, na.rm = TRUE))
      if (is.finite(support_bounds["lower"])) {
        x_min <- max(x_min, ceiling(support_bounds["lower"]))
      }
      if (is.finite(support_bounds["upper"])) {
        x_max <- min(x_max, floor(support_bounds["upper"]))
      }
    } else {
      x_min <- min(yy, na.rm = TRUE) - pad
      x_max <- max(yy, na.rm = TRUE) + pad
      if (is.finite(support_bounds["lower"])) {
        x_min <- max(x_min, support_bounds["lower"])
        if (x_min <= support_bounds["lower"]) {
          yy_inside <- yy[is.finite(yy) & yy > support_bounds["lower"]]
          if (length(yy_inside) > 0L) {
            range_floor <- support_bounds["lower"] + (x_max - support_bounds["lower"]) * 0.01
            x_min <- max(
              as.numeric(stats::quantile(yy_inside, probs = 0.025, type = 8, names = FALSE)),
              range_floor,
              support_bounds["lower"] + max(diff(range(yy_inside, na.rm = TRUE)), abs(support_bounds["lower"]), 1) * 1e-6
            )
          }
        }
      }
      if (is.finite(support_bounds["upper"])) {
        x_max <- min(x_max, support_bounds["upper"])
      }
    }
    if (!is.finite(x_min) || !is.finite(x_max) || x_max <= x_min) {
      x_min <- min(yy, na.rm = TRUE)
      x_max <- max(yy, na.rm = TRUE)
    }
    x_grid <- if (is_discrete) {
      seq.int(from = x_min, to = x_max, by = 1)
    } else {
      seq(x_min, x_max, length.out = grid_n)
    }
    density <- vapply(x_grid, function(x_value) {
      par_i <- .gl_subset_family_params(params, idx)
      d_value <- suppressWarnings(.gl_call_family_fun("d", family_name, rep(x_value, length(idx)), par_i))
      d_value[!is.finite(d_value)] <- NA_real_
      mean(d_value, na.rm = TRUE)
    }, numeric(1), USE.NAMES = FALSE)
    density[!is.finite(density)] <- NA_real_
    data.frame(
      response = x_grid,
      density = density,
      split_group = group_name,
      stringsAsFactors = FALSE
    )
  }

  if (is.null(group)) {
    return(build_one(seq_along(y)))
  }

  groups <- unique(group[!is.na(group)])
  out <- lapply(groups, function(g) build_one(which(group == g), group_name = g))
  do.call(rbind, out)
}
