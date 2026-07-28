.copula_v2_average_density_grid <- function(family_num, pair_data, grid_n = 35, max_pairs_overlay = 300) {
  grid <- seq(0.02, 0.98, length.out = grid_n)

  grid_df <- expand.grid(u1 = grid, u2 = grid)

  pair_data <- pair_data[is.finite(pair_data$theta_pair), , drop = FALSE]

  if (nrow(pair_data) == 0) {
    grid_df$density <- NA_real_

    return(grid_df)
  }

  if (nrow(pair_data) > max_pairs_overlay) {
    set.seed(1)

    pair_data <- pair_data[sample(seq_len(nrow(pair_data)), max_pairs_overlay), , drop = FALSE]
  }

  density_sum <- rep(0, nrow(grid_df))

  density_count <- 0L

  for (i in seq_len(nrow(pair_data))) {
    par <- pair_data$theta_pair[i]

    par2 <- pair_data$zeta_pair[i]

    density_i <- tryCatch(
      {
        if (is.finite(par2)) {
          .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = par2)
        } else {
          .copula_pdf(grid_df$u1, grid_df$u2, family = family_num, par = par, par2 = 0)
        }
      },
      error = function(e) rep(NA_real_, nrow(grid_df))
    )

    if (all(!is.finite(density_i))) next

    density_i[!is.finite(density_i)] <- 0

    density_sum <- density_sum + density_i

    density_count <- density_count + 1L
  }

  if (density_count == 0L) {
    grid_df$density <- NA_real_
  } else {
    grid_df$density <- density_sum / density_count
  }

  grid_df
}

.copula_v2_density_grid_for_plot <- function(family_num,
                                             pair_data_uniform,
                                             grid_n,
                                             max_pairs_overlay,
                                             transform,
                                             is_grouped) {
  if (is_grouped) {
    density_list <- lapply(split(pair_data_uniform, pair_data_uniform$split_group), function(x) {
      grid_i <- .copula_v2_average_density_grid(
        family_num = family_num,
        pair_data = x,
        grid_n = grid_n,
        max_pairs_overlay = max_pairs_overlay
      )
      grid_i$split_group <- as.character(x$split_group[1])
      grid_i
    })
    density_grid <- do.call(rbind, density_list)
  } else {
    density_grid <- .copula_v2_average_density_grid(
      family_num = family_num,
      pair_data = pair_data_uniform,
      grid_n = grid_n,
      max_pairs_overlay = max_pairs_overlay
    )
  }

  .plot_copula_transform_grid(density_grid, transform)
}

.copula_v2_axis_labels <- function(transform) {
  list(
    x = if (transform == "normal") {
      expression(Phi^-1 * (U[t]))
    } else {
      expression(U[t])
    },
    y = if (transform == "normal") {
      expression(Phi^-1 * (U[t + 1]))
    } else {
      expression(U[t + 1])
    }
  )
}

.copula_v2_empirical_density_grid <- function(pair_data, grid_n = 35, lims = NULL) {
  if (is.null(lims)) {
    x_rng <- range(pair_data$u1, na.rm = TRUE)

    y_rng <- range(pair_data$u2, na.rm = TRUE)

    x_pad <- max(0.001, 0.025 * diff(x_rng))

    y_pad <- max(0.001, 0.025 * diff(y_rng))

    lims <- c(x_rng[1] - x_pad, x_rng[2] + x_pad, y_rng[1] - y_pad, y_rng[2] + y_pad)
  }

  kde <- MASS::kde2d(pair_data$u1, pair_data$u2, n = grid_n, lims = lims)

  data.frame(
    u1 = rep(kde$x, each = length(kde$y)),
    u2 = rep(kde$y, times = length(kde$x)),
    density = as.vector(kde$z),
    stringsAsFactors = FALSE
  )
}
