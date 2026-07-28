.copula_v2_transform_data <- function(data, transform = "uniform") {
  # Transform uniform [0,1] data to normal scale or other scales

  if (transform == "normal") {
    # Clamp to avoid infinite values from qnorm at 0 or 1

    data$u1 <- stats::qnorm(.copula_v2_clamp01(data$u1))

    data$u2 <- stats::qnorm(.copula_v2_clamp01(data$u2))
  }

  data
}
