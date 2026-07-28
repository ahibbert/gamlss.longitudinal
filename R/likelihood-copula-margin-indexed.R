.calc_dlcopdpar_indexed <- function(
    row_id1,
    row_id2,
    dcdu1,
    dcdu2,
    copula_d,
    F_nd,
    n_obs,
    pair_complete = NULL) {
  row_id1 <- as.integer(row_id1)
  row_id2 <- as.integer(row_id2)
  n_pair <- length(row_id1)

  if (length(row_id2) != n_pair || length(dcdu1) != n_pair || length(dcdu2) != n_pair ||
    length(copula_d) != n_pair) {
    stop("Copula derivative inputs have inconsistent pair lengths.", call. = FALSE)
  }
  if (length(F_nd) != n_obs) {
    stop("F derivative length does not match the number of observations.", call. = FALSE)
  }

  if (is.null(pair_complete)) {
    pair_complete <- rep(TRUE, n_pair)
  } else {
    pair_complete <- as.logical(pair_complete)
    if (length(pair_complete) != n_pair) {
      stop("pair_complete length does not match copula pair length.", call. = FALSE)
    }
  }

  dlogc_row1 <- (as.numeric(dcdu1) * as.numeric(F_nd[row_id1])) / as.numeric(copula_d)
  dlogc_row2 <- (as.numeric(dcdu2) * as.numeric(F_nd[row_id2])) / as.numeric(copula_d)
  dlogc_row1[!pair_complete | !is.finite(dlogc_row1)] <- 0
  dlogc_row2[!pair_complete | !is.finite(dlogc_row2)] <- 0

  out <- numeric(n_obs)
  valid1 <- is.finite(row_id1) & row_id1 >= 1L & row_id1 <= n_obs
  valid2 <- is.finite(row_id2) & row_id2 >= 1L & row_id2 <= n_obs
  if (any(valid1)) {
    sum1 <- rowsum(dlogc_row1[valid1], row_id1[valid1], reorder = FALSE)
    out[as.integer(rownames(sum1))] <- out[as.integer(rownames(sum1))] + sum1[, 1]
  }
  if (any(valid2)) {
    sum2 <- rowsum(dlogc_row2[valid2], row_id2[valid2], reorder = FALSE)
    out[as.integer(rownames(sum2))] <- out[as.integer(rownames(sum2))] + sum2[, 1]
  }
  out
}
