.copula_v2_kendall_diagnostic <- function(pair_data, family_num) {
  if (nrow(pair_data) < 2) {
    return(data.frame())
  }

  emp_copula <- .copula_v2_empirical_cdf_2d(pair_data$u1, pair_data$u2)

  fit_copula <- .copula_v2_bicop_cdf(

    pair_data$u1,
    pair_data$u2,
    family_num = family_num,
    par = pair_data$theta_pair,
    par2 = pair_data$zeta_pair
  )

  ok <- is.finite(emp_copula) & is.finite(fit_copula)

  data.frame(
    empirical = sort(emp_copula[ok]),
    fitted = sort(fit_copula[ok]),
    stringsAsFactors = FALSE
  )
}

.copula_v2_empirical_cdf_2d <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out)

  x_ok <- as.numeric(x[ok])
  y_ok <- as.numeric(y[ok])
  n <- length(x_ok)
  y_rank <- match(y_ok, sort(unique(y_ok)))
  tree <- integer(max(y_rank))

  update_tree <- function(index) {
    while (index <= length(tree)) {
      tree[index] <<- tree[index] + 1L
      index <- index + bitwAnd(index, -index)
    }
  }
  query_tree <- function(index) {
    total <- 0L
    while (index > 0L) {
      total <- total + tree[index]
      index <- index - bitwAnd(index, -index)
    }
    total
  }

  order_x <- order(x_ok, y_ok)
  sorted_x <- x_ok[order_x]
  values <- numeric(n)
  start <- 1L
  while (start <= n) {
    end <- start
    while (end < n && sorted_x[end + 1L] == sorted_x[start]) end <- end + 1L
    group <- order_x[start:end]
    for (row in group) update_tree(y_rank[row])
    for (row in group) values[row] <- query_tree(y_rank[row]) / n
    start <- end + 1L
  }
  out[ok] <- values
  out
}

.copula_v2_tail_diagnostics <- function(pair_data, family_num, thresholds = c(0.05, 0.10, 0.20)) {
  rows <- lapply(thresholds, function(alpha) {
    c_lower <- .copula_v2_bicop_cdf(

      rep(alpha, nrow(pair_data)),
      rep(alpha, nrow(pair_data)),
      family_num = family_num,
      par = pair_data$theta_pair,
      par2 = pair_data$zeta_pair
    )

    upper_cut <- 1 - alpha

    c_upper_cut <- .copula_v2_bicop_cdf(

      rep(upper_cut, nrow(pair_data)),
      rep(upper_cut, nrow(pair_data)),
      family_num = family_num,
      par = pair_data$theta_pair,
      par2 = pair_data$zeta_pair
    )

    lower_fit <- mean(c_lower, na.rm = TRUE)

    upper_fit <- mean(1 - 2 * upper_cut + c_upper_cut, na.rm = TRUE)

    data.frame(
      threshold = alpha,
      tail = c("Lower", "Upper"),
      empirical = c(
        mean(pair_data$u1 <= alpha & pair_data$u2 <= alpha, na.rm = TRUE),
        mean(pair_data$u1 >= upper_cut & pair_data$u2 >= upper_cut, na.rm = TRUE)
      ),
      fitted = c(lower_fit, upper_fit),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.copula_v2_tail_long_data <- function(tail_df) {
  if (nrow(tail_df) == 0) {
    return(data.frame())
  }

  rbind(
    data.frame(
      threshold = tail_df$threshold,
      tail = tail_df$tail,
      source = "Empirical",
      probability = tail_df$empirical,
      stringsAsFactors = FALSE
    ),
    data.frame(
      threshold = tail_df$threshold,
      tail = tail_df$tail,
      source = "Fitted",
      probability = tail_df$fitted,
      stringsAsFactors = FALSE
    )
  )
}

.copula_v2_conditional_tail_diagnostics <- function(tail_df) {
  if (nrow(tail_df) == 0) {
    return(tail_df)
  }

  out <- tail_df

  out$empirical <- out$empirical / out$threshold

  out$fitted <- out$fitted / out$threshold

  out$empirical <- pmin(pmax(out$empirical, 0), 1)

  out$fitted <- pmin(pmax(out$fitted, 0), 1)

  out
}
