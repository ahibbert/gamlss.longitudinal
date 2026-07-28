#' Calculate the likelihood components for the joint model
#'
#' This function calculates the marginal and copula log likelihoods and components
#' for the joint model by organizing response data by margin and subject for
#' efficient pair-based copula calculations.
#'
#' @param response A numeric vector of response values.
#' @param response_margin A numeric vector indicating the margin (time) for each response.
#' @param response_subject A numeric vector indicating the subject for each response.
#' @return A list containing:
#' \item{log_lik}{A named vector with marginal, copula, and joint log-likelihoods.}
#' \item{margin_d}{A numeric vector of marginal densities.}
#' \item{copula_d}{A numeric vector of copula densities.}
#' \item{margin_p}{A numeric vector of marginal distribution function values.}
#' \item{Fx_1_2}{A matrix of marginal distribution function values for pairs of margins.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#' \item{margin_deriv}{A list of marginal derivatives.}
#' \item{order_copula}{A matrix indicating the order of margins and subjects for copula calculations.}
#'
#' @keywords internal
#' @noRd
build_copula_pair_cache <- function(response, response_margin, response_subject) {
  margin_names <- sort(unique(response_margin))
  num_margins <- length(margin_names)
  n_obs <- length(response)
  obs_response <- !is.na(response)

  base_df <- data.frame(
    row_id = seq_len(n_obs),
    time = response_margin,
    subject = response_subject,
    observed = obs_response,
    stringsAsFactors = FALSE
  )

  pair_df_all <- list()
  if (num_margins > 1) {
    for (i in seq_len(num_margins - 1)) {
      t1 <- margin_names[i]
      t2 <- margin_names[i + 1]

      left <- base_df[base_df$time == t1, c("row_id", "subject", "time", "observed")]
      right <- base_df[base_df$time == t2, c("row_id", "subject", "time", "observed")]
      names(left) <- c("row_id1", "subject", "time1", "observed1")
      names(right) <- c("row_id2", "subject", "time2", "observed2")

      pair_i <- merge(left, right, by = "subject", all = FALSE)
      if (nrow(pair_i) > 0) {
        pair_df_all[[length(pair_df_all) + 1]] <- pair_i
      }
    }
  }

  if (length(pair_df_all) == 0) {
    pair_df <- data.frame(
      subject = response_subject[0],
      row_id1 = integer(0),
      time1 = response_margin[0],
      observed1 = logical(0),
      row_id2 = integer(0),
      time2 = response_margin[0],
      observed2 = logical(0)
    )
  } else {
    pair_df <- do.call(rbind, pair_df_all)
  }

  order_copula <- as.matrix(pair_df[, c("time1", "subject", "time2", "subject")])
  colnames(order_copula) <- c("time1", "subject1", "time2", "subject2")

  observed_pair_base <- rep(FALSE, nrow(pair_df))
  if (nrow(pair_df) > 0) {
    observed_pair_base <- pair_df$observed1 & pair_df$observed2
  }

  Fx_1_2_template <- matrix(NA_real_, nrow = nrow(pair_df), ncol = 2)
  colnames(Fx_1_2_template) <- c("u1", "u2")

  theta_rows <- which(response_margin %in% margin_names[seq_len(max(1, num_margins - 1))])
  theta_index_map <- rep(NA_integer_, n_obs)
  theta_index_map[theta_rows] <- seq_along(theta_rows)

  cache <- list(
    row_id1 = pair_df$row_id1,
    row_id2 = pair_df$row_id2,
    Fx_1_2_template = Fx_1_2_template,
    order_copula = order_copula,
    observed_pair_base = observed_pair_base,
    theta_index_map = theta_index_map,
    margin_names = margin_names,
    num_margins = num_margins,
    n_obs = n_obs
  )
  cache
}
