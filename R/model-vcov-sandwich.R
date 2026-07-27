#' Resolve cluster labels for sandwich variance estimation
#'
#' @noRd
.gl_sandwich_cluster <- function(object, cluster = NULL) {
  response <- object$response
  if (is.null(cluster)) {
    cluster <- object$response_subject
  }

  if (length(cluster) != length(response)) {
    stop("'cluster' must have one value per observation.", call. = FALSE)
  }

  if (any(is.na(cluster))) {
    stop("'cluster' cannot contain missing values.", call. = FALSE)
  }

  as.character(cluster)
}

#' Compute cluster-level joint log-likelihood contributions
#'
#' @noRd
.gl_cluster_joint_loglik_contributions <- function(object,
                                                   par_cov = object$par,
                                                   par_s = object$par_s,
                                                   cluster = NULL) {
  cluster <- .gl_sandwich_cluster(object, cluster = cluster)
  response <- object$response
  response_margin <- object$response_margin
  response_subject <- object$response_subject
  margin_dist <- object$margin_dist
  copula_dist <- object$copula_dist
  copula_link <- get_copula_dist(copula_dist)$copula_link
  mm <- object$model_matrix

  eta_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)
  pair_cache <- build_copula_pair_cache(response, response_margin, response_subject)
  calc_lik <- calc_likelihood_minimal(
    eta_out$eta_inv,
    mm = mm$x,
    margin_dist = margin_dist,
    copula_dist = copula_dist,
    calc_d2 = FALSE,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject,
    pair_cache = pair_cache,
    calc_margin_deriv = FALSE
  )

  cluster_levels <- unique(cluster)
  out <- stats::setNames(rep(0, length(cluster_levels)), cluster_levels)

  margin_log <- rep(0, length(response))
  margin_ok <- !is.na(calc_lik$margin_d) & calc_lik$margin_d > 0
  margin_log[margin_ok] <- log(calc_lik$margin_d[margin_ok])
  margin_log[!is.finite(margin_log)] <- 0
  margin_sum <- rowsum(margin_log, cluster, reorder = FALSE)
  out[rownames(margin_sum)] <- out[rownames(margin_sum)] + margin_sum[, 1]

  row_id1 <- calc_lik$copula_row_id1
  row_id2 <- calc_lik$copula_row_id2
  if (length(row_id1) > 0L) {
    pair_cluster_1 <- cluster[row_id1]
    pair_cluster_2 <- cluster[row_id2]
    pair_ok <- calc_lik$pair_complete & pair_cluster_1 == pair_cluster_2

    if (any(calc_lik$pair_complete & pair_cluster_1 != pair_cluster_2)) {
      stop(
        "Copula pairs cross cluster boundaries; use a cluster variable that is constant within subject.",
        call. = FALSE
      )
    }

    copula_log <- rep(0, length(row_id1))
    copula_ok <- pair_ok & calc_lik$copula_d > 0
    copula_log[copula_ok] <- log(calc_lik$copula_d[copula_ok])
    copula_log[!is.finite(copula_log)] <- 0
    copula_sum <- rowsum(copula_log, pair_cluster_1, reorder = FALSE)
    out[rownames(copula_sum)] <- out[rownames(copula_sum)] + copula_sum[, 1]
  }

  attr(out, "joint_loglik") <- unname(calc_lik$log_lik["joint"])
  out
}

#' Compute finite-difference cluster scores for all fixed coefficients
#'
#' @noRd
.gl_cluster_score_matrix <- function(object,
                                     par_cov = object$par,
                                     par_s = object$par_s,
                                     cluster = NULL,
                                     h = 1e-5) {
  if (!is.numeric(h) || length(h) != 1L || !is.finite(h) || h <= 0) {
    stop("'sandwich_h' must be a single positive finite number.", call. = FALSE)
  }

  base <- .gl_cluster_joint_loglik_contributions(
    object,
    par_cov = par_cov,
    par_s = par_s,
    cluster = cluster
  )

  p <- length(par_cov)
  scores <- matrix(NA_real_, nrow = length(base), ncol = p)
  rownames(scores) <- names(base)
  colnames(scores) <- names(par_cov)

  for (j in seq_along(par_cov)) {
    step <- h * max(1, abs(par_cov[[j]]))
    plus <- par_cov
    minus <- par_cov
    plus[[j]] <- plus[[j]] + step
    minus[[j]] <- minus[[j]] - step

    plus_ll <- .gl_cluster_joint_loglik_contributions(
      object,
      par_cov = plus,
      par_s = par_s,
      cluster = cluster
    )
    minus_ll <- .gl_cluster_joint_loglik_contributions(
      object,
      par_cov = minus,
      par_s = par_s,
      cluster = cluster
    )

    scores[, j] <- (plus_ll[names(base)] - minus_ll[names(base)]) / (2 * step)
  }

  attr(scores, "base_contributions") <- base
  attr(scores, "joint_loglik") <- attr(base, "joint_loglik")
  scores
}

#' Small-sample correction for cluster-robust covariance
#'
#' @noRd
.gl_sandwich_correction <- function(n_obs, n_clusters, n_parameters, adjust = TRUE) {
  if (!isTRUE(adjust)) {
    return(1)
  }
  if (n_clusters <= 1L || n_obs <= n_parameters) {
    warning(
      "Cannot apply sandwich small-sample correction; using correction factor 1.",
      call. = FALSE
    )
    return(1)
  }
  (n_clusters / (n_clusters - 1)) * ((n_obs - 1) / (n_obs - n_parameters))
}

#' Assemble cluster-robust sandwich variance-covariance output
#'
#' @noRd
.gl_vcov_compute_sandwich <- function(object,
                                      par_cov,
                                      par_s,
                                      mm,
                                      margin_dist,
                                      response,
                                      response_margin,
                                      response_subject,
                                      cluster = NULL,
                                      score_h = 1e-5,
                                      bread_h = 1e-4,
                                      adjust = TRUE,
                                      bread_method = c("analytical", "numderiv", "analytical_only"),
                                      progress = interactive()) {
  bread_method <- match.arg(bread_method)

  object_eval <- object
  object_eval$par <- par_cov
  object_eval$par_s <- par_s

  copula_link <- get_copula_dist(object_eval$copula_dist)$copula_link
  eta_out <- calc_eta(par_cov, mm, margin_dist, copula_link, par_s = par_s)
  method_info <- .gl_vcov_apply_margin_preflight(
    bread_method,
    bread_method,
    margin_dist,
    eta_out$eta_inv
  )
  bread_method <- method_info$method
  calc_lik_out <- calc_likelihood_minimal(
    eta_out$eta_inv,
    mm = mm$x,
    margin_dist = margin_dist,
    copula_dist = object_eval$copula_dist,
    calc_d2 = FALSE,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject,
    calc_margin_deriv = FALSE
  )
  method_info <- .gl_vcov_apply_likelihood_preflight(
    bread_method,
    bread_method,
    calc_lik_out,
    response
  )
  bread_method <- method_info$method

  bread_path <- .gl_vcov_compute_primary(
    object = object_eval,
    par_cov = par_cov,
    mm = mm,
    margin_dist = margin_dist,
    response = response,
    response_margin = response_margin,
    response_subject = response_subject,
    method = bread_method,
    progress = progress,
    h = bread_h
  )

  H <- bread_path$hessian_nd
  if (is.null(H) || !is.matrix(H)) {
    stop("Sandwich vcov requires an observed Hessian bread matrix.", call. = FALSE)
  }

  scores <- .gl_cluster_score_matrix(
    object_eval,
    par_cov = par_cov,
    par_s = par_s,
    cluster = cluster,
    h = score_h
  )

  if (!is.null(colnames(H)) && all(colnames(scores) %in% colnames(H))) {
    H <- H[colnames(scores), colnames(scores), drop = FALSE]
  }

  bread_inv <- solve(-H)
  correction <- .gl_sandwich_correction(
    n_obs = length(response),
    n_clusters = nrow(scores),
    n_parameters = ncol(scores),
    adjust = adjust
  )
  meat <- crossprod(scores) * correction
  vc <- bread_inv %*% meat %*% bread_inv
  rownames(vc) <- colnames(vc) <- colnames(scores)
  se <- sqrt(pmax(0, diag(vc)))
  names(se) <- colnames(scores)

  eig <- tryCatch(eigen((H + t(H)) / 2, symmetric = TRUE, only.values = TRUE)$values, error = function(e) NA_real_)
  diagnostics <- list(
    condition_number = tryCatch(kappa(H), error = function(e) NA_real_),
    min_abs_eigen = if (any(is.finite(eig))) min(abs(eig[is.finite(eig)])) else NA_real_,
    max_abs_eigen = if (any(is.finite(eig))) max(abs(eig[is.finite(eig)])) else NA_real_,
    bread_method = bread_path$method_used,
    n_clusters = nrow(scores),
    n_observations = length(response),
    n_parameters = ncol(scores),
    score_step = score_h,
    bread_step = bread_h,
    small_sample_correction = correction,
    max_abs_cluster_score_sum = max(abs(colSums(scores)), na.rm = TRUE),
    joint_loglik_from_clusters = sum(attr(scores, "base_contributions")),
    joint_loglik = attr(scores, "joint_loglik")
  )

  list(
    vcov_final = vc,
    se_final = se,
    hessian_diagnostics = diagnostics,
    method_used = "sandwich_cluster",
    scores = scores
  )
}
