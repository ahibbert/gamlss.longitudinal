# The model-based joint approximation is inverse penalized information,
# not the frequentist A^-1 I A^-1 or a cluster sandwich covariance.
.gl_joint_solve <- function(H, P, q, raw_gradient, control,
                            source = "analytical", reference_hessian = NULL,
                            gradient_steps = NULL) {
  n <- length(q)
  valid_matrix <- function(x) is.matrix(x) && is.numeric(x) &&
    identical(dim(x), c(n, n)) && all(is.finite(x)) &&
    identical(rownames(x), names(q)) && identical(colnames(x), names(q))
  if (!valid_matrix(H) || !valid_matrix(P)) {
    .gl_joint_stop("joint_curvature_invalid", "Joint Hessian/penalty coordinates are malformed or nonfinite.")
  }
  scale <- max(1, max(abs(P)))
  if (max(abs(P - t(P))) > control$symmetry_tol * scale ||
      min(eigen((P + t(P)) / 2, symmetric = TRUE, only.values = TRUE)$values) < -control$rank_tol * scale) {
    .gl_joint_stop("joint_penalty_invalid", "The joint penalty must be symmetric positive semidefinite.")
  }
  hscale <- max(max(abs(H)), .Machine$double.eps)
  if (max(abs(H - t(H))) / hscale > control$symmetry_tol) {
    .gl_joint_stop("joint_hessian_asymmetric", "Joint likelihood curvature is not symmetric within tolerance.")
  }
  H <- (H + t(H)) / 2
  P <- (P + t(P)) / 2
  A <- -H + P
  eigen_A <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigen_A)) || min(eigen_A) <= 0) {
    .gl_joint_stop("joint_information_not_positive_definite",
      "Joint inverse-penalized inference requires positive definite information in the fitted constrained coordinates.",
      information_eigenvalues = eigen_A)
  }
  scaling <- 1 / sqrt(diag(A))
  scaled_A <- A * tcrossprod(scaling)
  rank_A <- qr(scaled_A, tol = control$rank_tol)$rank
  if (rank_A < n) {
    .gl_joint_stop("joint_information_rank_deficient", "Joint penalized information is numerically rank deficient.")
  }
  condition_A <- kappa(scaled_A, exact = TRUE)
  if (!is.finite(condition_A) || condition_A > control$condition_max) {
    .gl_joint_stop("joint_information_ill_conditioned", "Joint penalized information is too ill-conditioned for inference.")
  }
  if (!is.numeric(raw_gradient) || length(raw_gradient) != n ||
      any(!is.finite(raw_gradient)) || !identical(names(raw_gradient), names(q))) {
    .gl_joint_stop("joint_score_invalid", "The joint fitted likelihood score is nonfinite or misaligned.")
  }
  penalized_score <- raw_gradient - as.numeric(P %*% q)
  gradient <- list(gradient = penalized_score, steps = gradient_steps,
                  scaled_max = max(abs(penalized_score) / sqrt(diag(A))))
  # The common validator preserves the existing stationarity-warning policy.
  # Agreement is assessed separately on raw likelihood curvature so that a
  # large penalty cannot conceal an incorrect likelihood derivative.
  solve_control <- control
  solve_control$check_agreement <- FALSE
  result <- .gl_solve_hessian_vcov(-A, names(q), inference = solve_control,
                                  source = source, gradient = gradient)
  diagnostics <- result$hessian_diagnostics
  diagnostics$effective_thresholds$check_agreement <- control$check_agreement
  diagnostics$coefficients <- "joint"
  diagnostics$covariance_type <- "penalized_inverse"
  diagnostics$score_target <- "gradient of log likelihood minus penalty times coefficients"
  diagnostics$raw_likelihood_gradient <- raw_gradient
  diagnostics$penalized_gradient <- penalized_score
  diagnostics$penalty_rank <- qr(P, tol = control$rank_tol)$rank
  diagnostics$unpenalized_information_min_eigenvalue <-
    min(eigen(-H, symmetric = TRUE, only.values = TRUE)$values)
  diagnostics$likelihood_hessian_agreement_relative <- NA_real_
  if (!is.null(reference_hessian)) {
    if (!valid_matrix(reference_hessian)) {
      .gl_joint_stop("joint_reference_invalid", "Joint reference likelihood Hessian is malformed or nonfinite.")
    }
    hs <- 1 / sqrt(pmax(abs(diag(H)), max(hscale * sqrt(.Machine$double.eps), .Machine$double.eps)))
    delta <- (H - (reference_hessian + t(reference_hessian)) / 2) * tcrossprod(hs)
    base <- H * tcrossprod(hs)
    discrepancy <- sqrt(sum(delta^2)) / max(sqrt(sum(base^2)), .Machine$double.eps)
    diagnostics$likelihood_hessian_agreement_relative <- discrepancy
    diagnostics$hessian_agreement_relative <- discrepancy
    if (discrepancy > control$agreement_tol) {
      diagnostics$status <- "available_with_warnings"
      diagnostics$failure_codes <- unique(c(diagnostics$failure_codes, "hessian_method_disagreement"))
      diagnostics$warning_codes <- unique(c(diagnostics$warning_codes, "hessian_method_disagreement"))
      warning(.gl_inference_warning(diagnostics))
    }
  } else if (isTRUE(control$check_agreement) && identical(source, "analytical")) {
    .gl_joint_stop("joint_reference_missing", "Requested joint analytical curvature agreement was not checked.")
  }
  result$hessian_diagnostics <- diagnostics
  result$likelihood_hessian <- H
  result$penalty <- P
  result$penalized_information <- A
  result
}

.gl_joint_fd_gradient <- function(q, evaluate, step) {
  steps <- step * pmax(1, abs(q))
  out <- stats::setNames(numeric(length(q)), names(q))
  for (j in seq_along(q)) {
    plus <- minus <- q
    plus[j] <- plus[j] + steps[j]
    minus[j] <- minus[j] - steps[j]
    out[j] <- (evaluate(plus)$loglik - evaluate(minus)$loglik) / (2 * steps[j])
  }
  list(gradient = out, steps = steps)
}

.gl_joint_fd_hessian <- function(q, evaluate, step) {
  steps <- step * pmax(1, abs(q))
  H <- matrix(0, length(q), length(q), dimnames = list(names(q), names(q)))
  base <- evaluate(q)$loglik
  for (j in seq_along(q)) {
    plus <- minus <- q
    plus[j] <- plus[j] + steps[j]
    minus[j] <- minus[j] - steps[j]
    H[j, j] <- (evaluate(plus)$loglik + evaluate(minus)$loglik - 2 * base) / steps[j]^2
    if (j > 1L) for (k in seq_len(j - 1L)) {
      value <- 0
      for (a in c(-1, 1)) for (b in c(-1, 1)) {
        point <- q
        point[j] <- point[j] + a * steps[j]
        point[k] <- point[k] + b * steps[k]
        value <- value + a * b * evaluate(point)$loglik
      }
      H[j, k] <- H[k, j] <- value / (4 * steps[j] * steps[k])
    }
  }
  H
}
