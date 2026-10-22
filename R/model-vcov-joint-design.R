# Joint inference uses the fitted, constraint-absorbed CG coefficient coordinates.
# It never reconstructs a penalty from generic differences.
.gl_joint_stop <- function(code, message, ...) {
  stop(.gl_inference_unavailable(c(list(
    status = "unavailable", failure_codes = code, coefficients = "joint",
    message = message
  ), list(...))))
}

.gl_joint_valid_names <- function(x) {
  is.character(x) && !anyNA(x) && all(nzchar(x)) && !anyDuplicated(x)
}

.gl_joint_smooth_spec <- function(B, context) {
  spec <- attr(B, "gl_smooth_spec", exact = TRUE)
  if (!identical(attr(B, "gl_smooth_spec_version", exact = TRUE), 1L) ||
      !identical(attr(B, "gl_smooth_construction_count", exact = TRUE), 1L) ||
      !identical(attr(B, "gl_smooth_penalty_count", exact = TRUE), 1L) ||
      !is.list(spec) || length(spec$S) != 1L || !is.matrix(spec$X)) {
    .gl_joint_stop("joint_smooth_spec_unavailable", paste0(
      "Joint inference requires a retained single-component, single-penalty ",
      "fitted smooth specification: ", context,
      ". Legacy or unsupported smooth objects must be refitted or explicitly migrated and validated."
    ))
  }
  if (!identical(dim(spec$X), dim(B)) ||
      !isTRUE(all.equal(as.numeric(spec$X), as.numeric(B), tolerance = 1e-12))) {
    .gl_joint_stop("joint_smooth_basis_mismatch",
                  paste("Stored smooth basis differs from its fitted specification:", context))
  }
  spec
}

.gl_joint_design <- function(object, control) {
  mm <- object$model_matrix
  parameters <- names(mm$x)
  if (!is.list(mm$x) || !length(parameters) || !.gl_joint_valid_names(parameters) ||
      any(!parameters %in% c("mu", "sigma", "nu", "tau", "theta", "zeta"))) {
    .gl_joint_stop("joint_design_invalid", "Joint inference requires named distribution-parameter designs.")
  }
  fixed <- object$par
  if (!is.numeric(fixed) || !length(fixed) || any(!is.finite(fixed)) ||
      !.gl_joint_valid_names(names(fixed))) {
    .gl_joint_stop("joint_coefficients_invalid", "Joint inference requires finite, uniquely named fitted coefficients.")
  }
  nonempty <- function(x) names(x)[vapply(x, length, integer(1)) > 0L]
  if (length(setdiff(nonempty(object$par_s %||% list()), parameters)) ||
      length(setdiff(nonempty(mm$s %||% list()), parameters))) {
    .gl_joint_stop("joint_unmapped_smooth", "Smooth coefficients or designs have no matching parameter block.")
  }
  smooth <- setNames(rep(list(list()), length(parameters)), parameters)
  fixed_rows <- list()
  smooth_rows <- list()
  expected_fixed <- character()
  for (pn in parameters) {
    X <- as.matrix(mm$x[[pn]])
    if (!is.numeric(X) || any(!is.finite(X)) ||
        (ncol(X) > 0L && !.gl_joint_valid_names(colnames(X)))) {
      .gl_joint_stop("joint_design_invalid", paste("Invalid fixed design for", pn))
    }
    colnames(X) <- sub(paste0("^", pn, "\\."), "", colnames(X))
    fixed_names <- if (ncol(X)) paste(pn, colnames(X), sep = ".") else character()
    expected_fixed <- c(expected_fixed, fixed_names)
    if (length(fixed_names)) {
      fixed_rows[[pn]] <- data.frame(name = fixed_names, parameter = pn,
        kind = "fixed", term = colnames(X), basis_index = NA_integer_, stringsAsFactors = FALSE)
    }
    mm$x[[pn]] <- X
    blocks <- mm$s[[pn]] %||% list()
    fitted <- object$par_s[[pn]] %||% list()
    if (!is.list(fitted) && all(is.na(fitted)) && !length(blocks)) fitted <- list()
    if (!is.list(blocks) || !is.list(fitted) ||
        (length(blocks) && !.gl_joint_valid_names(names(blocks))) ||
        !setequal(names(blocks), names(fitted))) {
      .gl_joint_stop("joint_smooth_map_invalid", paste("Smooth design/coefficient map differs for", pn))
    }
    for (term in names(blocks)) {
      B <- blocks[[term]]
      b <- fitted[[term]]
      context <- paste(pn, term, sep = ":")
      if (!is.matrix(B) || !is.numeric(B) || !ncol(B) || nrow(B) != nrow(X) ||
          any(!is.finite(B)) || !is.numeric(b) || length(b) != ncol(B) ||
          any(!is.finite(b)) || !.gl_joint_valid_names(names(b))) {
        .gl_joint_stop("joint_smooth_map_invalid", paste("Invalid smooth coordinates:", context))
      }
      expected <- paste(pn, term, seq_len(ncol(B)), sep = ".")
      if (!setequal(names(b), expected)) {
        .gl_joint_stop("joint_smooth_names_invalid", paste("Unmatched smooth basis names:", context))
      }
      spec <- .gl_joint_smooth_spec(B, context)
      S <- attr(B, "penalty", exact = TRUE)
      lambda <- object$lambda_s[[pn]][[term]]
      if (!is.matrix(S) || !is.numeric(S) || !identical(dim(S), c(ncol(B), ncol(B))) ||
          any(!is.finite(S)) ||
          !isTRUE(all.equal(as.numeric(S), as.numeric(spec$S[[1L]]), tolerance = 1e-12)) ||
          !is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
        .gl_joint_stop("joint_penalty_invalid", paste("Missing or inconsistent fitted penalty/lambda:", context))
      }
      scale <- max(1, max(abs(S)))
      eigen_S <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
      if (max(abs(S - t(S))) > control$symmetry_tol * scale ||
          min(eigen_S) < -control$rank_tol * scale) {
        .gl_joint_stop("joint_penalty_invalid", paste("Smooth penalty is not symmetric positive semidefinite:", context))
      }
      smooth[[pn]][[term]] <- b[expected]
      smooth_rows[[length(smooth_rows) + 1L]] <- data.frame(
        name = expected, parameter = pn, kind = "smooth", term = term,
        basis_index = seq_len(ncol(B)), stringsAsFactors = FALSE)
    }
  }
  if (!.gl_joint_valid_names(expected_fixed) || !setequal(expected_fixed, names(fixed))) {
    .gl_joint_stop("joint_fixed_map_invalid", "Fixed design columns do not exactly match fitted coefficient names.")
  }
  augmented <- .gl_build_cg_model(mm, fixed, smooth)
  q <- augmented$beta
  if (!.gl_joint_valid_names(names(q))) {
    .gl_joint_stop("joint_duplicate_coordinate", "Joint coefficient names are not unique.")
  }
  # calc_eta currently selects groups by an unanchored name search. Refuse an
  # ambiguous existing name rather than silently evaluating the wrong group.
  for (pn in parameters) {
    selected <- startsWith(names(q), paste0(pn, "."))
    if (!identical(selected, grepl(pn, names(q), fixed = TRUE))) {
      .gl_joint_stop("joint_ambiguous_parameter_names", "A coefficient label overlaps a different parameter name.")
    }
    bare <- sub(paste0("^", pn, "\\."), "", names(q)[selected])
    augmented$mm$x[[pn]] <- augmented$mm$x[[pn]][, bare, drop = FALSE]
  }
  P <- .gl_build_cg_penalty(names(q), object$lambda_s, smooth, mm)
  if (!identical(dim(P), c(length(q), length(q))) || any(!is.finite(P))) {
    .gl_joint_stop("joint_penalty_invalid", "The assembled fitted penalty is nonfinite or nonconformable.")
  }
  map <- do.call(rbind, c(fixed_rows, smooth_rows))
  map <- map[match(names(q), map$name), , drop = FALSE]
  rownames(map) <- NULL
  map$index <- seq_along(q)
  evaluation <- object
  evaluation$par <- q
  evaluation$model_matrix <- augmented$mm
  evaluation$par_s <- setNames(rep(list(list()), length(parameters)), parameters)
  # Separate structured evaluation is ordered by the fixed design, not q order.
  structured <- object
  structured$model_matrix <- mm
  structured$par <- fixed[expected_fixed]
  structured$par_s <- smooth
  list(object = evaluation, structured = structured, q = q, penalty = P,
       map = map, fixed_names = names(fixed), smooth = smooth)
}

# New-data smooth contrasts are mapped through the fitted constrained basis.
# No call to smoothCon is made on new data and no legacy basis is guessed.
.gl_joint_smooth_prediction_matrix <- function(object, parameter, term, newdata) {
  B <- object$model_matrix$s[[parameter]][[term]]
  if (is.null(B)) .gl_joint_stop("joint_smooth_missing", "Requested fitted smooth is absent.")
  spec <- .gl_joint_smooth_spec(B, paste(parameter, term, sep = ":"))
  out <- mgcv::PredictMat(spec, newdata)
  if (!is.matrix(out) || ncol(out) != ncol(B) || any(!is.finite(out))) {
    .gl_joint_stop("joint_prediction_design_invalid", "Fitted-basis prediction returned an invalid matrix.")
  }
  colnames(out) <- paste(parameter, term, seq_len(ncol(B)), sep = ".")
  out
}
