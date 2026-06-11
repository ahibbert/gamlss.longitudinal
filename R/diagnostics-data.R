.gl_get_family_fun <- function(family_name, prefix) {

  fun_name <- paste0(prefix, family_name)

  if (exists(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE)) {

    return(get(fun_name, envir = asNamespace("gamlss.dist"), mode = "function", inherits = FALSE))

  }

  if (exists(fun_name, mode = "function")) {

    return(get(fun_name, mode = "function"))

  }

  stop(

    "Distribution function '", fun_name,

    "' is not available in gamlss.dist or the current session."

  )

}


.gl_call_family_fun <- function(prefix, family_name, x, params, extra_args = list()) {

  fun <- .gl_get_family_fun(family_name, prefix)

  arg_name <- switch(prefix, p = "q", d = "x", q = "p", "x")

  args <- c(stats::setNames(list(x), arg_name), extra_args, params)

  args <- args[names(args) %in% formalArgs(fun)]

  do.call(fun, args)

}


.gl_diag_data <- function(object) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")

  }


  copula_link <- get_copula_dist(object$copula_dist)$copula_link

  eta_out <- calc_eta(

    par_cov = object$par,

    mm = object$model_matrix,

    margin_dist = object$margin_dist,

    copula_link = copula_link,

    par_s = object$par_s

  )


  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]

  y <- object$response


  keep <- is.finite(y)

  for (par_name in names(margin_params)) {

    keep <- keep & is.finite(margin_params[[par_name]])

  }


  margin_params <- lapply(margin_params, function(x) x[keep])

  y <- y[keep]


  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]

  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(y, na.rm = TRUE), length(y))


  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)


  list(

    response = y,

    params = margin_params,

    mu_hat = as.numeric(mu_hat),

    sigma_hat = sigma_hat,

    family = object$margin_dist$family[1],

    subject = object$response_subject[keep],

    time = object$response_margin[keep]

  )

}


.gl_prepare_newdata_internal <- function(object, newdata, require_response = FALSE) {

  if (is.null(newdata)) return(NULL)


  if (is.null(object$formulas_int) || is.null(object$var_map)) {

    stop("newdata prediction requires a model fit created with stored formulas/variable map. Refit with the current package version.")

  }


  nd <- as.data.frame(newdata, stringsAsFactors = FALSE)


  # Translate user variable names to internal names used by model formulas.

  for (old_name in names(object$var_map)) {

    new_name <- object$var_map[[old_name]]

    if (old_name %in% names(nd) && !new_name %in% names(nd)) {

      names(nd)[names(nd) == old_name] <- new_name

    }

  }


  # Fitting keeps the original user time scale as `time_covariate` for

  # formulas, while `time` is used internally for ordering/pairing. Recreate

  # that column for prediction data supplied with either original or internal

  # names.

  if (!"time_covariate" %in% names(nd) && "time" %in% names(nd)) {

    nd$time_covariate <- nd$time

  }


  if (!"time" %in% names(nd) && "time" %in% names(object$model_matrix$x$mu)) {

    nd$time <- NA

  }

  if (!"subject" %in% names(nd)) {

    nd$subject <- seq_len(nrow(nd))

  }

  if (!"response" %in% names(nd)) {

    nd$response <- NA_real_

  }


  if (!is.null(object$dataset)) {

    factor_cols <- names(object$dataset)[vapply(object$dataset, is.factor, logical(1))]

    for (nm in intersect(factor_cols, names(nd))) {

      train_col <- object$dataset[[nm]]

      train_levels <- levels(train_col)

      nd_values <- as.character(nd[[nm]])

      unknown <- setdiff(unique(nd_values[!is.na(nd_values)]), train_levels)

      if (length(unknown) > 0L) {

        stop(

          "newdata column '", nm, "' contains level(s) not seen during fitting: ",

          paste(unknown, collapse = ", "),

          call. = FALSE

        )

      }

      nd[[nm]] <- factor(nd_values, levels = train_levels, ordered = is.ordered(train_col))

      if (is.ordered(train_col) && length(train_levels) > 1L) {

        contr <- contr.treatment(length(train_levels))

        colnames(contr) <- train_levels[-1]

        contrasts(nd[[nm]]) <- contr

      }

    }

  }


  if (require_response && all(is.na(nd$response))) {

    stop("newdata must include a response column (or mapped response variable) for this operation.")

  }


  nd

}


.gl_align_model_matrix_columns <- function(mm_use, mm_reference) {

  if (is.null(mm_use) || is.null(mm_reference)) {

    return(mm_use)

  }


  for (par_name in intersect(names(mm_reference$x), names(mm_use$x))) {

    ref_cols <- colnames(mm_reference$x[[par_name]])

    use_cols <- colnames(mm_use$x[[par_name]])

    missing_cols <- setdiff(ref_cols, use_cols)


    if (length(missing_cols) > 0L) {

      for (col_name in missing_cols) {

        mm_use$x[[par_name]][[col_name]] <- 0

      }

    }


    extra_cols <- setdiff(colnames(mm_use$x[[par_name]]), ref_cols)

    if (length(extra_cols) > 0L) {

      mm_use$x[[par_name]] <- mm_use$x[[par_name]][

        ,

        setdiff(colnames(mm_use$x[[par_name]]), extra_cols),

        drop = FALSE

      ]

    }


    mm_use$x[[par_name]] <- mm_use$x[[par_name]][, ref_cols, drop = FALSE]

  }


  mm_use

}


.gl_fitted_distribution <- function(object, newdata = NULL, require_response = TRUE) {

  if (!inherits(object, "gamlss.longitudinal")) {

    stop("Diagnostics are only available for objects of class 'gamlss.longitudinal'.")

  }


  copula_link <- get_copula_dist(object$copula_dist)$copula_link


  if (is.null(newdata)) {

    mm_use <- object$model_matrix

    response <- object$response

    response_margin <- object$response_margin

    response_subject <- object$response_subject

  } else {

    nd <- .gl_prepare_newdata_internal(object, newdata, require_response = require_response)


    mm_use <- do.call(

      create_model_matrices,

      list(

        mu.formula = object$formulas_int$mu,

        sigma.formula = object$formulas_int$sigma,

        nu.formula = object$formulas_int$nu,

        tau.formula = object$formulas_int$tau,

        theta.formula = object$formulas_int$theta,

        zeta.formula = object$formulas_int$zeta,

        margin.family = object$margin_dist,

        copula.family = object$copula_dist,

        copula.link = copula_link,

        dataset = nd,

        quiet_gamlss2 = TRUE,

        preserve_factor_levels = TRUE

      )

    )

    mm_use <- .gl_align_model_matrix_columns(mm_use, object$model_matrix)


    response <- nd$response

    response_margin <- nd$time

    response_subject <- nd$subject

  }


  eta_out <- calc_eta(

    par_cov = object$par,

    mm = mm_use,

    margin_dist = object$margin_dist,

    copula_link = copula_link,

    par_s = object$par_s

  )


  margin_params <- eta_out$eta_inv[names(object$margin_dist$parameters)]


  keep <- rep(TRUE, length(response))

  if (require_response) {

    keep <- is.finite(response)

  }

  for (par_name in names(margin_params)) {

    keep <- keep & is.finite(margin_params[[par_name]])

  }


  common_n <- min(

    length(response),

    length(response_margin),

    length(response_subject),

    if (length(margin_params) > 0) min(vapply(margin_params, length, integer(1))) else length(response)

  )

  if (!is.finite(common_n) || common_n < 0) {

    common_n <- 0L

  }


  response <- response[seq_len(common_n)]

  response_margin <- response_margin[seq_len(common_n)]

  response_subject <- response_subject[seq_len(common_n)]

  margin_params <- lapply(margin_params, function(x) x[seq_len(common_n)])

  keep <- keep[seq_len(common_n)]


  margin_params <- lapply(margin_params, function(x) x[keep])

  response <- response[keep]

  response_margin <- response_margin[keep]

  response_subject <- response_subject[keep]


  mu_hat <- if ("mu" %in% names(margin_params)) margin_params$mu else margin_params[[1]]

  sigma_hat <- if ("sigma" %in% names(margin_params)) margin_params$sigma else rep(stats::sd(response, na.rm = TRUE), length(response))

  sigma_hat <- pmax(as.numeric(sigma_hat), .Machine$double.eps)


  list(

    response = response,

    params = margin_params,

    mu_hat = as.numeric(mu_hat),

    sigma_hat = sigma_hat,

    family = object$margin_dist$family[1],

    subject = response_subject,

    time = response_margin,

    keep_mask = keep,

    keep_index = which(keep)

  )

}


.gl_pit <- function(object, randomize = FALSE) {

  diag_data <- .gl_fitted_distribution(object, newdata = NULL, require_response = TRUE)

  y <- diag_data$response

  params <- diag_data$params


  pit_upper <- .gl_call_family_fun("p", diag_data$family, y, params)

  pit <- pit_upper


  if (randomize && .is_discrete_margin(object$margin_dist)) {

    pit_lower <- .gl_call_family_fun("p", diag_data$family, y - 1, params)

    pit_lower <- pmin(pmax(as.numeric(pit_lower), 0), 1)

    pit_upper <- pmin(pmax(as.numeric(pit_upper), 0), 1)

    interval_width <- pmax(pit_upper - pit_lower, 0)

    pit <- pit_lower + stats::runif(length(pit_upper)) * interval_width

  }


  pit <- pmin(pmax(as.numeric(pit), 0), 1)

  list(diag = diag_data, pit = pit)

}


.gl_check_by_time <- function(by_time) {

  if (!is.logical(by_time) || length(by_time) != 1L || is.na(by_time)) {

    stop("'by_time' must be TRUE or FALSE.", call. = FALSE)

  }

}


.gl_resolve_diag_split <- function(by, diag_data, data = NULL, plot_name = "diagnostic") {

  if (is.null(by) || (is.character(by) && length(by) == 1L && !nzchar(by))) {

    return(NULL)

  }

  if (!is.character(by) || length(by) != 1L) {

    stop("'by' must be NULL or a single column name as a character string.", call. = FALSE)

  }


  if (by %in% c("time", "response_margin")) {

    return(as.factor(diag_data$time))

  }

  if (by %in% c("subject", "response_subject")) {

    return(as.factor(diag_data$subject))

  }


  if (is.null(data)) {

    stop("To split ", plot_name, " by '", by, "', provide data= containing that column.", call. = FALSE)

  }

  if (!is.data.frame(data)) {

    data <- as.data.frame(data, stringsAsFactors = FALSE)

  }

  if (!by %in% names(data)) {

    stop("Column '", by, "' not found in provided data.", call. = FALSE)

  }


  n_diag <- length(diag_data$response)

  if (nrow(data) == n_diag) {

    return(as.factor(data[[by]]))

  }


  keep_mask <- diag_data$keep_mask

  if (!is.null(keep_mask) && length(keep_mask) == nrow(data)) {

    vec <- data[[by]][keep_mask]

    if (length(vec) == n_diag) {

      return(as.factor(vec))

    }

  }


  keep_index <- diag_data$keep_index

  if (!is.null(keep_index) && length(keep_index) == n_diag && (n_diag == 0L || max(keep_index, na.rm = TRUE) <= nrow(data))) {

    return(as.factor(data[[by]][keep_index]))

  }


  stop(

    "Could not align data rows with ", plot_name, " rows for by='", by, "'. ",

    "Provide data with row count equal to either the diagnostic row count (", n_diag, ") or the original fit data rows.",

    call. = FALSE

  )

}


.gl_diag_split_info <- function(by_time, by, diag_data, data = NULL, plot_name = "diagnostic") {

  .gl_check_by_time(by_time)


  if (isTRUE(by_time) && is.null(by)) {

    by <- "time"

  } else if (isTRUE(by_time) && !is.null(by)) {

    warning("Both by_time and by were provided; using by='", by, "'.", call. = FALSE)

  }


  split_group <- .gl_resolve_diag_split(by, diag_data, data = data, plot_name = plot_name)

  list(by = by, group = split_group, split_by = !is.null(split_group))

}

