# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- coverage-runner.R ----

#' Run coverage simulation cases
#'
#' Runs a compact grid of coverage simulation cases and optionally writes CSV,
#' RDS, and LaTeX summary outputs for reviewer inspection.
#'
#' @param families Character vector of margin families to include. When `NULL`,
#'   the default coverage families are used.
#' @param copulas Character vector of copula families to include.
#' @param methods Character vector of fitting methods to include.
#' @param designs Character vector of simulation designs to include.
#' @param include_mixed Include mixed discrete/continuous cases when `TRUE`.

#' @param output_dir Directory for CSV/RDS outputs.

#' @param write_results Write result files when `TRUE`.

#' @param write_summary Write a LaTeX summary report when `TRUE` and

#'   `write_results` is also `TRUE`.

#' @param compile_summary_pdf Compile the LaTeX summary to PDF when `TRUE` and

#'   a LaTeX installation is available.

#' @param smooth_results_file Optional smooth smoke-test CSV to include in the

#'   generated summary report.

#' @param report_title Title for the generated summary report.

#' @param ... Passed to the grid runner, e.g. `n`, `times`, `max_outer_iter`,

#'   `dependence`, `missingness`, and `start_mode`.

#'

#' @return A data frame of per-fit results. The return value also carries

#'   `"parameter_results"` and `"runtime_summary"` attributes.

#' @export

run_coverage_simulations <- function(
    families = NULL,
    copulas = c("N", "C", "F", "G", "J", "t"),
    methods = .coverage_default_methods(),
    designs = "intercept",
    include_mixed = FALSE,
    output_dir = file.path("results", "coverage_simulations"),
    write_results = TRUE,
    write_summary = TRUE,
    compile_summary_pdf = FALSE,
    smooth_results_file = NULL,
    report_title = "Coverage Simulation Summary",
    ...) {
  grid <- .coverage_make_case_grid(
    families = families,
    copulas = copulas,
    methods = methods,
    designs = designs,
    include_mixed = include_mixed
  )

  results <- .coverage_run_grid(grid, ...)

  if (isTRUE(write_results)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

    parameter_results <- attr(results, "parameter_results")

    runtime_summary <- attr(results, "runtime_summary")

    results_file <- file.path(output_dir, paste0("coverage_results_", stamp, ".csv"))

    parameter_results_file <- file.path(output_dir, paste0("coverage_parameter_results_", stamp, ".csv"))

    runtime_summary_file <- file.path(output_dir, paste0("coverage_runtime_summary_", stamp, ".csv"))

    summary_file <- file.path(output_dir, paste0("coverage_summary_", stamp, ".tex"))

    utils::write.csv(results, results_file, row.names = FALSE)

    saveRDS(results, file.path(output_dir, paste0("coverage_results_", stamp, ".rds")))

    if (!is.null(parameter_results) && nrow(parameter_results) > 0L) {
      utils::write.csv(

        parameter_results,
        parameter_results_file,
        row.names = FALSE
      )
    }

    if (!is.null(runtime_summary) && nrow(runtime_summary) > 0L) {
      utils::write.csv(

        runtime_summary,
        runtime_summary_file,
        row.names = FALSE
      )
    }

    if (isTRUE(write_summary)) {
      write_coverage_summary_report(
        results = results,
        parameter_results = parameter_results,
        smooth_results_file = smooth_results_file,
        results_file = results_file,
        parameter_results_file = if (file.exists(parameter_results_file)) parameter_results_file else NULL,
        output_tex = summary_file,
        compile_pdf = compile_summary_pdf,
        title = report_title,
        run_label = paste0("coverage_results_", stamp)
      )
    }
  }

  results
}

# ---- coverage-runner-grid.R ----

#' @keywords internal

#' @noRd

.coverage_run_grid <- function(
    grid,
    n = 80,
    times = 1:3,
    seed = 1,
    max_outer_iter = 8,
    max_inner_iter = 8,
    max_elapsed_sec = 20,
    method_max_outer_iter = NULL,
    method_max_inner_iter = NULL,
    method_max_elapsed_sec = NULL,
    dependence = "moderate",
    missingness = "none",
    start_mode = "default") {
  row_bind_fill <- function(rows) {
    rows <- rows[!vapply(rows, is.null, logical(1))]

    if (length(rows) == 0L) {
      return(data.frame())
    }

    all_names <- unique(unlist(lapply(rows, names), use.names = FALSE))

    rows <- lapply(rows, function(x) {
      missing <- setdiff(all_names, names(x))

      for (nm in missing) {
        x[[nm]] <- NA
      }

      x[all_names]
    })

    do.call(rbind, rows)
  }

  rows <- vector("list", nrow(grid))

  parameter_rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    case_result <- .coverage_run_case(

      grid[i, , drop = FALSE],
      n = n,
      times = times,
      seed = seed + i,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      method_max_outer_iter = method_max_outer_iter,
      method_max_inner_iter = method_max_inner_iter,
      method_max_elapsed_sec = method_max_elapsed_sec,
      dependence = dependence,
      missingness = missingness,
      start_mode = start_mode
    )

    parameter_rows[[i]] <- attr(case_result, "parameter_results")

    rows[[i]] <- case_result
  }

  results <- .coverage_add_review_metrics(row_bind_fill(rows))

  parameter_rows <- parameter_rows[!vapply(parameter_rows, is.null, logical(1))]

  attr(results, "parameter_results") <- if (length(parameter_rows) == 0L) {
    data.frame()
  } else {
    row_bind_fill(parameter_rows)
  }

  attr(results, "runtime_summary") <- .coverage_runtime_summary(results)

  results
}


#' Run opt-in distribution/copula/method coverage simulations

#'

#' @param families Optional family names. Defaults to all supported non-mixed

#'   `gamlss.dist` families with `q`, `p`, and `d` functions.

#' @param copulas Copula codes.

#' @param methods Fit methods. Defaults to `"gamlss"`, `"rs_separate"`,

#'   `"rs_joint"`, and `"cg"`. Method `"gamlss2"` is also supported when the

#'   optional non-CRAN package is installed and explicitly requested. Standard

#'   comparator methods `"gee"`, `"glmm"`, `"gam"`, and `"glmmTMB"` are

#'   available for families that map to common mean-model families.

#' @param designs Simulation designs: `"intercept"`, `"covariate"`,

#'   `"scale"`, `"time_dependence"`, or `"smooth"`.

# ---- coverage-runner-case.R ----

.coverage_run_case <- function(
    case,
    n = 80,
    times = 1:3,
    seed = 1,
    max_outer_iter = 8,
    max_inner_iter = 8,
    max_elapsed_sec = 20,
    method_max_outer_iter = NULL,
    method_max_inner_iter = NULL,
    method_max_elapsed_sec = NULL,
    dependence = "moderate",
    missingness = "none",
    start_mode = c("default", "truth_adjacent")) {
  start_mode <- match.arg(start_mode)

  method_value <- as.character(case$method)

  pick_method_control <- function(default, override) {
    if (is.null(override)) {
      return(default)
    }

    if (is.list(override)) {
      value <- override[[method_value]]
    } else {
      value <- override[[method_value]]
    }

    if (is.null(value) || length(value) == 0L || is.na(value[[1L]])) default else value[[1L]]
  }

  max_outer_iter <- pick_method_control(max_outer_iter, method_max_outer_iter)

  max_inner_iter <- pick_method_control(max_inner_iter, method_max_inner_iter)

  max_elapsed_sec <- pick_method_control(max_elapsed_sec, method_max_elapsed_sec)

  dat <- .coverage_simulate_case(
    family = case$family,
    copula = case$copula,
    design = case$design,
    n = n,
    times = times,
    seed = seed,
    dependence = dependence
  )

  dat <- .coverage_apply_missingness(dat, missingness = missingness)

  start_from <- if (identical(start_mode, "truth_adjacent") && !case$method %in% c("gamlss", "gamlss2")) {
    .coverage_truth_adjacent_start(dat, case$family, case$copula, case$design)
  } else {
    NA
  }

  fit_row <- if (identical(case$method, "gamlss")) {
    .coverage_fit_gamlss(

      dat,
      case$family,
      case$copula,
      case$design,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec
    )
  } else if (identical(case$method, "gamlss2")) {
    .coverage_fit_gamlss2(dat, case$family, case$copula, case$design, max_elapsed_sec = max_elapsed_sec)
  } else if (case$method %in% c("gee", "glmm", "gam", "glmmTMB")) {
    .coverage_fit_standard_comparator(

      dat,
      family = case$family,
      copula = case$copula,
      design = case$design,
      method = case$method,
      max_elapsed_sec = max_elapsed_sec
    )
  } else {
    .coverage_fit_longitudinal(

      dat,
      family = case$family,
      copula = case$copula,
      method = case$method,
      design = case$design,
      max_outer_iter = max_outer_iter,
      max_inner_iter = max_inner_iter,
      max_elapsed_sec = max_elapsed_sec,
      start_from = start_from
    )
  }

  context <- data.frame(
    case_id = case$case_id %||% NA_integer_,
    family = case$family,
    copula = case$copula,
    design = case$design,
    n_subject = n,
    n_time = length(times),
    dependence = dependence,
    missingness = missingness,
    start_mode = start_mode,
    stringsAsFactors = FALSE
  )

  out <- cbind(

    context,
    fit_row,
    row.names = NULL
  )

  parameter_results <- attr(fit_row, "parameter_results")

  if (!is.null(parameter_results) && nrow(parameter_results) > 0L) {
    attr(out, "parameter_results") <- cbind(

      context[rep(1L, nrow(parameter_results)), , drop = FALSE],
      data.frame(
        method = fit_row$method,
        success = fit_row$success,
        elapsed_sec = fit_row$elapsed_sec,
        stringsAsFactors = FALSE
      ),
      parameter_results,
      row.names = NULL
    )
  }

  out
}

# ---- coverage-simulate-case.R ----

.coverage_simulate_case <- function(
    family,
    copula,
    design = .coverage_supported_designs(),
    n = 80,
    times = 1:3,
    seed = 1,
    dependence = "moderate") {
  design <- match.arg(design)

  margin_dist <- do.call(get(family, envir = asNamespace("gamlss.dist")), list())

  margin_params <- .coverage_default_margin_params(margin_dist)

  covariates <- NULL


  if (design %in% c("covariate", "scale", "smooth")) {
    covariates <- function(base) {
      x <- sim_rescale01(as.numeric(base$.sim_subject_index))

      data.frame(x = x - mean(x), stringsAsFactors = FALSE)
    }

    if (identical(design, "covariate") && "mu" %in% names(margin_params)) {
      base_mu <- margin_params$mu

      margin_params$mu <- function(data) {
        x <- sim_rescale01(as.numeric(data$.sim_subject_index))

        x <- x - mean(x)

        if (base_mu > 0) {
          pmax(base_mu * exp(0.15 * x), .Machine$double.eps)
        } else {
          base_mu + 0.15 * x
        }
      }
    }

    if (identical(design, "smooth")) {
      for (par_name in names(margin_params)) {
        linkfun <- margin_dist[[paste0(par_name, ".linkfun")]]

        linkinv <- margin_dist[[paste0(par_name, ".linkinv")]]

        if (!is.function(linkfun) || !is.function(linkinv)) next

        margin_params[[par_name]] <- .coverage_make_smooth_param(

          linkfun,
          linkinv,
          margin_params[[par_name]],
          amplitude = if (identical(par_name, "mu")) 0.22 else 0.14
        )
      }
    }

    if (identical(design, "scale") && "sigma" %in% names(margin_params)) {
      base_sigma <- margin_params$sigma

      margin_params$sigma <- function(data) {
        x <- sim_rescale01(as.numeric(data$.sim_subject_index))

        x <- x - mean(x)

        pmax(base_sigma * exp(0.45 * x), .Machine$double.eps)
      }
    }
  }


  copula_params <- if (identical(design, "time_dependence")) {
    .coverage_time_varying_copula_params(copula)
  } else if (identical(design, "smooth")) {
    copula_link <- get_copula_dist(copula)$copula_link

    base <- .coverage_copula_params(copula, dependence = dependence)

    if ("tau" %in% names(base)) {
      list(tau = function(edge_data) {
        eta <- stats::qlogis(base$tau) + .coverage_smooth_eta_component(edge_data, amplitude = 0.12)

        stats::plogis(eta)
      })
    } else {
      base$theta <- .coverage_make_smooth_param(copula_link$theta.linkfun, copula_link$theta.linkinv, base$theta, amplitude = 0.12)

      base
    }
  } else {
    .coverage_copula_params(copula, dependence = dependence)
  }


  simulate_longitudinal_dataset(
    n = n,
    times = times,
    margin_dist = margin_dist,
    copula_dist = copula,
    margin_params = margin_params,
    copula_params = copula_params,
    covariates = covariates,
    seed = seed,
    include_truth = TRUE,
    u_bounds = .coverage_simulation_u_bounds(family)
  )
}


#' @keywords internal

#' @noRd

.coverage_apply_missingness <- function(dat, missingness = c("none", "mcar", "drop_rows"), prop = 0.05) {
  missingness <- match.arg(missingness)

  if (identical(missingness, "none")) {
    return(dat)
  }


  n_subject <- length(unique(dat$subject))

  n_time <- length(unique(dat$time))

  n_target <- max(1L, floor(nrow(dat) * prop))


  if (identical(missingness, "mcar")) {
    eligible <- which(dat$time != min(dat$time))

    idx <- eligible[seq_len(min(length(eligible), n_target))]

    dat$response[idx] <- NA_real_

    return(dat)
  }


  subject_index <- as.integer(dat$subject)

  time_values <- sort(unique(dat$time))

  drop_idx <- subject_index <= max(1L, floor(n_subject * prop)) & dat$time == time_values[min(2L, n_time)]

  dat[!drop_idx, , drop = FALSE]
}


#' @keywords internal

#' @noRd

.coverage_fit_formulas <- function(design) {
  if (identical(design, "covariate")) {
    list(mu = response ~ x, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  } else if (identical(design, "scale")) {
    list(mu = response ~ 1, sigma = ~x, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  } else if (identical(design, "time_dependence")) {
    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~time, zeta = ~1)
  } else if (identical(design, "smooth")) {
    list(
      mu = response ~ s(x, bs = "ps", k = 6),
      sigma = ~ s(x, bs = "ps", k = 6),
      nu = ~ s(x, bs = "ps", k = 6),
      tau = ~ s(x, bs = "ps", k = 6),
      theta = ~ s(x, bs = "ps", k = 6),
      zeta = ~ s(x, bs = "ps", k = 6)
    )
  } else {
    list(mu = response ~ 1, sigma = ~1, nu = ~1, tau = ~1, theta = ~1, zeta = ~1)
  }
}


#' @keywords internal

#' @noRd
