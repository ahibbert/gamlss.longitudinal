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

  start_mode = c("default", "truth_adjacent")

) {

  start_mode <- match.arg(start_mode)

  method_value <- as.character(case$method)

  pick_method_control <- function(default, override) {

    if (is.null(override)) return(default)

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

  start_mode = "default"

) {

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


#' @keywords internal

#' @noRd

.coverage_margin_class <- function(method, gap_pct) {

  if (!is.finite(gap_pct)) return("missing")

  if (method %in% c("gamlss", "gamlss2")) return("reference")

  if (identical(method, "rs_separate")) {

    if (gap_pct <= 1) return("excellent")

    if (gap_pct <= 3) return("acceptable")

    return("review")

  }

  if (gap_pct <= 3) return("excellent")

  if (gap_pct <= 7.5) return("acceptable")

  "review"

}


#' @keywords internal

#' @noRd

.coverage_joint_class <- function(method, delta_pct) {

  if (identical(method, "gamlss2")) return("not_applicable")

  if (!is.finite(delta_pct)) return("missing")

  if (delta_pct >= -1) return("acceptable")

  "review"

}


#' @keywords internal

#' @noRd

.coverage_add_review_metrics <- function(results) {

  results$gamlss2_marginal_loglik <- NA_real_

  results$margin_gap_pct_vs_gamlss2 <- NA_real_

  results$reference_marginal_method <- NA_character_

  results$reference_marginal_loglik <- NA_real_

  results$margin_gap_pct_vs_reference <- NA_real_

  results$margin_review_class <- "missing"

  results$rs_separate_joint_loglik <- NA_real_

  results$joint_delta_vs_rs_separate <- NA_real_

  results$joint_delta_pct_vs_rs_separate <- NA_real_

  results$joint_review_class <- "missing"


  group_cols <- c("family", "copula", "design", "n_subject", "n_time", "dependence", "missingness", "start_mode")

  group_key <- interaction(results[group_cols], drop = TRUE, lex.order = TRUE)


  for (key in unique(group_key)) {

    idx <- which(group_key == key)

    group <- results[idx, , drop = FALSE]

    gamlss2_ll <- group$marginal_loglik[group$method == "gamlss2" & group$success][1]

    gamlss_ll <- group$marginal_loglik[group$method == "gamlss" & group$success][1]

    gamlss2_fit_method <- group$marginal_fit_method[group$method == "gamlss2" & group$success][1]

    gamlss_fit_method <- group$marginal_fit_method[group$method == "gamlss" & group$success][1]

    rs_joint_ll <- group$joint_loglik[group$method == "rs_separate" & group$success][1]


    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {

      results$gamlss2_marginal_loglik[idx] <- gamlss2_ll

      results$margin_gap_pct_vs_gamlss2[idx] <- 100 * (gamlss2_ll - results$marginal_loglik[idx]) / pmax(abs(gamlss2_ll), 1e-8)

    }

    reference_ll <- NA_real_

    reference_method <- NA_character_

    if (length(gamlss2_ll) == 1L && is.finite(gamlss2_ll)) {

      reference_ll <- gamlss2_ll

      reference_method <- if (length(gamlss2_fit_method) == 1L && !is.na(gamlss2_fit_method)) gamlss2_fit_method else "gamlss2"

    } else if (length(gamlss_ll) == 1L && is.finite(gamlss_ll)) {

      reference_ll <- gamlss_ll

      reference_method <- if (length(gamlss_fit_method) == 1L && !is.na(gamlss_fit_method)) gamlss_fit_method else "gamlss"

    }

    if (is.finite(reference_ll)) {

      results$reference_marginal_method[idx] <- reference_method

      results$reference_marginal_loglik[idx] <- reference_ll

      results$margin_gap_pct_vs_reference[idx] <- 100 * (reference_ll - results$marginal_loglik[idx]) / pmax(abs(reference_ll), 1e-8)

    }

    if (length(rs_joint_ll) == 1L && is.finite(rs_joint_ll)) {

      results$rs_separate_joint_loglik[idx] <- rs_joint_ll

      results$joint_delta_vs_rs_separate[idx] <- results$joint_loglik[idx] - rs_joint_ll

      results$joint_delta_pct_vs_rs_separate[idx] <- 100 * results$joint_delta_vs_rs_separate[idx] / pmax(abs(rs_joint_ll), 1e-8)

    }

  }


  results$margin_review_class <- vapply(seq_len(nrow(results)), function(i) {

    .coverage_margin_class(results$method[[i]], results$margin_gap_pct_vs_reference[[i]])

  }, character(1))

  results$joint_review_class <- vapply(seq_len(nrow(results)), function(i) {

    .coverage_joint_class(results$method[[i]], results$joint_delta_pct_vs_rs_separate[[i]])

  }, character(1))

  results

}


#' @keywords internal

#' @noRd

.coverage_runtime_summary <- function(results) {

  stats::aggregate(

    elapsed_sec ~ family + copula + method + design + success + failure_type,

    data = results,

    FUN = function(x) c(n = length(x), median = stats::median(x), max = max(x))

  )

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

#' @param include_mixed Logical; include mixed-support `gamlss.dist` families in

#'   the candidate family grid.

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

  ...

) {

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
