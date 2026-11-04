# Consolidated support module for reviewer navigation.
# Source blocks are grouped mechanically from the files named below.

# ---- coverage_report.R ----

#' Write a LaTeX summary report for coverage simulations
#'
#' @param prefix Optional coverage-run file prefix. When supplied, the function
#'   looks for `paste0(prefix, "_final_results.csv")` and
#'   `paste0(prefix, "_final_parameter_results.csv")`.
#' @param results Optional per-fit coverage results data frame.
#' @param parameter_results Optional eta-scale parameter results data frame.
#' @param smooth_results Optional smooth smoke-test results data frame.
#' @param results_file Optional CSV file containing per-fit coverage results.
#' @param parameter_results_file Optional CSV file containing eta-scale
#'   parameter results.
#' @param smooth_results_file Optional CSV file containing smooth smoke-test
#'   results.
#' @param output_tex Output `.tex` file. Defaults to
#'   `paste0(prefix, "_summary.tex")` when `prefix` is supplied.
#' @param compile_pdf Compile the report to PDF when a LaTeX installation is
#'   available.
#' @param title Report title.
#' @param run_label Optional run label used in the report date/scope text.
#' @param top_n Number of rows for "slowest" and "largest error" tables.
#'
#' @return Invisibly returns a list with `tex_file`, `pdf_file`, and `compiled`.
#' @export
write_coverage_summary_report <- function(
    prefix = NULL,
    results = NULL,
    parameter_results = NULL,
    smooth_results = NULL,
    results_file = NULL,
    parameter_results_file = NULL,
    smooth_results_file = NULL,
    output_tex = NULL,
    compile_pdf = FALSE,
    title = "Coverage Simulation Summary",
    run_label = NULL,
    top_n = 10) {
  inputs <- .coverage_report_prepare_inputs(
    prefix = prefix,
    results = results,
    parameter_results = parameter_results,
    smooth_results = smooth_results,
    results_file = results_file,
    parameter_results_file = parameter_results_file,
    smooth_results_file = smooth_results_file,
    output_tex = output_tex
  )
  prefix <- inputs$prefix
  results <- inputs$results
  parameter_results <- inputs$parameter_results
  smooth_results <- inputs$smooth_results
  results_file <- inputs$results_file
  parameter_results_file <- inputs$parameter_results_file
  output_tex <- inputs$output_tex

  dir.create(dirname(output_tex), recursive = TRUE, showWarnings = FALSE)

  tables <- .coverage_report_table_bundle(
    results = results,
    parameter_results = parameter_results,
    smooth_results = smooth_results,
    top_n = top_n
  )
  lines <- .coverage_report_document_lines(
    results = results,
    tables = tables,
    title = title,
    run_label = run_label,
    prefix = prefix,
    results_file = results_file,
    parameter_results_file = parameter_results_file
  )

  writeLines(lines, output_tex, useBytes = TRUE)

  pdf_file <- sub("\\.tex$", ".pdf", output_tex, ignore.case = TRUE)
  compiled <- FALSE
  if (isTRUE(compile_pdf)) {
    compiled <- .coverage_compile_latex(output_tex, pdf_file)
  }
  invisible(list(
    tex_file = normalizePath(output_tex, winslash = "/", mustWork = FALSE),
    pdf_file = normalizePath(pdf_file, winslash = "/", mustWork = FALSE),
    compiled = compiled
  ))
}

# ---- coverage-report-workflow.R ----

#' Resolve coverage report file inputs and defaults
#'
#' @noRd
.coverage_report_prepare_inputs <- function(
    prefix,
    results,
    parameter_results,
    smooth_results,
    results_file,
    parameter_results_file,
    smooth_results_file,
    output_tex) {
  if (!is.null(prefix)) {
    if (is.null(results_file)) {
      results_file <- paste0(prefix, "_final_results.csv")
    }
    if (is.null(parameter_results_file)) {
      parameter_results_file <- paste0(prefix, "_final_parameter_results.csv")
    }
    if (is.null(output_tex)) {
      output_tex <- paste0(prefix, "_summary.tex")
    }
  }
  if (is.null(results) && !is.null(results_file) && file.exists(results_file)) {
    results <- utils::read.csv(results_file, stringsAsFactors = FALSE)
  }
  if (is.null(parameter_results) &&
    !is.null(parameter_results_file) &&
    file.exists(parameter_results_file)) {
    parameter_results <- utils::read.csv(parameter_results_file, stringsAsFactors = FALSE)
  }
  if (is.null(smooth_results) &&
    !is.null(smooth_results_file) &&
    file.exists(smooth_results_file)) {
    smooth_results <- utils::read.csv(smooth_results_file, stringsAsFactors = FALSE)
  }
  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0L) {
    stop("Coverage results are required.", call. = FALSE)
  }
  if (is.null(output_tex)) {
    output_tex <- file.path(getwd(), paste0("coverage_summary_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tex"))
  }

  list(
    prefix = prefix,
    results = results,
    parameter_results = parameter_results,
    smooth_results = smooth_results,
    results_file = results_file,
    parameter_results_file = parameter_results_file,
    smooth_results_file = smooth_results_file,
    output_tex = output_tex
  )
}

#' Build coverage report table bundle
#'
#' @noRd
.coverage_report_table_bundle <- function(results, parameter_results, smooth_results, top_n) {
  fit_tables <- .coverage_report_fit_tables(results)
  runtime_tables <- .coverage_report_runtime_tables(results, top_n = top_n)
  eta_tables <- .coverage_report_eta_tables(parameter_results, top_n = top_n)
  likelihood_tables <- .coverage_report_likelihood_tables(results)

  c(
    fit_tables,
    runtime_tables,
    eta_tables,
    likelihood_tables,
    list(smooth_df = .coverage_report_smooth_table(smooth_results))
  )
}

# ---- coverage-report-document.R ----

#' Build the LaTeX lines for a coverage summary report
#'
#' @noRd
.coverage_report_document_lines <- function(
    results,
    tables,
    title,
    run_label,
    prefix,
    results_file,
    parameter_results_file) {
  copulas <- paste(sort(unique(results$copula)), collapse = ", ")
  methods <- paste(sort(unique(results$method)), collapse = ", ")
  designs <- paste(sort(unique(results$design)), collapse = ", ")
  n_families <- length(unique(results$family))
  n_cases <- nrow(results)
  total_success <- sum(results$success)
  total_failed <- nrow(results) - total_success
  converged <- if ("converged" %in% names(results)) as.logical(results$converged) else rep(NA, nrow(results))
  total_converged <- sum(results$success & !is.na(converged) & converged)
  total_nonconverged <- sum(results$success & !is.na(converged) & !converged)
  total_convergence_unknown <- sum(results$success & is.na(converged))
  wall_note <- "Wall time is not inferred from per-fit timings; per-fit elapsed times are summed by method below."
  label_text <- if (!is.null(run_label) && length(run_label) > 0L && !is.na(run_label[[1L]])) {
    run_label[[1L]]
  } else if (!is.null(prefix)) {
    basename(prefix)
  } else {
    format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  }

  lines <- c(
    "\\documentclass[11pt]{article}",
    "",
    "\\usepackage[margin=1in]{geometry}",
    "\\usepackage{booktabs}",
    "\\usepackage{array}",
    "\\usepackage{longtable}",
    "\\usepackage{hyperref}",
    "",
    paste0("\\title{", .coverage_tex_escape(title), "}"),
    "\\author{\\texttt{gamlss.longitudinal}}",
    paste0("\\date{", .coverage_tex_escape(label_text), "}"),
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\section*{Scope}",
    paste0(
      "This report was generated automatically from coverage simulation output. It includes ",
      n_families, " margin families, copula(s) ", .coverage_tex_texttt(copulas), ", design(s) ",
      .coverage_tex_texttt(designs), ", and method(s) ", .coverage_tex_texttt(methods), "."
    ),
    "",
    paste0(
      "Total fitted cases: ", n_cases,
      ". Completed without an execution error: ", total_success,
      ". Optimizer-converged: ", total_converged,
      ". Completed but not converged: ", total_nonconverged,
      ". Convergence status unavailable: ", total_convergence_unknown,
      ". Execution failures: ", total_failed, "."
    ),
    "",
    if (!is.null(results_file)) paste0("Results file: ", .coverage_tex_texttt(basename(results_file)), ".") else NULL,
    if (!is.null(parameter_results_file) && file.exists(parameter_results_file)) {
      paste0("Parameter recovery file: ", .coverage_tex_texttt(basename(parameter_results_file)), ".")
    } else {
      NULL
    },
    "",
    "\\section*{Execution Completion}",
    "A completed fit is not necessarily optimizer-converged; inspect the structured convergence output separately.",
    .coverage_tex_table(tables$fit_df, "Execution completion by method"),
    .coverage_tex_table(tables$by_copula, "Execution completion by copula and method"),
    "",
    "\\section*{Runtime}",
    wall_note,
    .coverage_tex_table(tables$runtime_df, "Runtime by method"),
    .coverage_tex_table(tables$family_runtime_df, "Slowest family/copula combinations by summed elapsed time"),
    .coverage_tex_table(tables$slow_case_df, "Slowest individual fits"),
    "",
    "\\section*{Parameter Recovery}",
    "Recovery is summarised on the eta scale for available parameter estimates.",
    if (!is.null(tables$eta_summary_df)) .coverage_tex_table(tables$eta_summary_df, "Eta-scale absolute error summary") else "No eta-scale parameter recovery rows were available.",
    if (!is.null(tables$eta_class_df)) .coverage_tex_table(tables$eta_class_df, "Eta-scale recovery classification") else NULL,
    if (!is.null(tables$eta_worst_df)) .coverage_tex_table(tables$eta_worst_df, "Largest non-reference eta-scale errors") else NULL,
    "",
    "\\section*{Likelihood Review}",
    "The percentage likelihood-gap review columns are useful screening diagnostics, but can be unstable when reference log-likelihoods are close to zero or positive. Prefer a per-observation absolute log-likelihood gap for future multi-copula review thresholds.",
    if (!is.null(tables$margin_df)) .coverage_tex_table(tables$margin_df, "Marginal likelihood review class versus reference") else NULL,
    if (!is.null(tables$joint_df)) .coverage_tex_table(tables$joint_df, "Joint likelihood review class versus rs\\_separate") else NULL,
    "",
    if (!is.null(tables$smooth_df)) "\\section*{Smooth-Term Smoke Layer}" else NULL,
    if (!is.null(tables$smooth_df)) "Smooth recovery rows are summarised when a smooth smoke-test result file is supplied." else NULL,
    if (!is.null(tables$smooth_df)) .coverage_tex_table(tables$smooth_df, "Smooth subset success and eta recovery") else NULL,
    "",
    "\\section*{Interpretation}",
    "A clean fit-success table means the supported family/copula/method combinations completed with finite likelihoods under the recorded simulation design. Runtime and recovery tables should be used to identify targeted optimisation or larger replicated recovery studies.",
    "",
    "\\appendix",
    "\\section{Detailed Review Rows}",
    "This appendix lists all rows requiring manual review under the current coverage thresholds.",
    "\\subsection{Parameter Recovery Concerns}",
    "Rows below are all eta-scale parameter recovery entries classified as \\texttt{concern}.",
    if (!is.null(tables$eta_concern_df)) .coverage_tex_longtable(tables$eta_concern_df, "All parameter recovery concern rows") else "No parameter recovery concern rows were recorded.",
    "\\subsection{Likelihood Review Rows}",
    "Rows below are all fits where either marginal or joint likelihood review class is \\texttt{review}.",
    if (!is.null(tables$likelihood_review_df)) .coverage_tex_longtable(tables$likelihood_review_df, "All likelihood review rows", font = "\\tiny") else "No likelihood review rows were recorded.",
    "",
    "\\end{document}"
  )
  unlist(lines, use.names = FALSE)
}

# ---- coverage-report-tex.R ----

#' Escape values for LaTeX coverage reports
#'
#' @noRd
.coverage_tex_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  repl <- c(
    "\\" = "\\textbackslash{}",
    "{" = "\\{",
    "}" = "\\}",
    "$" = "\\$",
    "&" = "\\&",
    "#" = "\\#",
    "_" = "\\_",
    "%" = "\\%",
    "~" = "\\textasciitilde{}",
    "^" = "\\textasciicircum{}"
  )
  for (pattern in names(repl)) {
    x <- gsub(pattern, repl[[pattern]], x, fixed = TRUE)
  }
  x
}

#' Wrap values in LaTeX texttt after escaping
#'
#' @noRd
.coverage_tex_texttt <- function(x) {
  paste0("\\texttt{", .coverage_tex_escape(x), "}")
}

#' Format numeric coverage report values
#'
#' @noRd
.coverage_format_number <- function(x, digits = 2) {
  out <- ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "--")
  as.character(out)
}

#' Build a compact LaTeX table
#'
#' @noRd
.coverage_tex_table <- function(df, caption, align = NULL) {
  if (is.null(df) || nrow(df) == 0L) {
    return(c(paste0("\\paragraph{", .coverage_tex_escape(caption), "} No rows available.")))
  }
  if (is.null(align)) align <- paste0("l", paste(rep("r", max(0L, ncol(df) - 1L)), collapse = ""))
  header <- paste(.coverage_tex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1L, function(row) paste(.coverage_tex_escape(row), collapse = " & "))
  c(
    "\\begin{table}[htbp]",
    "\\centering",
    paste0("\\caption{", .coverage_tex_escape(caption), "}"),
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "\\end{table}"
  )
}

#' Build a LaTeX longtable for detailed coverage review rows
#'
#' @noRd
.coverage_tex_longtable <- function(df, caption, align = NULL, font = "\\scriptsize") {
  if (is.null(df) || nrow(df) == 0L) {
    return(c(paste0("\\paragraph{", .coverage_tex_escape(caption), "} No rows available.")))
  }
  if (is.null(align)) align <- paste(rep("l", ncol(df)), collapse = "")
  header <- paste(.coverage_tex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1L, function(row) paste(.coverage_tex_escape(row), collapse = " & "))
  c(
    "\\begingroup",
    "\\setlength{\\tabcolsep}{2pt}",
    "\\begin{center}",
    font,
    paste0("\\begin{longtable}{", align, "}"),
    paste0("\\caption{", .coverage_tex_escape(caption), "}\\\\"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    "\\endfirsthead",
    paste0("\\multicolumn{", ncol(df), "}{l}{\\scriptsize\\emph{", .coverage_tex_escape(caption), " continued}}\\\\"),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    "\\endhead",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{longtable}",
    "\\end{center}",
    "\\endgroup"
  )
}

# ---- coverage-report-latex.R ----

#' @keywords internal
#' @noRd
.coverage_compile_latex <- function(tex_file, pdf_file = sub("\\.tex$", ".pdf", tex_file, ignore.case = TRUE)) {
  tex_file <- normalizePath(tex_file, winslash = "/", mustWork = TRUE)
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(dirname(tex_file))
  tex_name <- basename(tex_file)

  ok <- FALSE
  if (requireNamespace("tinytex", quietly = TRUE)) {
    ok <- tryCatch(
      {
        tinytex::latexmk(tex_name, clean = FALSE)
        file.exists(basename(pdf_file))
      },
      error = function(e) FALSE
    )
  }
  if (!ok) {
    engine <- Sys.which("pdflatex")
    if (nzchar(engine)) {
      ok <- tryCatch(
        {
          system2(engine, c("-interaction=nonstopmode", "-halt-on-error", tex_name), stdout = TRUE, stderr = TRUE)
          system2(engine, c("-interaction=nonstopmode", "-halt-on-error", tex_name), stdout = TRUE, stderr = TRUE)
          file.exists(basename(pdf_file))
        },
        error = function(e) FALSE
      )
    }
  }
  isTRUE(ok)
}

# ---- coverage-report-eta-tables.R ----

#' Build coverage report eta recovery tables
#'
#' @noRd
.coverage_report_eta_tables <- function(parameter_results, top_n) {
  eta_summary_df <- NULL
  eta_class_df <- NULL
  eta_worst_df <- NULL
  eta_concern_df <- NULL

  if (!is.null(parameter_results) && nrow(parameter_results) > 0L && "abs_eta_error" %in% names(parameter_results)) {
    eta <- parameter_results[is.finite(parameter_results$abs_eta_error), , drop = FALSE]
    if (nrow(eta) > 0L) {
      eta_summary <- do.call(rbind, lapply(split(eta$abs_eta_error, eta$method), function(v) {
        data.frame(
          N = length(v),
          Median = stats::median(v),
          P75 = as.numeric(stats::quantile(v, 0.75)),
          P90 = as.numeric(stats::quantile(v, 0.9)),
          P95 = as.numeric(stats::quantile(v, 0.95)),
          Max = max(v),
          stringsAsFactors = FALSE
        )
      }))
      eta_summary$Method <- rownames(eta_summary)
      rownames(eta_summary) <- NULL
      eta_summary_df <- data.frame(
        Method = eta_summary$Method,
        N = eta_summary$N,
        Median = .coverage_format_number(eta_summary$Median, 3),
        `75th pct.` = .coverage_format_number(eta_summary$P75, 3),
        `90th pct.` = .coverage_format_number(eta_summary$P90, 3),
        `95th pct.` = .coverage_format_number(eta_summary$P95, 3),
        Max = .coverage_format_number(eta_summary$Max, 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      eta_summary_df <- eta_summary_df[order(eta_summary_df$Method), , drop = FALSE]

      if ("eta_error_class" %in% names(eta)) {
        class_tab <- as.data.frame.matrix(table(eta$method, eta$eta_error_class))
        eta_class_df <- data.frame(Method = rownames(class_tab), class_tab, row.names = NULL, check.names = FALSE)
      }
      eta_worst <- eta[!eta$method %in% c("gamlss2", "gamlss"), , drop = FALSE]
      eta_worst <- head(eta_worst[order(-eta_worst$abs_eta_error), , drop = FALSE], top_n)
      eta_worst_df <- data.frame(
        Family = eta_worst$family,
        Copula = eta_worst$copula,
        Method = eta_worst$method,
        Parameter = eta_worst$parameter,
        `Estimate eta` = .coverage_format_number(eta_worst$estimate_eta, 3),
        `True eta` = .coverage_format_number(eta_worst$true_eta, 3),
        `Abs. error` = .coverage_format_number(eta_worst$abs_eta_error, 3),
        Class = eta_worst$eta_error_class,
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
      if ("eta_error_class" %in% names(eta)) {
        eta_concern <- eta[eta$eta_error_class == "concern", , drop = FALSE]
        eta_concern <- eta_concern[order(-eta_concern$abs_eta_error), , drop = FALSE]
        if (nrow(eta_concern) > 0L) {
          eta_concern_df <- data.frame(
            Family = eta_concern$family,
            Copula = eta_concern$copula,
            Method = eta_concern$method,
            Parameter = eta_concern$parameter,
            Term = if ("term" %in% names(eta_concern)) eta_concern$term else "",
            `Estimate eta` = .coverage_format_number(eta_concern$estimate_eta, 3),
            `True eta` = .coverage_format_number(eta_concern$true_eta, 3),
            `Abs. error` = .coverage_format_number(eta_concern$abs_eta_error, 3),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  list(
    eta_summary_df = eta_summary_df,
    eta_class_df = eta_class_df,
    eta_worst_df = eta_worst_df,
    eta_concern_df = eta_concern_df
  )
}

# ---- coverage-report-likelihood-tables.R ----

#' Build coverage report likelihood review tables
#'
#' @noRd
.coverage_report_likelihood_tables <- function(results) {
  margin_df <- NULL
  if ("margin_review_class" %in% names(results)) {
    margin_tab <- as.data.frame.matrix(table(results$method, results$margin_review_class))
    margin_df <- data.frame(Method = rownames(margin_tab), margin_tab, row.names = NULL, check.names = FALSE)
  }
  joint_df <- NULL
  if ("joint_review_class" %in% names(results)) {
    joint_tab <- as.data.frame.matrix(table(results$method, results$joint_review_class))
    joint_df <- data.frame(Method = rownames(joint_tab), joint_tab, row.names = NULL, check.names = FALSE)
  }
  likelihood_review_df <- NULL
  likelihood_review <- rep(FALSE, nrow(results))
  if ("margin_review_class" %in% names(results)) {
    likelihood_review <- likelihood_review | results$margin_review_class == "review"
  }
  if ("joint_review_class" %in% names(results)) {
    likelihood_review <- likelihood_review | results$joint_review_class == "review"
  }
  likelihood_review <- results[likelihood_review %in% TRUE, , drop = FALSE]
  if (nrow(likelihood_review) > 0L) {
    likelihood_review <- likelihood_review[order(
      likelihood_review$family,
      likelihood_review$copula,
      likelihood_review$method
    ), , drop = FALSE]
    likelihood_review_df <- data.frame(
      Family = likelihood_review$family,
      Copula = likelihood_review$copula,
      Method = likelihood_review$method,
      `M class` = likelihood_review$margin_review_class,
      `M gap %` = .coverage_format_number(likelihood_review$margin_gap_pct_vs_reference, 2),
      `J class` = likelihood_review$joint_review_class,
      `J delta %` = .coverage_format_number(likelihood_review$joint_delta_pct_vs_rs_separate, 2),
      `Marg LL` = .coverage_format_number(likelihood_review$marginal_loglik, 2),
      `Joint LL` = .coverage_format_number(likelihood_review$joint_loglik, 2),
      Sec = .coverage_format_number(likelihood_review$elapsed_sec, 2),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  list(
    margin_df = margin_df,
    joint_df = joint_df,
    likelihood_review_df = likelihood_review_df
  )
}

# ---- coverage-report-fit-runtime-tables.R ----

#' Build coverage report fit-success tables
#'
#' @noRd
.coverage_report_fit_tables <- function(results) {
  success_tab <- as.data.frame.matrix(table(results$method, results$success))
  if (!"TRUE" %in% names(success_tab)) success_tab$`TRUE` <- 0L
  if (!"FALSE" %in% names(success_tab)) success_tab$`FALSE` <- 0L
  fit_df <- data.frame(
    Method = rownames(success_tab),
    Successful = as.integer(success_tab$`TRUE`),
    Failed = as.integer(success_tab$`FALSE`),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  fit_df <- fit_df[order(fit_df$Method), , drop = FALSE]

  by_copula <- stats::aggregate(
    success ~ copula + method,
    data = results,
    FUN = function(x) c(successful = sum(x %in% TRUE), total = length(x))
  )
  by_copula <- data.frame(
    Copula = by_copula$copula,
    Method = by_copula$method,
    Successful = by_copula$success[, "successful"],
    Total = by_copula$success[, "total"],
    stringsAsFactors = FALSE
  )
  by_copula$Failed <- by_copula$Total - by_copula$Successful
  by_copula <- by_copula[order(by_copula$Copula, by_copula$Method), c("Copula", "Method", "Successful", "Failed", "Total")]

  list(fit_df = fit_df, by_copula = by_copula)
}

#' Build coverage report runtime tables
#'
#' @noRd
.coverage_report_runtime_tables <- function(results, top_n) {
  runtime <- stats::aggregate(elapsed_sec ~ method, results, sum)
  runtime$elapsed_min <- runtime$elapsed_sec / 60
  runtime$median_sec <- vapply(split(results$elapsed_sec, results$method), stats::median, numeric(1))[runtime$method]
  runtime$p90_sec <- vapply(split(results$elapsed_sec, results$method), function(x) as.numeric(stats::quantile(x, 0.9)), numeric(1))[runtime$method]
  runtime$max_sec <- vapply(split(results$elapsed_sec, results$method), max, numeric(1))[runtime$method]
  runtime_df <- data.frame(
    Method = runtime$method,
    `Total seconds` = .coverage_format_number(runtime$elapsed_sec, 2),
    `Total minutes` = .coverage_format_number(runtime$elapsed_min, 2),
    `Median seconds` = .coverage_format_number(runtime$median_sec, 2),
    `90th pct. seconds` = .coverage_format_number(runtime$p90_sec, 2),
    `Max seconds` = .coverage_format_number(runtime$max_sec, 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  runtime_df <- runtime_df[order(runtime_df$Method), , drop = FALSE]

  family_group <- c("family", if ("copula" %in% names(results)) "copula" else character(0))
  family_runtime <- stats::aggregate(
    elapsed_sec ~ .,
    data = results[c(family_group, "elapsed_sec")],
    sum
  )
  family_runtime <- family_runtime[order(-family_runtime$elapsed_sec), , drop = FALSE]
  family_runtime <- head(family_runtime, top_n)
  family_runtime_df <- data.frame(
    Family = family_runtime$family,
    Copula = if ("copula" %in% names(family_runtime)) family_runtime$copula else "",
    `Total seconds` = .coverage_format_number(family_runtime$elapsed_sec, 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  slow_cases <- results[order(-results$elapsed_sec), , drop = FALSE]
  slow_cases <- head(slow_cases, top_n)
  slow_case_df <- data.frame(
    Family = slow_cases$family,
    Copula = slow_cases$copula,
    Method = slow_cases$method,
    Seconds = .coverage_format_number(slow_cases$elapsed_sec, 2),
    Success = as.character(slow_cases$success),
    Converged = if ("converged" %in% names(slow_cases)) as.character(slow_cases$converged) else "",
    stringsAsFactors = FALSE
  )

  list(
    runtime_df = runtime_df,
    family_runtime_df = family_runtime_df,
    slow_case_df = slow_case_df
  )
}

# ---- coverage-report-smooth-tables.R ----

#' Build coverage report smooth smoke-test table
#'
#' @noRd
.coverage_report_smooth_table <- function(smooth_results) {
  smooth_df <- NULL
  if (!is.null(smooth_results) && nrow(smooth_results) > 0L && "smooth_eta_rmse" %in% names(smooth_results)) {
    smooth_split <- split(smooth_results, smooth_results$method)
    smooth_df <- do.call(rbind, lapply(names(smooth_split), function(method) {
      dat <- smooth_split[[method]]
      data.frame(
        Method = method,
        `Successful fits` = sum(dat$success),
        `Median eta RMSE` = .coverage_format_number(stats::median(dat$smooth_eta_rmse, na.rm = TRUE), 3),
        `Max eta RMSE` = .coverage_format_number(max(dat$smooth_eta_rmse, na.rm = TRUE), 3),
        `Max abs. eta error` = .coverage_format_number(max(dat$smooth_eta_max_abs_error, na.rm = TRUE), 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  }
  smooth_df
}
