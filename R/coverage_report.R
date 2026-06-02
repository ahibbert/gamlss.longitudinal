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
  top_n = 10
) {
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
  dir.create(dirname(output_tex), recursive = TRUE, showWarnings = FALSE)

  esc <- function(x) {
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
  tt <- function(x) paste0("\\texttt{", esc(x), "}")
  fmt_num <- function(x, digits = 2) {
    out <- ifelse(is.finite(x), format(round(x, digits), nsmall = digits, trim = TRUE), "--")
    as.character(out)
  }
  tex_table <- function(df, caption, align = NULL) {
    if (is.null(df) || nrow(df) == 0L) {
      return(c("\\paragraph{" %+% esc(caption) %+% "} No rows available."))
    }
    if (is.null(align)) align <- paste0("l", paste(rep("r", max(0L, ncol(df) - 1L)), collapse = ""))
    header <- paste(esc(names(df)), collapse = " & ")
    body <- apply(df, 1L, function(row) paste(esc(row), collapse = " & "))
    c(
      "\\begin{table}[htbp]",
      "\\centering",
      paste0("\\caption{", esc(caption), "}"),
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
  tex_longtable <- function(df, caption, align = NULL, font = "\\scriptsize") {
    if (is.null(df) || nrow(df) == 0L) {
      return(c("\\paragraph{" %+% esc(caption) %+% "} No rows available."))
    }
    if (is.null(align)) align <- paste(rep("l", ncol(df)), collapse = "")
    header <- paste(esc(names(df)), collapse = " & ")
    body <- apply(df, 1L, function(row) paste(esc(row), collapse = " & "))
    c(
      "\\begingroup",
      "\\setlength{\\tabcolsep}{2pt}",
      "\\begin{center}",
      font,
      paste0("\\begin{longtable}{", align, "}"),
      paste0("\\caption{", esc(caption), "}\\\\"),
      "\\toprule",
      paste0(header, " \\\\"),
      "\\midrule",
      "\\endfirsthead",
      paste0("\\multicolumn{", ncol(df), "}{l}{\\scriptsize\\emph{", esc(caption), " continued}}\\\\"),
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
  `%+%` <- function(a, b) paste0(a, b)

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

  runtime <- stats::aggregate(elapsed_sec ~ method, results, sum)
  runtime$elapsed_min <- runtime$elapsed_sec / 60
  runtime$median_sec <- vapply(split(results$elapsed_sec, results$method), stats::median, numeric(1))[runtime$method]
  runtime$p90_sec <- vapply(split(results$elapsed_sec, results$method), function(x) as.numeric(stats::quantile(x, 0.9)), numeric(1))[runtime$method]
  runtime$max_sec <- vapply(split(results$elapsed_sec, results$method), max, numeric(1))[runtime$method]
  runtime_df <- data.frame(
    Method = runtime$method,
    `Total seconds` = fmt_num(runtime$elapsed_sec, 2),
    `Total minutes` = fmt_num(runtime$elapsed_min, 2),
    `Median seconds` = fmt_num(runtime$median_sec, 2),
    `90th pct. seconds` = fmt_num(runtime$p90_sec, 2),
    `Max seconds` = fmt_num(runtime$max_sec, 2),
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
    `Total seconds` = fmt_num(family_runtime$elapsed_sec, 2),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  slow_cases <- results[order(-results$elapsed_sec), , drop = FALSE]
  slow_cases <- head(slow_cases, top_n)
  slow_case_df <- data.frame(
    Family = slow_cases$family,
    Copula = slow_cases$copula,
    Method = slow_cases$method,
    Seconds = fmt_num(slow_cases$elapsed_sec, 2),
    Success = as.character(slow_cases$success),
    Converged = if ("converged" %in% names(slow_cases)) as.character(slow_cases$converged) else "",
    stringsAsFactors = FALSE
  )

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
        Median = fmt_num(eta_summary$Median, 3),
        `75th pct.` = fmt_num(eta_summary$P75, 3),
        `90th pct.` = fmt_num(eta_summary$P90, 3),
        `95th pct.` = fmt_num(eta_summary$P95, 3),
        Max = fmt_num(eta_summary$Max, 3),
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
        `Estimate eta` = fmt_num(eta_worst$estimate_eta, 3),
        `True eta` = fmt_num(eta_worst$true_eta, 3),
        `Abs. error` = fmt_num(eta_worst$abs_eta_error, 3),
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
            `Estimate eta` = fmt_num(eta_concern$estimate_eta, 3),
            `True eta` = fmt_num(eta_concern$true_eta, 3),
            `Abs. error` = fmt_num(eta_concern$abs_eta_error, 3),
            check.names = FALSE,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

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
      `M gap %` = fmt_num(likelihood_review$margin_gap_pct_vs_reference, 2),
      `J class` = likelihood_review$joint_review_class,
      `J delta %` = fmt_num(likelihood_review$joint_delta_pct_vs_rs_separate, 2),
      `Marg LL` = fmt_num(likelihood_review$marginal_loglik, 2),
      `Joint LL` = fmt_num(likelihood_review$joint_loglik, 2),
      Sec = fmt_num(likelihood_review$elapsed_sec, 2),
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }

  smooth_df <- NULL
  if (!is.null(smooth_results) && nrow(smooth_results) > 0L && "smooth_eta_rmse" %in% names(smooth_results)) {
    smooth_split <- split(smooth_results, smooth_results$method)
    smooth_df <- do.call(rbind, lapply(names(smooth_split), function(method) {
      dat <- smooth_split[[method]]
      data.frame(
        Method = method,
        `Successful fits` = sum(dat$success),
        `Median eta RMSE` = fmt_num(stats::median(dat$smooth_eta_rmse, na.rm = TRUE), 3),
        `Max eta RMSE` = fmt_num(max(dat$smooth_eta_rmse, na.rm = TRUE), 3),
        `Max abs. eta error` = fmt_num(max(dat$smooth_eta_max_abs_error, na.rm = TRUE), 3),
        check.names = FALSE,
        stringsAsFactors = FALSE
      )
    }))
  }

  copulas <- paste(sort(unique(results$copula)), collapse = ", ")
  methods <- paste(sort(unique(results$method)), collapse = ", ")
  designs <- paste(sort(unique(results$design)), collapse = ", ")
  n_families <- length(unique(results$family))
  n_cases <- nrow(results)
  total_success <- sum(results$success)
  total_failed <- nrow(results) - total_success
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
    paste0("\\title{", esc(title), "}"),
    "\\author{\\texttt{gamlss.longitudinal}}",
    paste0("\\date{", esc(label_text), "}"),
    "",
    "\\begin{document}",
    "\\maketitle",
    "",
    "\\section*{Scope}",
    paste0("This report was generated automatically from coverage simulation output. It includes ",
           n_families, " margin families, copula(s) ", tt(copulas), ", design(s) ",
           tt(designs), ", and method(s) ", tt(methods), "."),
    "",
    paste0("Total fitted cases: ", n_cases, ". Successful: ", total_success,
           ". Failed: ", total_failed, "."),
    "",
    if (!is.null(results_file)) paste0("Results file: ", tt(basename(results_file)), ".") else NULL,
    if (!is.null(parameter_results_file) && file.exists(parameter_results_file)) {
      paste0("Parameter recovery file: ", tt(basename(parameter_results_file)), ".")
    } else NULL,
    "",
    "\\section*{Fit Success}",
    tex_table(fit_df, "Fit success by method"),
    tex_table(by_copula, "Fit success by copula and method"),
    "",
    "\\section*{Runtime}",
    wall_note,
    tex_table(runtime_df, "Runtime by method"),
    tex_table(family_runtime_df, "Slowest family/copula combinations by summed elapsed time"),
    tex_table(slow_case_df, "Slowest individual fits"),
    "",
    "\\section*{Parameter Recovery}",
    "Recovery is summarised on the eta scale for available parameter estimates.",
    if (!is.null(eta_summary_df)) tex_table(eta_summary_df, "Eta-scale absolute error summary") else "No eta-scale parameter recovery rows were available.",
    if (!is.null(eta_class_df)) tex_table(eta_class_df, "Eta-scale recovery classification") else NULL,
    if (!is.null(eta_worst_df)) tex_table(eta_worst_df, "Largest non-reference eta-scale errors") else NULL,
    "",
    "\\section*{Likelihood Review}",
    "The percentage likelihood-gap review columns are useful screening diagnostics, but can be unstable when reference log-likelihoods are close to zero or positive. Prefer a per-observation absolute log-likelihood gap for future multi-copula review thresholds.",
    if (!is.null(margin_df)) tex_table(margin_df, "Marginal likelihood review class versus reference") else NULL,
    if (!is.null(joint_df)) tex_table(joint_df, "Joint likelihood review class versus rs\\_separate") else NULL,
    "",
    if (!is.null(smooth_df)) "\\section*{Smooth-Term Smoke Layer}" else NULL,
    if (!is.null(smooth_df)) "Smooth recovery rows are summarised when a smooth smoke-test result file is supplied." else NULL,
    if (!is.null(smooth_df)) tex_table(smooth_df, "Smooth subset success and eta recovery") else NULL,
    "",
    "\\section*{Interpretation}",
    "A clean fit-success table means the supported family/copula/method combinations completed with finite likelihoods under the recorded simulation design. Runtime and recovery tables should be used to identify targeted optimisation or larger replicated recovery studies.",
    "",
    "\\appendix",
    "\\section{Detailed Review Rows}",
    "This appendix lists all rows requiring manual review under the current coverage thresholds.",
    "\\subsection{Parameter Recovery Concerns}",
    "Rows below are all eta-scale parameter recovery entries classified as \\texttt{concern}.",
    if (!is.null(eta_concern_df)) tex_longtable(eta_concern_df, "All parameter recovery concern rows") else "No parameter recovery concern rows were recorded.",
    "\\subsection{Likelihood Review Rows}",
    "Rows below are all fits where either marginal or joint likelihood review class is \\texttt{review}.",
    if (!is.null(likelihood_review_df)) tex_longtable(likelihood_review_df, "All likelihood review rows", font = "\\tiny") else "No likelihood review rows were recorded.",
    "",
    "\\end{document}"
  )
  lines <- unlist(lines, use.names = FALSE)
  writeLines(lines, output_tex, useBytes = TRUE)

  pdf_file <- sub("\\.tex$", ".pdf", output_tex, ignore.case = TRUE)
  compiled <- FALSE
  if (isTRUE(compile_pdf)) {
    compiled <- .coverage_compile_latex(output_tex, pdf_file)
  }
  invisible(list(tex_file = normalizePath(output_tex, winslash = "/", mustWork = FALSE),
                 pdf_file = normalizePath(pdf_file, winslash = "/", mustWork = FALSE),
                 compiled = compiled))
}

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
    ok <- tryCatch({
      tinytex::latexmk(tex_name, clean = FALSE)
      file.exists(basename(pdf_file))
    }, error = function(e) FALSE)
  }
  if (!ok) {
    engine <- Sys.which("pdflatex")
    if (nzchar(engine)) {
      ok <- tryCatch({
        system2(engine, c("-interaction=nonstopmode", "-halt-on-error", tex_name), stdout = TRUE, stderr = TRUE)
        system2(engine, c("-interaction=nonstopmode", "-halt-on-error", tex_name), stdout = TRUE, stderr = TRUE)
        file.exists(basename(pdf_file))
      }, error = function(e) FALSE)
    }
  }
  isTRUE(ok)
}
