mvt_suite_current_file <- function() {
  frame_files <- vapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value) || !nzchar(value)) "" else value
  }, character(1))
  frame_files <- frame_files[nzchar(frame_files)]
  if (length(frame_files) > 0L) {
    return(tail(frame_files, 1L))
  }
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) > 0L) {
    return(sub("^--file=", "", file_arg[[1L]]))
  }
  ""
}

mvt_suite_setup_path <- function() {
  current_file <- mvt_suite_current_file()
  candidates <- c(
    if (nzchar(current_file)) file.path(dirname(current_file), "00-multivariate-setup.R") else character(),
    file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R")
  )
  candidates <- normalizePath(candidates, winslash = "/", mustWork = FALSE)
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) "" else candidates[[1L]]
}

setup_path <- mvt_suite_setup_path()
if (is.na(setup_path) || !nzchar(setup_path)) {
  stop("Could not locate 00-multivariate-setup.R for the publication suite.", call. = FALSE)
}
setup_repo_hint <- normalizePath(file.path(dirname(setup_path), "..", "..", ".."), winslash = "/", mustWork = FALSE)
setup_old_wd <- getwd()
tryCatch({
  if (dir.exists(setup_repo_hint) && file.exists(file.path(setup_repo_hint, "DESCRIPTION"))) {
    setwd(setup_repo_hint)
  }
  source(setup_path)
}, finally = {
  setwd(setup_old_wd)
})

suite_run_enabled <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE", FALSE)
suite_strict_readiness <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUITE_STRICT_READINESS", TRUE)
suite_preflight_only <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUITE_PREFLIGHT_ONLY", FALSE)
suite_id <- mvt_env("GAMLSS_LONGITUDINAL_MVT_SUITE_ID", paste0("publication_suite_", mvt_timestamp()))
suite_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_SUITE_DIR", file.path(mvt_output_root, suite_id))
suite_roles <- mvt_env_vector(
  "GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES",
  c("pilot", "main_core", "appendix")
)
if (mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUITE_INCLUDE_SPECIAL", FALSE)) {
  suite_roles <- c(suite_roles, "special_gamlss")
}
if (mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_SUITE_INCLUDE_GLMM_SENSITIVITY", FALSE)) {
  suite_roles <- c(suite_roles, "glmm_sensitivity")
}
suite_roles <- unique(suite_roles)

mvt_suite_csv <- function(x) paste(x, collapse = ",")

mvt_suite_env_default <- function(role, suffix, default) {
  mvt_env(paste0("GAMLSS_LONGITUDINAL_MVT_SUITE_", toupper(role), "_", suffix), default)
}

mvt_suite_script <- function(name) file.path(mvt_script_dir, name)

mvt_suite_role_specs <- function(suite_dir) {
  default_methods <- mvt_suite_csv(mvt_default_comparators())
  list(
    pilot = list(
      role = "pilot",
      script = mvt_suite_script("01-run-pilot-grid.R"),
      run_dir = file.path(suite_dir, "pilot"),
      env = c(
        GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5,t20",
        GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,gamma,binomial",
        GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_covariate_dependent_adjacent",
        GAMLSS_LONGITUDINAL_MVT_REPS = mvt_suite_env_default("pilot", "REPS", "5"),
        GAMLSS_LONGITUDINAL_MVT_COMPARATORS = default_methods,
        GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "false",
        GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "",
        GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "5"
      )
    ),
    main_core = list(
      role = "main_core",
      script = mvt_suite_script("02-run-main-grid.R"),
      run_dir = file.path(suite_dir, "main_core"),
      env = c(
        GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "core",
        GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t20",
        GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial",
        GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_time_varying_adjacent,native_covariate_dependent_adjacent",
        GAMLSS_LONGITUDINAL_MVT_REPS = mvt_suite_env_default("main_core", "REPS", "100"),
        GAMLSS_LONGITUDINAL_MVT_COMPARATORS = default_methods,
        GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "false",
        GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "",
        GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "10",
        GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = mvt_env("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", "30")
      )
    ),
    appendix = list(
      role = "appendix",
      script = mvt_suite_script("02-run-main-grid.R"),
      run_dir = file.path(suite_dir, "appendix"),
      env = c(
        GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "appendix",
        GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5,t20,t50",
        GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial",
        GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_time_varying_adjacent,native_covariate_dependent_adjacent",
        GAMLSS_LONGITUDINAL_MVT_REPS = mvt_suite_env_default("appendix", "REPS", "100"),
        GAMLSS_LONGITUDINAL_MVT_COMPARATORS = default_methods,
        GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "false",
        GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "",
        GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "10",
        GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC = mvt_env("GAMLSS_LONGITUDINAL_MVT_GEE_UNSTRUCTURED_TIMEOUT_SEC", "30")
      )
    ),
    special_gamlss = list(
      role = "special_gamlss",
      script = mvt_suite_script("02-run-main-grid.R"),
      run_dir = file.path(suite_dir, "special_gamlss"),
      env = c(
        GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "core",
        GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t20",
        GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gg_continuous",
        GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_ar1,native_covariate_dependent_adjacent",
        GAMLSS_LONGITUDINAL_MVT_REPS = mvt_suite_env_default("special_gamlss", "REPS", "100"),
        GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "gamlss.longitudinal",
        GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "true",
        GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "",
        GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "10"
      )
    ),
    glmm_sensitivity = list(
      role = "glmm_sensitivity",
      script = mvt_suite_script("02-run-main-grid.R"),
      run_dir = file.path(suite_dir, "glmm_sensitivity"),
      env = c(
        GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE = "appendix",
        GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS = "t5,t20,t50",
        GAMLSS_LONGITUDINAL_MVT_FAMILIES = "gaussian,poisson,gamma,binomial",
        GAMLSS_LONGITUDINAL_MVT_DEPENDENCE = "external_exchangeable,external_ar1,native_time_varying_adjacent,native_covariate_dependent_adjacent",
        GAMLSS_LONGITUDINAL_MVT_REPS = mvt_suite_env_default("glmm_sensitivity", "REPS", "100"),
        GAMLSS_LONGITUDINAL_MVT_COMPARATORS = "glm,glmm,glmm_slope,gamCopula_markov,gamCopula_vine,gamlss.longitudinal",
        GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL = "false",
        GAMLSS_LONGITUDINAL_MVT_N_SUBJECT = "",
        GAMLSS_LONGITUDINAL_MVT_CHECKPOINT_EVERY = "10"
      )
    )
  )
}

mvt_suite_with_env <- function(vars, expr) {
  vars <- stats::setNames(as.character(vars), names(vars))
  old <- Sys.getenv(names(vars), unset = NA_character_)
  on.exit({
    for (i in seq_along(old)) {
      if (is.na(old[[i]])) {
        Sys.unsetenv(names(old)[[i]])
      } else {
        do.call(Sys.setenv, stats::setNames(as.list(old[[i]]), names(old)[[i]]))
      }
    }
  }, add = TRUE)
  do.call(Sys.setenv, as.list(vars))
  eval.parent(substitute(expr))
}

mvt_suite_plan <- function(specs) {
  rows <- lapply(specs, function(spec) {
    env <- spec$env
    estimated_grid <- mvt_suite_grid(spec, quiet = TRUE)
    data.frame(
      role = spec$role,
      run_dir = normalizePath(spec$run_dir, winslash = "/", mustWork = FALSE),
      script = normalizePath(spec$script, winslash = "/", mustWork = FALSE),
      reps = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_REPS"),
      estimated_cases = nrow(estimated_grid),
      max_total_rows = if (nrow(estimated_grid) > 0L) max(estimated_grid$total_rows, na.rm = TRUE) else NA_integer_,
      timepoints = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS"),
      families = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_FAMILIES"),
      dependence = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_DEPENDENCE"),
      comparators = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_COMPARATORS"),
      main_scope = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_MAIN_SCOPE"),
      include_special = mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL", "false"),
      stringsAsFactors = FALSE
    )
  })
  mvt_bind_rows_fill(rows)
}

mvt_suite_env_value <- function(env, name, default = "") {
  if (name %in% names(env)) env[[name]] else default
}

mvt_suite_grid <- function(spec, quiet = FALSE) {
  env <- spec$env
  tryCatch(
    mvt_suite_with_env(env, {
      mvt_expand_grid(
        time_names = mvt_split_csv_value(mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_TIMEPOINTS")),
        family_names = mvt_split_csv_value(mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_FAMILIES")),
        dependence_names = mvt_split_csv_value(mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_DEPENDENCE")),
        reps = as.integer(mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_REPS", "1")),
        include_special = identical(tolower(mvt_suite_env_value(env, "GAMLSS_LONGITUDINAL_MVT_INCLUDE_SPECIAL", "false")), "true"),
        include_appendix = TRUE
      )
    }),
    error = function(e) {
      if (!quiet) warning("Could not expand suite grid for role ", spec$role, ": ", conditionMessage(e), call. = FALSE)
      data.frame()
    }
  )
}

mvt_write_suite_plan <- function(plan, suite_dir, run_enabled) {
  dir.create(suite_dir, recursive = TRUE, showWarnings = FALSE)
  mvt_write_csv(plan, file.path(suite_dir, "publication_suite_plan.csv"))
  lines <- c(
    "# Publication Simulation Suite Plan",
    "",
    paste("Suite directory:", normalizePath(suite_dir, winslash = "/", mustWork = FALSE)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Run enabled:", run_enabled),
    "",
    "This script writes the plan by default. To execute the full suite, set:",
    "",
    "```powershell",
    '$env:GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE = "true"',
    'Rscript "paper/R/09-simulation-multivariate-longitudinal/10-run-publication-suite.R"',
    "```",
    "",
    "## Roles",
    ""
  )
  role_lines <- paste0(
    "- ",
    plan$role,
    ": reps=",
    plan$reps,
    ", cases=",
    plan$estimated_cases,
    ", max_rows=",
    plan$max_total_rows,
    ", T=",
    plan$timepoints,
    ", families=",
    plan$families,
    ", dependence=",
    plan$dependence
  )
  writeLines(c(lines, role_lines), file.path(suite_dir, "publication_suite_plan.md"), useBytes = TRUE)
  invisible(plan)
}

mvt_suite_artifacts <- function(plan, suite_dir) {
  rows <- lapply(seq_len(nrow(plan)), function(i) {
    run_dir <- plan$run_dir[[i]]
    data.frame(
      role = plan$role[[i]],
      run_dir = run_dir,
      run_dir_exists = dir.exists(run_dir),
      review_audit = file.exists(file.path(run_dir, "review_audit.md")),
      review_bundle = file.exists(file.path(run_dir, "review_bundle", "README.md")),
      fit_status = file.exists(file.path(run_dir, "fit_status_by_rep.csv")),
      paper_tables = dir.exists(file.path(run_dir, "paper_tables")),
      figures = dir.exists(file.path(run_dir, "figures")),
      artifact_manifest = file.exists(file.path(run_dir, "artifact_manifest.csv")),
      stringsAsFactors = FALSE
    )
  })
  artifacts <- mvt_bind_rows_fill(rows)
  mvt_write_csv(artifacts, file.path(suite_dir, "publication_suite_artifacts.csv"))
  artifacts
}

mvt_suite_readiness_env <- function(all_specs, roles) {
  role_to_env <- c(
    pilot = "GAMLSS_LONGITUDINAL_MVT_PILOT_RUN_DIR",
    main_core = "GAMLSS_LONGITUDINAL_MVT_MAIN_RUN_DIR",
    appendix = "GAMLSS_LONGITUDINAL_MVT_APPENDIX_RUN_DIR",
    special_gamlss = "GAMLSS_LONGITUDINAL_MVT_SPECIAL_RUN_DIR",
    glmm_sensitivity = "GAMLSS_LONGITUDINAL_MVT_GLMM_SENSITIVITY_RUN_DIR"
  )
  env <- stats::setNames(rep("", length(role_to_env)), unname(role_to_env))
  for (role in intersect(roles, names(role_to_env))) {
    env[[role_to_env[[role]]]] <- all_specs[[role]]$run_dir
  }
  env
}

mvt_write_suite_index <- function(plan, suite_dir, run_enabled, strict_readiness, readiness = NULL) {
  artifacts <- mvt_suite_artifacts(plan, suite_dir)
  readiness_path <- file.path(suite_dir, "publication_readiness_audit.md")
  preflight_path <- file.path(suite_dir, "publication_suite_preflight.md")
  readiness_line <- if (is.null(readiness)) {
    "Not run."
  } else {
    paste("Required evidence ready:", readiness$ready)
  }
  artifact_lines <- if (nrow(artifacts) == 0L) {
    "No role artifacts found."
  } else {
    paste0(
      "- ",
      artifacts$role,
      ": run_dir=",
      artifacts$run_dir_exists,
      ", audit=",
      artifacts$review_audit,
      ", bundle=",
      artifacts$review_bundle,
      ", tables=",
      artifacts$paper_tables,
      ", figures=",
      artifacts$figures
    )
  }
  lines <- c(
    "# Multivariate Publication Suite",
    "",
    paste("Suite directory:", normalizePath(suite_dir, winslash = "/", mustWork = FALSE)),
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Run enabled:", run_enabled),
    paste("Strict publication readiness:", strict_readiness),
    paste("Readiness:", readiness_line),
    "",
    "## Start Here",
    "",
    "- `publication_suite_plan.md`",
    "- `publication_suite_artifacts.csv`",
    if (file.exists(preflight_path)) "- `publication_suite_preflight.md`" else character(),
    if (file.exists(readiness_path)) "- `publication_readiness_audit.md`" else character(),
    "",
    "## Role Artifacts",
    "",
    artifact_lines
  )
  writeLines(lines, file.path(suite_dir, "README.md"), useBytes = TRUE)
  invisible(file.path(suite_dir, "README.md"))
}

mvt_suite_write_preflight_role <- function(spec, seed_base = 20260818L) {
  role_env <- c(
    spec$env,
    GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = spec$run_dir,
    GAMLSS_LONGITUDINAL_MVT_RUN_DIR = spec$run_dir
  )
  mvt_suite_with_env(role_env, {
    dir.create(spec$run_dir, recursive = TRUE, showWarnings = FALSE)
    grid <- mvt_suite_grid(spec)
    resume <- mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_RESUME", TRUE)
    mvt_write_csv(grid, file.path(spec$run_dir, "scenario_grid.csv"))
    preflight <- mvt_write_preflight(grid, spec$run_dir, require_gamcopula = TRUE, resume = resume)
    mvt_write_run_metadata(spec$run_dir, grid, seed_base = seed_base, require_gamcopula = TRUE, resume = resume)
    mvt_write_artifact_manifest(spec$run_dir)
    preflight
  })
}

mvt_suite_write_preflight <- function(specs, suite_dir, plan) {
  seed_base <- mvt_env_int("GAMLSS_LONGITUDINAL_MVT_SEED", 20260818L)
  rows <- lapply(specs, function(spec) {
    preflight <- mvt_suite_write_preflight_role(spec, seed_base = seed_base)
    data.frame(
      role = spec$role,
      checks = nrow(preflight),
      failures = sum(preflight$status == "fail"),
      warnings = sum(preflight$status == "warn"),
      run_dir = normalizePath(spec$run_dir, winslash = "/", mustWork = FALSE),
      stringsAsFactors = FALSE
    )
  })
  summary <- mvt_bind_rows_fill(rows)
  mvt_write_csv(summary, file.path(suite_dir, "publication_suite_preflight.csv"))
  lines <- c(
    "# Publication Suite Preflight",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste("Roles:", paste(plan$role, collapse = ", ")),
    "",
    "## Summary",
    "",
    paste0("- ", summary$role, ": failures=", summary$failures, ", warnings=", summary$warnings, ", checks=", summary$checks)
  )
  writeLines(lines, file.path(suite_dir, "publication_suite_preflight.md"), useBytes = TRUE)
  mvt_write_suite_index(plan, suite_dir, run_enabled = FALSE, strict_readiness = suite_strict_readiness)
  invisible(summary)
}

all_specs <- mvt_suite_role_specs(suite_dir)
unknown_roles <- setdiff(suite_roles, names(all_specs))
if (length(unknown_roles) > 0L) {
  stop("Unknown publication suite role(s): ", paste(unknown_roles, collapse = ", "), call. = FALSE)
}
specs <- all_specs[suite_roles]
plan <- mvt_suite_plan(specs)
mvt_write_suite_plan(plan, suite_dir, suite_run_enabled)
mvt_write_suite_index(plan, suite_dir, suite_run_enabled, suite_strict_readiness)

message("Publication suite plan written to: ", suite_dir)
if (isTRUE(suite_preflight_only)) {
  preflight_summary <- mvt_suite_write_preflight(specs, suite_dir, plan)
  message("Publication suite preflight written to: ", suite_dir)
  if (any(preflight_summary$failures > 0L)) {
    stop("Publication suite preflight found failing role checks. See ", file.path(suite_dir, "publication_suite_preflight.md"), call. = FALSE)
  }
} else if (!isTRUE(suite_run_enabled)) {
  message("Dry run only. Set GAMLSS_LONGITUDINAL_MVT_RUN_PUBLICATION_SUITE=true to execute the suite.")
} else {
  required_roles <- c("pilot", "main_core")
  missing_required_roles <- setdiff(required_roles, suite_roles)
  if (isTRUE(suite_strict_readiness) && length(missing_required_roles) > 0L) {
    stop(
      "Strict publication readiness requires suite role(s): ",
      paste(missing_required_roles, collapse = ", "),
      ". Add them to GAMLSS_LONGITUDINAL_MVT_SUITE_ROLES or set ",
      "GAMLSS_LONGITUDINAL_MVT_SUITE_STRICT_READINESS=false for a partial run.",
      call. = FALSE
    )
  }

  source_used <- mvt_load_package()
  mvt_require_namespaces(
    c("gamlss", "gamlss.dist", "mvtnorm", "VineCopula", "gamCopula", "geepack", "lme4", "callr", "ggplot2"),
    strict = TRUE
  )
  message("Package source: ", source_used)

  post_scripts <- mvt_suite_script(c(
    "04-add-variogram-scores.R",
    "05-write-paper-tables.R",
    "06-make-diagnostic-plots.R",
    "07-review-audit.R",
    "09-make-review-bundle.R"
  ))

  for (spec in specs) {
    role_env <- c(
      spec$env,
      GAMLSS_LONGITUDINAL_MVT_OUTPUT_DIR = spec$run_dir,
      GAMLSS_LONGITUDINAL_MVT_RUN_DIR = spec$run_dir
    )
    message("Starting publication suite role: ", spec$role)
    mvt_suite_with_env(role_env, {
      workflow_old_wd <- getwd()
      tryCatch({
        setwd(mvt_repo_root)
        source(spec$script, local = new.env(parent = globalenv()))
        for (post_script in post_scripts) {
          source(post_script, local = new.env(parent = globalenv()))
        }
      }, finally = {
        setwd(workflow_old_wd)
      })
    })
    message("Completed publication suite role: ", spec$role)
    mvt_write_suite_index(plan, suite_dir, suite_run_enabled, suite_strict_readiness)
  }

  readiness <- NULL
  if (isTRUE(suite_strict_readiness)) {
    readiness_env <- mvt_suite_readiness_env(all_specs, suite_roles)
    mvt_suite_with_env(readiness_env, {
      readiness <<- mvt_publication_readiness_audit(suite_dir)
      if (!isTRUE(readiness$ready)) {
        stop("Publication readiness audit did not pass. See ", file.path(suite_dir, "publication_readiness_audit.md"), call. = FALSE)
      }
    })
  } else {
    message("Strict publication readiness skipped for partial suite run.")
  }

  mvt_write_suite_index(plan, suite_dir, suite_run_enabled, suite_strict_readiness, readiness = readiness)

  if (isTRUE(suite_strict_readiness)) {
    message("Publication suite complete and readiness audit passed: ", suite_dir)
  } else {
    message("Publication suite complete: ", suite_dir)
  }
}
