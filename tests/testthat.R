library(testthat)
library(gamlss.longitudinal)

# Some Windows runners inherit the POSIX-only C.UTF-8 locale. Removing it here
# prevents testthat from emitting one locale warning per expectation.
if (identical(.Platform$OS.type, "windows") &&
    grepl("^C\\.UTF-?8$", Sys.getenv("LC_ALL"), ignore.case = TRUE)) {
  Sys.unsetenv("LC_ALL")
}

# Retain the CRAN reporter's conditions even if tests fail. Older testthat
# versions store warnings among problems; newer versions keep them separately.
check_reporter <- testthat::CheckReporter$new()
tryCatch(
  test_check("gamlss.longitudinal", reporter = check_reporter),
  finally = {
    problems <- check_reporter$problems$as_list()
    warnings <- if (is.null(check_reporter$warnings)) {
      Filter(function(x) inherits(x, "expectation_warning"), problems)
    } else {
      check_reporter$warnings$as_list()
    }
    saveRDS(
      list(
        warnings = warnings,
        failures = Filter(function(x) !inherits(x, "expectation_warning"), problems),
        skips = check_reporter$skips$as_list(),
        passed = check_reporter$n_ok,
        session_info = sessionInfo(),
        locale = Sys.getlocale()
      ),
      "testthat-conditions.rds",
      version = 2
    )
  }
)
