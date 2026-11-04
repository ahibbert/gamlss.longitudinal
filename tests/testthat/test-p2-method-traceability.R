traceability_split <- function(x) {
  trimws(strsplit(x, ";", fixed = TRUE)[[1]])
}

traceability_function_defs <- function(path) {
  exprs <- parse(path, keep.source = FALSE)
  defs <- character()

  for (expr in exprs) {
    if (!is.call(expr)) next
    op <- as.character(expr[[1]])[1]
    if (!op %in% c("<-", "=") || length(expr) < 3L || !is.call(expr[[3]])) next
    if (!identical(as.character(expr[[3]][[1]])[1], "function")) next
    defs <- c(defs, as.character(expr[[2]])[1])
  }

  defs
}

traceability_package_root <- function() {
  candidates <- unique(normalizePath(c(
    getwd(),
    file.path(getwd(), ".."),
    file.path(getwd(), "..", ".."),
    testthat::test_path("..", "..")
  ), winslash = "/", mustWork = FALSE))

  for (candidate in candidates) {
    if (file.exists(file.path(candidate, "inst", "standards", "method-traceability.csv"))) {
      return(candidate)
    }
  }

  getwd()
}

test_that("method traceability table references existing source and test files", {
  root <- traceability_package_root()
  traceability_path <- file.path(root, "inst", "standards", "method-traceability.csv")
  testthat::skip_if_not(
    file.exists(traceability_path) && dir.exists(file.path(root, "R")),
    "method traceability source-tree audit requires package source files"
  )

  traceability <- utils::read.csv(traceability_path, stringsAsFactors = FALSE)
  expect_true(all(c("area", "implementation_entry_point", "primary_tests") %in% names(traceability)))
  expect_equal(anyDuplicated(traceability$area), 0L)
  expect_true(all(nzchar(traceability$area)))
  expect_true(all(nzchar(traceability$implementation_entry_point)))
  expect_true(all(nzchar(traceability$primary_tests)))

  source_entries <- unlist(lapply(traceability$implementation_entry_point, traceability_split), use.names = FALSE)
  source_entries <- source_entries[nzchar(source_entries)]

  source_files <- character()
  exact_entries <- source_entries[grepl("::", source_entries, fixed = TRUE)]
  wildcard_entries <- source_entries[grepl("\\*", source_entries)]
  file_only_entries <- setdiff(source_entries[!grepl("::", source_entries, fixed = TRUE)], wildcard_entries)

  for (entry in exact_entries) {
    parts <- strsplit(entry, "::", fixed = TRUE)[[1]]
    source_files <- c(source_files, parts[1])
    expect_true(file.exists(file.path(root, parts[1])), info = entry)
  }

  for (entry in file_only_entries) {
    source_files <- c(source_files, entry)
    expect_true(file.exists(file.path(root, entry)), info = entry)
  }

  for (entry in wildcard_entries) {
    matches <- Sys.glob(file.path(root, entry))
    expect_true(length(matches) > 0L, info = entry)
    source_files <- c(source_files, matches)
  }

  test_entries <- unlist(lapply(traceability$primary_tests, traceability_split), use.names = FALSE)
  test_entries <- test_entries[nzchar(test_entries)]
  for (entry in test_entries) {
    expect_true(file.exists(file.path(root, entry)), info = entry)
  }

  expect_true(any(grepl("^R/", source_files)))
  expect_true(any(grepl("^tests/testthat/", test_entries)))

  covered_source_files <- normalizePath(
    ifelse(grepl("^R/", source_files), file.path(root, source_files), source_files),
    winslash = "/",
    mustWork = FALSE
  )
  package_source_files <- normalizePath(
    list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE),
    winslash = "/",
    mustWork = FALSE
  )
  uncovered_source_files <- setdiff(package_source_files, unique(covered_source_files))

  expect_equal(
    basename(uncovered_source_files),
    character(),
    info = paste(
      "Unmapped R source files:",
      paste(basename(uncovered_source_files), collapse = ", ")
    )
  )
})

test_that("method traceability table references defined R entry point functions", {
  root <- traceability_package_root()
  traceability_path <- file.path(root, "inst", "standards", "method-traceability.csv")
  testthat::skip_if_not(
    file.exists(traceability_path) && dir.exists(file.path(root, "R")),
    "method traceability source-tree audit requires package source files"
  )
  traceability <- utils::read.csv(traceability_path, stringsAsFactors = FALSE)
  source_entries <- unlist(lapply(traceability$implementation_entry_point, traceability_split), use.names = FALSE)
  exact_entries <- source_entries[grepl("::", source_entries, fixed = TRUE)]

  for (entry in exact_entries) {
    parts <- strsplit(entry, "::", fixed = TRUE)[[1]]
    file <- file.path(root, parts[1])
    fn <- sub("\\(\\)$", "", parts[2])
    expect_true(file.exists(file), info = entry)
    expect_true(fn %in% traceability_function_defs(file), info = entry)
  }
})
