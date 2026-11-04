test_that("withdrawn application identifiers are absent from the public repository", {
  root <- normalizePath(file.path("..", ".."), winslash = "/", mustWork = TRUE)
  skip_if_not(file.exists(file.path(root, "DESCRIPTION")), "repository checkout unavailable")

  prohibited <- c(
    paste0("li", "pid"),
    paste0("prava", "statin"),
    paste0("chole", "sterol")
  )
  pattern <- paste(prohibited, collapse = "|")
  relative <- suppressWarnings(system2("git", c("-C", shQuote(root), "ls-files"), stdout = TRUE, stderr = FALSE))
  skip_if(is.null(relative) || !length(relative) || !is.null(attr(relative, "status")), "Git index unavailable")
  paths <- file.path(root, relative)
  present <- file.exists(paths)
  relative <- relative[present]
  paths <- paths[present]
  expect_false(any(grepl(pattern, relative, ignore.case = TRUE, perl = TRUE)))

  text_paths <- paths[grepl("[.](R|Rd|Rmd|qmd|md|tex|csv|tsv|txt|ya?ml|json)$", paths, ignore.case = TRUE)]
  text_paths <- text_paths[!grepl("/(paper/)?data/", text_paths)]
  hits <- vapply(text_paths, function(path) {
    content <- tryCatch(readLines(path, warn = FALSE, encoding = "UTF-8"), error = function(error) character())
    any(grepl(pattern, content, ignore.case = TRUE, perl = TRUE))
  }, logical(1))
  expect_false(any(hits), info = paste(relative[match(text_paths[hits], paths)], collapse = ", "))
})
