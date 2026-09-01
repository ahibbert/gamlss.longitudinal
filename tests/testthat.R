library(testthat)
library(gamlss.longitudinal)

# Some Windows runners inherit the POSIX-only C.UTF-8 locale. Removing it here
# prevents testthat from emitting one locale warning per expectation.
if (identical(.Platform$OS.type, "windows") &&
    grepl("^C\\.UTF-?8$", Sys.getenv("LC_ALL"), ignore.case = TRUE)) {
  Sys.unsetenv("LC_ALL")
}

test_check("gamlss.longitudinal")
