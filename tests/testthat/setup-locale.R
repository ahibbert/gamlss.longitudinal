if (identical(.Platform$OS.type, "windows")) {
  unsupported_utf8_locale <- function(x) {
    grepl("^C\\.UTF-?8$", x, ignore.case = TRUE)
  }

  if (unsupported_utf8_locale(Sys.getenv("LC_ALL"))) {
    Sys.unsetenv("LC_ALL")
  }
  if (unsupported_utf8_locale(Sys.getenv("LANG"))) {
    Sys.unsetenv("LANG")
  }
}
