#!/usr/bin/env Rscript

# Inventory public gamlss.dist family constructors without assuming that every
# exported uppercase function is a distribution family.

if (!requireNamespace("gamlss.dist", quietly = TRUE)) {
  stop("Install gamlss.dist before running the family audit.")
}

namespace <- asNamespace("gamlss.dist")
exports <- getNamespaceExports("gamlss.dist")
constructors <- sort(exports[grepl("^[A-Z][A-Z0-9.]*$", exports)])

instantiate <- function(name) {
  fun <- get(name, envir = namespace, inherits = FALSE)
  if (!is.function(fun)) return(NULL)
  old <- setwd(tempdir())
  on.exit(setwd(old), add = TRUE)
  suppressMessages(withCallingHandlers(
    tryCatch(do.call(fun, list()), error = function(e) NULL),
    warning = function(w) invokeRestart("muffleWarning")
  ))
}

formal_summary <- function(fun_name, parameters, response_arg) {
  if (!exists(fun_name, envir = namespace, mode = "function", inherits = FALSE)) {
    return(list(exists = FALSE, extra = "", required = ""))
  }
  fml <- formals(get(fun_name, envir = namespace, inherits = FALSE))
  standard <- c(response_arg, parameters, "log", "log.p", "lower.tail", "...")
  extra <- setdiff(names(fml), standard)
  required <- extra[vapply(fml[extra], identical, logical(1), quote(expr = ))]
  list(
    exists = TRUE,
    extra = paste(extra, collapse = ","),
    required = paste(required, collapse = ",")
  )
}

rows <- lapply(constructors, function(constructor) {
  object <- instantiate(constructor)
  if (!inherits(object, "gamlss.family")) return(NULL)
  family <- as.character(object$family[[1L]])
  parameters <- names(object$parameters)
  density <- formal_summary(paste0("d", family), parameters, "x")
  probability <- formal_summary(paste0("p", family), parameters, "q")
  quantile <- formal_summary(paste0("q", family), parameters, "p")
  type <- tolower(paste(as.character(object$type), collapse = " "))
  proposed <- if (family %in% c("MN3", "MN4", "MN5", "MULTIN")) {
    "exclude_multinomial"
  } else if (family %in% c("LG", "ZALG")) {
    "exclude_logarithmic_series"
  } else if (grepl("mixed", type, fixed = TRUE)) {
    "exclude_mixed_measure"
  } else if (family %in% c("BI", "BB", "DBI", "ZABB", "ZABI", "ZIBB", "ZIBI")) {
    "requires_binomial_denominator"
  } else if (grepl("continuous", type, fixed = TRUE)) {
    "candidate_continuous"
  } else if (grepl("discrete", type, fixed = TRUE)) {
    "candidate_discrete"
  } else {
    "manual_review"
  }
  data.frame(
    gamlss_dist_version = as.character(utils::packageVersion("gamlss.dist")),
    constructor = constructor,
    family = family,
    description = paste(as.character(object$family[-1L]), collapse = "; "),
    type = paste(as.character(object$type), collapse = "; "),
    parameters = paste(parameters, collapse = ","),
    linked_parameters = paste(
      parameters[vapply(parameters, function(parameter) {
        all(vapply(
          paste0(parameter, c(".linkfun", ".linkinv", ".dr")),
          function(field) is.function(object[[field]]), logical(1)
        ))
      }, logical(1))],
      collapse = ","
    ),
    has_density = density$exists,
    has_probability = probability$exists,
    has_quantile = quantile$exists,
    density_extra_args = density$extra,
    probability_extra_args = probability$extra,
    quantile_extra_args = quantile$extra,
    required_extra_args = paste(
      unique(c(density$required, probability$required, quantile$required)[
        nzchar(c(density$required, probability$required, quantile$required))
      ]),
      collapse = ","
    ),
    proposed_route = proposed,
    stringsAsFactors = FALSE
  )
})

audit <- do.call(rbind, Filter(Negate(is.null), rows))
audit$is_alias <- duplicated(audit$family) | duplicated(audit$family, fromLast = TRUE)
audit <- audit[order(audit$proposed_route, audit$family, audit$constructor), ]
rownames(audit) <- NULL

output <- file.path(
  "paper", "reviews", "2026-09-08-family-capability", "gamlss-family-audit.csv"
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(audit, output, row.names = FALSE, na = "")

cat("Wrote", nrow(audit), "constructor rows for", length(unique(audit$family)),
    "canonical families to", output, "\n")
print(as.data.frame(table(audit$proposed_route)), row.names = FALSE)
