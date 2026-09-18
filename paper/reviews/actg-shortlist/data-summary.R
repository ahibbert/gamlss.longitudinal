#!/usr/bin/env Rscript

# Descriptive summaries of the ACTG cd4 series used by fit-final-model.R:
# per-time moments/zero-inflation, cross-time correlation (Pearson + Kendall's
# tau), and empirical pairwise tail-dependence matrices to sanity-check the
# copula screening in fit-final-model.R (Student t vs Gumbel).

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools with install.packages('devtools').")
}
review_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))))
if (!length(review_dir) || is.na(review_dir)) review_dir <- getwd()
package_root <- normalizePath(file.path(review_dir, "..", "..", ".."))
devtools::load_all(package_root, quiet = TRUE)

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "gamlss.longitudinal", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

actg_data <- load_package_data("actg_cd4")
actg_data$week_factor <- factor(actg_data$week)

# ---- Per-time moment summaries ---------------------------------------------

sample_skewness <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m3 <- mean((x - m)^3)
  m3 / m2^1.5
}

sample_excess_kurtosis <- function(x) {
  x <- x[is.finite(x)]
  m <- mean(x)
  m2 <- mean((x - m)^2)
  m4 <- mean((x - m)^4)
  m4 / m2^2 - 3
}

time_summary <- do.call(rbind, lapply(split(actg_data, actg_data$week_factor), function(d) {
  cd4 <- d$cd4[is.finite(d$cd4)]
  data.frame(
    week = d$week[1L],
    n = length(cd4),
    n_zero = sum(cd4 == 0),
    zero_prop = mean(cd4 == 0),
    mean = mean(cd4),
    variance = stats::var(cd4),
    skewness = sample_skewness(cd4),
    excess_kurtosis = sample_excess_kurtosis(cd4)
  )
}))
time_summary <- time_summary[order(time_summary$week), ]
rownames(time_summary) <- NULL

cat("\n=== Per-time summary (mean/variance/skewness/kurtosis/zero-inflation) ===\n")
print(time_summary, digits = 4)

# ---- Cross-time correlation matrices ---------------------------------------
# Wide subject x week matrix; pairwise-complete correlations across weeks.

wide_cd4 <- reshape(
  actg_data[, c("subject", "week", "cd4")],
  idvar = "subject",
  timevar = "week",
  direction = "wide"
)
wide_cd4 <- wide_cd4[, -1L, drop = FALSE]
colnames(wide_cd4) <- sub("^cd4\\.", "week_", colnames(wide_cd4))

pearson_cor <- stats::cor(wide_cd4, use = "pairwise.complete.obs", method = "pearson")
kendall_cor <- stats::cor(wide_cd4, use = "pairwise.complete.obs", method = "kendall")

cat("\n=== Pearson correlation across weeks ===\n")
print(round(pearson_cor, 3))

cat("\n=== Kendall's tau across weeks ===\n")
print(round(kendall_cor, 3))

# ---- Empirical pairwise tail-dependence matrices ---------------------------
# Rank-transform each week's observed values to pseudo-uniform margins, then
# estimate upper/lower tail dependence as the empirical joint-exceedance
# probability at a threshold, conditional on both margins exceeding it.

pseudo_uniform <- function(x) {
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  out[ok] <- rank(x[ok], ties.method = "average") / (sum(ok) + 1)
  out
}

pseudo_obs <- as.data.frame(lapply(wide_cd4, pseudo_uniform))

empirical_tail_dependence <- function(u1, u2, q = 0.05, tail = c("upper", "lower")) {
  tail <- match.arg(tail)
  keep <- is.finite(u1) & is.finite(u2)
  u1 <- u1[keep]; u2 <- u2[keep]
  if (length(u1) < 2L) return(NA_real_)
  if (tail == "upper") {
    mean(u1 > 1 - q & u2 > 1 - q) / q
  } else {
    mean(u1 <= q & u2 <= q) / q
  }
}

weeks <- colnames(pseudo_obs)
upper_tail_matrix <- matrix(NA_real_, length(weeks), length(weeks), dimnames = list(weeks, weeks))
lower_tail_matrix <- upper_tail_matrix

for (i in seq_along(weeks)) {
  for (j in seq_along(weeks)) {
    if (i == j) {
      upper_tail_matrix[i, j] <- lower_tail_matrix[i, j] <- 1
      next
    }
    upper_tail_matrix[i, j] <- empirical_tail_dependence(pseudo_obs[[i]], pseudo_obs[[j]], tail = "upper")
    lower_tail_matrix[i, j] <- empirical_tail_dependence(pseudo_obs[[i]], pseudo_obs[[j]], tail = "lower")
  }
}

cat("\n=== Empirical upper-tail dependence (lambda_U, q = 0.05) ===\n")
print(round(upper_tail_matrix, 3))

cat("\n=== Empirical lower-tail dependence (lambda_L, q = 0.05) ===\n")
print(round(lower_tail_matrix, 3))

invisible(list(
  time_summary = time_summary,
  pearson_cor = pearson_cor,
  kendall_cor = kendall_cor,
  upper_tail_matrix = upper_tail_matrix,
  lower_tail_matrix = lower_tail_matrix
))
