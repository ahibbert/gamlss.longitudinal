# Fit-free: only refit/simulation boundaries are mocked; the real bootstrap loop,
# condition capture and callback ordering execute. No models are fitted or read.
bootstrap_retention_fixture <- function(refit, original = FALSE) {
  ns <- asNamespace("gamlss.longitudinal")
  e <- new.env(parent = ns)
  if (original) {
    sys.source(testthat::test_path("fixtures", "bootstrap-before-retention.R"), envir = e)
  } else {
    for (name in c("bootstrap_inference", ".gl_bootstrap_capture_attempt", ".gl_bootstrap_call_attempt_callback")) {
      fun <- get(name, ns, inherits = FALSE); environment(fun) <- e; assign(name, fun, e)
    }
  }
  e$.gl_require_converged_fit <- function(...) invisible(NULL)
  e$.gl_validate_bootstrap_args <- function(object, R, fit_args) list(R = as.integer(R), fit_args = fit_args)
  e$.resolve_coefficient_terms <- function(terms, names, ...) seq_along(names)
  e$.fixed_term_coefficient_names <- function(object) list()
  e$.gl_bootstrap_summary <- function(estimates, terms_use, boot_coef, level) boot_coef
  e$.gl_bootstrap_result <- function(...) list(...)
  e$simulate <- function(object, nsim, ...) as.data.frame(matrix(seq_len(3L * nsim), 3L, nsim))
  e$calls <- 0L
  e$.gl_bootstrap_refit <- function(object, response, fit_args) {
    e$calls <- e$calls + 1L; refit(e$calls, response)
  }
  e$.gl_bootstrap_refit_dataset <- function(object, dat, fit_args) {
    e$calls <- e$calls + 1L; refit(e$calls, dat$response)
  }
  e$.gl_cluster_bootstrap_dataset <- function(object) {
    dat <- object$dataset; dat$response <- stats::runif(3); dat
  }
  e$object <- list(coefficients = c(mu = 2), response = c(1, NA, 3),
    dataset = data.frame(subject = c(1, 1, 2), time = c(1, 2, 1), response = c(1, NA, 3)))
  e
}

bootstrap_retention_native <- function(value = 1, converged = TRUE, name = "mu") {
  list(coefficients = setNames(value, name), convergence = list(converged = converged),
    complete_native_payload = list(untouched = "retained"))
}

test_that("all native outcomes reach callback before eligibility filtering", {
  original_error <- structure(list(message = "native failure", call = quote(native_refit()),
    private_detail = "preserved"), class = c("fixture_native_error", "error", "condition"))
  native <- list(bootstrap_retention_native(), bootstrap_retention_native(converged = FALSE),
    bootstrap_retention_native(name = "other"))
  refit <- function(i, response) if (i == 4L) stop(original_error) else native[[i]]
  e <- bootstrap_retention_fixture(refit); records <- list()
  out <- e$bootstrap_inference(e$object, R = 4, seed = 71, keep_fits = TRUE,
    attempt_callback = function(record) records[[length(records) + 1L]] <<- record)
  expect_identical(e$calls, 4L)
  expect_identical(vapply(records, `[[`, integer(1), "attempt"), 1:4)
  for (i in 1:3) expect_identical(records[[i]]$fit, native[[i]])
  expect_identical(records[[4]]$error, original_error)
  expect_null(records[[4]]$fit)
  expect_identical(records[[4]]$conditions[[1]]$condition, original_error)
  expect_equal(out$boot_coef[, 1], c(1, NA, NA, NA))
  expect_identical(out$fits[[1]], native[[1]])
  expect_true(all(vapply(out$fits[2:4], is.null, logical(1))))
  expect_equal(records[[1]]$data$response, c(1, NA, 3))
  expect_identical(records[[1]]$data$subject, e$object$dataset$subject)
  expect_identical(records[[1]]$seed, 71)
  expect_identical(records[[1]]$coefficient_names, "mu")
})

test_that("duplicate warnings and messages retain order and original conditions", {
  original_warning <- structure(list(message = "duplicate", call = quote(native_warn()), tag = 7),
    class = c("fixture_warning", "warning", "condition"))
  refit <- function(i, response) {
    warning(original_warning); message("native message"); warning(original_warning)
    bootstrap_retention_native()
  }
  e <- bootstrap_retention_fixture(refit); record <- NULL; outward <- character(); warnings <- character()
  withCallingHandlers(e$bootstrap_inference(e$object, R = 1,
    attempt_callback = function(x) record <<- x),
    message = function(m) { outward <<- c(outward, conditionMessage(m)); invokeRestart("muffleMessage") },
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") })
  expect_length(warnings, 0)
  expect_length(outward, 1)
  expect_identical(vapply(record$conditions, `[[`, character(1), "type"), c("warning", "message", "warning"))
  expect_identical(vapply(record$conditions, `[[`, integer(1), "sequence"), 1:3)
  expect_identical(record$conditions[[1]]$condition, original_warning)
  expect_identical(record$conditions[[3]]$condition, original_warning)
})

test_that("callback storage failure aborts before the next refit and restores RNG", {
  e <- bootstrap_retention_fixture(function(i, response) bootstrap_retention_native())
  failure <- structure(list(message = "storage verification failed", call = NULL),
    class = c("fixture_storage_error", "error", "condition"))
  set.seed(19); before <- .Random.seed
  seen <- NULL
  caught <- tryCatch(e$bootstrap_inference(e$object, R = 4, seed = 71,
    attempt_callback = function(record) { seen <<- record; stats::runif(4); stop(failure) }), error = identity)
  expect_identical(caught, failure)
  expect_identical(e$calls, 1L)
  expect_identical(seen$fit, bootstrap_retention_native())
  expect_identical(.Random.seed, before)
})

test_that("NULL callback matches the frozen previous loop and keep_fits semantics", {
  refit <- function(i, response) {
    warning("same suppressed warning")
    if (i == 3L) stop("same native failure")
    bootstrap_retention_native(i, converged = i != 2L)
  }
  old <- bootstrap_retention_fixture(refit, original = TRUE)
  current <- bootstrap_retention_fixture(refit)
  expect_identical(current$bootstrap_inference(current$object, R = 4, seed = 91, keep_fits = TRUE),
    old$bootstrap_inference(old$object, R = 4, seed = 91, keep_fits = TRUE))
  current <- bootstrap_retention_fixture(refit); observed <- bootstrap_retention_fixture(refit)
  expect_identical(current$bootstrap_inference(current$object, R = 4, seed = 91),
    observed$bootstrap_inference(observed$object, R = 4, seed = 91, attempt_callback = function(record) NULL))
  expect_identical(formals(current$bootstrap_inference)$R, 100)
})

test_that("cluster callback retains sampled data without perturbing subsequent draws", {
  refit <- function(i, response) bootstrap_retention_native(sum(response))
  a <- bootstrap_retention_fixture(refit); b <- bootstrap_retention_fixture(refit); records <- list()
  plain <- a$bootstrap_inference(a$object, R = 3, seed = 29, type = "cluster")
  observed <- b$bootstrap_inference(b$object, R = 3, seed = 29, type = "cluster",
    attempt_callback = function(record) {
      records[[length(records) + 1L]] <<- record; stats::runif(100)
    })
  expect_identical(plain, observed)
  expect_length(records, 3)
  expect_identical(colnames(observed$boot_coef), "mu")
  for (i in 1:3) {
    expect_identical(records[[i]]$type, "cluster")
    expect_equal(sum(records[[i]]$data$response), unname(observed$boot_coef[i, 1]))
    expect_identical(records[[i]]$rng_before_refit, records[[i]]$rng_after_refit)
  }
  expect_error(b$bootstrap_inference(b$object, R = 1, attempt_callback = TRUE), "NULL or a function")
  expect_identical(b$calls, 3L)
})
