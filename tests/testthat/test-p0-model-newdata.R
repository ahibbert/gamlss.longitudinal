make_newdata_policy_object <- function() {
  list(
    formulas_int = list(mu = response ~ 1),
    var_map = list(y = "response", visit = "time", id = "subject"),
    model_matrix = list(x = list(mu = list(time = 1))),
    dataset = data.frame(
      group = factor(c("control", "active"), levels = c("control", "active"), ordered = TRUE),
      stringsAsFactors = FALSE
    )
  )
}

test_that("newdata name translation and defaults preserve internal prediction columns", {
  object <- make_newdata_policy_object()
  nd <- data.frame(
    y = c(1, 2),
    visit = c("t1", "t2"),
    id = c(10, 10),
    group = c("active", "control"),
    stringsAsFactors = FALSE
  )

  out <- gamlss.longitudinal:::.gl_prepare_newdata_internal(object, nd, require_response = TRUE)

  expect_true(all(c("response", "time", "subject", "time_covariate") %in% names(out)))
  expect_equal(out$response, c(1, 2))
  expect_equal(out$time, c("t1", "t2"))
  expect_equal(out$subject, c(10, 10))
  expect_equal(out$time_covariate, out$time)
  expect_true(is.factor(out$group))
  expect_true(is.ordered(out$group))
  expect_equal(levels(out$group), c("control", "active"))
})

test_that("newdata helpers reject unseen factor levels and missing required responses", {
  object <- make_newdata_policy_object()
  nd_bad <- data.frame(
    visit = "t1",
    id = 1,
    group = "other",
    stringsAsFactors = FALSE
  )

  expect_error(
    gamlss.longitudinal:::.gl_prepare_newdata_internal(object, nd_bad, require_response = FALSE),
    "newdata column 'group' contains level(s) not seen during fitting: other",
    fixed = TRUE
  )

  nd_no_response <- data.frame(
    visit = "t1",
    id = 1,
    group = "active",
    stringsAsFactors = FALSE
  )

  expect_error(
    gamlss.longitudinal:::.gl_prepare_newdata_internal(object, nd_no_response, require_response = TRUE),
    "newdata must include a response column (or mapped response variable) for this operation.",
    fixed = TRUE
  )
})

test_that("model-matrix alignment preserves reference column order and fills missing columns", {
  mm_use <- list(
    x = list(
      mu = data.frame(intercept = c(1, 1), extra = c(9, 9), check.names = FALSE),
      sigma = data.frame(intercept = c(1, 1), x = c(2, 3), check.names = FALSE)
    )
  )
  mm_reference <- list(
    x = list(
      mu = data.frame(intercept = c(1, 1), x = c(0, 0), check.names = FALSE),
      sigma = data.frame(x = c(0, 0), intercept = c(1, 1), check.names = FALSE)
    )
  )

  out <- gamlss.longitudinal:::.gl_align_model_matrix_columns(mm_use, mm_reference)

  expect_equal(colnames(out$x$mu), c("intercept", "x"))
  expect_equal(out$x$mu$x, c(0, 0))
  expect_equal(colnames(out$x$sigma), c("x", "intercept"))
  expect_equal(out$x$sigma$x, c(2, 3))
})
