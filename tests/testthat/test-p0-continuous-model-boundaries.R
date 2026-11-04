test_that("all families and backends use actual likelihood probability clamps", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (backend in c("native","vinecopula"))
  for (family in c("N","C","F","G","J","t"))
  for (residual in c(-6.6, -7.5)) {
    options(gamlss.longitudinal.copula_backend=backend)
    object <- curvature_model_fixture(family)
    object$reference_backend <- backend
    object$response[2L] <- 0.3 + 0.7 * object$dataset$x[2L] + exp(-0.2) * residual
    object$dataset$response <- object$response
    object$dataset_original <- object$dataset
    reference <- curvature_reference_hessian(object)
    actual <- calc_analytical_hessian(object, progress = FALSE)
    error <- sqrt(sum((actual-reference)^2))/sqrt(sum(reference^2))
    if (exists("validation_details", inherits = TRUE)) {
      validation_details$hessians[[paste(backend,family,residual,sep="/")]] <-
        list(backend=backend,family=family,residual=residual,relative_error=error,actual=actual,reference=reference)
    }
    expect_true(error < 2e-4, info=paste(backend,family,"standardized residual",residual,"relative error",error))
  }
})

test_that("native full Hessians keep the actual Clayton parameter", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  object <- curvature_model_fixture("C")
  object$par[["theta.intercept"]] <- log(30)
  reference <- curvature_reference_hessian(object)
  actual <- calc_analytical_hessian(object, progress = FALSE)
  error <- sqrt(sum((actual-reference)^2))/sqrt(sum(reference^2))
  if (exists("validation_details", inherits = TRUE)) {
    validation_details$hessians[["native-Clayton-30"]] <-
      list(relative_error=error,actual=actual,reference=reference)
  }
  expect_lt(error, 2e-4)
})

test_that("native full Hessians respect saturated correlation and df parameters", {
  old <- options(gamlss.longitudinal.copula_backend = "native")
  on.exit(options(old), add = TRUE)
  for (family in c("N","t")) {
    object <- curvature_model_fixture(family)
    target <- if (family=="N") "theta.intercept" else "zeta.intercept"
    object$par[[target]] <- if (family=="N") atanh(0.9999995) else log(1e-8)
    if (family=="N") {
      # Keep included densities representable near perfect correlation.
      object$response <- rep(object$response[seq(1L,length(object$response),3L)],each=3L)
      object$dataset$response <- object$response
      object$dataset_original <- object$dataset
    }
    reference <- curvature_reference_hessian(object)
    actual <- calc_analytical_hessian(object, progress = FALSE)
    index <- match(target,names(object$par))
    expect_equal(unname(reference[index,]),rep(0,length(object$par)),tolerance=1e-10)
    expect_equal(unname(actual[index,]),rep(0,length(object$par)),tolerance=1e-10)
  }
})

test_that("delegated out-of-domain Clayton parameters are never shifted", {
  old <- options(gamlss.longitudinal.copula_backend="vinecopula")
  on.exit(options(old),add=TRUE)
  cache <- list(row_id1=1L,row_id2=2L,theta_index_map=c(1L,NA_integer_),observed_pair_base=TRUE)
  eta <- list(mu=c(0,0),theta=30,margin_p_cache=c(0.2,0.4))
  input <- .copula_hessian_pair_inputs(eta,cache,"C")
  expect_identical(input$par1_e,30)
  expect_error(.copula_hessian_pair_derivatives(input$u1,input$u2,"C",
    input$par1_e,input$par2_e,input$pair_ok,FALSE),"Clayton")
})
