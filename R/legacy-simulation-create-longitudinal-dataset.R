#' Legacy simulation and research-era helpers
#'
#' These helpers are retained for compatibility and internal historical support.
#' They are not part of the core reviewer path for fitting, prediction, or diagnostics.
#'
#' @noRd
NULL

#' @keywords internal
#' @noRd
create_longitudinal_dataset <- function(response, covariates, labels = NA) {
  num_time_points <- ncol(response)
  if (num_time_points <= 1) {
    print("Not enough time points")
  }

  dataset <- matrix(data = NA, ncol = 2 + length(covariates), nrow = 0)
  subject <- as.factor(seq(1:nrow(response)))

  for (t in 1:ncol(response)) {
    dataset_temp <- cbind(subject, t, response[, t])

    for (i in 1:length(covariates)) {
      if (ncol(covariates[[i]]) == 1) {
        covariate_for_time <- covariates[[i]]
      } else {
        covariate_for_time <- covariates[[i]][, t]
      }
      dataset_temp <- cbind(dataset_temp, covariate_for_time)
    }

    ### Add dataset temp to full table
    dataset <- rbind(dataset, dataset_temp)
  }

  if (!all(is.na(labels))) {
    colnames(dataset) <- labels
  }

  dataset <- dataset[order(dataset$time, dataset$subject), ] ### NOTE THIS WILL BREAK GLMM
  rownames(dataset) <- 1:nrow(dataset)

  return(dataset)
}
