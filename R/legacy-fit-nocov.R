#' @keywords internal

#' @noRd

fit_jointreg_nocov <- function(
    input_par, margin_dist, copula_dist, data,
    use_dlcopdpar = TRUE, verbose = TRUE, plot_results = TRUE,
    crit_lik_change = 0.05, start_step_size = .5, step_adjustment = .9, max_steps = 5,
    true_val = NA) {
  log_lik_history <- matrix(ncol = 3 + 2, nrow = 0)

  par_history <- matrix(ncol = length(input_par) + 2, nrow = 0)

  ### Run fit for separate and joint optimisation

  copula_deriv <- if (use_dlcopdpar == TRUE) {
    1
  } else {
    0
  }

  ### CORE ITERATION

  change <- 1
  log_lik_start <- 0
  log_lik_change <- 1000
  run_counter <- 1
  step_size <- start_step_size

  while (abs(log_lik_change) > crit_lik_change) {
    step_size <- step_size * (step_adjustment^min(max_steps, run_counter))

    par_history <- rbind(par_history, c(copula_deriv, run_counter, input_par))

    # Run optimisation

    outer_optim_output <- optim_outer(par = input_par, dataset, margin_dist, copula_dist, use_dlcopdpar = use_dlcopdpar, verbose = FALSE, step_size = step_size)

    # Capture outputs

    input_par <- outer_optim_output$par_end

    change <- sum(outer_optim_output$par_change)

    # print(outer_optim_output$log_lik)

    log_lik <- outer_optim_output$log_lik["joint"]

    log_lik_change <- log_lik - log_lik_start

    log_lik_start <- log_lik

    # Capture changes in parameters

    log_lik_history <- rbind(log_lik_history, c(copula_deriv, run_counter, outer_optim_output$log_lik))

    run_counter <- run_counter + 1
  }

  par_history <- rbind(par_history, c(copula_deriv, run_counter, input_par))

  outer_optim_output <- optim_outer(par = input_par, dataset, margin_dist, copula_dist, use_dlcopdpar = use_dlcopdpar, verbose = FALSE, step_size = step_size)

  log_lik <- outer_optim_output$log_lik["joint"]

  log_lik_change <- log_lik - log_lik_start

  log_lik_start <- log_lik

  log_lik_history <- rbind(log_lik_history, c(copula_deriv, run_counter, outer_optim_output$log_lik))

  colnames(log_lik_history)[1:2] <- colnames(par_history)[1:2] <- c("use_dlcopdpar", "run_counter")

  # Plot likelihood and parameters

  if (plot_results == TRUE) {
    plot.new()

    par_count <- round(sqrt((ncol(par_history) + 1)), 0) + 1

    par(mfrow = c(par_count, par_count))

    for (i in colnames(log_lik_history)[3:5]) {
      plot(log_lik_history[, i], xlab = "Iteration", ylab = "LogLik", main = i, type = "l", col = "blue", xlim = c(1, max(log_lik_history[, "run_counter"])))

      # lines(log_lik_history[log_lik_history[,"use_dlcopdpar"]==0,i],xlab="LogLik",ylab="Iteration",main=i,type = "l",col="red",xlim=c(1,max(log_lik_history[,"run_counter"])),ylim=range(log_lik_history[,i]))

      # legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    }

    for (i in 1:(ncol(par_history) - 2)) {
      # lines(par_nodlcop[,i+2],col="red",type="l")

      if (!all(is.na(true_val))) {
        plot(par_history[, i + 2], col = "blue", type = "l", main = colnames(par_history)[i + 2], ylab = "Parameter estimate", ylim = range(c(par_history[, i + 2], true_val[i])))

        abline(h = true_val[i])
      } else {
        plot(par_history[, i + 2], col = "blue", type = "l", main = colnames(par_history)[i + 2], ylab = "Parameter estimate")
      }

      # legend("bottomright",c("Joint","Separate"), lwd=c(5,2), col=c("blue","red"))
    }
  }

  return_list <- list(par_history, log_lik_history)

  names(return_list) <- c("par_history", "log_lik_history")

  return(return_list)
}
