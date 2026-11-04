calc_deriv_copula_wrt_margin <- function(input, margin_par, par_name, calc_d2 = FALSE) {
  # Calculate copula derivative with respect to marginal parameters
  # input=copula_merged
  num_col <- function(nm) .copula_margin_derivative_numeric_column(input, nm)

  if (calc_d2 == FALSE) {
    return(.calc_deriv_copula_wrt_margin_d1(input, margin_par, par_name))
  } else {
    d2lcopdpar2 <- matrix(0, nrow = nrow(input), ncol = length(margin_par))
    i <- 1
    for (inner_par_name in margin_par) {
      if (inner_par_name == par_name) {
        # Take parameters from input for clarity
        dc_tplus_du_t <- num_col("dcdu1")
        dc_tplus_du_tplus <- num_col("dcdu2")
        # l_t=input[,paste(paste("dld",inner_par_name,sep=""),".x",sep="")]
        # l_t_plus=input[,paste(paste("dld",inner_par_name,sep=""),".y",sep="")]
        # x_t=input[,"response.x"]
        # x_t_plus=input[,"response.y"]
        # f_t=input[,"margin_d.x"]
        # f_t_plus=input[,"margin_d.y"]
        c_tplus <- num_col("copula_d")
        mu_t <- num_col("mu.x")
        mu_t_plus <- num_col("mu.y")

        F_nd_t <- num_col("F_nd.x")
        F_nd_t_plus <- num_col("F_nd.y")

        du_t_dmu <- F_nd_t
        du_t_plus_dmu <- F_nd_t_plus

        dc_plus_dt_dmu <- dc_tplus_du_t * du_t_dmu
        dc_plus_dt_plus_dmu <- dc_tplus_du_tplus * du_t_plus_dmu
        dc_plus_dt_dmu[is.nan(dc_plus_dt_dmu)] <- 0
        dc_plus_dt_plus_dmu[is.nan(dc_plus_dt_plus_dmu)] <- 0
        dcdmu_tplus <- ((dc_plus_dt_dmu + dc_plus_dt_plus_dmu) / c_tplus)
        dcdmu_tplus[is.nan(dcdmu_tplus) | is.na(dcdmu_tplus)] <- 0

        # dlcopdpar[,i]=dcdmu_tplus

        ####### NOW FOR SECOND DERIVATIVE OF COPULA TERM

        F_nd2 <- num_col("F_nd2.x")
        F_nd2_plus <- num_col("F_nd2.y")

        d2u_t_dmu2 <- F_nd2
        d2u_t_plus_dmu2 <- F_nd2_plus

        d2cdu_t2 <- num_col("d2cdu12")
        d2cdu_t_plus2 <- num_col("d2cdu22")
        d2cdu_t2[is.nan(d2cdu_t2)] <- 0
        d2cdu_t_plus2[is.nan(d2cdu_t_plus2)] <- 0

        d2cdmu2 <- d2cdu_t2 * du_t_dmu^2 +
          dc_tplus_du_t * d2u_t_dmu2 +
          d2cdu_t_plus2 * du_t_plus_dmu^2 +
          dc_tplus_du_tplus * d2u_t_plus_dmu2

        d2lcdmu2 <- as.matrix((d2cdmu2 * c_tplus - (dcdmu_tplus^2)) / (c_tplus^2))
        d2lcdmu2 <- num_col("c_nd2")

        d2lcopdpar2[, i] <- d2lcdmu2
      }
      i <- i + 1
    }
    colnames(d2lcopdpar2) <- paste("d2lcopd", margin_par, sep = "")

    par_d2lcopdpar <- d2lcopdpar2[, paste("d2lcopd", margin_par, sep = "")]
    merged_d2lcopdpar <- merge(cbind(input[, c("time1", "time2", "subject1", "subject2")], par_d2lcopdpar),
      cbind(input[, c("time1", "time2", "subject1", "subject2")], par_d2lcopdpar),
      by.x = c("time2", "subject2"), by.y = c("time1", "subject1"), all = TRUE
    )
    merged_d2lcopdpar[is.na(merged_d2lcopdpar)] <- 0

    x_comp <- grepl("d2lcopd", colnames(merged_d2lcopdpar)) & grepl(".x", colnames(merged_d2lcopdpar))
    y_comp <- grepl("d2lcopd", colnames(merged_d2lcopdpar)) & grepl(".y", colnames(merged_d2lcopdpar))

    d2_cop <- 0.5 * (merged_d2lcopdpar[, x_comp] + merged_d2lcopdpar[, y_comp])

    # plot(d2lcopdpar2[,paste("d2lcopd",par_name,sep="")],input[,"c_nd2"],main="d2",ylab="numerical")

    return(d2_cop)
  }
}
