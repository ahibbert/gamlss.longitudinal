#region Setup and Libraries
source("R/common_functions.R");source("R/link_functions.R"); library("VineCopula");library("moments"); 
library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix"); library("MASS")
library(gamlss2); library(gamlss)
set.seed(100)
#endregion

#region Simulation Configuration
########### Generate Dataset ###########
n=500; d=4

# Missingness configuration:
# - "increasing_time": p_miss(i) = i/(T+1) by ordered time index i.
# - "mar": missing completely at random across all rows at mar_missing_rate.
# mar_missing_rate controls the overall target missingness level for both modes.
missingness_mode = "mar" # or "mar"
mar_missing_rate = 0.1

#copula_dist="t"; margin_dist=GG(); mu=2; sigma=-2; nu=1; theta=1; tau=2; zeta=2.5; simOption=10
copula_dist="t"; margin_dist=BCPE(); mu=1; sigma=0.3; nu=2; theta=0; tau=0.5; zeta=0.5; simOption=10

# Parameterization control for simulation starts:
# - TRUE: intercept values above are interpreted on eta (linear predictor) scale.
# - FALSE: values are interpreted on natural parameter scale.
margin_inputs_on_eta_scale <- TRUE
copula_inputs_on_eta_scale <- TRUE
#endregion

#region Parameter Guards and Starts
# Reproducibility guard: validate simulation inputs explicitly so stale
# workspace variables cannot leak into loadDataset() calls.
if (!exists("simOption", inherits = FALSE) || length(simOption) != 1 || !is.numeric(simOption) || !is.finite(simOption)) {
  stop("simOption must be defined as one finite numeric scalar before simulation.")
}

required_margin_par <- names(margin_dist$parameters)
if (length(required_margin_par) == 0) {
  stop("margin_dist does not expose any parameters; cannot build par.margin.")
}

margin_start <- vapply(required_margin_par, function(pn) {
  if (!exists(pn, envir = .GlobalEnv, inherits = FALSE)) {
    stop("Missing required margin parameter variable: '", pn, "'.")
  }
  val <- get(pn, envir = .GlobalEnv, inherits = FALSE)
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
    stop("Margin parameter '", pn, "' must be one finite numeric scalar.")
  }
  if (isTRUE(margin_inputs_on_eta_scale)) {
    linkinv_name <- paste0(pn, ".linkinv")
    if (!(linkinv_name %in% names(margin_dist)) || !is.function(margin_dist[[linkinv_name]])) {
      stop("Missing link-inverse function for margin parameter '", pn, "'.")
    }
    natural_val <- suppressWarnings(margin_dist[[linkinv_name]](as.numeric(val)))
    if (!is.numeric(natural_val) || length(natural_val) != 1 || !is.finite(natural_val)) {
      stop(
        "Margin eta parameter '", pn, "' = ", as.numeric(val),
        " does not map to a finite natural parameter for family ", margin_dist$family[1],
        "."
      )
    }
    return(as.numeric(natural_val))
  }

  linkfun_name <- paste0(pn, ".linkfun")
  if (linkfun_name %in% names(margin_dist) && is.function(margin_dist[[linkfun_name]])) {
    eta_val <- suppressWarnings(margin_dist[[linkfun_name]](as.numeric(val)))
    if (!is.numeric(eta_val) || length(eta_val) != 1 || !is.finite(eta_val)) {
      stop(
        "Margin parameter '", pn, "' = ", as.numeric(val),
        " is outside the link-function domain for family ", margin_dist$family[1],
        ". Provide a valid starting value."
      )
    }
  }
  as.numeric(val)
}, numeric(1))

copula_spec <- get_copula_dist(copula_dist)
required_copula_par <- copula_spec$parameters
if (length(required_copula_par) == 0) {
  stop("copula_dist does not expose any parameters; cannot build par.copula.")
}

copula_start <- vapply(required_copula_par, function(pn) {
  if (!exists(pn, envir = .GlobalEnv, inherits = FALSE)) {
    stop("Missing required copula parameter variable: '", pn, "'.")
  }
  val <- get(pn, envir = .GlobalEnv, inherits = FALSE)
  if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
    stop("Copula parameter '", pn, "' must be one finite numeric scalar.")
  }
  if (isTRUE(copula_inputs_on_eta_scale)) {
    linkinv_name <- paste0(pn, ".linkinv")
    if (!(linkinv_name %in% names(copula_spec$copula_link)) || !is.function(copula_spec$copula_link[[linkinv_name]])) {
      stop("Missing link-inverse function for copula parameter '", pn, "'.")
    }
    natural_val <- suppressWarnings(copula_spec$copula_link[[linkinv_name]](as.numeric(val)))
    if (!is.numeric(natural_val) || length(natural_val) != 1 || !is.finite(natural_val)) {
      stop(
        "Copula eta parameter '", pn, "' = ", as.numeric(val),
        " does not map to a finite natural parameter for copula ", copula_dist,
        "."
      )
    }
    return(as.numeric(natural_val))
  }
  linkfun_name <- paste0(pn, ".linkfun")
  if (linkfun_name %in% names(copula_spec$copula_link) && is.function(copula_spec$copula_link[[linkfun_name]])) {
    eta_val <- suppressWarnings(copula_spec$copula_link[[linkfun_name]](as.numeric(val)))
    if (!is.numeric(eta_val) || length(eta_val) != 1 || !is.finite(eta_val)) {
      stop(
        "Copula parameter '", pn, "' = ", as.numeric(val),
        " is outside the link-function domain for copula ", copula_dist,
        ". Provide a valid starting value."
      )
    }
  }
  as.numeric(val)
}, numeric(1))
#endregion

#region Covariate Effects
# USE THIS WITH SIMOPTION 10
covariates_input=list( mu.time=.1   ,sigma.time=.2   ,nu.time=0    ,tau.time=0   ,theta.time=0.5  ,zeta.time=1
                        ,mu.age=0.5    ,sigma.age=0.5     ,nu.age=0     ,tau.age=0    ,theta.age=0    ,zeta.age=0
                        ,mu.gender=1 ,sigma.gender=1  ,nu.gender=0  ,tau.gender=0 ,theta.gender=.2 ,zeta.gender=0)
#endregion

#region Generate Dataset
#########Generate dataset

if (exists("dataset", inherits = FALSE)) rm(dataset)
# Safety fallback: if this section is run independently (without running the
# guard block above), rebuild validated parameter vectors on the fly.
if (!exists("margin_start", inherits = FALSE) || !exists("copula_start", inherits = FALSE)) {
  required_margin_par <- names(margin_dist$parameters)
  margin_start <- vapply(required_margin_par, function(pn) {
    if (!exists(pn, envir = .GlobalEnv, inherits = FALSE)) {
      stop("Missing required margin parameter variable: '", pn, "'.")
    }
    val <- get(pn, envir = .GlobalEnv, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop("Margin parameter '", pn, "' must be one finite numeric scalar.")
    }
    if (isTRUE(margin_inputs_on_eta_scale)) {
      linkinv_name <- paste0(pn, ".linkinv")
      if (!(linkinv_name %in% names(margin_dist)) || !is.function(margin_dist[[linkinv_name]])) {
        stop("Missing link-inverse function for margin parameter '", pn, "'.")
      }
      natural_val <- suppressWarnings(margin_dist[[linkinv_name]](as.numeric(val)))
      if (!is.numeric(natural_val) || length(natural_val) != 1 || !is.finite(natural_val)) {
        stop(
          "Margin eta parameter '", pn, "' = ", as.numeric(val),
          " does not map to a finite natural parameter for family ", margin_dist$family[1],
          "."
        )
      }
      return(as.numeric(natural_val))
    }

    linkfun_name <- paste0(pn, ".linkfun")
    if (linkfun_name %in% names(margin_dist) && is.function(margin_dist[[linkfun_name]])) {
      eta_val <- suppressWarnings(margin_dist[[linkfun_name]](as.numeric(val)))
      if (!is.numeric(eta_val) || length(eta_val) != 1 || !is.finite(eta_val)) {
        stop(
          "Margin parameter '", pn, "' = ", as.numeric(val),
          " is outside the link-function domain for family ", margin_dist$family[1],
          ". Provide a valid starting value."
        )
      }
    }
    as.numeric(val)
  }, numeric(1))

  copula_spec <- get_copula_dist(copula_dist)
  required_copula_par <- copula_spec$parameters
  copula_start <- vapply(required_copula_par, function(pn) {
    if (!exists(pn, envir = .GlobalEnv, inherits = FALSE)) {
      stop("Missing required copula parameter variable: '", pn, "'.")
    }
    val <- get(pn, envir = .GlobalEnv, inherits = FALSE)
    if (!is.numeric(val) || length(val) != 1 || !is.finite(val)) {
      stop("Copula parameter '", pn, "' must be one finite numeric scalar.")
    }
    if (isTRUE(copula_inputs_on_eta_scale)) {
      linkinv_name <- paste0(pn, ".linkinv")
      if (!(linkinv_name %in% names(copula_spec$copula_link)) || !is.function(copula_spec$copula_link[[linkinv_name]])) {
        stop("Missing link-inverse function for copula parameter '", pn, "'.")
      }
      natural_val <- suppressWarnings(copula_spec$copula_link[[linkinv_name]](as.numeric(val)))
      if (!is.numeric(natural_val) || length(natural_val) != 1 || !is.finite(natural_val)) {
        stop(
          "Copula eta parameter '", pn, "' = ", as.numeric(val),
          " does not map to a finite natural parameter for copula ", copula_dist,
          "."
        )
      }
      return(as.numeric(natural_val))
    }
    linkfun_name <- paste0(pn, ".linkfun")
    if (linkfun_name %in% names(copula_spec$copula_link) && is.function(copula_spec$copula_link[[linkfun_name]])) {
      eta_val <- suppressWarnings(copula_spec$copula_link[[linkfun_name]](as.numeric(val)))
      if (!is.numeric(eta_val) || length(eta_val) != 1 || !is.finite(eta_val)) {
        stop(
          "Copula parameter '", pn, "' = ", as.numeric(val),
          " is outside the link-function domain for copula ", copula_dist,
          ". Provide a valid starting value."
        )
      }
    }
    as.numeric(val)
  }, numeric(1))
}
dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                    , par.margin=margin_start, par.copula=copula_start,covariates_input=covariates_input)
#endregion

#region Inject Missingness
# Inject missingness into response using selected mode.
if (missingness_mode == "increasing_time") {
  time_points_missing = sort(unique(dataset$time))
  T_missing = length(time_points_missing)

  # Base increasing probabilities by time rank.
  base_p_by_time = seq_along(time_points_missing) / (T_missing + 1)

  # Scale to target mar_missing_rate on average while preserving the trend.
  n_by_time = sapply(time_points_missing, function(t_val) sum(dataset$time == t_val))
  base_mean_p = sum(base_p_by_time * n_by_time) / sum(n_by_time)
  scale_factor = ifelse(base_mean_p > 0, mar_missing_rate / base_mean_p, 1)
  p_by_time = pmin(base_p_by_time * scale_factor, 1)

  for (i in seq_along(time_points_missing)) {
    t_val = time_points_missing[i]
    p_miss = p_by_time[i]
    idx_t = which(dataset$time == t_val)
    miss_flags = runif(length(idx_t)) < p_miss
    dataset$response[idx_t[miss_flags]] = NA
  }
} else if (missingness_mode == "mar") {
  miss_flags = runif(nrow(dataset)) < mar_missing_rate
  dataset$response[miss_flags] = NA
} else {
  stop("Invalid missingness_mode. Use 'increasing_time' or 'mar'.")
}
#endregion

#region Exploratory Data Plots
plotDist(dataset, margin_dist, offdiag_scale = "pseudo")

########### PLOTTING ###########
#Group dataset by age categories in buckets of 10 years then calculate the correlation in each bucket and plot
library(dplyr)
# Create age groups with extended range and filter out any remaining NAs
age_range = range(dataset$age, na.rm=TRUE)
age_breaks = seq(floor(age_range[1]/10)*10, ceiling(age_range[2]/10)*10 + 10, by=10)
dataset=dataset%>%
  mutate(age_group=cut(age, breaks=age_breaks, include.lowest=TRUE)) %>%
  filter(!is.na(age_group))

# Calculate correlations between consecutive time points
time_points = sort(unique(dataset$time))
age_corrs_list = list()

for(i in 1:(length(time_points)-1)) {
  t1 = time_points[i]
  t2 = time_points[i+1]
  temp_corrs = dataset %>% 
    group_by(age_group) %>% 
    summarise(cor = cor(response[time==t1], response[time==t2], use = "pairwise.complete.obs"), .groups='drop') %>%
    mutate(time_pair = paste0("T", t1, " vs T", t2)) %>%
    filter(!is.na(cor))
  age_corrs_list[[i]] = temp_corrs
}

age_corrs = bind_rows(age_corrs_list) %>% filter(!is.na(age_group))

p1 = ggplot(age_corrs, aes(x=age_group, y=cor, color=time_pair, group=time_pair)) + 
  geom_point() + 
  geom_line() + 
  labs(title="Correlation by Age Group", 
       x="Age Group", 
       y="Correlation",
       color="Time Comparison") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate mean response by age group and time
age_response = dataset %>%
  group_by(age_group, time) %>%
  summarise(mean_response = mean(response, na.rm = TRUE), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p2 = ggplot(age_response, aes(x=age_group, y=mean_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Mean Response by Age Group",
       x="Age Group",
       y="Mean Response",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate standard deviation by age group and time
age_sd = dataset %>%
  group_by(age_group, time) %>%
  summarise(sd_response = sd(response, na.rm = TRUE), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p3 = ggplot(age_sd, aes(x=age_group, y=sd_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="SD of Response by Age Group",
       x="Age Group",
       y="Standard Deviation",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate kurtosis by age group and time
library(moments)
age_kurtosis = dataset %>%
  group_by(age_group, time) %>%
  summarise(kurtosis_response = kurtosis(response, na.rm = TRUE), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p4 = ggplot(age_kurtosis, aes(x=age_group, y=kurtosis_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Kurtosis of Response by Age Group",
       x="Age Group",
       y="Kurtosis",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Calculate skewness by age group and time
age_skewness = dataset %>%
  group_by(age_group, time) %>%
  summarise(skewness_response = skewness(response, na.rm = TRUE), .groups='drop') %>%
  mutate(time = as.factor(time)) %>%
  filter(!is.na(age_group))

p5 = ggplot(age_skewness, aes(x=age_group, y=skewness_response, color=time, group=time)) +
  geom_point() +
  geom_line() +
  labs(title="Skewness of Response by Age Group",
       x="Age Group",
       y="Skewness",
       color="Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Combine all five plots
library(ggpubr)
ggarrange(p1, p2, p3, p4, p5, ncol=2, nrow=3, common.legend=FALSE)
#endregion

#region Model Fit
########## FIT ###########
source("R/common_functions.R")
mu_formula="random_name ~ time_of_observation_random_name + s(age_new_name) + gender"
sigma_formula="~ time_of_observation_random_name + s(age_new_name) + gender"
nu_formula="~ 1"
tau_formula="~ 1"
theta_formula="~ time_of_observation_random_name + gender"
zeta_formula="~ 1"

### FOR TESTING DATASETS WITH NON STANDARD NAMING
if("age_group" %in% names(dataset)) dataset$age_group = NULL
colnames(dataset)=c("person","time_of_observation_random_name","random_name","age_new_name","year","gender")
data_in=dataset
data_in$time_of_observation_random_name = factor(data_in$time_of_observation_random_name)
data_in$gender = factor(data_in$gender)
rm(dataset)

options(scipen=999)
source("R/common_functions.R")
source("R/link_functions.R")
fit=gamlss.longitudinal(dataset = data_in
                   , margin_dist = margin_dist
                   , copula_dist = copula_dist
                   , time_var = "time_of_observation_random_name"
                   , subject_var = "person"
                   , mu.formula = mu_formula
                   , sigma.formula = sigma_formula
                   , nu.formula = nu_formula
                   , tau.formula = tau_formula
                   , theta.formula = theta_formula
                   , zeta.formula = zeta_formula
                   , verbose = 3
                   , compute_vcov = TRUE
                   , include_dlcopdpar=TRUE
                   , vcov_method="analytical"
                   , method="RS"
)
#endregion

#region Model Diagnostics
#################### PLOT METHOD ####################
source("R/common_functions.R")
source("R/diagnostics_topmodels.R")
source("R (testing)/plot_copula_v2.R")
summary(fit)
plot(fit)
plot.terms(fit,  data = data_in)
wormplot(fit, by="time_of_observation_random_name", data=data_in)
wormplot(fit, by="gender", data=data_in)
plot.copula(fit, data=data_in)
#endregion

#region Hessian Method Comparison
# Compare the analytical Hessian against the numerical ground truth.
# Run the model fit block first, then this block.
source("R/common_functions.R")
source("R/link_functions.R")
source("R/analytical_hessian.R")

# Optional diagnostics toggles for this section
run_hessian_plots <- TRUE
hessian_block_tol <- 0.05

cat("=== Hessian Method Comparison ===\n\n")

# --- Numerical (source of truth) ---
cat("Computing NUMERICAL Hessian...\n")
t_num <- system.time(
  vcov_num <- vcov(fit, method = "numderiv", progress = TRUE)
)
cat(sprintf("Numerical time: %.1f sec\n\n", t_num["elapsed"]))

# --- Analytical ---
cat("Computing ANALYTICAL Hessian...\n")
source("R/analytical_hessian.R")
t_ana <- system.time(
  vcov_ana <- vcov(fit, method = "analytical", progress = TRUE)
)
cat(sprintf("Analytical time: %.1f sec\n\n", t_ana["elapsed"]))

cat(sprintf("Speed-up: %.1fx\n\n", t_num["elapsed"] / t_ana["elapsed"]))

# --- Comparison ---
H_num <- -solve(vcov_num$vcov$overall)
H_ana <- -solve(vcov_ana$vcov$overall)

se_num <- sqrt(diag(vcov_num$vcov$overall))
se_ana <- sqrt(diag(vcov_ana$vcov$overall))

comparison_df <- data.frame(
  parameter = names(fit$par),
  se_numerical  = se_num,
  se_analytical = se_ana,
  se_rel_diff_pct = round(100 * (se_ana - se_num) / abs(se_num), 2)
)
cat("SE comparison (numerical vs analytical):\n")
print(comparison_df, row.names = FALSE)

# Hessian element-wise comparison
hess_rel_err <- abs(H_ana - H_num) / (abs(H_num) + 1e-12)
cat(sprintf("\nHessian element-wise relative error:\n"))
cat(sprintf("  Max:    %.4f\n", max(hess_rel_err, na.rm = TRUE)))
cat(sprintf("  Median: %.6f\n", median(hess_rel_err, na.rm = TRUE)))
cat(sprintf("  Mean:   %.6f\n", mean(hess_rel_err, na.rm = TRUE)))

# Block-level drift checks (more stable than element-wise max near tiny denominators)
block_rel_frob <- function(A, B) {
  num <- sqrt(sum((A - B)^2, na.rm = TRUE))
  den <- sqrt(sum(B^2, na.rm = TRUE))
  num / (den + 1e-12)
}

theta_rows <- grepl("theta", rownames(H_ana))
mu_rows    <- grepl("^mu\\.", rownames(H_ana))
sigma_rows <- grepl("^sigma\\.", rownames(H_ana))
nu_rows    <- grepl("^nu\\.", rownames(H_ana))

block_checks <- data.frame(
  block = c("theta-theta", "mu-mu", "sigma-sigma", "nu-nu", "mu-sigma"),
  rel_frob = c(
    block_rel_frob(H_ana[theta_rows, theta_rows, drop = FALSE], H_num[theta_rows, theta_rows, drop = FALSE]),
    block_rel_frob(H_ana[mu_rows, mu_rows, drop = FALSE], H_num[mu_rows, mu_rows, drop = FALSE]),
    block_rel_frob(H_ana[sigma_rows, sigma_rows, drop = FALSE], H_num[sigma_rows, sigma_rows, drop = FALSE]),
    block_rel_frob(H_ana[nu_rows, nu_rows, drop = FALSE], H_num[nu_rows, nu_rows, drop = FALSE]),
    block_rel_frob(H_ana[mu_rows, sigma_rows, drop = FALSE], H_num[mu_rows, sigma_rows, drop = FALSE])
  )
)
cat("\nBlock drift check (relative Frobenius norm):\n")
print(block_checks, row.names = FALSE)

bad_blocks <- block_checks$block[is.finite(block_checks$rel_frob) & block_checks$rel_frob > hessian_block_tol]
if (length(bad_blocks) > 0) {
  warning(
    "Hessian block drift above tolerance (",
    hessian_block_tol,
    "): ",
    paste(bad_blocks, collapse = ", "),
    call. = FALSE
  )
}

# Visual: scatter plot of Hessian diagonal
if (isTRUE(run_hessian_plots)) {
  par(mfrow = c(1, 2))
  diag_num <- diag(H_num); diag_ana <- diag(H_ana)
  plot(diag_num, diag_ana,
       xlab = "Numerical H diag", ylab = "Analytical H diag",
       main = "Hessian diagonal: Numerical vs Analytical",
       pch = 19, col = "steelblue")
  abline(0, 1, col = "red", lty = 2)

  # Off-diagonal scatter (upper triangle)
  idx_upper <- which(upper.tri(H_num))
  plot(H_num[idx_upper], H_ana[idx_upper],
       xlab = "Numerical H off-diag", ylab = "Analytical H off-diag",
       main = "Hessian off-diagonal: Numerical vs Analytical",
       pch = 19, col = "darkorange", cex = 0.5)
  abline(0, 1, col = "red", lty = 2)
  par(mfrow = c(1, 1))

  # Full matrix heatmap of relative error
  image(hess_rel_err, main = "Hessian relative error |H_ana - H_num| / |H_num|",
        xlab = "Parameter index", ylab = "Parameter index",
        col = heat.colors(50))
}

cat("\n=== Hessian Diagonal Ratio (analytical / numerical) ===\n")
diag_df <- data.frame(
  par = rownames(H_num),
  se_num  = round(sqrt(diag(vcov_num$vcov$overall)), 5),
  se_ana  = round(sqrt(diag(vcov_ana$vcov$overall)), 5),
  H_ratio = round(diag(H_ana) / diag(H_num), 3)
)
print(diag_df, row.names = FALSE)

cat("\n=== Off-diagonal block checks ===\n")

cat("Theta x theta block:\n")
print(round(H_ana[theta_rows, theta_rows], 2))
cat("Numerical:\n")
print(round(H_num[theta_rows, theta_rows], 2))

cat("Mu x sigma cross block (first 4x4):\n")
nr <- min(4, sum(mu_rows)); nc <- min(4, sum(sigma_rows))
print(round(H_ana[mu_rows, sigma_rows][1:nr, 1:nc], 4))
cat("Numerical:\n")
print(round(H_num[mu_rows, sigma_rows][1:nr, 1:nc], 4))
#endregion

#region Debug sigma copula contribution
run_sigma_debug <- FALSE
if (isTRUE(run_sigma_debug)) {
source("R/common_functions.R")
source("R/link_functions.R")
source("R/analytical_hessian.R")

mm2 <- fit$model_matrix
copula_link2 <- get_copula_dist(fit$copula_dist)$copula_link
eta_out2 <- calc_eta(fit$par, mm2, fit$margin_dist, copula_link2, fit$par_s)
pair_cache2 <- build_copula_pair_cache(fit$response, fit$response_margin, fit$response_subject)
calc_lik2 <- calc_likelihood_minimal(
  eta_inv = eta_out2$eta_inv, mm = mm2$x,
  margin_dist = fit$margin_dist, copula_dist = fit$copula_dist,
  calc_d2 = FALSE, response = fit$response,
  response_margin = fit$response_margin, response_subject = fit$response_subject,
  pair_cache = pair_cache2
)
eta_inv2 <- eta_out2$eta_inv
eta_inv2[["margin_p_cache"]] <- calc_lik2$margin_p
dF2 <- .calc_dFdpar(eta_inv2, mm2, fit$margin_dist, fit$response, h = 1e-4)
d2F2 <- .calc_d2Fdpar2(eta_inv2, mm2, fit$margin_dist, fit$response, h = 1e-4)
d2Fx2 <- .calc_d2Fdpar_cross(eta_inv2, mm2, fit$margin_dist, fit$response, h = 1e-4)
cop2 <- suppressWarnings(.calc_copula_hessian_contributions(
  eta_inv2, pair_cache2, fit$copula_dist, dF2, d2F2, d2Fx2
))

# For a single pair, verify
k <- which(cop2$pair_ok)[1]
i1 <- cop2$row_id1[k]; i2 <- cop2$row_id2[k]
cat("Pair", k, ": obs i1=", i1, "i2=", i2, "\n")

# Check copula family number
fam_num <- as.numeric(VineCopula::BiCopName(fit$copula_dist))
cat("Copula fam_num:", fam_num, "\n")
th_row <- pair_cache2$theta_index_map[i1]
cat("th_row:", th_row, "\n")
theta_val <- copula_link2$theta.linkinv(as.numeric(eta_out2$eta_inv$theta[th_row]))
cat("theta_val:", theta_val, "\n")

u1_v <- calc_lik2$margin_p[i1]; u2_v <- calc_lik2$margin_p[i2]
s1_0 <- eta_inv2$sigma[i1]; s2_0 <- eta_inv2$sigma[i2]
cat("sigma i1:", s1_0, "sigma i2:", s2_0, "\n")

h_s <- 1e-4
pfun <- eval(parse(text = paste0("p", fit$margin_dist$family[1])))
y1 <- fit$response[i1]; y2 <- fit$response[i2]
mu1 <- eta_inv2$mu[i1]; mu2 <- eta_inv2$mu[i2]
nu1 <- eta_inv2$nu[i1]; nu2 <- eta_inv2$nu[i2]

logc_s <- function(s1, s2) {
  u1 <- pmax(pmin(pfun(y1, mu=mu1, sigma=s1, nu=nu1), 1-1e-7), 1e-7)
  u2 <- pmax(pmin(pfun(y2, mu=mu2, sigma=s2, nu=nu2), 1-1e-7), 1e-7)
  log(VineCopula::BiCopPDF(u1, u2, fam_num, theta_val))
}

true_d2_i1 <- (logc_s(s1_0+h_s, s2_0) - 2*logc_s(s1_0, s2_0) + logc_s(s1_0-h_s, s2_0)) / h_s^2
true_d2_i2 <- (logc_s(s1_0, s2_0+h_s) - 2*logc_s(s1_0, s2_0) + logc_s(s1_0, s2_0-h_s)) / h_s^2
true_cross  <- (logc_s(s1_0+h_s, s2_0+h_s) - logc_s(s1_0+h_s, s2_0-h_s) -
                logc_s(s1_0-h_s, s2_0+h_s) + logc_s(s1_0-h_s, s2_0-h_s)) / (4*h_s^2)

cat("\nFD d2 log c / d s_i1^2:", true_d2_i1, "\n")
cat("FD d2 log c / d s_i2^2:", true_d2_i2, "\n")
cat("FD d2 log c / (d s_i1 d s_i2):", true_cross, "\n")

cat("\nFormula cop_d2l_margin sigma[i1]:", cop2$cop_d2l_margin$sigma$sigma[i1], "\n")
cat("Formula cross_pair[k]:", cop2$cross_pair_contribs$sigma$sigma[k], "\n")

# TRUE theta: eta_inv$theta is already linkinv-transformed
true_theta <- as.numeric(eta_out2$eta_inv$theta[pair_cache2$theta_index_map[i1]])
cat("True natural theta:", true_theta, "\n")

u1_v <- calc_lik2$margin_p[i1]; u2_v <- calc_lik2$margin_p[i2]
cat("u1:", u1_v, "u2:", u2_v, "\n")
s1_0 <- eta_inv2$sigma[i1]; s2_0 <- eta_inv2$sigma[i2]

h_s <- 1e-4
pfun <- eval(parse(text = paste0("p", fit$margin_dist$family[1])))
y1 <- fit$response[i1]; y2 <- fit$response[i2]
mu1 <- eta_inv2$mu[i1]; mu2 <- eta_inv2$mu[i2]
nu1 <- eta_inv2$nu[i1]; nu2 <- eta_inv2$nu[i2]

# FD using TRUE theta
logc_s2 <- function(s1, s2) {
  uu1 <- pmax(pmin(pfun(y1, mu=mu1, sigma=s1, nu=nu1), 1-1e-7), 1e-7)
  uu2 <- pmax(pmin(pfun(y2, mu=mu2, sigma=s2, nu=nu2), 1-1e-7), 1e-7)
  log(VineCopula::BiCopPDF(uu1, uu2, fam_num, true_theta))
}
true_d2_i1_corrected  <- (logc_s2(s1_0+h_s, s2_0) - 2*logc_s2(s1_0, s2_0) + logc_s2(s1_0-h_s, s2_0)) / h_s^2
true_d2_i2_corrected  <- (logc_s2(s1_0, s2_0+h_s) - 2*logc_s2(s1_0, s2_0) + logc_s2(s1_0, s2_0-h_s)) / h_s^2
true_cross_corrected  <- (logc_s2(s1_0+h_s, s2_0+h_s) - logc_s2(s1_0+h_s, s2_0-h_s) -
                          logc_s2(s1_0-h_s, s2_0+h_s) + logc_s2(s1_0-h_s, s2_0-h_s)) / (4*h_s^2)

cat("\nFD d2 log c / d s_i1^2 (correct theta):", true_d2_i1_corrected, "\n")
cat("FD d2 log c / d s_i2^2 (correct theta):", true_d2_i2_corrected, "\n")
cat("FD d2 log c / (d s_i1 d s_i2) (correct theta):", true_cross_corrected, "\n")
cat("\ncop_d2l_margin sigma[i1]:", cop2$cop_d2l_margin$sigma$sigma[i1], "\n")
cat("cross_pair_contribs sigma[k]:", cop2$cross_pair_contribs$sigma$sigma[k], "\n")
cat("Error diag_i1:", cop2$cop_d2l_margin$sigma$sigma[i1] - true_d2_i1_corrected, "\n")
cat("Error cross:", cop2$cross_pair_contribs$sigma$sigma[k] - true_cross_corrected, "\n")
} else {
  cat("Debug sigma copula contribution skipped (set run_sigma_debug <- TRUE to run).\n")
}
#endregion

