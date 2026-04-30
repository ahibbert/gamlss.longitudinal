#region Setup and Libraries
source("R/common_functions.R");source("R/link_functions.R"); library("VineCopula");library("moments"); 
library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix"); library("MASS")
#library(gamlss.longitudinal); 
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

copula_dist="N"; margin_dist=GG(); mu=1.56; sigma=-2.1; nu=.7; theta=0.8; tau=NA; zeta=NA; simOption=10

# Parameterization control for simulation starts:
# - TRUE: mu/sigma/nu/tau values above are interpreted on eta (linear predictor) scale.
# - FALSE: values are interpreted on natural parameter scale.
margin_inputs_on_eta_scale <- TRUE
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
covariates_input=list( mu.time=.1   ,sigma.time=.2   ,nu.time=0    ,tau.time=0   ,theta.time=.2  ,zeta.time=0
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
                   , verbose = 1
                   , compute_vcov = TRUE
                   , include_dlcopdpar=TRUE
)
#endregion

#region Model Diagnostics
#################### PLOT METHOD ####################
source("R/common_functions.R")
source("R/diagnostics_topmodels.R")
source("R (testing)/plot_copula_v2.R")
summary(fit)
plot.terms(fit,  data = data_in)
plot(fit)
plot(fit, time_stratified = TRUE)
plot.copula(fit, contour_bins=5, time_stratified = TRUE, plot2_cuts=10)
#endregion

