source("R/common_functions.R");source("R/link_functions.R"); library("VineCopula");library("moments"); 
library("ggplot2"); library("latex2exp"); library("ggpubr"); library("Matrix"); library("MASS")
#library(gamlss.longitudinal); 
library(gamlss2); library(gamlss)
set.seed(100)
#########DATASET
n=500; d=4

# Missingness configuration:
# - "increasing_time": p_miss(i) = i/(T+1) by ordered time index i.
# - "mar": missing completely at random across all rows at mar_missing_rate.
# mar_missing_rate controls the overall target missingness level for both modes.
missingness_mode = "mar" # or "mar"
mar_missing_rate = 0.1

copula_dist="N"; margin_dist=BCPEo(); mu=1; sigma=0.5;nu=-1; tau=1; theta=-0.5; zeta=NA; simOption=10;
copula_dist="C"; margin_dist=BCPEo(); mu=1; sigma=0.5;nu=-1; tau=1; theta=-2; zeta=NA; simOption=10;


# USE THIS WITH SIMOPTION 10
covariates_input=list( mu.time=0.1   ,sigma.time=0.1   ,nu.time=1    ,tau.time=0.1   ,theta.time=1  ,zeta.time=0
                        ,mu.age=1    ,sigma.age=0.5     ,nu.age=0     ,tau.age=0    ,theta.age=0    ,zeta.age=0
                        ,mu.gender=0.1 ,sigma.gender=0.1  ,nu.gender=0  ,tau.gender=0 ,theta.gender=0 ,zeta.gender=0)

#########Generate dataset

rm(dataset)
dataset=loadDataset(simOption=simOption, n=n,d=d, copula_dist=copula_dist, margin_dist=margin_dist
                    , par.margin=c(mu,sigma,nu,tau), par.copula=c(theta=theta),covariates_input=covariates_input)

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

plotDist(dataset, margin_dist, offdiag_scale = "pseudo")

#########PLOTTING#############
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

########## FIT ###########
source("R/common_functions.R")
mu_formula="random_name ~ time_of_observation_random_name + s(age_new_name,bs='ps') + gender"
sigma_formula="~ time_of_observation_random_name + gender + s(age_new_name,bs='ps')"
nu_formula="~ time_of_observation_random_name"
tau_formula="~ time_of_observation_random_name"
theta_formula="~ time_of_observation_random_name"
zeta_formula="~ time_of_observation_random_name"

### FOR TESTING DATASETS WITH NON STANDARD NAMING
if("age_group" %in% names(dataset)) dataset$age_group = NULL
colnames(dataset)=c("person","time_of_observation_random_name","random_name","age_new_name","year","gender")
data_in=dataset; data_in$gender=as.factor(data_in$gender); data_in$time_of_observation_random_name=as.factor(data_in$time_of_observation_random_name)
rm(dataset)

source("R/common_functions.R")
fit=gamlss.longitudinal(dataset=data_in
                   , margin_dist=margin_dist
                   , copula_dist=copula_dist
                   , time_var="time_of_observation_random_name"
                   , subject_var="person"
                   , mu.formula = mu_formula
                   , sigma.formula = sigma_formula
                   , nu.formula = nu_formula
                   , tau.formula = tau_formula
                   , theta.formula=theta_formula
                   , zeta.formula=zeta_formula
                   , include_dlcopdpar=FALSE
                   #, verbose=1, plot_results=FALSE,  true_val=par_to_eta(input_par,copula_dist,margin_dist)
                   #, use_Rcpp=FALSE, start_step_size=0.5, step_adjustment = 0.5, inner_stop_crit=.1, outer_stop_crit=.1
                   #, lambda_start = 5
                   #, lambda_penalty_K = 2 #Optimising for AIC
)

#vcov_fit=vcov.gamlss.longitudinal(fit, numderiv=TRUE)

#################### PLOT METHOD ####################
source("R/common_functions.R")
source("R/diagnostics_topmodels.R")
source("R (testing)/plot_copula_v2.R")
summary(fit)
plot(fit)
plot(fit, time_stratified = TRUE)
plot.copula(fit, contour_bins=5, time_stratified = TRUE, plot2_cuts=10)
plot.copula_contour_compare(fit, time_stratified = TRUE, transform="normal", diff_scale_limit=.1)
plot.terms(fit,  data = data_in)
