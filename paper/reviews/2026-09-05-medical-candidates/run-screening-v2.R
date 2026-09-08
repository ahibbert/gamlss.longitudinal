#!/usr/bin/env Rscript

# Local exploratory comparison of the package-shipped public medical datasets.
# This deliberately remains outside the manuscript targets graph.

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) {
  normalizePath(sub("^--file=", "", script_argument[[1L]]), mustWork = FALSE)
} else {
  file.path(getwd(), "paper", "reviews", "2026-09-05-medical-candidates", "run-screening-v2.R")
}
package_root <- normalizePath(
  file.path(dirname(script_path), "..", "..", ".."),
  mustWork = FALSE
)
if (!file.exists(file.path(package_root, "DESCRIPTION"))) {
  stop("Could not locate the gamlss.longitudinal package root from this script.")
}

if (!requireNamespace("devtools", quietly = TRUE)) {
  stop("Install devtools with install.packages('devtools').")
}

setwd(package_root)
devtools::load_all(quiet = TRUE)

output_dir <- file.path("paper", "reviews", "2026-09-05-medical-candidates", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_package_data <- function(name) {
  environment <- new.env(parent = emptyenv())
  utils::data(list = name, package = "gamlss.longitudinal", envir = environment)
  get(name, envir = environment, inherits = FALSE)
}

read_tlc <- function() {
  data <- load_package_data("tlc_lead")
  data$week_factor <- factor(data$week)
  data$week_placebo <- data$week * as.integer(data$treatment == "placebo")
  data
}

read_actg <- function() {
  data <- load_package_data("actg_cd4")
  data$week_factor <- factor(data$week)
  data$baseline_log_z <- as.numeric(scale(log1p(data$baseline_cd4)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

read_skin <- function() {
  data <- load_package_data("skin_cancer")
  data$year_factor <- factor(data$year)
  data$exposure_z <- as.numeric(scale(log1p(data$exposure_count)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

read_epilepsy <- function() {
  data <- load_package_data("epilepsy_counts")
  data$period_factor <- factor(data$period)
  data$baseline_log_z <- as.numeric(scale(log1p(data$baseline_seizures)))
  data$age_z <- as.numeric(scale(data$age))
  data
}

tlc_data = read_tlc()
actg_data = read_actg()
skin_data = read_skin()
epi_data = read_epilepsy()

head(tlc_data)
head(actg_data)
head(skin_data)
head(epi_data)

sample_patients <- function(data, id_var, n = 10, seed = 1) {
  set.seed(seed)
  ids <- sample(unique(data[[id_var]]), size = n)
  data[data[[id_var]] %in% ids, ]
}

actg_patient_sample <- sample_patients(actg_data, id_var = "subject", n = 300)
print(actg_patient_sample)


#tlc_mar   = select_margin(response_var="lead"   ,data=tlc_data  ,time_intercepts=TRUE ,time_var="week", type="realAll")
actg_mar  = select_margin(response_var="cd4"    ,data=actg_patient_sample ,time_intercepts=TRUE ,time_var="week_factor", type="counts")
#skin_mar  = select_margin(response_var="cancers",data=skin_data ,time_intercepts=TRUE ,time_var="year_factor", type="counts")
#epi_mar   = select_margin(response_var="seizures",data=epi_data,time_intercepts=TRUE,time_var="period_factor", type="counts")

plot_margin_fit(response_var="cd4"    ,data=actg_patient_sample ,time_intercepts=TRUE ,time_var="week_factor",
margin_dist=NBI())
plot_margin_fit(response_var="cd4"    ,data=actg_patient_sample ,time_intercepts=TRUE ,time_var="week_factor",
margin_dist=DEL(),by_time=TRUE)

#plot_dist(response_var="lead"   ,data=tlc_data  ,time_var="week", subject_var="subject")
plot_dist(response_var="cd4"    ,data=actg_patient_sample ,time_var="week_factor",subject_var="subject")
#plot_dist(response_var="cancers",data=skin_data ,time_var="year_factor")
#plot_dist(response_var="seizures",data=epi_data,time_var="period_factor",subject_var="subject")

select_copula(response_var="cd4", data=actg_patient_sample, time_intercepts=TRUE, time_var="week_factor",
margin_dist=DEL(), subject_var="subject")

fit_actg=gamlss_longitudinal(data  = actg_patient_sample
                    , time_var  = "week_factor"
                    , subject_var = "subject"
                    , margin_dist = DEL()
                    , copula_dist= "N"
                    , mu.formula = cd4 ~ treatment + s(age) + sex + baseline_cd4 + week_factor + treatment * week_factor
                    , missingness = "segment"
)
