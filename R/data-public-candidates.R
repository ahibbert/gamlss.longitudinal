#' US patent and R&D longitudinal panel
#'
#' An analysis-ready, balanced panel of annual patent counts and R&D measures
#' for 346 US firms observed from 1970 through 1979. Original firm identifiers
#' have been replaced by sequential integers, and the outcome-derived `sumpat`
#' field has been omitted to prevent predictor leakage.
#'
#' @format A data frame with 3,460 rows and 7 variables:
#' \describe{
#'   \item{firm}{Sequential firm identifier.}
#'   \item{year}{Calendar year, 1970--1979.}
#'   \item{patents}{Patent applications in that year that were eventually granted.}
#'   \item{rd}{Annual R&D spending in 1972 dollars.}
#'   \item{scientific}{Whether the firm belongs to the scientific sector.}
#'   \item{capital_1972}{Book value of capital in 1972.}
#'   \item{industry_code}{Two-digit applied-R&D industry classification; missing for 10 firms.}
#' }
#' @source Derived from `pglm::PatentsRDUS`, pglm 0.2-4, distributed under
#'   GPL (>= 2). Hall, Griliches, and Hausman (1986), "Patents and R&D: Is
#'   There a Lag?" See \url{https://stat.ethz.ch/CRAN/web/packages/pglm/refman/pglm.html}.
#' @seealso [gamlss_longitudinal()]
"patents_panel"

#' Mayo Clinic PBC sequential prothrombin measurements
#'
#' Public longitudinal prothrombin measurements for the 312 randomized
#' participants in the Mayo Clinic primary biliary cirrhosis trial. Original
#' patient identifiers have been normalized to sequential integers. Baseline
#' covariates are repeated to support distributional and dependence modelling.
#'
#' @format A data frame with 1,945 rows and 13 variables:
#' \describe{
#'   \item{subject}{Sequential participant identifier.}
#'   \item{visit}{Order of the recorded visit within participant.}
#'   \item{day}{Days since enrollment.}
#'   \item{years}{Years since enrollment, computed as `day / 365.25`.}
#'   \item{prothrombin}{Standardized blood clotting time.}
#'   \item{stage}{Histological disease stage at the current record.}
#'   \item{baseline_stage}{Histological stage at enrollment.}
#'   \item{baseline_bilirubin}{Serum bilirubin at enrollment.}
#'   \item{baseline_age}{Age at enrollment in years.}
#'   \item{sex}{Recorded sex.}
#'   \item{treatment}{Randomized treatment arm.}
#'   \item{followup_days}{Follow-up time to endpoint or censoring.}
#'   \item{endpoint}{Censoring, transplant, or death.}
#' }
#' @source Derived from `survival::pbcseq`, survival 3.8-11, distributed under
#'   LGPL (>= 2), copyright Mayo Foundation for Medical Education and Research.
#'   See \url{https://stat.ethz.ch/R-manual/R-devel/library/survival/html/pbcseq.html}.
#' @seealso [gamlss_longitudinal()]
"pbc_prothrombin"

#' Daily step counts in Vietnamese adolescents
#'
#' Seven daily pedometer measurements for 475 sixth-grade students in Ho Chi
#' Minh City, with participant-level BMI and Physical Activity Questionnaire
#' for Older Children (PAQ-C) score. Original study identifiers and raw
#' anthropometric/questionnaire items are omitted. Five invalid or visibly
#' fractional step-count cells are represented as missing and identified by
#' `step_status`.
#'
#' @format A data frame with 3,325 rows and 9 variables:
#' \describe{
#'   \item{subject}{Sequential participant identifier.}
#'   \item{day}{Day index, Monday = 1 through Sunday = 7.}
#'   \item{day_name}{Ordered day-of-week factor.}
#'   \item{steps}{Daily integer step count; five invalid/fractional cells are `NA`.}
#'   \item{step_status}{Whether the source cell was observed, out of range, or fractional.}
#'   \item{age}{Age in years.}
#'   \item{sex}{Boy or girl, inferred from the coding and published group totals.}
#'   \item{bmi}{BMI from the means of duplicate weight and height measurements.}
#'   \item{paqc}{PAQ-C score reconstructed according to the article's scoring rule.}
#' }
#' @source Truong, Huynh, and To (2025), "Comparison between self-reported and
#'   pedometer-measured physical activity in Vietnamese adolescents."
#'   \doi{10.1371/journal.pgph.0004725}. Dataset
#'   \doi{10.1371/journal.pgph.0004725.s002}, CC BY 4.0. The data were modified
#'   as documented in `data-raw/build-public-candidate-data.R`.
#' @seealso [gamlss_longitudinal()], [check_missingness()]
"vietnam_steps"
