#' Treatment of Lead-Exposed Children trial
#'
#' Four blood lead measurements from a randomized, placebo-controlled trial of
#' succimer in 100 children with baseline blood lead concentrations of 20--44
#' micrograms per decilitre. The data are a documented subsample of the TLC
#' trial and form a complete, balanced longitudinal panel.
#'
#' @format A data frame with 400 rows and 4 variables:
#' \describe{
#'   \item{subject}{Sequential participant identifier, 1--100.}
#'   \item{treatment}{Randomized treatment arm: `succimer` or `placebo`.}
#'   \item{week}{Scheduled measurement week: 0, 1, 4, or 6.}
#'   \item{lead}{Blood lead concentration in micrograms per decilitre.}
#' }
#' @source Derived from `ALA::tlc`, ALA 1.0, distributed under GPL-3.
#'   Treatment of Lead-exposed Children Trial Group (2000), "Safety and
#'   Efficacy of Succimer in Toddlers with Blood Lead Levels of 20--44
#'   micrograms/dL", *Pediatric Research*, 48, 593--599. See
#'   \url{https://r-forge.r-project.org/projects/applong/}.
#' @seealso [gamlss_longitudinal()]
"tlc_lead"

#' ACTG 193A CD4 cell counts
#'
#' Four post-baseline CD4 measurements for 1,177 participants in AIDS Clinical
#' Trial Group Study 193A. Participants with advanced immune suppression were
#' randomized to one of four reverse-transcriptase-inhibitor regimens. The
#' source stores CD4 measurements with one added; this version restores the
#' original non-negative counts. The scheduled panel is balanced, with 1,356
#' post-baseline responses missing as in the source.
#'
#' @format A data frame with 4,708 rows and 7 variables:
#' \describe{
#'   \item{subject}{Participant identifier.}
#'   \item{week}{Nominal visit week: 8, 16, 24, or 32.}
#'   \item{treatment}{Randomized treatment regimen.}
#'   \item{age}{Age in years.}
#'   \item{sex}{Recorded sex: `female` or `male`.}
#'   \item{baseline_cd4}{Baseline CD4 cell count.}
#'   \item{cd4}{Post-baseline CD4 cell count; 1,356 values are missing.}
#' }
#' @source Derived from `bcmixed::aidscd4`, bcmixed 0.1.6, distributed under
#'   GPL (>= 2). Henry et al. (1998), "A Randomized, Controlled, Double-Blind
#'   Study Comparing the Survival Benefit of Four Different Reverse
#'   Transcriptase Inhibitor Therapies", *Journal of Acquired Immune
#'   Deficiency Syndromes and Human Retrovirology*, 19, 339--349.
#'   \doi{10.1097/00042560-199812010-00004}.
#' @seealso [gamlss_longitudinal()]
"actg_cd4"

#' Skin Cancer Prevention Study
#'
#' Annual counts of new non-melanoma skin cancers for 1,683 participants in a
#' randomized trial of beta-carotene supplementation. Participants contribute
#' between one and five annual observations.
#'
#' @format A data frame with 7,081 rows and 9 variables:
#' \describe{
#'   \item{subject}{Sequential participant identifier.}
#'   \item{center}{Study center.}
#'   \item{age}{Age in years at study entry.}
#'   \item{skin_type}{Whether the participant's skin burns rather than tans.}
#'   \item{sex}{Recorded sex: `female` or `male`.}
#'   \item{exposure_count}{Number of previous skin cancers at study entry.}
#'   \item{cancers}{Number of new skin cancers during the year.}
#'   \item{treatment}{Randomized arm: `placebo` or `beta_carotene`.}
#'   \item{year}{Follow-up year, 1--5.}
#' }
#' @source Derived from `ALA::skin`, ALA 1.0, distributed under GPL-3.
#'   Greenberg et al. (1989), "A Clinical Trial of Beta Carotene to Prevent
#'   Basal-Cell and Squamous-Cell Cancers of the Skin", *New England Journal
#'   of Medicine*, 323, 789--795. See
#'   \url{https://rdrr.io/rforge/ALA/man/skin.html}.
#' @seealso [gamlss_longitudinal()]
"skin_cancer"

#' Seizure counts for participants with epilepsy
#'
#' Four two-week seizure counts for 59 participants randomized to progabide or
#' placebo. Baseline eight-week seizure counts and age are repeated for each
#' participant.
#'
#' @format A data frame with 236 rows and 7 variables:
#' \describe{
#'   \item{subject}{Participant identifier.}
#'   \item{period}{Two-week follow-up period, 1--4.}
#'   \item{seizures}{Seizure count during the period.}
#'   \item{treatment}{Randomized arm: `placebo` or `progabide`.}
#'   \item{baseline_seizures}{Seizure count during the eight-week baseline.}
#'   \item{age}{Age in years.}
#'   \item{fourth_visit}{Whether the row is the fourth follow-up period.}
#' }
#' @source Derived from `MASS::epil`, MASS 7.3-66, distributed under GPL-2 or
#'   GPL-3. Thall and Vail (1990), "Some Covariance Models for Longitudinal
#'   Count Data with Overdispersion", *Biometrics*, 46, 657--671. See
#'   \url{https://stat.ethz.ch/R-manual/R-devel/library/MASS/html/epil.html}.
#' @seealso [gamlss_longitudinal()]
"epilepsy_counts"

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
