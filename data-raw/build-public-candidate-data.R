# Build the public datasets used by the paper analyses.
#
# This script is intentionally not run during package installation. It records
# the exact upstream versions, checksums, and transformations used to create the
# compressed objects in data/. See inst/DATA-LICENSES.md before updating an
# upstream source.

dir.create("data", showWarnings = FALSE)

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("The digest package is required to verify source archives.", call. = FALSE)
}

assert_true <- function(x, message) {
  if (!isTRUE(x)) stop(message, call. = FALSE)
}

download_first <- function(urls, destination) {
  errors <- character()
  for (url in urls) {
    ok <- tryCatch({
      utils::download.file(url, destination, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) {
      errors <<- c(errors, paste(url, conditionMessage(e), sep = ": "))
      FALSE
    })
    if (ok) return(invisible(destination))
  }
  stop("Unable to download upstream source:\n", paste(errors, collapse = "\n"), call. = FALSE)
}

temporary_root <- tempfile("gamlss-longitudinal-public-data-")
dir.create(temporary_root)
on.exit(unlink(temporary_root, recursive = TRUE, force = TRUE), add = TRUE)

# -------------------------------------------------------------------------
# Treatment of Lead-Exposed Children trial, ALA 1.0, GPL-3
# -------------------------------------------------------------------------

ala_archive <- file.path(temporary_root, "ALA_1.0.tar.gz")
download_first(
  "https://r-forge.r-project.org/src/contrib/ALA_1.0.tar.gz",
  ala_archive
)
ala_sha256 <- digest::digest(ala_archive, algo = "sha256", file = TRUE)
assert_true(
  identical(ala_sha256, "ac1b0afe0beda5eefbbf286ec4ddd87b854a1bf40bd813bf4c5c20b5bb1ab058"),
  paste("Unexpected ALA source checksum:", ala_sha256)
)

ala_extract <- file.path(temporary_root, "ALA-source")
dir.create(ala_extract)
utils::untar(ala_archive, exdir = ala_extract)
ala_environment <- new.env(parent = emptyenv())
load(file.path(ala_extract, "ALA", "data", "tlc.rda"), envir = ala_environment)
tlc_source <- ala_environment$tlc

assert_true(nrow(tlc_source) == 400L, "ALA::tlc row count changed.")
assert_true(length(unique(tlc_source$id)) == 100L, "ALA::tlc subject count changed.")
assert_true(
  identical(sort(unique(as.numeric(tlc_source$week))), c(0, 1, 4, 6)),
  "ALA::tlc visit schedule changed."
)

tlc_lead <- data.frame(
  subject = as.integer(as.character(tlc_source$id)),
  treatment = factor(
    as.character(tlc_source$treatment),
    levels = c("succimer", "placebo")
  ),
  week = as.integer(tlc_source$week),
  lead = as.numeric(tlc_source$lead)
)
tlc_lead <- tlc_lead[order(tlc_lead$subject, tlc_lead$week), ]
rownames(tlc_lead) <- NULL

assert_true(all(table(tlc_lead$subject) == 4L), "TLC panel is no longer balanced.")
assert_true(!anyNA(tlc_lead), "Unexpected missing TLC values.")
assert_true(all(tlc_lead$lead > 0), "TLC blood lead measurements must be positive.")

# -------------------------------------------------------------------------
# Skin Cancer Prevention Study, ALA 1.0, GPL-3
# -------------------------------------------------------------------------

load(file.path(ala_extract, "ALA", "data", "skin.rda"), envir = ala_environment)
skin_source <- ala_environment$skin

assert_true(nrow(skin_source) == 7081L, "ALA::skin row count changed.")
assert_true(length(unique(skin_source$id)) == 1683L, "ALA::skin subject count changed.")

skin_cancer <- data.frame(
  subject = match(as.character(skin_source$id), unique(as.character(skin_source$id))),
  center = factor(as.character(skin_source$center)),
  age = as.integer(skin_source$age),
  skin_type = factor(
    as.character(skin_source$skin),
    levels = c("otherwise", "burns")
  ),
  sex = factor(
    as.character(skin_source$gender),
    levels = c("F", "M"),
    labels = c("female", "male")
  ),
  exposure_count = as.integer(skin_source$nExposure),
  cancers = as.integer(skin_source$nCancers),
  treatment = factor(
    as.character(skin_source$treatment),
    levels = c("placebo", "B-carotene"),
    labels = c("placebo", "beta_carotene")
  ),
  year = as.integer(skin_source$year)
)
skin_cancer <- skin_cancer[order(skin_cancer$subject, skin_cancer$year), ]
rownames(skin_cancer) <- NULL

assert_true(!anyNA(skin_cancer), "Unexpected missing ALA::skin values.")
assert_true(all(skin_cancer$cancers >= 0L), "Skin cancer counts must be non-negative.")
assert_true(
  all(skin_cancer$cancers == floor(skin_cancer$cancers)),
  "Skin cancer outcomes must be integers."
)

# -------------------------------------------------------------------------
# ACTG 193A CD4 trial, bcmixed 0.1.6, GPL >= 2
# -------------------------------------------------------------------------

bcmixed_archive <- file.path(temporary_root, "bcmixed_0.1.6.tar.gz")
download_first(
  c(
    "https://cran.r-project.org/src/contrib/bcmixed_0.1.6.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/bcmixed/bcmixed_0.1.6.tar.gz"
  ),
  bcmixed_archive
)
bcmixed_sha256 <- digest::digest(bcmixed_archive, algo = "sha256", file = TRUE)
assert_true(
  identical(bcmixed_sha256, "876d85c20e73eb6dacee32d7b165a31f386f8fec103371d1a1b4c60150c4e7ec"),
  paste("Unexpected bcmixed source checksum:", bcmixed_sha256)
)

bcmixed_extract <- file.path(temporary_root, "bcmixed-source")
dir.create(bcmixed_extract)
utils::untar(bcmixed_archive, exdir = bcmixed_extract)
bcmixed_environment <- new.env(parent = emptyenv())
load(file.path(bcmixed_extract, "bcmixed", "data", "aidscd4.rda"), envir = bcmixed_environment)
actg_source <- bcmixed_environment$aidscd4

assert_true(nrow(actg_source) == 4708L, "bcmixed::aidscd4 row count changed.")
assert_true(length(unique(actg_source$id)) == 1177L, "bcmixed::aidscd4 subject count changed.")

actg_cd4 <- data.frame(
  subject = as.integer(actg_source$id),
  week = as.integer(actg_source$weekc),
  treatment = factor(
    actg_source$treatment,
    levels = 1:4,
    labels = c(
      "alternating_zidovudine_didanosine",
      "zidovudine_zalcitabine",
      "zidovudine_didanosine",
      "zidovudine_didanosine_nevirapine"
    )
  ),
  age = as.numeric(actg_source$age),
  sex = factor(actg_source$sex, levels = c(0, 1), labels = c("female", "male")),
  baseline_cd4 = as.numeric(actg_source$cd4.bl) - 1,
  cd4 = as.integer(actg_source$cd4) - 1L
)
actg_cd4 <- actg_cd4[order(actg_cd4$subject, actg_cd4$week), ]
rownames(actg_cd4) <- NULL

assert_true(all(table(actg_cd4$subject) == 4L), "ACTG CD4 panel is no longer balanced.")
assert_true(sum(is.na(actg_cd4$cd4)) == 1356L, "ACTG CD4 missing-response pattern changed.")
assert_true(
  !anyNA(actg_cd4[setdiff(names(actg_cd4), "cd4")]),
  "Unexpected missing ACTG covariate values."
)
assert_true(all(actg_cd4$cd4[!is.na(actg_cd4$cd4)] >= 0L), "ACTG CD4 counts must be non-negative.")

# -------------------------------------------------------------------------
# Epilepsy seizure trial, MASS 7.3-66, GPL-2 | GPL-3
# -------------------------------------------------------------------------

mass_archive <- file.path(temporary_root, "MASS_7.3-66.tar.gz")
download_first(
  c(
    "https://cran.r-project.org/src/contrib/MASS_7.3-66.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/MASS/MASS_7.3-66.tar.gz"
  ),
  mass_archive
)
mass_sha256 <- digest::digest(mass_archive, algo = "sha256", file = TRUE)
assert_true(
  identical(mass_sha256, "59c619c70bfca653f1c74aef880e69405daefd0134e716a1862ecc5334b4e155"),
  paste("Unexpected MASS source checksum:", mass_sha256)
)

mass_extract <- file.path(temporary_root, "MASS-source")
dir.create(mass_extract)
utils::untar(mass_archive, exdir = mass_extract)
mass_environment <- new.env(parent = emptyenv())
load(file.path(mass_extract, "MASS", "data", "epil.rda"), envir = mass_environment)
epilepsy_source <- mass_environment$epil

assert_true(nrow(epilepsy_source) == 236L, "MASS::epil row count changed.")
assert_true(length(unique(epilepsy_source$subject)) == 59L, "MASS::epil subject count changed.")

epilepsy_counts <- data.frame(
  subject = as.integer(epilepsy_source$subject),
  period = as.integer(epilepsy_source$period),
  seizures = as.integer(epilepsy_source$y),
  treatment = factor(epilepsy_source$trt, levels = c("placebo", "progabide")),
  baseline_seizures = as.integer(epilepsy_source$base),
  age = as.integer(epilepsy_source$age),
  fourth_visit = as.logical(epilepsy_source$V4)
)
epilepsy_counts <- epilepsy_counts[order(epilepsy_counts$subject, epilepsy_counts$period), ]
rownames(epilepsy_counts) <- NULL

assert_true(!anyNA(epilepsy_counts), "Unexpected missing epilepsy values.")
assert_true(all(table(epilepsy_counts$subject) == 4L), "Epilepsy panel is no longer balanced.")
assert_true(all(epilepsy_counts$seizures >= 0L), "Seizure counts must be non-negative.")

# -------------------------------------------------------------------------
# PatentsRDUS, pglm 0.2-4, GPL >= 2
# -------------------------------------------------------------------------

pglm_archive <- file.path(temporary_root, "pglm_0.2-4.tar.gz")
download_first(
  c(
    "https://cran.r-project.org/src/contrib/pglm_0.2-4.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/pglm/pglm_0.2-4.tar.gz"
  ),
  pglm_archive
)
pglm_sha256 <- digest::digest(pglm_archive, algo = "sha256", file = TRUE)
assert_true(
  identical(pglm_sha256, "caf085e7f5d693efdb06a96057c6e1d9e132a9525b0bd6a595804e2b853659e3"),
  paste("Unexpected pglm source checksum:", pglm_sha256)
)

pglm_extract <- file.path(temporary_root, "pglm-source")
dir.create(pglm_extract)
utils::untar(pglm_archive, exdir = pglm_extract)
pglm_environment <- new.env(parent = emptyenv())
load(file.path(pglm_extract, "pglm", "data", "PatentsRDUS.rda"), envir = pglm_environment)
patents_source <- pglm_environment$PatentsRDUS

assert_true(nrow(patents_source) == 3460L, "PatentsRDUS row count changed.")
assert_true(length(unique(patents_source$cusip)) == 346L, "PatentsRDUS firm count changed.")

firm_levels <- sort(unique(patents_source$cusip))
patents_panel <- data.frame(
  firm = match(patents_source$cusip, firm_levels),
  year = as.integer(as.character(patents_source$year)),
  patents = as.integer(patents_source$patents),
  rd = as.numeric(patents_source$rd),
  scientific = factor(patents_source$scisect, levels = c("no", "yes")),
  capital_1972 = as.numeric(patents_source$capital72),
  industry_code = as.integer(as.character(patents_source$ardssic))
)
patents_panel <- patents_panel[order(patents_panel$firm, patents_panel$year), ]
rownames(patents_panel) <- NULL

assert_true(all(table(patents_panel$firm) == 10L), "Patents panel is no longer balanced.")
assert_true(
  all(patents_panel$patents >= 0L & patents_panel$patents == floor(patents_panel$patents)),
  "Patent outcomes must be non-negative integers."
)

# -------------------------------------------------------------------------
# Mayo PBC sequential data, survival 3.8-11, LGPL >= 2
# -------------------------------------------------------------------------

survival_archive <- file.path(temporary_root, "survival_3.8-11.tar.gz")
download_first(
  c(
    "https://cran.r-project.org/src/contrib/survival_3.8-11.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/survival/survival_3.8-11.tar.gz"
  ),
  survival_archive
)
survival_sha256 <- digest::digest(survival_archive, algo = "sha256", file = TRUE)
assert_true(
  identical(survival_sha256, "4a87aea323d477e142c36601509f3771f150b441df6db75537ea1777e6546888"),
  paste("Unexpected survival source checksum:", survival_sha256)
)

survival_extract <- file.path(temporary_root, "survival-source")
dir.create(survival_extract)
utils::untar(survival_archive, exdir = survival_extract)
pbc_environment <- new.env(parent = emptyenv())
load(file.path(survival_extract, "survival", "data", "pbc.rda"), envir = pbc_environment)
pbc_source <- pbc_environment$pbcseq

assert_true(nrow(pbc_source) == 1945L, "pbcseq row count changed.")
assert_true(length(unique(pbc_source$id)) == 312L, "pbcseq subject count changed.")

pbc_source <- pbc_source[order(pbc_source$id, pbc_source$day), ]
pbc_source$visit <- ave(pbc_source$day, pbc_source$id, FUN = seq_along)
first_record <- !duplicated(pbc_source$id)
baseline <- pbc_source[first_record, c("id", "age", "sex", "trt", "stage", "bili")]
names(baseline) <- c(
  "id", "baseline_age", "baseline_sex", "baseline_treatment",
  "baseline_stage", "baseline_bilirubin"
)
pbc_source <- merge(pbc_source, baseline, by = "id", all.x = TRUE, sort = FALSE)
pbc_source <- pbc_source[order(pbc_source$id, pbc_source$day), ]
patient_levels <- sort(unique(pbc_source$id))

pbc_prothrombin <- data.frame(
  subject = match(pbc_source$id, patient_levels),
  visit = as.integer(pbc_source$visit),
  day = as.integer(pbc_source$day),
  years = as.numeric(pbc_source$day / 365.25),
  prothrombin = as.numeric(pbc_source$protime),
  stage = ordered(pbc_source$stage, levels = 1:4, labels = paste("stage", 1:4)),
  baseline_stage = ordered(
    pbc_source$baseline_stage, levels = 1:4, labels = paste("stage", 1:4)
  ),
  baseline_bilirubin = as.numeric(pbc_source$baseline_bilirubin),
  baseline_age = as.numeric(pbc_source$baseline_age),
  sex = factor(pbc_source$baseline_sex, levels = c("m", "f"), labels = c("male", "female")),
  treatment = factor(
    pbc_source$baseline_treatment,
    # survival 3.8-11 stores pbcseq treatment as 1/0. This mapping was
    # cross-checked against the 1/2 coding in the baseline pbc object.
    levels = c(1, 0),
    labels = c("D-penicillamine", "placebo")
  ),
  followup_days = as.integer(pbc_source$futime),
  endpoint = factor(
    pbc_source$status,
    levels = c(0, 1, 2),
    labels = c("censored", "transplant", "death")
  )
)
rownames(pbc_prothrombin) <- NULL

assert_true(!anyNA(pbc_prothrombin$prothrombin), "Unexpected missing prothrombin values.")
assert_true(all(pbc_prothrombin$prothrombin > 0), "Prothrombin values must be positive.")

# -------------------------------------------------------------------------
# Vietnamese adolescent pedometer data, CC BY 4.0
# -------------------------------------------------------------------------

assert_true(requireNamespace("readxl", quietly = TRUE), "The readxl package is required.")
steps_workbook <- file.path(temporary_root, "pgph.0004725.s002.xlsx")
download_first("https://ndownloader.figshare.com/files/55180168", steps_workbook)
steps_md5 <- unname(tools::md5sum(steps_workbook))
assert_true(
  identical(steps_md5, "722cb5aa627ff9273712c485e130800b"),
  paste("Unexpected PLOS workbook checksum:", steps_md5)
)

steps_source <- as.data.frame(readxl::read_excel(
  steps_workbook,
  sheet = "Second sample (N=475)",
  skip = 2,
  col_names = FALSE,
  .name_repair = "minimal"
))
for (column in seq_len(ncol(steps_source))) {
  steps_source[[column]] <- as.numeric(steps_source[[column]])
}

names(steps_source)[1:14] <- c(
  "source_id", "age", "sex", "weight_1", "weight_2", "height_1", "height_2",
  "monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday"
)
day_names <- c("monday", "tuesday", "wednesday", "thursday", "friday", "saturday", "sunday")
steps_source$bmi <- rowMeans(steps_source[c("weight_1", "weight_2")]) /
  (rowMeans(steps_source[c("height_1", "height_2")]) / 100)^2
steps_source$paqc <- rowMeans(cbind(
  q1 = rowMeans(steps_source[15:26]),
  steps_source[27:33],
  q9 = rowMeans(steps_source[34:40])
))

assert_true(nrow(steps_source) == 475L, "PLOS participant count changed.")
assert_true(
  identical(as.integer(table(steps_source$sex)), c(241L, 234L)),
  "PLOS sex coding changed."
)

vietnam_steps <- do.call(rbind, lapply(seq_along(day_names), function(day_index) {
  original_steps <- steps_source[[day_names[[day_index]]]]
  out_of_range <- original_steps < 1000 | original_steps > 30000
  fractional <- abs(original_steps - round(original_steps)) > sqrt(.Machine$double.eps)
  status <- ifelse(out_of_range, "out_of_range", ifelse(fractional, "fractional", "observed"))
  cleaned_steps <- original_steps
  cleaned_steps[out_of_range | fractional] <- NA_real_

  data.frame(
    subject = seq_len(nrow(steps_source)),
    day = day_index,
    day_name = ordered(day_names[[day_index]], levels = day_names),
    steps = cleaned_steps,
    step_status = factor(status, levels = c("observed", "out_of_range", "fractional")),
    age = steps_source$age,
    sex = factor(steps_source$sex, levels = c(1, 2), labels = c("boy", "girl")),
    bmi = steps_source$bmi,
    paqc = steps_source$paqc
  )
}))
vietnam_steps <- vietnam_steps[order(vietnam_steps$subject, vietnam_steps$day), ]
rownames(vietnam_steps) <- NULL

assert_true(nrow(vietnam_steps) == 3325L, "Vietnam step panel size changed.")
assert_true(sum(is.na(vietnam_steps$steps)) == 5L, "Expected five excluded step cells.")
assert_true(
  all(vietnam_steps$steps[!is.na(vietnam_steps$steps)] ==
        floor(vietnam_steps$steps[!is.na(vietnam_steps$steps)])),
  "Cleaned step outcomes must be integers."
)

save(tlc_lead, file = file.path("data", "tlc_lead.rda"), compress = "xz", version = 3)
save(actg_cd4, file = file.path("data", "actg_cd4.rda"), compress = "xz", version = 3)
save(skin_cancer, file = file.path("data", "skin_cancer.rda"), compress = "xz", version = 3)
save(epilepsy_counts, file = file.path("data", "epilepsy_counts.rda"), compress = "xz", version = 3)
save(patents_panel, file = file.path("data", "patents_panel.rda"), compress = "xz", version = 3)
save(pbc_prothrombin, file = file.path("data", "pbc_prothrombin.rda"), compress = "xz", version = 3)
save(vietnam_steps, file = file.path("data", "vietnam_steps.rda"), compress = "xz", version = 3)

message("Wrote package data assets to data/.")
