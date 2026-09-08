# Build the three public candidate datasets used by the paper analyses.
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

save(patents_panel, file = file.path("data", "patents_panel.rda"), compress = "xz", version = 3)
save(pbc_prothrombin, file = file.path("data", "pbc_prothrombin.rda"), compress = "xz", version = 3)
save(vietnam_steps, file = file.path("data", "vietnam_steps.rda"), compress = "xz", version = 3)

message("Wrote package data assets to data/.")
