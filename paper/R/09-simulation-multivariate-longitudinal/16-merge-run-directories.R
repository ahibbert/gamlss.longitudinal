source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()

source_dirs <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_MERGE_SOURCE_DIRS", character())
if (length(source_dirs) == 0L) {
  stop("Set GAMLSS_LONGITUDINAL_MVT_MERGE_SOURCE_DIRS to a comma-separated list of run directories.", call. = FALSE)
}

target_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_MERGE_TARGET_DIR", "")
if (!nzchar(target_dir)) {
  stop("Set GAMLSS_LONGITUDINAL_MVT_MERGE_TARGET_DIR to the merged output directory.", call. = FALSE)
}

if (dir.exists(target_dir) && !mvt_env_flag("GAMLSS_LONGITUDINAL_MVT_MERGE_OVERWRITE", FALSE)) {
  stop(
    "Merge target already exists: ",
    target_dir,
    ". Set GAMLSS_LONGITUDINAL_MVT_MERGE_OVERWRITE=true to replace it.",
    call. = FALSE
  )
}
if (dir.exists(target_dir)) unlink(target_dir, recursive = TRUE, force = TRUE)

mvt_merge_run_shards(source_dirs, target_dir)

metadata <- mvt_read_optional_csv(file.path(target_dir, "run_metadata.csv"))
metadata <- mvt_bind_rows_fill(
  metadata,
  data.frame(
    name = c("merged_at", "merged_source_dirs", "merged_note"),
    value = c(
      format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
      paste(normalizePath(source_dirs, winslash = "/", mustWork = FALSE), collapse = ";"),
      mvt_env("GAMLSS_LONGITUDINAL_MVT_MERGE_NOTE", "Merged run directory for review or publication evidence.")
    ),
    stringsAsFactors = FALSE
  )
)
mvt_write_csv(metadata, file.path(target_dir, "run_metadata.csv"))
mvt_summarise_results(target_dir)
mvt_write_artifact_manifest(target_dir)

message("Merged run directory written to: ", normalizePath(target_dir, winslash = "/", mustWork = FALSE))
