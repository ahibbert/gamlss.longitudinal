source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

source_dirs <- mvt_env_vector("GAMLSS_LONGITUDINAL_MVT_SHARD_DIRS", character())
if (length(source_dirs) == 0L) {
  stop(
    "Set GAMLSS_LONGITUDINAL_MVT_SHARD_DIRS to a comma-separated list of shard directories.",
    call. = FALSE
  )
}

target_dir <- mvt_env("GAMLSS_LONGITUDINAL_MVT_MERGED_DIR", "")
if (!nzchar(target_dir)) {
  stop("Set GAMLSS_LONGITUDINAL_MVT_MERGED_DIR to the merged output directory.", call. = FALSE)
}

message("Merging ", length(source_dirs), " shard(s) into: ", target_dir)
mvt_merge_run_shards(source_dirs, target_dir)
message("Merged run written to: ", target_dir)
