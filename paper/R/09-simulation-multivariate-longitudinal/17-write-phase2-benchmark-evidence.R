source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

run_dir <- mvt_read_run_dir()
evidence <- mvt_write_phase2_benchmark_evidence(run_dir)

message("Phase 2 benchmark evidence written to: ", run_dir)
message(
  "Audit checks: ",
  sum(evidence$audit$status == "pass"), " pass, ",
  sum(evidence$audit$status == "fail"), " fail."
)
