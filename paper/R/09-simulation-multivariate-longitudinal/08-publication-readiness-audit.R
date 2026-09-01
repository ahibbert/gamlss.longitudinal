source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

mvt_load_package()

readiness <- mvt_publication_readiness_audit(mvt_output_root)

message("Publication readiness audit written to: ", readiness$csv)
message("Publication readiness report written to: ", readiness$md)
message("Required evidence ready: ", readiness$ready)
