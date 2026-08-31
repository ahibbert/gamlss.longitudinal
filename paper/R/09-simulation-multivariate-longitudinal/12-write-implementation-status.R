source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

status <- mvt_write_implementation_status(mvt_output_root)

message("Implementation status written to: ", status$md)
message("Publication review ready: ", status$ready)
