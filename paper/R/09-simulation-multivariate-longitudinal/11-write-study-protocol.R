source(file.path("paper", "R", "09-simulation-multivariate-longitudinal", "00-multivariate-setup.R"))

protocol <- mvt_write_study_protocol(mvt_output_root)

message("Study protocol written to: ", protocol$md)
message("Study protocol tables written to: ", file.path(mvt_output_root, "study_protocol"))
