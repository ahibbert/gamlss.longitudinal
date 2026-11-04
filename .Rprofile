if (Sys.getenv("RSTUDIO_PANDOC") == "") {
	quarto_exe <- Sys.which("quarto")
	candidate_dirs <- character()

	if (nzchar(quarto_exe)) {
		candidate_dirs <- c(candidate_dirs, file.path(dirname(quarto_exe), "tools"))
	}

	candidate_dirs <- c(
		candidate_dirs,
		"C:/Program Files/Quarto/bin/tools",
		file.path(Sys.getenv("LOCALAPPDATA"), "Programs", "Quarto", "bin", "tools")
	)

	pandoc_dir <- candidate_dirs[vapply(candidate_dirs, dir.exists, logical(1))][1]

	if (!is.na(pandoc_dir) && nzchar(pandoc_dir)) {
		Sys.setenv(RSTUDIO_PANDOC = pandoc_dir)
		options(rmarkdown.pandoc.dir = pandoc_dir)
	}
}

