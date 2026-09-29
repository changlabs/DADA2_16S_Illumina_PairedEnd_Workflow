# Resolve the workflow's mutable data and result roots.
#
# Normal runs keep the historical defaults (data/ and results/). A caller can
# isolate another dataset by setting DADA2_DATA_DIR and DADA2_RESULTS_DIR to
# absolute paths or paths relative to the project root. This is used by the
# bundled tutorial example so its inputs and outputs never mix with user data.

resolve_workflow_root <- function(environment_variable,
                                  default_directory,
                                  project_root = here::here()) {
  configured_directory <- trimws(Sys.getenv(environment_variable, unset = ""))

  if (!nzchar(configured_directory)) {
    configured_directory <- default_directory
  }

  configured_directory <- path.expand(configured_directory)
  is_absolute <- grepl(
    "^(?:/|[A-Za-z]:[/\\\\]|\\\\\\\\)",
    configured_directory,
    perl = TRUE
  )

  if (!is_absolute) {
    configured_directory <- file.path(project_root, configured_directory)
  }

  normalizePath(configured_directory, winslash = "/", mustWork = FALSE)
}

workflow_data_dir <- function(..., project_root = here::here()) {
  file.path(
    resolve_workflow_root("DADA2_DATA_DIR", "data", project_root),
    ...
  )
}

workflow_results_dir <- function(..., project_root = here::here()) {
  file.path(
    resolve_workflow_root("DADA2_RESULTS_DIR", "results", project_root),
    ...
  )
}

workflow_path_configuration <- function(project_root = here::here()) {
  data.frame(
    Setting = c("Data root", "Results root"),
    Path = c(
      workflow_data_dir(project_root = project_root),
      workflow_results_dir(project_root = project_root)
    ),
    stringsAsFactors = FALSE
  )
}
