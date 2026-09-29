#!/usr/bin/env Rscript

# Run the repository's core workflow against the bundled DADA2 tutorial data.
# Inputs, generated outputs, and rendered reports stay below example/ so the
# normal data/fastq and results/ trees are never read or modified.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Run this file with Rscript: Rscript example/run_example.R", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

data_root <- file.path(project_root, "example", "data")
run_results <- file.path(project_root, "example", "run_results")
reference_results <- file.path(project_root, "example", "reference_results")
configuration_workbook <- file.path(
  project_root, "example", "config", "dada2_filter_parameters.xlsx"
)

Sys.setenv(
  DADA2_DATA_DIR = data_root,
  DADA2_RESULTS_DIR = run_results,
  DADA2_TAXONOMY_DATABASE = "GTDB",
  DADA2_FORWARD_PRIMER = "GTGCCAGCMGCCGCGGTAA",
  DADA2_REVERSE_PRIMER = "GGACTACHVGGGTWTCTAAT",
  DADA2_TRIM_READ_THROUGH = "false",
  DADA2_CUTADAPT_TIMES = "1"
)

required_packages <- c("rmarkdown", "openxlsx", "dada2", "phyloseq")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install the workflow dependencies first; missing: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# rmarkdown needs Pandoc. RStudio normally configures it automatically; when
# running from a plain terminal, reuse Quarto's bundled Pandoc if available.
if (!rmarkdown::pandoc_available()) {
  quarto_executable <- Sys.which("quarto")
  if (nzchar(quarto_executable)) {
    quarto_executable <- normalizePath(quarto_executable, mustWork = TRUE)
    pandoc_candidates <- list.files(
      file.path(dirname(quarto_executable), "tools"),
      pattern = "^pandoc$",
      recursive = TRUE,
      full.names = TRUE
    )
    pandoc_candidates <- pandoc_candidates[file.access(pandoc_candidates, 1L) == 0L]
    if (length(pandoc_candidates)) {
      Sys.setenv(RSTUDIO_PANDOC = dirname(pandoc_candidates[[1]]))
    }
  }
}
if (!rmarkdown::pandoc_available()) {
  stop("Pandoc is required. Run from RStudio or install Quarto/Pandoc.", call. = FALSE)
}

gtdb_files <- file.path(
  project_root,
  "tools", "trainsets", "GTDB",
  c(
    "GTDB_bac120_arc53_ssu_r220_genus.fa.gz",
    "GTDB_bac120_arc53_ssu_r220_species.fa.gz"
  )
)
if (!all(file.exists(gtdb_files))) {
  stop(
    "The GTDB trainset is missing. Run setup/download_reference_databases.R ",
    "before the example (the example itself uses only GTDB).",
    call. = FALSE
  )
}

cutadapt_executable <- file.path(
  project_root, "tools", "cutadapt", "venv", "bin", "cutadapt"
)
if (file.access(cutadapt_executable, 1L) != 0L) {
  stop(
    "Cutadapt is missing. Run setup/install_required_tools.R before the example.",
    call. = FALSE
  )
}

fastq_files <- sort(list.files(
  file.path(data_root, "fastq"),
  pattern = "[.]fastq[.]gz$",
  full.names = TRUE
))
if (length(fastq_files) != 40L) {
  stop("Expected 40 bundled FASTQ files but found ", length(fastq_files), call. = FALSE)
}

parameter_directory <- file.path(run_results, "4_dada2_parameter_selection")
dir.create(parameter_directory, recursive = TRUE, showWarnings = FALSE)
if (!file.copy(
  configuration_workbook,
  file.path(parameter_directory, basename(configuration_workbook)),
  overwrite = TRUE
)) {
  stop("Could not stage the example DADA2 parameter workbook.", call. = FALSE)
}

report_directory <- file.path(run_results, "reports")
dir.create(report_directory, recursive = TRUE, showWarnings = FALSE)

render_step <- function(filename) {
  message("\nRendering ", filename, " ...")
  rmarkdown::render(
    input = file.path(project_root, "R", "notebooks", filename),
    output_format = "html_document",
    output_dir = report_directory,
    envir = new.env(parent = globalenv()),
    clean = TRUE,
    quiet = FALSE
  )
}

render_step("1_data_integrity_check.Rmd")
render_step("3_cutadapt_primer_trimming.Rmd")
render_step("5_dada2_pipeline.Rmd")
render_step("9_phyloseq_object.Rmd")

message(
  "\nExample run complete.\n",
  "Generated results: ", run_results, "\n",
  "Bundled reference results: ", reference_results, "\n",
  "Compare the two trees without touching data/fastq or results/."
)
