#!/usr/bin/env Rscript

# Run the repository's complete executable workflow against the bundled DADA2
# tutorial data. Step 4 is represented by the validated, preconfigured parameter
# workbook because its Shiny app is intentionally interactive.
# Inputs, generated outputs, and rendered reports stay below example/ so the
# normal data/fastq and results/ trees are never read or modified.

script_argument <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("Run this file with Rscript: Rscript example/run_example.R", call. = FALSE)
}

script_path <- normalizePath(sub("^--file=", "", script_argument[[1]]), mustWork = TRUE)
project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
setwd(project_root)

run_arguments <- commandArgs(trailingOnly = TRUE)
if (length(run_arguments)) {
  stop(
    "This runner does not accept arguments. Use: Rscript example/run_example.R",
    call. = FALSE
  )
}

data_root <- file.path(project_root, "example", "data")
run_results <- file.path(project_root, "example", "run_results")
reference_results <- file.path(project_root, "example", "reference_results")
configuration_workbook <- file.path(
  data_root, "dada2_filter_parameters.xlsx"
)

Sys.setenv(
  DADA2_DATA_DIR = data_root,
  DADA2_RESULTS_DIR = run_results,
  DADA2_PARAMETER_FILE = configuration_workbook,
  DADA2_TAXONOMY_DATABASE = "BOTH",
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

# Validate the prepared Step 4 handoff before running any notebooks. Step 5
# reads this example-data artifact directly through DADA2_PARAMETER_FILE; its
# generic fallback values must never be used for the tutorial dataset.
expected_step4_parameters <- c(
  truncation_length_forward = 240,
  truncation_length_reverse = 160,
  max_expected_errors_forward = 2,
  max_expected_errors_reverse = 2,
  amplicon_min_length = 250,
  amplicon_max_length = 256
)

configuration_sheets <- openxlsx::getSheetNames(configuration_workbook)
if (!all(c("Info", "Parameters") %in% configuration_sheets)) {
  stop(
    "The example Step 4 workbook must contain Info and Parameters sheets: ",
    configuration_workbook,
    call. = FALSE
  )
}

step4_parameter_table <- openxlsx::read.xlsx(
  configuration_workbook,
  sheet = "Parameters",
  check.names = FALSE
)
if (!all(c("Parameter", "Value") %in% names(step4_parameter_table))) {
  stop("The example Step 4 Parameters sheet must contain Parameter and Value columns.", call. = FALSE)
}

parameter_rows <- match(names(expected_step4_parameters), step4_parameter_table$Parameter)
if (anyNA(parameter_rows)) {
  stop(
    "The example Step 4 workbook is missing required parameter(s): ",
    paste(names(expected_step4_parameters)[is.na(parameter_rows)], collapse = ", "),
    call. = FALSE
  )
}
observed_step4_parameters <- suppressWarnings(as.numeric(
  step4_parameter_table$Value[parameter_rows]
))
names(observed_step4_parameters) <- names(expected_step4_parameters)
if (anyNA(observed_step4_parameters) ||
    !identical(unname(observed_step4_parameters),
               unname(as.numeric(expected_step4_parameters)))) {
  stop(
    "The example Step 4 workbook does not contain the validated V4/2x250 settings.\n",
    "Expected: ",
    paste(names(expected_step4_parameters), expected_step4_parameters,
          sep = "=", collapse = ", "),
    call. = FALSE
  )
}

step4_info <- openxlsx::read.xlsx(
  configuration_workbook,
  sheet = "Info",
  check.names = FALSE
)
if (!all(c("Parameter", "Value") %in% names(step4_info))) {
  stop("The example Step 4 Info sheet must contain Parameter and Value columns.", call. = FALSE)
}
step4_info_values <- setNames(as.character(step4_info$Value), step4_info$Parameter)
if (!identical(unname(step4_info_values[["sequencing_platform"]]),
               "Illumina MiSeq (2x250)") ||
    !identical(unname(step4_info_values[["target_region"]]), "16S V4")) {
  stop(
    "The example Step 4 workbook must identify Illumina MiSeq (2x250) and 16S V4.",
    call. = FALSE
  )
}

message(
  "Validated example Step 4 settings: Illumina MiSeq 2x250, 16S V4; ",
  "truncLen 240/160; maxEE 2/2; amplicon 250-256 bp."
)

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

taxonomy_files <- c(
  file.path(
    project_root, "tools", "trainsets", "SILVA",
    c(
      "silva_nr99_v138.2_toGenus_trainset.fa.gz",
      "silva_v138.2_assignSpecies.fa.gz"
    )
  ),
  file.path(
    project_root, "tools", "trainsets", "GTDB",
    c(
      "GTDB_bac120_arc53_ssu_r220_genus.fa.gz",
      "GTDB_bac120_arc53_ssu_r220_species.fa.gz"
    )
  )
)
if (!all(file.exists(taxonomy_files))) {
  stop(
    "The complete SILVA and GTDB trainsets are required. Run ",
    "setup/download_reference_databases.R before the example. Missing:\n  - ",
    paste(taxonomy_files[!file.exists(taxonomy_files)], collapse = "\n  - "),
    call. = FALSE
  )
}

required_executables <- c(
  Cutadapt = file.path(project_root, "tools", "cutadapt", "venv", "bin", "cutadapt"),
  FastQC = file.path(project_root, "tools", "FastQC", "fastqc"),
  MultiQC = file.path(project_root, "tools", "multiqc", "venv", "bin", "multiqc"),
  FastTree = file.path(project_root, "tools", "fasttree", "FastTree")
)
missing_executables <- names(required_executables)[
  file.access(required_executables, 1L) != 0L
]
if (length(missing_executables)) {
  stop(
    "Required project-local tool(s) missing or not executable: ",
    paste(missing_executables, collapse = ", "), ".\n",
    "Run setup/install_required_tools.R before the example.",
    call. = FALSE
  )
}

if (!nzchar(Sys.which("conda"))) {
  stop(
    "Conda is required for Step 7 but was not found on PATH. Install PICRUSt2 ",
    "as documented in setup/install_picrust2.sh before running the example.",
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

message(
  "Step 5 will read the validated Step 4 workbook directly from example data: ",
  configuration_workbook
)

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
render_step("2_fastqc_quality_reports.Rmd")
render_step("3_cutadapt_primer_trimming.Rmd")
render_step("5_dada2_pipeline.Rmd")
render_step("6_phylogenetic_tree.Rmd")
message(
  "\nThe bundled cell counts used by Steps 7-8 are synthetic test values and ",
  "must not be interpreted as measurements from the source study."
)
render_step("7_copy_number_correction.Rmd")
render_step("8_microbial_load_correction.Rmd")
render_step("9_phyloseq_object.Rmd")

required_dual_taxonomy_outputs <- c(
  file.path(run_results, "5_dada2_pipeline", c(
    "silva_taxonomy_table.csv",
    "gtdb_taxonomy_table.csv"
  )),
  file.path(run_results, "6_phylogenetic_tree", c(
    "phylogenetic_tree_SILVA_labeled.nwk",
    "phylogenetic_tree_GTDB_labeled.nwk",
    "phylogenetic_tree_SILVA.pdf",
    "phylogenetic_tree_GTDB.pdf"
  )),
  file.path(run_results, "9_phyloseq_object", "SILVA", "phyloseq_objects", c(
    "phyloseq_object_silva_raw_counts.RData",
    "phyloseq_object_silva_copy_number_corrected.RData",
    "phyloseq_object_silva_microbial_load_corrected.RData"
  )),
  file.path(run_results, "9_phyloseq_object", "GTDB", "phyloseq_objects", c(
    "phyloseq_object_gtdb_raw_counts.RData",
    "phyloseq_object_gtdb_copy_number_corrected.RData",
    "phyloseq_object_gtdb_microbial_load_corrected.RData"
  ))
)
missing_dual_taxonomy_outputs <- required_dual_taxonomy_outputs[
  !file.exists(required_dual_taxonomy_outputs)
]
if (length(missing_dual_taxonomy_outputs)) {
  stop(
    "The example run did not produce the complete SILVA + GTDB deliverable set:\n  - ",
    paste(missing_dual_taxonomy_outputs, collapse = "\n  - "),
    call. = FALSE
  )
}

message(
  "\nExample run complete.\n",
  "Completed: Steps 1-3, preconfigured Step 4 handoff, and Steps 5-9.\n",
  "Taxonomy: SILVA and GTDB, including both labelled trees and all six ",
  "database-by-abundance phyloseq objects.\n",
  "Generated results: ", run_results, "\n",
  "Bundled reference results: ", reference_results, "\n",
  "Compare the two trees without touching data/fastq or results/."
)
