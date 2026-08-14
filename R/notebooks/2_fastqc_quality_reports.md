Step 2: FastQC Quality Reports
================

- [Introduction](#introduction)
  - [Purpose](#purpose)
  - [Prerequisites](#prerequisites)
  - [What This Notebook Does](#what-this-notebook-does)
  - [Expected Input](#expected-input)
  - [Expected Output](#expected-output)
  - [FastQC Quality Modules](#fastqc-quality-modules)
- [Environment Setup](#environment-setup)
  - [Load Required Packages](#load-packages)
- [Configuration](#configuration)
  - [Define Path Parameters](#define-paths)
  - [Define FastQC Parameters](#define-fastqc-params)
- [Initialize Output Directory](#initialize-output-directory)
- [File Discovery](#file-discovery)
  - [Scan for FASTQ Files](#scan-files)
  - [Create File Inventory](#file-inventory)
- [FastQC Analysis](#fastqc-analysis)
  - [Run FastQC](#run-fastqc)
  - [Verify Generated Reports](#verify-reports)
  - [Run MultiQC](#run-multiqc)
- [Basic Statistics Extraction](#basic-statistics-extraction)
  - [Parse Basic Statistics](#parse-basic-stats)
  - [Build Basic Statistics Table](#build-stats-table)
- [QC Status Overview](#qc-status-overview)
  - [Summary Statistics](#summary-statistics)
  - [Save to Excel](#save-excel)
  - [Document and Export Column Dictionary](#column-dictionary)
  - [Cleanup](#cleanup)
  - [HTML Report Links](#report-links)
- [Output File Summary](#output-file-summary)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [Methods](#methods)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Understanding QC Module Status](#understanding-qc-module-status)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [Per Base Sequence Quality](#per-base-sequence-quality)
    - [Per Base Sequence Content](#per-base-sequence-content)
    - [Sequence Duplication Levels](#sequence-duplication-levels)
    - [Adapter Content](#adapter-content)
    - [Overrepresented Sequences](#overrepresented-sequences)

<style type="text/css">
/* Custom styling for improved readability */
body {
  font-size: 16px;
  line-height: 1.6;
}
&#10;h1, h2, h3 {
  color: #2c3e50;
  margin-top: 1.5em;
}
&#10;/* Code and code-output font size is pinned explicitly (rather than left to
   inherit from body) so that raising the prose size above does not also
   enlarge code blocks, inline code, or printed console output -- those stay
   at this project's original 14px regardless of any future prose-size
   changes. */
code {
  background-color: #F2F2F2;
  color: #2c3e50;
  padding: 2px 6px;
  border-radius: 3px;
  font-size: 14px;
}
&#10;pre, pre code {
  font-size: 14px;
}
&#10;.alert-info {
  background-color: #2c3e50;
  color: #ffffff;
  border-left: 4px solid #f39c12;
  padding: 12px;
  margin: 15px 0;
}
&#10;.alert-warning {
  background-color: #2c3e50;
  color: #ffffff;
  border-left: 4px solid #e74c3c;
  padding: 12px;
  margin: 15px 0;
}
&#10;.alert-success {
  background-color: #d4edda;
  color: #155724;
  border-left: 4px solid #28a745;
  padding: 12px;
  margin: 15px 0;
}
&#10;/* Table cells wrap long content (long inline code, sentences) instead of
   overflowing or forcing horizontal scroll. */
table {
  width: 100%;
}
&#10;th, td {
  white-space: normal;
  word-wrap: break-word;
  overflow-wrap: break-word;
}
&#10;</style>

# Introduction

## Purpose

This notebook is **Step 2** of the 16S rRNA sequencing data processing
pipeline. It generates comprehensive quality control reports for FASTQ
files using
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), a
widely-used tool for assessing the quality of high-throughput sequencing
data. Reviewing the FastQC and [MultiQC](https://multiqc.info/)
summaries produced here – before any trimming or filtering happens in
later steps – makes it possible to catch problems (unexpectedly low
quality scores, adapter/primer contamination, or a single sample far out
of line with the rest of the run) while they are still cheap to act on,
rather than discovering them only after DADA2’s error-learning step
further downstream.

## Prerequisites

Before running this notebook, ensure that:

1.  FASTQ files (`.fastq`, `.fq`, and their gzip-compressed variants)
    are present in the [data/fastq/](../../data/fastq/) input directory
2.  FastQC is installed in the [tools/FastQC/](../../tools/FastQC/) R
    project directory by running the
    [setup/install_required_tools.R](../../setup/install_required_tools.R)
    script
3.  MultiQC is installed in the [tools/multiqc/](../../tools/multiqc/) R
    project directory by running the
    [setup/install_required_tools.R](../../setup/install_required_tools.R)
    script
4.  **Step 1 (Data Integrity Check)** has been completed successfully
    (**optional**)

## What This Notebook Does

1.  **Discovers FASTQ Files**: Scans the input directory for `.fastq`,
    `.fq`, and gzip-compressed variants
2.  **Runs FastQC**: Generates HTML quality reports for each file using
    parallel processing
3.  **Runs MultiQC**: Aggregates all individual FastQC reports into a
    single interactive HTML summary
4.  **Extracts Basic Statistics**: Creates a summary table with key
    metrics per sample
5.  **Exports Results**: Saves the aggregated statistics to Excel for
    documentation

## Expected Input

- **Location**: FASTQ files should be placed in the
  [data/fastq/](../../data/fastq/) directory relative to your [project
  root](../../)
- **Format**: Standard FASTQ format (4-line records per read). `.fastq`,
  `.fq`, and their gzip-compressed variants are allowed,
  case-insensitively.
- **Naming Convention**: FastQC analyzes every matching file
  independently. For the paired summary, this notebook recognizes a
  terminal Illumina `R1`/`R2` token, preserves underscores in biological
  sample identifiers, and retains a full collision-safe `Pair_Key`
  shared by the two mates.

## Expected Output

- Individual HTML reports for each FASTQ file in
  [results/2_fastqc_quality_reports/FastQC/quality_reports/](../../results/2_fastqc_quality_reports/FastQC/quality_reports/)
- A single aggregated MultiQC HTML report
  ([multiqc_report.html](../../results/2_fastqc_quality_reports/MultiQC/multiqc_report.html))
  in
  [results/2_fastqc_quality_reports/MultiQC/](../../results/2_fastqc_quality_reports/MultiQC/)
- Aggregated summary Excel file with basic statistics in
  [results/2_fastqc_quality_reports/FastQC/fastqc_basic_statistics.xlsx](../../results/2_fastqc_quality_reports/FastQC/fastqc_basic_statistics.xlsx)
- Interactive summary tables in this notebook

## FastQC Quality Modules

FastQC evaluates sequencing data across 12 quality modules:

| Module | Description |
|:---|:---|
| **Basic Statistics** | Total reads, sequence length, GC content |
| **Per Base Sequence Quality** | Quality scores at each position |
| **Per Sequence Quality Scores** | Distribution of average quality scores |
| **Per Tile Sequence Quality** | Systematic errors by flow cell tile |
| **Per Base Sequence Content** | Base composition (A, T, C, G) at each position |
| **Per Sequence GC Content** | GC content distribution |
| **Per Base N Content** | Proportion of ambiguous bases (N) |
| **Sequence Length Distribution** | Fragment size distribution |
| **Sequence Duplication Levels** | PCR duplicate assessment |
| **Overrepresented Sequences** | Unusually frequent sequences |
| **Adapter Content** | Presence of adapter contamination |
| **K-mer Content** | Biases in short nucleotide sequences |

Detailed information can be found in the [FastQC
Documentation](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).

# Environment Setup

## Load Required Packages

The chunk below loads every R package this notebook depends on –
including
[fastqcr](https://cran.r-project.org/web/packages/fastqcr/readme/README.html)
for parsing FastQC archives and `data.table`/`DT` for building and
displaying the results tables below – plus the shared Excel,
column-dictionary, and output-link helpers sourced from
[R/functions/](../functions).

``` r
# data.table: High-performance data manipulation
# Provides fast operations for handling large QC result tables
library(data.table)

# dplyr: Data manipulation grammar
# Used for data transformation and summarization operations
library(dplyr)

# DT: Interactive tables in R Markdown
library(DT)

# fastqcr: Parse the ZIP reports produced by the FastQC executable
# Enables running FastQC and parsing results directly from R
library(fastqcr)

# fs: Cross-platform filesystem operations
# Consistent functions for file and directory manipulation
library(fs)

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# openxlsx: Excel file creation and manipulation
# Read and write Excel files without Java dependencies
library(openxlsx)

# parallel: Parallel processing support
# Built-in R package for multi-core processing
library(parallel)

# stringr: Consistent string manipulation
# Intuitive string processing functions
library(stringr)

# knitr: Dynamic report generation
# Used for global chunk options (knitr::opts_chunk$set() in the setup chunk above)
library(knitr)

# Source our custom Excel utility function from the project's function library
# This function handles creating or appending sheets to Excel workbooks
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source our column-dictionary utility function from the project's function library
# This function documents every column of a data frame exported to Excel
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source our output-links utility function from the project's function library
# This function renders clickable Markdown links to a chunk's output file(s)
source(here("R", "functions", "render_output_links_function.R"))

# Source the custom render_output_tree utility function from the project's
# function library. This function prints this notebook's entire output
# folder as a clickable directory tree (used in "Output File Summary" below),
# scanned live from disk at knit time so it always matches what was actually
# produced on this run.
source(here("R", "functions", "render_output_tree_function.R"))
```

------------------------------------------------------------------------

# Configuration

## Define Path Parameters

Set up the input and output directory paths. Using `here()` ensures
paths are relative to the [project root](../../), making the script
portable across different systems.

``` r
# Input folder containing raw FASTQ files to be analyzed
# This should be the same directory used in Step 1 (Data Integrity Check)
fastq_input_folder <- here("data", "fastq")

# Base results folder for all pipeline outputs
results_folder <- here("results")

# Base FastQC folder for Step 2 -- holds the aggregated statistics workbook
# directly, plus the nested quality_reports/ subfolder of individual
# per-file HTML/ZIP reports
# Numbered prefix maintains pipeline organization
fastqc_folder <- here(results_folder, "2_fastqc_quality_reports", "FastQC")

# Specific output folder for FastQC's individual per-file reports (Step 2),
# nested inside fastqc_folder
output_folder <- here(fastqc_folder, "quality_reports")

# Specific output folder for the MultiQC aggregated report (Step 2)
multiqc_output_folder <- here(results_folder, "2_fastqc_quality_reports", "MultiQC")

# Output Excel filename for aggregated statistics
output_excel_filename <- "fastqc_basic_statistics.xlsx"

# Name of the worksheet for FastQC basic statistics
output_excel_sheet <- "FastQC_Basic_Statistics"

# Construct the full path to the output Excel file -- written directly
# inside fastqc_folder (not inside output_folder/quality_reports), so it
# sits alongside quality_reports/ rather than among the individual reports
output_excel_path <- here(fastqc_folder, output_excel_filename)
```

## Define FastQC Parameters

Configure parameters for FastQC execution, including the FASTQ
file-matching pattern, how many CPU threads to use for parallel
processing, and the paths to the project-local FastQC and MultiQC
executables installed via the Prerequisites step above.

``` r
# Regular expression to match FASTQ file extensions
# Matches case-insensitively: .fastq, .fastq.gz, .fq, .fq.gz
fastq_extensions_regex <- "(?i)\\.(fastq|fq)(\\.gz)?$"

# Determine the number of CPU threads for parallel processing
# Reserve 2 cores for system operations to maintain responsiveness
# Use at least 1 thread even on single-core systems
detected_cores <- suppressWarnings(detectCores())
available_cores <- if (length(detected_cores) == 1L && is.finite(detected_cores) && detected_cores >= 1) {
  as.integer(detected_cores)
} else {
  1L
}
nr_threads <- min(4L, max(1L, available_cores - 2L))

# Path to FastQC executable (project-local installation)
fastqc_path <- here("tools", "FastQC", "fastqc")

# Path to MultiQC executable (project-local installation)
# MultiQC aggregates individual FastQC reports into one interactive summary
multiqc_path <- here("tools", "multiqc", "venv", "bin", "multiqc")

# Fail before creating outputs if an input or required executable is missing.
if (!dir_exists(fastq_input_folder)) {
  stop("Raw FASTQ input folder not found: ", fastq_input_folder,
       "\nComplete Step 1 and verify the project paths before running Step 2.")
}
tool_paths <- c(FastQC = fastqc_path, MultiQC = multiqc_path)
for (tool_name in names(tool_paths)) {
  tool_path <- tool_paths[[tool_name]]
  if (!file_exists(tool_path) || file.access(tool_path, mode = 1) != 0) {
    stop(tool_name, " executable not found or not executable: ", tool_path,
         "\nRun setup/install_required_tools.R.")
  }

  # Execute a harmless command as well: an executable bit alone does not catch
  # a relocated virtual environment with a stale interpreter shebang.
  version_result <- suppressWarnings(system2(
    tool_path, "--version", stdout = TRUE, stderr = TRUE
  ))
  version_status <- attr(version_result, "status")
  if (!is.null(version_status) && version_status != 0L) {
    stop(tool_name, " exists but cannot run:\n", paste(version_result, collapse = "\n"),
         "\nRe-run setup/install_required_tools.R.")
  }
}

# Display configuration
cat("FastQC / MultiQC Configuration:\n",
    "- Input folder:", fastq_input_folder, "\n",
    "- FastQC output folder:", output_folder, "\n",
    "- MultiQC output folder:", multiqc_output_folder, "\n",
    "- Available CPU cores:", available_cores, "\n",
    "- Threads for FastQC:", nr_threads, "\n",
    "- FastQC path:", fastqc_path, "\n",
    "- MultiQC path:", multiqc_path, "\n")
```

------------------------------------------------------------------------

# Initialize Output Directory

Create the output folder structure for FastQC reports.

``` r
# Create the main results folder if it doesn't exist
# This is shared across all pipeline steps
if (!dir_exists(results_folder)) {
  dir_create(results_folder, recurse = TRUE)
  cat("Created results directory:", results_folder, "\n")
} else {
  cat("Using existing results directory:", results_folder, "\n")
}

# Create the FastQC-specific output folder
# recurse = TRUE allows creation of nested directories
if (!dir_exists(output_folder)) {
  dir_create(output_folder, recurse = TRUE)
  cat("Created FastQC output directory:", output_folder, "\n")
} else {
  cat("Using existing FastQC output directory:", output_folder, "\n",
      "  Note: Existing reports may be overwritten.\n")
}

# Create the MultiQC-specific output folder
if (!dir_exists(multiqc_output_folder)) {
  dir_create(multiqc_output_folder, recurse = TRUE)
  cat("Created MultiQC output directory:", multiqc_output_folder, "\n")
} else {
  cat("Using existing MultiQC output directory:", multiqc_output_folder, "\n",
      "  Note: Existing report may be overwritten.\n")
}
```

------------------------------------------------------------------------

# File Discovery

## Scan for FASTQ Files

Discover all FASTQ files in the input directory.

``` r
# List all files in the input folder matching the FASTQ extension pattern
# recurse = FALSE: Only search the top-level directory
# type = "file": Only return files, not directories
# regexp: Filter using our FASTQ extension regex
all_fastq_paths <- dir_ls(
  path = fastq_input_folder,
  recurse = FALSE,
  type = "file",
  regexp = fastq_extensions_regex
)

# Validate that we found FASTQ files
if (length(all_fastq_paths) == 0) {
  stop("No FASTQ files found in: ", fastq_input_folder)
}

# Report discovery results
cat(
  "File discovery summary:\n",
  "- FASTQ files found:", length(all_fastq_paths), "\n",
  "- Input directory:", fastq_input_folder, "\n"
)
```

## Create File Inventory

Build a summary table of all discovered FASTQ files.

``` r
# Create a data.table with file information
file_inventory <- data.table(
  # Sequential file number for easy reference
  Nr = seq_along(all_fastq_paths),
  # Extract just the filename without directory path
  Filename = path_file(all_fastq_paths),
  # Calculate file size in megabytes for reference
  Size_MB = round(file_size(all_fastq_paths) / 1e6, 2),
  # Store full path for FastQC processing
  Full_Path = as.character(all_fastq_paths)
)

# Display the file inventory as an interactive, scrollable table rather
# than a single static block -- this list grows one row per FASTQ file, so
# a plain static table would dump the entire run's file list into the HTML
# report. DT::datatable() is used this way throughout the rest of this
# repository (see e.g. Steps 1, 3, 4, 5, 8), including the other tables
# further down this notebook; scrollY + scrollCollapse + paging = FALSE
# gives it a fixed-height scroll box instead of Prev/Next page controls.
datatable(
  file_inventory[, .(Nr, Filename, Size_MB)],
  colnames = c("#", "Filename", "Size (MB)"),
  options = list(pageLength = 10, scrollX = TRUE,
                 scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
  rownames = FALSE
)
```

------------------------------------------------------------------------

# FastQC Analysis

## Run FastQC

Execute FastQC on all discovered FASTQ files using parallel processing.

<div class="alert alert-info">

**Note**: This step may take several minutes depending on the number and
size of FASTQ files. Progress messages will be displayed during
execution!

</div>

``` r
# Remove only FastQC artefacts from an earlier run. Otherwise files from a
# different input set could be included in MultiQC and the final report links.
previous_fastqc_artifacts <- dir_ls(
  output_folder,
  type = "file",
  regexp = "_fastqc\\.(html|zip)$",
  fail = FALSE
)
if (length(previous_fastqc_artifacts) > 0L) {
  file_delete(previous_fastqc_artifacts)
}

# Record start time for performance tracking
start_time <- Sys.time()

# Invoke FastQC with the exact FASTQ inventory discovered above. The fastqcr
# directory wrapper passes every file in the directory to FastQC, including
# README files; using system2() here prevents non-FASTQ files from being
# mistaken for sequence data while fastqcr remains responsible for parsing
# the generated archives below.
fastqc_result <- system2(
  fastqc_path,
  args = c(
    "--threads", as.character(nr_threads),
    "--outdir", shQuote(as.character(output_folder)),
    shQuote(as.character(all_fastq_paths))
  ),
  stdout = TRUE,
  stderr = TRUE
)
fastqc_exit_status <- attr(fastqc_result, "status")
if (!is.null(fastqc_exit_status) && fastqc_exit_status != 0) {
  stop("FastQC failed with exit status ", fastqc_exit_status, ":\n",
       paste(tail(fastqc_result, 40), collapse = "\n"))
}

# Calculate elapsed time
end_time <- Sys.time()
elapsed_time <- difftime(end_time, start_time, units = "mins")

cat("\nFastQC analysis complete.\n",
    "Elapsed time:", round(as.numeric(elapsed_time), 2), "minutes\n")

# Link to the folder of individual FastQC HTML reports and ZIP archives
# just written above -- dozens of files are produced per run, so the whole
# folder is linked rather than every individual file.
render_output_links(output_folder, labels = "FastQC HTML reports and ZIP archives (folder)")
```

## Verify Generated Reports

Confirm that FastQC reports were successfully generated.

``` r
# List all generated HTML report files
fastqc_html_files <- dir_ls(
  path = output_folder,
  regexp = "_fastqc\\.html$"
) %>% sort()

# List all generated ZIP archive files (contain raw data)
fastqc_zip_files <- dir_ls(
  path = output_folder,
  regexp = "_fastqc\\.zip$"
) %>% sort()

# Both files are required at this point: HTML is the durable report and ZIP
# is parsed below before being removed during Cleanup.
if (length(fastqc_html_files) != length(all_fastq_paths) ||
    length(fastqc_zip_files) != length(all_fastq_paths)) {
  stop(
    "FastQC output mismatch: expected ", length(all_fastq_paths),
    " HTML and ZIP reports, found ", length(fastqc_html_files), " HTML and ",
    length(fastqc_zip_files), " ZIP files."
  )
}
```

## Run MultiQC

Aggregate all individual FastQC reports into a single interactive HTML
summary using [MultiQC](https://multiqc.info/).

<div class="alert alert-info">

**Note**: MultiQC scans the FastQC output folder automatically and
produces a single `multiqc_report.html` that allows cross-sample
comparison of all quality metrics in one place.

</div>

``` r
# Run MultiQC to aggregate all FastQC reports into one interactive HTML report.
# MultiQC scans output_folder for FastQC result files and combines them
# automatically — no need to specify individual files.
multiqc_start <- Sys.time()

result <- system2(
  multiqc_path,
  args = c(
    shQuote(as.character(output_folder)), # Directory to search for FastQC results
    "--outdir", shQuote(as.character(multiqc_output_folder)),
    "--filename", "multiqc_report",     # Output filename prefix (no extension needed)
    "--force"                           # Overwrite any existing report without prompting
  ),
  stdout = TRUE,
  stderr = TRUE
)

multiqc_elapsed <- difftime(Sys.time(), multiqc_start, units = "secs")

# Check exit status and report outcome
exit_status <- attr(result, "status")
multiqc_report_path <- here(multiqc_output_folder, "multiqc_report.html")
if ((!is.null(exit_status) && exit_status != 0) || !file_exists(multiqc_report_path)) {
  stop("MultiQC did not produce its expected report. Exit status: ",
       if (is.null(exit_status)) 0 else exit_status, "\n",
       paste(result, collapse = "\n"))
} else {
  cat("\nMultiQC report generated.\n",
      "  Elapsed time:", round(as.numeric(multiqc_elapsed), 1), "seconds\n",
      "  Report:", multiqc_report_path, "\n")

  # Link to the aggregated MultiQC report just written above -- only shown
  # when this branch actually ran, i.e. MultiQC succeeded.
  render_output_links(multiqc_report_path, labels = "MultiQC aggregated report (HTML)")
}
```

Detailed information can be found in the [MultiQC
Documentation](https://docs.seqera.io/multiqc).

------------------------------------------------------------------------

# Basic Statistics Extraction

## Parse Basic Statistics

Extract the “Basic Statistics” module data which contains key metrics
for each sample.

``` r
# Read detailed statistics from each FastQC zip file
# qc_read parses the full FastQC output including all metrics
# We'll use parallel processing to speed up reading multiple files

# Define a safe function to read QC stats from a single file
safe_qc_read <- function(zip_path) {
  tryCatch({
    list(path = zip_path, result = qc_read(zip_path), error = NA_character_)
  }, error = function(e) {
    list(path = zip_path, result = NULL, error = conditionMessage(e))
  })
}

# Read all QC results (this extracts detailed statistics)
cat("Extracting detailed statistics from", length(fastqc_zip_files), "reports...\n")

# Use lapply to read each zip file
# For very large datasets, consider using parallel::mclapply
qc_parse_results <- lapply(fastqc_zip_files, safe_qc_read)

# Never produce a silently incomplete workbook. Preserve every ZIP on failure
# and report the exact archives that require attention.
failed_qc_parses <- qc_parse_results[vapply(
  qc_parse_results, function(x) !is.na(x$error), logical(1)
)]
if (length(failed_qc_parses) > 0L) {
  stop(
    "Could not parse ", length(failed_qc_parses), " FastQC archive(s):\n",
    paste(vapply(failed_qc_parses, function(x) {
      paste0("  - ", basename(x$path), ": ", x$error)
    }, character(1)), collapse = "\n")
  )
}
qc_details_list <- lapply(qc_parse_results, `[[`, "result")
```

## Build Basic Statistics Table

Construct a comprehensive table of basic statistics for all samples.

``` r
# Extract basic statistics from each QC result
# Each qc_read result contains a $basic_statistics data frame
extract_basic_stats <- function(qc_result) {
  # Get the basic statistics component
  basic_stats <- qc_result$basic_statistics

  # Pivot from long to wide format for easier viewing
  # Original format: Measure | Value
  # Target format: One column per measure
  stats_wide <- data.table(
    Measure = basic_stats$Measure,
    Value = basic_stats$Value
  ) %>%
    # Transpose to wide format
    dcast(. ~ Measure, value.var = "Value") %>%
    # Remove the dummy grouping column
    select(-.)

  return(stats_wide)
}

# Apply extraction to all QC results and combine into single table
basic_stats_list <- lapply(qc_details_list, extract_basic_stats)

# Combine all rows into a single data.table
basic_stats_dt <- rbindlist(basic_stats_list, fill = TRUE)

# Clean up column names for better readability
# Replace spaces and special characters with underscores
original_names <- names(basic_stats_dt)
clean_names <- str_replace_all(original_names, "[^[:alnum:]]", "_")
clean_names <- str_replace_all(clean_names, "_+", "_")
clean_names <- str_remove(clean_names, "^_|_$")
setnames(basic_stats_dt, original_names, clean_names)

# Extract the sample ID and read direction from each row's Filename. Filename
# here is the raw FASTQ file FastQC actually analysed, as FastQC itself
# recorded it (e.g. "F3D0_S188_L001_R1_001.fastq") -- NOT the
# "*_fastqc.zip"/".html" report filename, so no "_fastqc" suffix is ever
# present here to strip.
#
# Derive direction only from the terminal Illumina read token, then retain a
# full pairing key. Unlike truncating at the first underscore, this preserves
# biological sample names that themselves contain underscores.
fastq_stem <- str_remove(basic_stats_dt$Filename, regex(fastq_extensions_regex))
direction_match <- str_match(fastq_stem, "_(R[12])(_[^_]*)?$")
basic_stats_dt[, Read_Direction := fifelse(
  direction_match[, 2] == "R1", "Forward",
  fifelse(direction_match[, 2] == "R2", "Reverse", NA_character_)
)]
basic_stats_dt[, Pair_Key := ifelse(
  is.na(direction_match[, 2]), fastq_stem,
  paste0(substr(fastq_stem, 1L, nchar(fastq_stem) - nchar(direction_match[, 1])),
         "_RX", fifelse(is.na(direction_match[, 3]), "", direction_match[, 3]))
)]
basic_stats_dt[, SampleID := str_remove(Pair_Key, "_S[0-9]+_L[0-9]+_RX.*$")]

duplicate_stat_rows <- basic_stats_dt[, .N, by = .(Pair_Key, Read_Direction)][N > 1L]
if (nrow(duplicate_stat_rows) > 0L) {
  stop("Multiple FastQC rows resolve to the same pairing key and direction:\n",
       paste0("  - ", duplicate_stat_rows$Pair_Key, " / ",
              duplicate_stat_rows$Read_Direction, collapse = "\n"))
}

# Reorder columns: SampleID first, Read_Direction immediately after it, then
# every other column in its existing order.
setcolorder(basic_stats_dt, c("SampleID", "Pair_Key", "Read_Direction",
                               setdiff(names(basic_stats_dt), c("SampleID", "Pair_Key", "Read_Direction"))))

# Convert FastQC's numeric-looking Basic Statistics fields from character to
# numeric. dcast() above pivots FastQC's raw Measure/Value text pairs
# straight into columns, so every field -- including these -- starts out as
# character regardless of its actual content; left as character, these
# would sort/filter as text rather than numbers in the exported Excel sheet
# (e.g. "9" > "80") and silently feed the wrong type into any downstream
# numeric calculation. Sequence_length is deliberately left as character,
# since it can be a range (e.g. "35-151") rather than a single number.

# Total_Sequences, Sequences_flagged_as_poor_quality, and GC are plain
# numbers in FastQC's own output (confirmed directly against this project's
# FastQC HTML reports, e.g. "Total Sequences</td><td>7793"); commas are
# stripped defensively before conversion in case a FastQC version ever adds
# thousands separators (e.g. "1,234,567"), even though none was observed.
plain_numeric_basic_statistics_columns <- c("Total_Sequences", "Sequences_flagged_as_poor_quality", "GC")
for (numeric_column in intersect(plain_numeric_basic_statistics_columns, names(basic_stats_dt))) {
  set(basic_stats_dt, j = numeric_column,
      value = as.numeric(gsub(",", "", basic_stats_dt[[numeric_column]])))
}

# Total_Bases needs its own parser rather than the plain conversion above:
# FastQC (confirmed against this project's own FastQC HTML reports, e.g.
# "Total Bases</td><td>1.9 Mbp" and "Total Bases</td><td>798.9 kbp" -- note
# the inconsistent capitalisation of the unit prefix between samples) writes
# this field as a human-readable string with a Kbp/Mbp/Gbp unit suffix
# rather than a raw base count. as.numeric() on a string like "1.9 Mbp"
# returns NA, which is why this column was coming out empty. This parses the
# numeric magnitude and unit prefix separately and rescales everything to a
# plain base-pair count, so Total_Bases ends up directly comparable across
# samples regardless of which unit FastQC chose to display it in for any
# given sample.
parse_fastqc_total_bases <- function(total_bases_text) {
  cleaned_text <- trimws(gsub(",", "", total_bases_text))

  # A small minority of FastQC versions report a plain integer (no unit
  # suffix, e.g. "1900000") -- handle that case first, since it needs no
  # unit rescaling.
  plain_number <- suppressWarnings(as.numeric(cleaned_text))
  if (!is.na(plain_number)) {
    return(plain_number)
  }

  # Otherwise expect "<number><optional space><optional K/M/G prefix>bp",
  # matched case-insensitively since FastQC does not capitalise consistently
  # (observed both "Mbp" and "kbp" within the same run's reports).
  matched_groups <- regmatches(
    cleaned_text,
    regexec("^([0-9.]+)\\s*([KMG]?)bp$", cleaned_text, ignore.case = TRUE)
  )[[1]]
  if (length(matched_groups) != 3) {
    return(NA_real_)
  }

  magnitude      <- as.numeric(matched_groups[2])
  unit_prefix    <- toupper(matched_groups[3])
  unit_multiplier <- switch(unit_prefix, "K" = 1e3, "M" = 1e6, "G" = 1e9, 1)

  magnitude * unit_multiplier
}

if ("Total_Bases" %in% names(basic_stats_dt)) {
  set(basic_stats_dt, j = "Total_Bases",
      value = vapply(basic_stats_dt$Total_Bases, parse_fastqc_total_bases, numeric(1)))
}

# Show a preview of the table
basic_stats_dt %>%
  select(-Filename) %>%  # Remove redundant filename column
  datatable(options = list(pageLength = 10, scrollX = TRUE,
                            scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
            rownames = FALSE,
            caption = "FastQC basic statistics per sample")
```

# QC Status Overview

## Summary Statistics

Calculate overall summary statistics across all samples.

``` r
# Calculate summary statistics from the basic stats table. Total_Sequences,
# GC, and Sequences_flagged_as_poor_quality are already numeric (converted
# in "Build Basic Statistics Table" above), so no per-use as.numeric()/comma
# stripping is needed here anymore.
total_reads <- sum(basic_stats_dt$Total_Sequences, na.rm = TRUE)
avg_gc <- mean(basic_stats_dt$GC, na.rm = TRUE)
total_poor_quality <- sum(basic_stats_dt$Sequences_flagged_as_poor_quality, na.rm = TRUE)

# Display summary
cat(
  "Overall QC summary:\n",
  "Sequencing statistics:\n",
  "- Biological sample/pair keys:", uniqueN(basic_stats_dt[!is.na(Read_Direction)]$Pair_Key), "\n",
  "- FASTQ files/read directions:", nrow(basic_stats_dt), "\n",
  "- Files with unclassified direction:", sum(is.na(basic_stats_dt$Read_Direction)), "\n",
  "- Total reads:", format(total_reads, big.mark = ","), "\n",
  "- Average GC content:", round(avg_gc, 1), "%\n",
  "- Total poor quality sequences:", format(total_poor_quality, big.mark = ","), "\n"
)
```

## Save to Excel

Export the basic statistics to the shared pipeline Excel file.

``` r
# Verify the output directory exists
excel_dir <- dirname(output_excel_path)
if (!dir_exists(excel_dir)) {
  dir_create(excel_dir, recurse = TRUE)
}

# This notebook owns the complete workbook. Rebuild it so obsolete sheets from
# an older schema cannot survive a rerun.
if (file_exists(output_excel_path)) {
  file_delete(output_excel_path)
}

# Add the basic statistics to this step's workbook. Re-running the notebook
# replaces the existing sheet instead of failing on a duplicate sheet name.
add_sheet_to_excel(
  workbook_path = output_excel_path,
  sheet_name = output_excel_sheet,
  data = basic_stats_dt,
  rownames = FALSE,
  overwrite = TRUE
)
```

## Document and Export Column Dictionary

Build a Sheet / Column / Explanation dictionary for every sheet written
to `output_excel_path`, so the workbook documents its own columns
without a hand-maintained description list that can drift out of sync
with the real exported data.

## Cleanup

This chunk deletes the per-file ZIP archives FastQC produced, now that
`qc_read()` has already extracted their statistics into `basic_stats_dt`
above – keeping only the HTML reports (and the Excel workbook already
written) as this notebook’s lasting output.

``` r
# Remove the intermediate ZIP archives; their contents were already
# extracted above, so keeping them around would only duplicate data already
# captured in the Excel export and HTML reports.
file.remove(fastqc_zip_files)
```

------------------------------------------------------------------------

## HTML Report Links

The following FastQC HTML reports are available in the [FastQC reports
folder](../../results/2_fastqc_quality_reports/FastQC/quality_reports/),
arranged in a table with each sample’s forward report on the left and
its reverse report on the right.

``` r
# Extract the same collision-safe SampleID, Pair_Key, and terminal direction
# used by the Basic Statistics table, so biological identifiers containing
# underscores are not collapsed when reports are aligned.
report_stem <- str_remove(basename(fastqc_html_files), "_fastqc\\.html$")
report_direction_match <- str_match(report_stem, "_(R[12])(_[^_]*)?$")
report_pair_key <- ifelse(
  is.na(report_direction_match[, 2]), report_stem,
  paste0(substr(report_stem, 1L, nchar(report_stem) - nchar(report_direction_match[, 1])),
         "_RX", fifelse(is.na(report_direction_match[, 3]), "", report_direction_match[, 3]))
)
report_links_dt <- data.table(
  SampleID = str_remove(report_pair_key, "_S[0-9]+_L[0-9]+_RX.*$"),
  Pair_Key = report_pair_key,
  Read_Direction = fifelse(
    report_direction_match[, 2] == "R1", "Forward",
    fifelse(report_direction_match[, 2] == "R2", "Reverse", NA_character_)
  ),
  Report_Path = fastqc_html_files
)

# Retain unclassified reports in their own visible table.
unclassified_reports <- report_links_dt[is.na(Read_Direction)]
if (nrow(unclassified_reports) > 0) {
  warning(
    "Could not determine forward/reverse direction for ", nrow(unclassified_reports),
    " report(s) from their terminal filename token; they are listed separately:\n",
    paste("  -", basename(unclassified_reports$Report_Path), collapse = "\n")
  )
}
report_links_dt <- report_links_dt[!is.na(Read_Direction)]

duplicate_report_keys <- report_links_dt[, .N, by = .(Pair_Key, Read_Direction)][N > 1L]
if (nrow(duplicate_report_keys) > 0L) {
  stop("Multiple HTML reports resolve to the same pairing key and direction:\n",
       paste0("  - ", duplicate_report_keys$Pair_Key, " / ",
              duplicate_report_keys$Read_Direction, collapse = "\n"))
}

# Reshape from one row per file to one row per sample, with separate
# Forward/Reverse path columns, so each sample's pair of reports lines up
# in the same row of the table below.
forward_links_dt <- report_links_dt[Read_Direction == "Forward", .(SampleID, Pair_Key, Forward_Path = Report_Path)]
reverse_links_dt  <- report_links_dt[Read_Direction == "Reverse", .(SampleID, Pair_Key, Reverse_Path = Report_Path)]
paired_links_dt   <- merge(forward_links_dt, reverse_links_dt,
                           by = c("SampleID", "Pair_Key"), all = TRUE)
setorder(paired_links_dt, SampleID, Pair_Key)

# Build a real HTML <a> tag for each column, using the report's own
# filename as the link text (e.g. "F3D0_S188_L001_R1_001_fastqc.html") --
# datatable(escape = FALSE) below leaves this markup intact instead of
# escaping it to literal "<a...>" text, which is what actually makes these
# clickable links. as.character() defensively coerces away fs_path's class
# first, since fs_path's own print/format methods are meant for console
# display, not for embedding in an HTML string. A missing report (an
# unpaired forward or reverse file) becomes plain "(not available)" text
# instead of a broken link.
encode_relative_href <- function(path) {
  relative <- as.character(fs::path_rel(path, start = here("R", "notebooks")))
  paste(vapply(strsplit(relative, "/", fixed = TRUE)[[1]],
               utils::URLencode, character(1), reserved = TRUE), collapse = "/")
}
build_report_link <- function(path) {
  if (is.na(path)) return("<em>(not available)</em>")
  paste0('<a href="', htmltools::htmlEscape(encode_relative_href(path), attribute = TRUE), '">',
         htmltools::htmlEscape(basename(path)), '</a>')
}
paired_links_dt[, Forward_Link := vapply(Forward_Path, build_report_link, character(1))]
paired_links_dt[, Reverse_Link := vapply(Reverse_Path, build_report_link, character(1))]

# Display as an interactive, scrollable HTML table -- the same datatable()
# presentation used for the Basic Statistics table above, so every table in
# this notebook looks and behaves the same way. escape = FALSE is required
# so the <a> tags above render as real clickable links instead of literal,
# escaped HTML text.
paired_links_dt[, .(Sample = SampleID, `Pair Key` = Pair_Key,
                    `Forward Report` = Forward_Link, `Reverse Report` = Reverse_Link)] %>%
  datatable(options = list(pageLength = 10, scrollX = TRUE,
                            scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
            rownames = FALSE,
            escape = FALSE,
            caption = "Paired FastQC HTML reports per sample")

if (nrow(unclassified_reports) > 0L) {
  unclassified_reports[, Report := vapply(Report_Path, build_report_link, character(1))]
  print(datatable(
    unclassified_reports[, .(`Unclassified Report` = Report)],
    rownames = FALSE, escape = FALSE,
    caption = "FastQC reports without a terminal R1/R2 token"
  ))
}
```

------------------------------------------------------------------------

# Output File Summary

The tree below lists every file this notebook has written to its own
output folder,
[results/2_fastqc_quality_reports/](../../results/2_fastqc_quality_reports/),
as a clickable, portable link (relative to this notebook’s own location)
with a short description – built live from what is actually on disk at
knit time via the project’s shared
[render_output_tree_function.R](../functions/render_output_tree_function.R)
helper, so it always matches this run’s real output rather than a
hand-maintained list.

------------------------------------------------------------------------

# Recommended Next Step

The FastQC and MultiQC reports above are diagnostic only – nothing
downstream reads them programmatically – so once you have reviewed them
(paying particular attention to per-base quality, adapter content, and
any sample that stands out from the rest), proceed to [Step 3 — Cutadapt
Primer Trimming](3_cutadapt_primer_trimming.md), which lists this
notebook as an optional but recommended prerequisite. Reviewing quality
here first lets you set informed trimming parameters in Step 3 and gives
you a pre-trimming baseline to compare its own before/after primer
counts against.

------------------------------------------------------------------------

# Session Information

Record the R environment for reproducibility.

------------------------------------------------------------------------

# References

## Methods

- Andrews S (2010). FastQC: A Quality Control Tool for High Throughput
  Sequence Data. Babraham Bioinformatics.
  <https://www.bioinformatics.babraham.ac.uk/projects/fastqc/>
- Ewels P, Magnusson M, Lundin S, Käller M (2016). MultiQC: summarize
  analysis results for multiple tools and samples in a single report.
  *Bioinformatics*, 32(19), 3047-3048.
  <https://doi.org/10.1093/bioinformatics/btw354> ([MultiQC
  GitHub](https://github.com/MultiQC/MultiQC))
- [fastqcr R
  Package](https://cran.r-project.org/web/packages/fastqcr/readme/README.html)
  — parser used for the FastQC ZIP output in this notebook.

## Related

- [Step 1 — Data Integrity Check](1_data_integrity_check.md) — optional
  prerequisite; confirms FASTQ file integrity before quality is assessed
  here.
- [Step 3 — Cutadapt Primer Trimming](3_cutadapt_primer_trimming.md) —
  this notebook’s recommended next step.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Understanding QC Module Status

FastQC assigns one of three status levels to each module:

| Status   | Symbol | Interpretation                            |
|:---------|:-------|:------------------------------------------|
| **PASS** | ✓      | Results fall within normal/expected range |
| **WARN** | ⚠      | Results are unusual but may be acceptable |
| **FAIL** | ✗      | Results indicate a potential problem      |

## Common Issues and Solutions

### Per Base Sequence Quality

**Warning/Fail**: Quality scores drop significantly toward the end of
reads.

**Action**: Consider trimming low-quality bases from read ends using
tools like [Cutadapt](https://cutadapt.readthedocs.io/en/stable/)

### Per Base Sequence Content

**Warning/Fail**: Non-uniform base composition at read starts.

**Context**: This is common and often expected in RNA-Seq data due to
random hexamer priming bias. Usually not a concern for 16S data.

### Sequence Duplication Levels

**Warning/Fail**: High levels of duplicate sequences.

**Context**: Some duplication is expected in amplicon sequencing (16S).
Very high levels may indicate PCR over-amplification.

### Adapter Content

**Warning/Fail**: Adapter sequences detected in reads.

**Action**: Remove adapters using trimming tools. This commonly occurs
when insert size is shorter than read length.

### Overrepresented Sequences

**Warning**: Specific sequences appear more frequently than expected.

**Action**: Check if overrepresented sequences match adapters, primers,
or known contaminants using
[BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi).
