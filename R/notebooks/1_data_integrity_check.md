Step 1: FASTQ Data Integrity Check
================

- [Introduction](#introduction)
  - [Purpose](#purpose)
  - [Prerequisites](#prerequisites)
  - [What This Notebook Does](#what-this-notebook-does)
  - [Expected Input](#expected-input)
  - [Expected Output](#expected-output)
- [Environment Setup](#environment-setup)
  - [Load Required Packages](#load-packages)
- [Configuration](#configuration)
  - [Define Path Parameters](#define-paths)
  - [Define File Pattern Parameters](#define-patterns)
- [Define Helper Functions](#define-helper-functions)
  - [Sample ID Extraction](#func-sample-id)
  - [Pair Key Derivation](#func-pair-key)
  - [Complete FASTQ Inspection](#func-fastq-inspection)
  - [File Size Calculation](#func-file-size)
  - [Output Table Rendering](#func-render-table)
- [File Discovery](#file-discovery)
  - [Scan for FASTQ Files](#scan-files)
  - [Create File Inventory Table](#create-inventory)
- [Pair Classification](#pair-classification)
  - [Identify Read Direction](#classify-direction)
  - [Separate and Validate](#separate-validate)
- [Pair Matching](#pair-matching)
  - [Generate Pairing Keys](#generate-keys)
  - [Merge Paired Files](#merge-pairs)
  - [Report Pairing Statistics](#pairing-stats)
- [Data Integrity Analysis](#data-integrity-analysis)
  - [Compute Quality Metrics](#compute-metrics)
  - [Verify Read Count Concordance](#verify-concordance)
- [Results Summary](#results-summary)
  - [Generate Summary Statistics](#summary-stats)
  - [Preview Results Table](#preview-results)
- [Export Results](#export-results)
  - [Save to Excel](#save-excel)
  - [Document and Export Column Dictionary](#column-dictionary)
- [Output File Summary](#output-file-summary)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [Methods](#methods)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [No FASTQ Files Found](#no-fastq-files-found)
    - [Mismatched Read Counts](#mismatched-read-counts)
    - [Ambiguous File Classification](#ambiguous-file-classification)
    - [Invalid FASTQ Files](#invalid-fastq-files)

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
&#10;  background-color: #2c3e50;
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
</style>

# Introduction

## Purpose

This notebook is **Step 1** of the 16S rRNA sequencing workflow. It
performs an optional but recommended structural integrity check of
paired-end [FASTQ](https://doi.org/10.1093/nar/gkp1137) files before
downstream processing. It confirms that files are readable and
parseable, matches forward and reverse files by name, counts complete
FASTQ records, and compares paired read counts. It does not replace
checksum verification against values supplied by a sequencing facility
or data repository.

## Prerequisites

Before running this notebook, ensure that your raw FASTQ files
(`.fastq`, `.fq`, and their gzip-compressed variants) are present in the
[data/fastq/](../../data/fastq/) input directory.

## What This Notebook Does

1.  **File Discovery**: Scans the input folder for `.fastq`, `.fq`, and
    gzip-compressed variants
2.  **Pair Matching**: Identifies and matches forward and reverse reads
    using the configured filename tokens.
3.  **Integrity Validation**: For each file, the notebook:
    - Parses the complete file and counts its FASTQ records
    - Measures file size in megabytes
    - Reports whether the complete file can be parsed as FASTQ
4.  **Concordance Check**: Verifies that paired files contain the same
    number of reads. Equal counts support correct pairing but do not
    prove sample identity.
5.  **Report Generation**: Exports a summary table to an Excel file for
    documentation and review

## Expected Input

- **Location**: FASTQ files should be placed in the
  [data/fastq/](../../data/fastq/) directory relative to your [project
  root](../../)
- **Naming Convention**: The default configuration expects:
  - Forward reads: `{sample_id}_L001_R1_001` followed by `.fastq`,
    `.fq`, or a gzip-compressed variant
  - Reverse reads: `{sample_id}_L001_R2_001` followed by `.fastq`,
    `.fq`, or a gzip-compressed variant
  - Edit `forward_token` and `reverse_token` in the configuration
    section when your files use another consistent convention.
- **Format**: Standard FASTQ format (4-line records per read)

## Expected Output

- An Excel file
  ([data_integrity_check.xlsx](../../results/1_data_integrity_check/data_integrity_check.xlsx))
  saved to
  [results/1_data_integrity_check/](../../results/1_data_integrity_check/)
- Console output with pairing statistics and any warnings

------------------------------------------------------------------------

# Environment Setup

## Load Required Packages

The chunk below loads every R package this notebook depends on, plus
four helper functions shared across the pipeline: Excel writing,
column-dictionary documentation, clickable output links, and the
output-tree summary (each sourced from [R/functions/](../functions/)).

``` r
# data.table: High-performance data manipulation package
# Provides fast operations on large datasets using reference semantics
library(data.table)

# DT: Interactive data tables in R Markdown
library(DT)

# fs: Cross-platform filesystem operations
# Offers consistent, intuitive functions for file and directory manipulation
library(fs)

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# openxlsx: Excel file creation and manipulation
# Allows reading and writing Excel files without requiring Java dependencies
library(openxlsx)

# ShortRead: Bioconductor package for FASTQ processing
# Provides efficient tools for reading, validating, and summarizing FASTQ files
library(ShortRead)

# stringr: Consistent string manipulation functions
# Part of the tidyverse, offers intuitive string processing functions
library(stringr)

# Source the custom Excel utility function from the project's function library
# This function handles creating or appending sheets to Excel workbooks
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source the custom column-dictionary utility function from the project's
# function library. This function documents every column of a sheet being
# exported, so a trailing Column_Dictionary sheet can be added to the workbook
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source the custom render_output_links utility function from the project's
# function library. This function prints clickable Markdown links to output
# files/folders written by later chunks, for reproducibility across machines
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
# Input folder containing raw FASTQ files to be checked
# The here() function constructs the path relative to the project root
fastq_input_folder <- here("data", "fastq")

# Base results folder for all pipeline outputs
results_folder <- here("results")

# Output folder where the integrity check results will be saved
# Results are organized in a numbered folder for pipeline step tracking
output_folder <- fs::path(results_folder, "1_data_integrity_check")

# Create the output directory and any necessary parent directories
# If the directory already exists, this operation has no effect (safe to re-run)
dir_create(output_folder, recurse = TRUE)

# Name of the Excel file that will contain the integrity report
output_excel_filename <- "data_integrity_check.xlsx"

# Name of the worksheet within the Excel file
output_excel_sheet <- "Data_Integrity_Check"

# Construct the full path to the output Excel file
output_excel_path <- fs::path(output_folder, output_excel_filename)
```

## Define File Pattern Parameters

Configure the regular expression and tokens used to identify and pair
FASTQ files.

``` r
# Regular expression to match FASTQ file extensions
# Matches, case-insensitively: .fastq, .fastq.gz, .fq, .fq.gz
# The \\. escapes the literal dot character
# The (fastq) captures the extension name
# The (\\.gz)? optionally matches gzip compression suffix
fastq_extensions_regex <- "(?i)\\.(fastq|fq)(\\.gz)?$"

# Token identifying forward (Read 1) files in paired-end sequencing
# Common conventions include: _1, _R1, .1, _read1
forward_token <- "_L001_R1_001"

# Token identifying reverse (Read 2) files in paired-end sequencing
# Common conventions include: _2, _R2, .2, _read2
reverse_token <- "_L001_R2_001"
```

# Define Helper Functions

## Sample ID Extraction

This function extracts the display identifier used consistently by this
workflow. By convention, it uses the first underscore-delimited filename
component; change this configuration if underscores are part of your
biological sample IDs or change the SampleIDs.

``` r
# Extract Sample ID from FASTQ Filename
#
# Parses a FASTQ filename and returns the sample identifier portion.
# Uses the workflow's first-component naming convention.
#
# Arguments:
#   filename_basename - Character string of the filename (without path)
# Returns:
#   Character string containing the sample ID
#
# Examples:
#   extract_sample_id("Sample001_1.fastq.gz")  # Returns "Sample001"
#   extract_sample_id("WT-Control_S1_L001_R1_001.fastq.gz")  # Returns "WT-Control"
extract_sample_id <- function(filename_basename) {
  # Split the filename at underscores and return the first component.
  str_split(filename_basename, "_", simplify = TRUE)[, 1]
}
```

## Pair Key Derivation

This function creates a matching key to pair forward and reverse reads.
It replaces the read direction token with a generic placeholder,
allowing forward and reverse files to share the same key.

``` r
# Derive Pairing Key from FASTQ Filename
#
# Creates a standardized key for matching forward and reverse read files.
# Replaces the direction token (_1 or _2) with a placeholder (_RX_).
#
# Arguments:
#   filename_basename - Character string of the filename (without path)
#   token              - Character string of the token to replace (_1 or _2)
# Returns:
#   Character string containing the pairing key
#
# Examples:
#   derive_pair_key("Sample001_1.fastq.gz", "_1")  # Returns "Sample001_RX_.fastq.gz"
#   derive_pair_key("Sample001_2.fastq.gz", "_2")  # Returns "Sample001_RX_.fastq.gz"
derive_pair_key <- function(filename_basename, token) {
  # Compression/extension is deliberately excluded: mates remain a pair when
  # one is compressed and the other is not.
  stem <- str_remove(filename_basename, regex(fastq_extensions_regex))
  ifelse(
    endsWith(stem, token),
    paste0(substr(stem, 1L, nchar(stem) - nchar(token)), "_RX_"),
    stem
  )
}
```

## Complete FASTQ Inspection

This function parses an entire FASTQ file with
[ShortRead](https://doi.org/10.1093/bioinformatics/btp450). A successful
parse supplies the record count and marks the file as valid; malformed
records, truncated gzip streams, empty files, missing files, and other
read errors are reported as invalid.

``` r
# Inspect One FASTQ File
#
# Arguments:
#   fastq_path - Full path to one FASTQ file
#
# Returns:
#   A list with:
#   - reads: number of complete FASTQ records, or NA when parsing fails
#   - valid: TRUE only when the non-empty file parses completely
#   - issue: parsing or file-access problem, or NA for a valid file
inspect_fastq_file <- function(fastq_path) {
  invalid_result <- function(issue) {
    list(reads = NA_real_, valid = FALSE, issue = issue)
  }

  if (length(fastq_path) != 1 || is.na(fastq_path) || !file_exists(fastq_path)) {
    return(invalid_result("File is missing"))
  }

  size_bytes <- tryCatch(
    as.numeric(file_info(fastq_path)$size),
    error = function(e) NA_real_
  )
  if (!is.finite(size_bytes) || size_bytes <= 0) {
    return(invalid_result("File is empty or its size could not be read"))
  }

  tryCatch({
    counts <- countFastq(fastq_path)
    reads <- as.numeric(counts$records[[1]])

    if (!is.finite(reads)) {
      invalid_result("FASTQ record count was unavailable")
    } else {
      list(reads = reads, valid = TRUE, issue = NA_character_)
    }
  }, error = function(e) {
    invalid_result(conditionMessage(e))
  })
}
```

## File Size Calculation

Calculate file size in megabytes with error handling for missing files.

``` r
# Safely Calculate File Size in Megabytes
#
# Arguments:
#   file_path - Full path to one file
#
# Returns:
#   File size in decimal megabytes, or NA if the file is unavailable
safe_file_size_mb <- function(file_path) {
  if (length(file_path) != 1 || is.na(file_path) || !file_exists(file_path)) {
    return(NA_real_)
  }

  tryCatch(
    as.numeric(file_info(file_path)$size) / 1e6,
    error = function(e) NA_real_
  )
}
```

## Output Table Rendering

Render interactive tables in HTML and ordinary Markdown tables for
non-HTML output formats.

``` r
render_notebook_table <- function(data, caption = NULL) {
  output_target <- knitr::opts_knit$get("rmarkdown.pandoc.to")

  if (identical(output_target, "html")) {
    datatable(
      data,
      caption = caption,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        scrollY = "400px",
        scrollCollapse = TRUE,
        paging = FALSE
      )
    )
  } else {
    knitr::kable(data, caption = caption)
  }
}
```

------------------------------------------------------------------------

# File Discovery

## Scan for FASTQ Files

Scan the input directory for all FASTQ files matching our given
extension patterns for forward and reverse reads.

``` r
# Fail clearly when the configured input directory does not exist.
if (!dir_exists(fastq_input_folder)) {
  stop("FASTQ input directory does not exist: ", fastq_input_folder)
}

# List all files in the input folder matching the FASTQ extension pattern
# recurse = FALSE: Only search the top-level directory (no subdirectories)
# type = "file": Only return files, not directories or symlinks
# regexp: Filter files using our FASTQ extension regular expression
all_fastq_paths <- dir_ls(
  path = fastq_input_folder,
  recurse = FALSE,
  type = "file",
  regexp = fastq_extensions_regex
)

# Validate that we found at least one FASTQ file
# Stop execution with an informative error if none found
if (length(all_fastq_paths) == 0) {
  stop("No FASTQ files found in: ", fastq_input_folder)
}

# Report the number of files found
cat("Found", length(all_fastq_paths), "FASTQ files in", fastq_input_folder, "\n")
```

## Create File Inventory Table

Build a data.table containing file paths and metadata for all discovered
FASTQ files.

``` r
# Create a data.table with file path information
# data.table provides fast operations for our subsequent processing
all_fastq_dt <- data.table(
  # Store just the filename (no directory) for pattern matching
  filename = path_file(all_fastq_paths),
  # Store the absolute (full) path for reliable file access
  absolute_path = path_abs(all_fastq_paths)
)

# Preview the discovered files.
render_notebook_table(all_fastq_dt)
```

------------------------------------------------------------------------

# Pair Classification

## Identify Read Direction

Classify each file as forward (Read 1) or reverse (Read 2) based on the
filename tokens.

``` r
# Require the configured token at the end of the filename stem, immediately
# before an accepted FASTQ extension. This prevents tokens embedded in a sample
# name from being mistaken for read-direction markers.
# Note: data.table's `:=` updates the table by reference and returns it
# invisibly, so these two assignments do not print anything on their own --
# the fully-assembled table is instead previewed once, later in this
# notebook, rather than being re-rendered as a full interactive table after
# every single column is added.
all_fastq_dt[, filename_stem := str_remove(filename, regex(fastq_extensions_regex))]
all_fastq_dt[, is_forward := endsWith(filename_stem, forward_token)]
all_fastq_dt[, is_reverse := endsWith(filename_stem, reverse_token)]

# Summary of classification
cat("Classification summary:\n",
    "- Forward reads:", sum(all_fastq_dt$is_forward), "\n",
    "- Reverse reads:", sum(all_fastq_dt$is_reverse), "\n")
```

## Separate and Validate

Separate files into forward and reverse sets, and check for ambiguous
classifications.

``` r
# Extract forward reads: a filename must contain only the forward token.
forward_dt <- all_fastq_dt[is_forward == TRUE & is_reverse == FALSE]

# Extract reverse reads: a filename must contain only the reverse token.
reverse_dt <- all_fastq_dt[is_reverse == TRUE & is_forward == FALSE]

# Identify problematic files that match both or neither pattern
# is_forward == is_reverse captures both cases:
#   - TRUE == TRUE: file has both tokens (ambiguous)
#   - FALSE == FALSE: file has neither token (unpaired)
ambiguous_dt <- all_fastq_dt[is_forward == is_reverse]

# Warn about ambiguous files if any are found
if (nrow(ambiguous_dt) > 0) {
  warning(
    "Found ", nrow(ambiguous_dt), " file(s) that match both or neither pattern:\n",
    paste("  -", ambiguous_dt$filename, collapse = "\n")
  )
}

# Report file counts
cat("\nFile classification results:\n",
    "- Forward read files:", nrow(forward_dt), "\n",
    "- Reverse read files:", nrow(reverse_dt), "\n",
    "- Ambiguous files:", nrow(ambiguous_dt), "\n")
```

# Pair Matching

## Generate Pairing Keys

Create standardized keys that allow matching of forward and reverse
files.

``` r
# Replace the configured forward token with a common placeholder.
forward_dt[, pair_key := derive_pair_key(filename, forward_token)]

# Apply the same replacement to reverse files so mates share a pair_key.
reverse_dt[, pair_key := derive_pair_key(filename, reverse_token)]

# A many-to-many merge would manufacture pair combinations. Reject duplicate
# keys explicitly instead of allowing an ambiguous report.
duplicate_forward_keys <- unique(forward_dt$pair_key[duplicated(forward_dt$pair_key)])
duplicate_reverse_keys <- unique(reverse_dt$pair_key[duplicated(reverse_dt$pair_key)])
if (length(duplicate_forward_keys) > 0 || length(duplicate_reverse_keys) > 0) {
  stop(
    "Multiple FASTQ files resolve to the same pairing key: ",
    paste(unique(c(duplicate_forward_keys, duplicate_reverse_keys)), collapse = ", ")
  )
}
```

## Merge Paired Files

Join forward and reverse file tables to create a unified paired dataset.

``` r
# Perform a full outer join on the pair_key column
# all = TRUE ensures all files are kept (even unpaired ones)
# suffixes distinguish columns from each source table
paired_dt <- merge(
  forward_dt,
  reverse_dt,
  by = "pair_key",
  all = TRUE,               # Keep all records (full outer join)
  suffixes = c("_fwd", "_rev")
)

# Extract sample ID from whichever filename is available
# fifelse is data.table's fast conditional function
# Priority: use forward filename if available, otherwise use reverse
paired_dt[, sample := fifelse(
  !is.na(filename_fwd),
  extract_sample_id(filename_fwd),
  fifelse(!is.na(filename_rev), extract_sample_id(filename_rev), NA_character_))]

duplicate_sample_ids <- unique(paired_dt$sample[duplicated(paired_dt$sample)])
if (length(duplicate_sample_ids) > 0) {
  warning(
    "Shortened SampleID values are not unique: ",
    paste(duplicate_sample_ids, collapse = ", "),
    ". Full pairing keys are retained in the report."
  )
}
```

## Report Pairing Statistics

Summarize the pairing results, identifying any orphaned files.

``` r
# Count fully paired samples (both forward and reverse present)
n_paired <- sum(!is.na(paired_dt$filename_fwd) & !is.na(paired_dt$filename_rev))

# Count samples with only forward read (missing reverse)
n_fwd_only <- sum(!is.na(paired_dt$filename_fwd) & is.na(paired_dt$filename_rev))

# Count samples with only reverse read (missing forward)
n_rev_only <- sum(is.na(paired_dt$filename_fwd) & !is.na(paired_dt$filename_rev))

# Display pairing summary
cat(
  "Pairing summary:\n",
  "- Complete pairs:", n_paired, "\n",
  "- Forward only:", n_fwd_only, "\n",
  "- Reverse only:", n_rev_only, "\n"
)
```

------------------------------------------------------------------------

# Data Integrity Analysis

## Compute Quality Metrics

For each file, parse the complete FASTQ once, collect its record count,
calculate its size, and record whether parsing succeeded.

<div class="alert alert-info">

**Note**: This step may take several minutes for large datasets because
every FASTQ record is parsed.

</div>

``` r
# Parse each existing file once; countFastq() scans the complete file and raises
# an error when it encounters malformed FASTQ records.
forward_inspections <- lapply(paired_dt$absolute_path_fwd, inspect_fastq_file)
reverse_inspections <- lapply(paired_dt$absolute_path_rev, inspect_fastq_file)

# Add all derived columns in one by-reference update so the intermediate table
# is not repeatedly printed.
paired_dt[, `:=`(
  reads_fwd = vapply(forward_inspections, `[[`, numeric(1), "reads"),
  reads_rev = vapply(reverse_inspections, `[[`, numeric(1), "reads"),
  size_mb_fwd = vapply(paired_dt$absolute_path_fwd, safe_file_size_mb, numeric(1)),
  size_mb_rev = vapply(paired_dt$absolute_path_rev, safe_file_size_mb, numeric(1)),
  valid_fwd = vapply(forward_inspections, `[[`, logical(1), "valid"),
  valid_rev = vapply(reverse_inspections, `[[`, logical(1), "valid"),
  issue_fwd = vapply(forward_inspections, `[[`, character(1), "issue"),
  issue_rev = vapply(reverse_inspections, `[[`, character(1), "issue")
)]

# A missing mate is represented as NA rather than an invalid file, because no
# file was available to inspect in that direction.
paired_dt[is.na(absolute_path_fwd), `:=`(valid_fwd = NA, issue_fwd = "Forward mate missing")]
paired_dt[is.na(absolute_path_rev), `:=`(valid_rev = NA, issue_rev = "Reverse mate missing")]

# Report parse failures for files that were present.
invalid_files <- rbind(
  paired_dt[!is.na(absolute_path_fwd) & valid_fwd == FALSE,
            .(sample, direction = "Forward", file = filename_fwd, issue = issue_fwd)],
  paired_dt[!is.na(absolute_path_rev) & valid_rev == FALSE,
            .(sample, direction = "Reverse", file = filename_rev, issue = issue_rev)]
)
if (nrow(invalid_files) > 0) {
  warning(
    "Found ", nrow(invalid_files), " FASTQ file(s) that could not be parsed completely:\n",
    paste0("  - ", invalid_files$file, ": ", invalid_files$issue, collapse = "\n")
  )
}

cat("FASTQ inspection complete.\n")
```

## Verify Read Count Concordance

Check that paired files have matching read counts—a critical quality
control step.

``` r
# Compare read counts between paired files
# Only compare when both values are available (not NA)
# Equal counts are expected for a complete pair, although they do not prove
# that the files came from the same biological sample.
paired_dt[, reads_match := fifelse(
  !is.na(reads_fwd) & !is.na(reads_rev),  # Only compare if both counts exist
  reads_fwd == reads_rev,                 # TRUE if counts match
  NA)]                                    # NA if either count is missing

# One explicit status makes every row actionable in the workbook.
paired_dt[, integrity_status := fifelse(
  is.na(filename_fwd), "FAIL: reverse file has no forward mate",
  fifelse(
    is.na(filename_rev), "FAIL: forward file has no reverse mate",
    fifelse(
      valid_fwd == FALSE | valid_rev == FALSE, "FAIL: FASTQ parsing failed",
      fifelse(reads_match == FALSE, "FAIL: paired read counts differ", "PASS")
    )
  ))]

# Identify samples with mismatched read counts
mismatched <- paired_dt[reads_match == FALSE]

# Only render the mismatch table (and warn) when there is something to show --
# an empty interactive table adds noise without conveying anything a reader
# needs, since "no mismatches" is already good news on its own.
if (nrow(mismatched) > 0) {
  warning(
    "Found ", nrow(mismatched), " paired sample(s) with mismatched read counts:\n",
    paste("  -", mismatched$sample, collapse = "\n")
  )
  render_notebook_table(
    mismatched,
    caption = "Samples with mismatched forward/reverse read counts"
  )
} else if (sum(!is.na(paired_dt$reads_match)) == 0L) {
  cat("No paired files had two available read counts to compare.\n")
} else {
  cat("No mismatched read counts -- all paired samples have concordant forward/reverse counts.\n")
}
```

------------------------------------------------------------------------

# Results Summary

## Generate Summary Statistics

Display an overview of the integrity check results.

## Preview Results Table

View the detailed results for each sample.

``` r
# Create a formatted export table with user-friendly column names
# Select and rename columns for the final report
export_dt <- paired_dt[, .(
  SampleID = sample,
  `Pair Key` = pair_key,
  `Forward File` = filename_fwd,
  `Reads (Forward)` = reads_fwd,
  `Size (MB, Forward)` = round(size_mb_fwd, 2),
  `Forward Issue` = issue_fwd,
  `Reverse File` = filename_rev,
  `Reads (Reverse)` = reads_rev,
  `Size (MB, Reverse)` = round(size_mb_rev, 2),
  `Reverse Issue` = issue_rev,
  `Valid (Forward)` = valid_fwd,
  `Valid (Reverse)` = valid_rev,
  `Reads Match` = reads_match,
  `Integrity Status` = integrity_status,
  `Overall Status` = overall_integrity_status
)]

# Sort by sample name for consistent ordering
setorder(export_dt, SampleID)

# Display the results table.
render_notebook_table(export_dt)

# Every discovered input file, including files that matched both/neither token,
# is retained in a separate inventory rather than disappearing after warnings.
file_inventory_dt <- copy(all_fastq_dt)[, .(
  Filename = filename,
  `Absolute Path` = absolute_path,
  Classification = fifelse(
    is_forward & !is_reverse, "Forward",
    fifelse(is_reverse & !is_forward, "Reverse",
            fifelse(is_forward & is_reverse, "Ambiguous: both tokens", "Unclassified: neither token"))
  )
)]
setorder(file_inventory_dt, Filename)
render_notebook_table(file_inventory_dt, caption = "Inventory of every discovered FASTQ file")
```

------------------------------------------------------------------------

# Export Results

## Save to Excel

Export the integrity check results to an Excel file for documentation
and sharing, using the project’s shared Excel-writing helper,
[add_sheet_to_excel_function.R](../functions/add_sheet_to_excel_function.R).

``` r
# This notebook owns the complete workbook. Start each run clean so sheets from
# an older schema cannot survive unnoticed or fall outside the dictionary.
if (file_exists(output_excel_path)) {
  file_delete(output_excel_path)
}

# Add the results as a sheet in the Excel workbook
# Replace this notebook's sheet when the workbook already exists, making the
# notebook safe to rerun without manual file deletion.
add_sheet_to_excel(
  workbook_path = output_excel_path,
  sheet_name = output_excel_sheet,
  data = export_dt,
  rownames = FALSE,
  overwrite = TRUE
)

add_sheet_to_excel(
  workbook_path = output_excel_path,
  sheet_name = "File_Inventory",
  data = file_inventory_dt,
  rownames = FALSE,
  overwrite = FALSE
)
```

## Document and Export Column Dictionary

Build a trailing `Column_Dictionary` sheet that documents every column
of every sheet already written to `output_excel_path`, using the
project’s shared
[build_column_dictionary_function.R](../functions/build_column_dictionary_function.R)
helper, then append it to the workbook as the final sheet.

------------------------------------------------------------------------

# Output File Summary

The tree below lists every file this notebook has written to its own
output folder,
[`results/1_data_integrity_check/`](../../results/1_data_integrity_check/),
as a clickable, portable link (relative to this notebook’s own location)
with a short description – built live from what is actually on disk at
knit time via the project’s shared
[`render_output_tree_function.R`](../functions/render_output_tree_function.R)
helper, so it always matches this run’s real output rather than a
hand-maintained list.

------------------------------------------------------------------------

# Recommended Next Step

This notebook’s integrity check is optional but recommended before
continuing. Proceed to [Step 2 — FastQC Quality
Reports](2_fastqc_quality_reports.md), which generates FastQC and
MultiQC quality reports for the same FASTQ files in
[data/fastq/](../../data/fastq/); Step 2 re-derives its own file
discovery and sample identification directly from that folder, so it can
be run whether or not this notebook was completed first.

------------------------------------------------------------------------

# Session Information

Record the R environment for reproducibility.

------------------------------------------------------------------------

# References

## Methods

- Cock PJA, Fields CJ, Goto N, Heuer ML, Rice PM (2010). The Sanger
  FASTQ file format for sequences with quality scores, and the
  Solexa/Illumina FASTQ variants. *Nucleic Acids Research*,
  38(6):1767-1771. <https://doi.org/10.1093/nar/gkp1137>
- Morgan M, Anders S, Lawrence M, Aboyoun P, Pagès H, Gentleman R
  (2009). ShortRead: a bioconductor package for input, quality
  assessment and exploration of high-throughput sequence data.
  *Bioinformatics*, 25(19):2607-2608. (provides `countFastq()`, used for
  complete-file parsing and record counts in this notebook)
  <https://doi.org/10.1093/bioinformatics/btp450>
- Rivest R (1992). The MD5 Message-Digest Algorithm. RFC 1321, Internet
  Engineering Task Force. <https://www.rfc-editor.org/rfc/rfc1321>
- Illumina, Inc. FASTQ Files Explained — background on the Illumina
  FASTQ naming convention (e.g. `_L001_R1_001`) this notebook’s
  file-pairing logic is built around.
  <https://support.illumina.com/bulletins/2016/04/fastq-files-explained.html>

## Related

- [Step 2 — FastQC Quality Reports](2_fastqc_quality_reports.md) — this
  notebook’s recommended next step.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### No FASTQ Files Found

**Error**: `stop("No FASTQ files found in: ...")`

**Solutions**:

- Verify the `fastq_input_folder` path is correct
- Check that files have `.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz`
  extensions
- Ensure files are in the directory root (not in subdirectories, as
  `recurse = FALSE`)

### Mismatched Read Counts

**Warning**: `"Found X paired sample(s) with mismatched read counts"`

**Possible causes**:

- Incomplete file transfer or download
- Files were processed differently (e.g., one was trimmed)
- Incorrectly paired files (wrong samples matched together)

**Actions**:

- Re-download or re-transfer the affected files
- Verify file integrity using [MD5
  checksums](https://www.rfc-editor.org/rfc/rfc1321)
- Check file naming conventions

### Ambiguous File Classification

**Warning**: `"Found X file(s) that match both or neither pattern"`

**Possible causes**:

- Non-standard naming convention
- Single-end reads (not paired)
- Filenames contain both configured direction tokens, or neither token

**Actions**:

- Rename files to follow expected convention
- Modify `forward_token` and `reverse_token` to match your naming scheme
- Exclude single-end reads from this workflow

### Invalid FASTQ Files

If `Valid (Forward)` or `Valid (Reverse)` is `FALSE`:

**Possible causes**:

- Corrupted file (incomplete download)
- Empty file
- Wrong file format (not FASTQ)

**Actions**:

- Re-download the file
- Inspect the reported parsing error and the first FASTQ records
- Ensure gzip compression is not corrupted
