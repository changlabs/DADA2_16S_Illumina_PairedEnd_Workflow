Step 7: 16S rRNA Gene Copy Number Correction
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
  - [Adjust Run Parameters](#adjust-run-parameters)
  - [Define Path Parameters](#define-paths)
- [Prepare Input for PICRUSt2](#prepare-input-for-picrust2)
  - [Convert ASV Sequences to FASTA](#convert-to-fasta)
- [Verify PICRUSt2 Installation](#verify-picrust2-installation)
- [Run Copy Number Prediction](#run-copy-number-prediction)
  - [Define Helper Function](#func-run-verify)
  - [Phylogenetic Placement](#run-place-seqs)
  - [Predict 16S rRNA Copy Number](#run-hsp)
- [Data Import](#data-import)
  - [Import Predicted Copy Numbers](#import-copy-numbers)
  - [Import Step 5 ASV Count Table](#import-raw-abundance)
- [Apply Copy Number Correction](#apply-copy-number-correction)
  - [Identify ASVs Excluded From Correction](#identify-excluded)
  - [Compute the Corrected Abundance Table](#compute-correction)
  - [Export the Corrected Abundance Table (CSV)](#export-csv)
- [Results Summary: Before vs. After](#results-summary-before-vs-after)
  - [Per-Sample Comparison](#per-sample-comparison)
  - [Per-ASV Comparison](#per-asv-comparison)
  - [Prediction Completeness](#excluded-asvs-table)
- [Export to Excel](#export-to-excel)
  - [Write Data and Comparison Sheets](#export-comparison-excel)
  - [Document and Export Column Dictionary](#column-dictionary)
  - [Preserve Raw Prediction File](#cleanup-raw-file)
- [Output File Summary](#output-file-summary)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [PICRUSt2](#picrust2)
  - [Copy Number Correction](#copy-number-correction)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [conda / `picrust2` Environment
      Issues](#conda--picrust2-environment-issues)
    - [`hsp.py` Fails With an Rcpp / Architecture-Mismatch
      Error](#hsppy-fails-with-an-rcpp--architecture-mismatch-error)
    - [`hsp.py`’s Log File Is Empty](#hsppys-log-file-is-empty)
    - [Sequence Placement Failures / Very Few ASVs
      Retained](#sequence-placement-failures--very-few-asvs-retained)

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
&#10;/* Style for parameter display */
.param-box {
  font-family: monospace;
  background-color: #f0f0f0;
  padding: 2px 6px;
  border-radius: 3px;
  font-size: 14px;
}
</style>

# Introduction

## Purpose

This notebook is **Step 7** of the 16S rRNA sequencing data processing
pipeline, and is **entirely optional**. It performs 16S rRNA gene
copy-number correction on the ASV table produced in [Step
5](5_dada2_pipeline.md), using
[PICRUSt2](https://github.com/picrust/picrust2/wiki)’s
phylogenetic-placement-based hidden-state prediction — the same method
PICRUSt2 uses internally to normalize abundances before predicting
community function, made available here as a standalone step for users
who only need copy-number-corrected ASV abundances.

<div class="alert alert-info">

**Why you might want this step**: [Step 8 (Microbial Load
Correction)](8_microbial_load_correction.md), which converts relative
ASV abundances into cell-count-scaled Quantitative Microbiome Profiles,
specifically requires a **copy-number-corrected** abundance table as its
input (Vandeputte et al. 2017). If you have flow-cytometry cell-count
data for your samples and plan to run [Step
8](8_microbial_load_correction.md), run this step first. [Step 9
(Phyloseq Object)](9_phyloseq_object.md) also picks up this step’s
output automatically if present, building an extra phyloseq object from
copy-number-corrected counts alongside its raw-counts object. If you do
not need copy-number correction for any of that, skip this notebook
entirely — it does not block [Step 6](6_phylogenetic_tree.md)
(Phylogenetic Tree) or [Step 9 (Phyloseq Object)](9_phyloseq_object.md),
both of which run directly off [Step 5’s](5_dada2_pipeline.md) output.

</div>

The 16S rRNA gene is present in a variable number of copies per genome
across bacterial and archaeal taxa — anywhere from a single copy to more
than fifteen, depending on the lineage. Because amplicon sequencing
counts 16S rRNA gene *molecules*, not *cells* or *genomes*, raw ASV read
counts are systematically biased toward high-copy-number taxa: two taxa
present at identical cell abundance in a sample can yield very different
read counts purely because of this copy-number difference. Copy-number
correction estimates each ASV’s most likely genomic 16S copy number —
via hidden-state prediction from its position in a reference phylogeny
built from genomes with known, sequenced copy numbers — and divides its
raw read count by that estimate, producing an abundance table that more
closely approximates relative cell (genome) abundance rather than
relative 16S-molecule abundance ([Kembel et
al. 2012](https://doi.org/10.1371/journal.pcbi.1002743); [Langille et
al. 2013](https://doi.org/10.1038/nbt.2676)).

## Prerequisites

Before running this notebook, ensure that:

1.  **Step 5 (DADA2 Pipeline)** has been completed successfully,
    producing
    [asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv)
    and
    [asv_sequences.csv](../../results/5_dada2_pipeline/asv_sequences.csv)
    in [results/5_dada2_pipeline/](../../results/5_dada2_pipeline/).
2.  **PICRUSt2 is installed** in a dedicated conda environment — run
    [install_picrust2.sh](../../setup/install_picrust2.sh) from the
    repository root once per machine. This notebook expects the
    environment to be named `picrust2` (the default created by that
    script).
3.  Required R packages are installed — run
    [setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
    once per R environment if you have not already done so.

## What This Notebook Does

1.  **Converts Step 5’s ASV sequences to FASTA**:
    [asv_sequences.csv](../../results/5_dada2_pipeline/asv_sequences.csv)
    is written out as a FASTA file for PICRUSt2’s placement tool.
2.  **Verifies the PICRUSt2 installation**: confirms the conda
    environment and required tools exist.
3.  **Phylogenetically places** the study ASVs into PICRUSt2’s reference
    tree (`place_seqs.py`).
4.  **Predicts each ASV’s 16S rRNA gene copy number** (and NSTI, as a
    free byproduct) via hidden-state prediction (`hsp.py -i 16S -n`).
5.  **Applies copy-number correction**: divides each ASV’s raw
    per-sample counts by its predicted copy number.
6.  **Exports the corrected table** as a standalone `.csv` file, in the
    same sample-rows / ASV-columns layout as Step 5’s
    `asv_count_table.csv`, so it can be dropped into any downstream
    analysis that expects that format — including [Step 8 (Microbial
    Load Correction)](8_microbial_load_correction.md) and [Step 9
    (Phyloseq Object)](9_phyloseq_object.md).
7.  **Builds a before/after Excel comparison workbook**: the raw
    predicted copy numbers from `hsp.py`, per-sample and per-ASV
    read-count and relative-abundance shifts caused by the correction,
    and a trailing `Column_Dictionary` sheet documenting every sheet.

## Expected Input

- [results/5_dada2_pipeline/asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv)
  — sample x ASV abundance table from [Step 5](5_dada2_pipeline.md)
  (`SampleID` column followed by one column per ASV ID).
- [results/5_dada2_pipeline/asv_sequences.csv](../../results/5_dada2_pipeline/asv_sequences.csv)
  — ASV ID -\> representative sequence mapping from Step 5 (`ASV_ID`,
  `Sequence` columns).

## Expected Output

- [results/7_copy_number_correction/asv_sequences.fasta](../../results/7_copy_number_correction/asv_sequences.fasta)
  — FASTA-format representative sequences, converted from Step 5’s
  `asv_sequences.csv` for PICRUSt2.
- [results/7_copy_number_correction/placed_seqs.tre](../../results/7_copy_number_correction/placed_seqs.tre)
  — the phylogenetic placement tree produced by `place_seqs.py`.
- [results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv](../../results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv)
  — the corrected ASV count table, in the same sample x ASV layout as
  Step 5’s `asv_count_table.csv` (this notebook’s primary deliverable,
  and [Step 8](8_microbial_load_correction.md)’s expected input; [Step
  9](9_phyloseq_object.md) also picks it up automatically if present).
- [results/7_copy_number_correction/copy_number_correction_summary.xlsx](../../results/7_copy_number_correction/copy_number_correction_summary.xlsx)
  — provenance and before/after comparison workbook: `Run_Provenance`,
  `Predicted_16S_Copy_Numbers`, `Per_Sample_Comparison`,
  `Per_ASV_Comparison`, and a trailing `Column_Dictionary` sheet.
- [results/7_copy_number_correction/16S_marker_predicted_and_nsti.tsv](../../results/7_copy_number_correction/16S_marker_predicted_and_nsti.tsv)
  — original machine-readable PICRUSt2 copy-number and NSTI output.
- [results/7_copy_number_correction/logs/](../../results/7_copy_number_correction/logs/)
  — console logs for `place_seqs.py` and `hsp.py`.

<div class="alert alert-info">

The raw PICRUSt2 prediction TSV is retained alongside the workbook so
the original machine-readable prediction, including NSTI values, remains
independently auditable and reusable.

</div>

------------------------------------------------------------------------

# Environment Setup

<div class="alert alert-info">

**Before running this notebook**: install the required R packages by
running
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R),
and install PICRUSt2 itself by running
[install_picrust2.sh](../../setup/install_picrust2.sh) from the
repository root once per machine. Neither needs to be repeated before
every run of this notebook.

</div>

## Load Required Packages

The chunk below loads every R package this notebook depends on, plus
three helper functions shared across the whole pipeline: Excel writing,
column-dictionary documentation, and clickable output links (each
sourced from [R/functions/](../functions/)).

``` r
# NOTE: This chunk assumes the packages below are already installed. If any
# library() call fails with "there is no package called ...", run
# setup/install_R_dependencies.R once to install every package this
# workflow needs, then re-run this notebook.

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# openxlsx: Excel file creation and manipulation
library(openxlsx)

# readr: Fast, consistent delimited-file import/export
library(readr)

# dplyr: Data manipulation grammar (select/rename/filter), used here for
# tidy renaming and validation of PICRUSt2's raw hsp.py output
library(dplyr)

# Biostrings: Biological string manipulation (Bioconductor)
# Used here to write Step 5's ASV sequences out as a FASTA file for PICRUSt2
library(Biostrings)

# jsonlite: Parse conda's machine-readable environment and package metadata
library(jsonlite)

# Source the custom Excel utility function from the project's function library
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source the custom column-dictionary utility function from the project's
# function library. This function documents every column of a sheet being
# exported, so a trailing "Column_Dictionary" sheet can be added to the
# workbook (see the Document and Export Column Dictionary section below).
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source the reproducibility-link utility function from the project's
# function library. This function renders one or more output file paths as
# clickable Markdown links, so this notebook can document exactly where each
# chunk's output landed (see the render_output_links() calls throughout).
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

## Adjust Run Parameters

Set the number of CPU threads PICRUSt2 should use, whether to keep
`place_seqs.py`’s intermediate alignment/placement files afterward, and
confirm the conda environment/executable this notebook will run PICRUSt2
through. `num_threads` is the setting most users need to change to match
their own machine; the conda settings are resolved automatically below
and rarely need manual editing.

``` r
# ------------------------------------------------------------------
# Number of CPU threads to use for PICRUSt2's placement and hidden-state
# prediction steps (their own -p argument). As general guidance, leave 1-2
# cores free on shared/interactive systems.
# ------------------------------------------------------------------
num_threads <- 4  # <-- EDIT THIS to match your available CPU cores

if (length(num_threads) != 1L || !is.numeric(num_threads) ||
    !is.finite(num_threads) || num_threads < 1 || num_threads != round(num_threads)) {
  stop("num_threads must be one positive whole number.")
}
num_threads <- as.integer(num_threads)

# PICRUSt2's conventional maximum NSTI for 16S predictions. Predictions above
# this value are retained but called out explicitly for scientific review.
max_nsti <- 2.0

# Retire prior final deliverables as soon as a new Linux run starts. This
# prevents Steps 8 and 9 from consuming an older correction if conda or either
# PICRUSt2 command fails during the current run.
previous_final_deliverables <- here(
  "results", "7_copy_number_correction",
  c("copy_number_corrected_asv_count_table.csv",
    "copy_number_correction_summary.xlsx")
)
unlink(previous_final_deliverables[file.exists(previous_final_deliverables)])

detected_cores <- suppressWarnings(parallel::detectCores())
available_cores <- if (length(detected_cores) == 1L && is.finite(detected_cores) && detected_cores >= 1) as.integer(detected_cores) else 1L
cat(
  "CPU Configuration:\n",
  "- Threads requested (num_threads):", num_threads, "\n",
  "- Cores detected on this machine:", available_cores, "\n"
)
if (num_threads > available_cores) {
  warning(
    "num_threads (", num_threads, ") exceeds the ", available_cores,
    " cores detected on this machine. PICRUSt2 will still run, but consider ",
    "lowering num_threads to avoid oversubscribing the CPU."
  )
}

# ------------------------------------------------------------------
# Whether to keep or remove place_seqs.py's intermediate alignment/placement
# working files afterward. Defaults to FALSE to save disk space by default;
# set to TRUE if you need to keep the intermediate alignment/placement files
# to troubleshoot a placement issue. Applied manually here (via unlink()
# after a verified-successful run) because place_seqs.py run standalone has
# no built-in --remove_intermediate flag of its own.
# ------------------------------------------------------------------
save_intermediate_files <- FALSE  # <-- set to TRUE to keep intermediate/ for troubleshooting

# ------------------------------------------------------------------
# Name of the conda environment created by install_picrust2.sh
# ------------------------------------------------------------------
conda_env_name <- "picrust2"

# ------------------------------------------------------------------
# Path to the conda executable to use for every command in this notebook.
# ------------------------------------------------------------------
# The Linux installer requires conda on PATH, so use that same explicit
# installation here. This avoids launching an interactive login shell and
# accidentally selecting a different conda installation from shell startup
# output.
#
# Manual override (only needed if automatic detection still picks the wrong
# installation): conda_executable <- "/full/path/to/your/conda"
resolve_conda_executable <- function() {
  path_lookup <- unname(Sys.which("conda"))
  if (nzchar(path_lookup)) path_lookup else NA_character_
}

conda_executable <- resolve_conda_executable()

if (is.na(conda_executable)) {
  stop(
    "Could not automatically determine which conda executable to use.\n",
    "Set it manually in this Configuration section, e.g.:\n",
    "  conda_executable <- \"/full/path/to/your/conda\""
  )
}

cat(
  "Conda Configuration:\n",
  "- Resolved conda executable:", conda_executable, "\n",
  "- (If this is not the installation containing your 'picrust2' environment,",
  " override conda_executable manually above.)\n"
)
```

## Define Path Parameters

This chunk resolves every input and output path relative to the [project
root](../../) via `here()`, creates this step’s own output folder (and
its [logs/](../../results/7_copy_number_correction/logs/) subfolder) if
they do not already exist, and fails fast with an informative error if
either of Step 5’s required input files cannot be found – so a missing
prerequisite is caught immediately, before any computation or PICRUSt2
call begins.

``` r
# Base results folder for all pipeline outputs
results_folder <- here("results")

# Folder containing Step 5's outputs (this step's inputs)
step5_output_folder <- here(results_folder, "5_dada2_pipeline")

# Specific output folder for this step (Step 7: Copy Number Correction)
# Numbered prefix maintains pipeline organization
output_folder <- here(results_folder, "7_copy_number_correction")

# Working directory for place_seqs.py's intermediate alignment/placement files
intermediate_folder <- here(output_folder, "intermediate")

# Folder for console logs from place_seqs.py and hsp.py
log_folder <- here(output_folder, "logs")
checkpoints_folder <- here(output_folder, "checkpoints")

dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(log_folder, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoints_folder, recursive = TRUE, showWarnings = FALSE)

# Path to Step 5's ASV count table (sample x ASV) and sequence mapping
input_asv_count_table_path <- here(step5_output_folder, "asv_count_table.csv")
input_asv_sequences_path   <- here(step5_output_folder, "asv_sequences.csv")

# FASTA file generated from asv_sequences.csv for place_seqs.py's -s argument
asv_fasta_path <- here(output_folder, "asv_sequences.fasta")

# place_seqs.py output: the study sequences placed onto PICRUSt2's reference tree
placed_tree_path <- here(output_folder, "placed_seqs.tre")

# hsp.py output: predicted 16S copy number and NSTI for every placed ASV.
# Retained as the original machine-readable prediction provenance.
marker_prediction_path <- here(output_folder, "16S_marker_predicted_and_nsti.tsv")

# This notebook's primary deliverable: the corrected ASV count table, kept in
# the same sample x ASV layout as Step 5's asv_count_table.csv
output_csv_path <- here(output_folder, "copy_number_corrected_asv_count_table.csv")

# Before/after comparison workbook
output_excel_path <- here(output_folder, "copy_number_correction_summary.xlsx")

# Log file paths
place_seqs_log_path <- here(log_folder, "place_seqs.log")
hsp_log_path         <- here(log_folder, "hsp_16S.log")

# Confirm the Step 5 inputs actually exist before going any further
if (!file.exists(input_asv_count_table_path)) {
  stop(
    "Could not find Step 5's ASV count table: ", input_asv_count_table_path, "\n",
    "Run Step 5 (5_dada2_pipeline.md) first."
  )
}
if (!file.exists(input_asv_sequences_path)) {
  stop(
    "Could not find Step 5's ASV sequences file: ", input_asv_sequences_path, "\n",
    "Run Step 5 (5_dada2_pipeline.md) first."
  )
}

cat(
  "Path Configuration:\n",
  "- Input ASV count table:", input_asv_count_table_path, "\n",
  "- Input ASV sequences:", input_asv_sequences_path, "\n",
  "- Placement tree output:", placed_tree_path, "\n",
  "- Log folder:", log_folder, "\n"
)
```

------------------------------------------------------------------------

# Prepare Input for PICRUSt2

## Convert ASV Sequences to FASTA

`place_seqs.py` requires a FASTA file of representative ASV sequences;
Step 5 exports them as a two-column CSV instead (`ASV_ID`, `Sequence`),
so this section converts one to the other.

``` r
if (file.exists(asv_fasta_path)) unlink(asv_fasta_path)

asv_sequences_raw <- read_csv(
  input_asv_sequences_path, show_col_types = FALSE, name_repair = "minimal"
)

required_sequence_columns <- c("ASV_ID", "Sequence")
missing_sequence_columns <- setdiff(required_sequence_columns, colnames(asv_sequences_raw))
if (length(missing_sequence_columns) > 0) {
  stop(
    "Step 5's asv_sequences.csv is missing expected column(s): ",
    paste(missing_sequence_columns, collapse = ", "),
    ". Found columns: ", paste(colnames(asv_sequences_raw), collapse = ", "), "."
  )
}
if (nrow(asv_sequences_raw) == 0L || anyNA(asv_sequences_raw$ASV_ID) ||
    any(asv_sequences_raw$ASV_ID == "") || anyDuplicated(asv_sequences_raw$ASV_ID)) {
  stop("Step 5's asv_sequences.csv must contain at least one non-missing, unique ASV_ID.")
}
if (any(grepl("[[:space:]]", asv_sequences_raw$ASV_ID))) {
  stop("ASV_ID values must not contain whitespace because they are used as FASTA identifiers.")
}
if (anyNA(asv_sequences_raw$Sequence) || any(asv_sequences_raw$Sequence == "") ||
    any(!grepl("^[ACGTRYSWKMBDHVN]+$", toupper(asv_sequences_raw$Sequence)))) {
  stop("Step 5's asv_sequences.csv contains missing, empty, or non-IUPAC DNA sequences.")
}

# Validate the paired Step 5 count table before running either external tool.
# Minimal name repair preserves duplicate headers so they can be rejected
# rather than silently renamed by readr.
asv_count_table_raw <- read_csv(
  input_asv_count_table_path, show_col_types = FALSE, name_repair = "minimal"
)
if (anyDuplicated(names(asv_count_table_raw))) {
  stop("Step 5's ASV count table contains duplicated column names: ",
       paste(unique(names(asv_count_table_raw)[duplicated(names(asv_count_table_raw))]),
             collapse = ", "))
}
if (!"SampleID" %in% names(asv_count_table_raw)) {
  stop("Step 5's asv_count_table.csv is missing the expected SampleID column.")
}
if (nrow(asv_count_table_raw) == 0L || anyNA(asv_count_table_raw$SampleID) ||
    any(!nzchar(trimws(asv_count_table_raw$SampleID))) ||
    anyDuplicated(asv_count_table_raw$SampleID)) {
  stop("Step 5's ASV count table must contain non-missing, unique SampleID values.")
}
count_asv_ids <- setdiff(names(asv_count_table_raw), "SampleID")
if (!setequal(count_asv_ids, asv_sequences_raw$ASV_ID)) {
  stop("Step 5's ASV sequence and count files contain different ASV_ID sets. Rerun Step 5.")
}
asv_count_table_raw <- asv_count_table_raw[, c("SampleID", asv_sequences_raw$ASV_ID)]
asv_count_matrix <- as.matrix(asv_count_table_raw[, -1, drop = FALSE])
rownames(asv_count_matrix) <- asv_count_table_raw$SampleID
if (!is.numeric(asv_count_matrix) || any(!is.finite(asv_count_matrix)) ||
    any(asv_count_matrix < 0) || any(asv_count_matrix != round(asv_count_matrix))) {
  stop("Step 5's ASV count table must contain finite, non-negative whole-number counts.")
}
if (sum(asv_count_matrix) <= 0) {
  stop("Step 5's ASV count table contains no reads; copy-number correction is undefined.")
}

# Build a DNAStringSet keyed by ASV_ID and write it out as FASTA
asv_sequence_set <- Biostrings::DNAStringSet(asv_sequences_raw$Sequence)
names(asv_sequence_set) <- asv_sequences_raw$ASV_ID
Biostrings::writeXStringSet(asv_sequence_set, filepath = asv_fasta_path, format = "fasta")

cat(
  "Wrote", length(asv_sequence_set), "ASV sequences to FASTA:\n", asv_fasta_path, "\n"
)
```

------------------------------------------------------------------------

# Verify PICRUSt2 Installation

Confirm that conda is available and that the `picrust2` environment
(created by `install_picrust2.sh`) exists, before attempting to run
anything. The environment is resolved by its full directory path
(`-p <prefix>`) rather than by name (`-n <name>`), since
`conda run -n <name>` can silently resolve to the wrong installation on
machines with more than one conda base install.

``` r
conda_version_output <- tryCatch(
  system2(conda_executable, args = "--version", stdout = TRUE, stderr = TRUE),
  error = function(e) NULL
)

if (is.null(conda_version_output) || !is.null(attr(conda_version_output, "status"))) {
  stop(
    "Could not run '", conda_executable, "'. PICRUSt2 is installed via conda — ",
    "run install_picrust2.sh from the repository root first."
  )
}

# Ask conda directly for a machine-readable environment list (--json), rather
# than pattern-matching its human-readable table.
conda_env_list_json <- system2(
  conda_executable, args = c("env", "list", "--json"),
  stdout = TRUE, stderr = FALSE
)
conda_env_info <- tryCatch(
  jsonlite::fromJSON(paste(conda_env_list_json, collapse = "\n")),
  error = function(e) stop("Could not parse 'conda env list --json': ", conditionMessage(e))
)
conda_env_paths <- as.character(conda_env_info$envs)
conda_env_names <- basename(conda_env_paths)

if (!(conda_env_name %in% conda_env_names)) {
  stop(
    "No conda environment named '", conda_env_name, "' was found by '", conda_executable, "'.\n",
    "Environments visible to this conda installation: ",
    if (length(conda_env_names) > 0) paste(conda_env_names, collapse = ", ") else "(none found)", "\n\n",
    "Run install_picrust2.sh from the repository root to create the environment, ",
    "or see the Troubleshooting appendix if you believe it already exists."
  )
}

# Resolve the environment's full directory (its "prefix"). Every `conda run`
# call below uses `-p <prefix>` rather than `-n <name>`, for the same reason
# noted above.
matching_env_prefixes <- conda_env_paths[basename(conda_env_paths) == conda_env_name]
if (length(matching_env_prefixes) != 1L) {
  stop("Expected exactly one conda environment named '", conda_env_name,
       "', but found ", length(matching_env_prefixes), ": ",
       paste(matching_env_prefixes, collapse = ", "))
}
conda_env_prefix <- matching_env_prefixes[[1]]

# Confirm place_seqs.py and hsp.py are actually available inside the
# environment before attempting to run either.
check_tool_available <- function(tool_name) {
  result <- system2(
    conda_executable,
    args = c("run", "-p", shQuote(conda_env_prefix), "--no-capture-output", tool_name, "--help"),
    stdout = TRUE, stderr = TRUE
  )
  status <- attr(result, "status")
  list(available = is.null(status) || status == 0, output = result)
}

place_seqs_availability <- check_tool_available("place_seqs.py")
hsp_availability         <- check_tool_available("hsp.py")

# Record the installed PICRUSt2 package build, which also identifies the
# bundled reference-data release used by place_seqs.py and hsp.py.
picrust_package_json <- system2(
  conda_executable,
  args = c("list", "-p", shQuote(conda_env_prefix), "picrust2", "--json"),
  stdout = TRUE, stderr = FALSE
)
picrust_package_info <- tryCatch(
  jsonlite::fromJSON(paste(picrust_package_json, collapse = "\n")),
  error = function(e) stop("Could not parse PICRUSt2 package metadata: ", conditionMessage(e))
)
if (nrow(picrust_package_info) != 1L) {
  stop("Expected exactly one installed picrust2 package record in: ", conda_env_prefix)
}
picrust2_version <- as.character(picrust_package_info$version[[1]])
picrust2_build <- as.character(picrust_package_info$build_string[[1]])
```

------------------------------------------------------------------------

# Run Copy Number Prediction

## Define Helper Function

Defines a single helper used for both commands below: run a
conda-wrapped command with its output logged to file, then confirm it
actually succeeded by checking for its expected output file rather than
trusting the reported exit code alone (`conda run` does not always
reliably propagate a real failure’s exit code).

``` r
# Run a Conda-Wrapped Command and Verify It Produced an Expected Output File
#
# Arguments:
#   step_label            - Character scalar label used in messages (e.g. "Phylogenetic
#                            placement").
#   command_args          - Character vector of arguments passed to the conda executable.
#   log_path              - Character scalar path to write this command's combined
#                            stdout/stderr to.
#   expected_output_path  - Character scalar path to the file this command should produce
#                            if it actually succeeded.
# Returns:
#   A list with `exit_code`, `duration_minutes`, and `succeeded` (logical).
run_and_verify_step <- function(step_label, command_args, log_path, expected_output_path) {

  # A current command must create the artefact used to judge its success.
  # Remove only this step's previous output and log before invoking it.
  if (file.exists(expected_output_path)) unlink(expected_output_path)
  if (file.exists(log_path)) unlink(log_path)

  cat(
    "Running", step_label, "...\n",
    "Command: conda", paste(command_args, collapse = " "), "\n\n"
  )

  # ---------------------------------------------------------------------------
  # hsp.py internally shells out to a bundled Rscript (castor_nsti.R) to do
  # the actual copy-number / NSTI calculation. On machines that also have a
  # system-wide RStudio/R installation, R-specific environment variables
  # inherited from the parent RStudio R session (R_HOME, R_LIBS, ...) can
  # leak into this child process. When that happens, the nested Rscript
  # resolves R packages (in particular Rcpp) against the system-wide R
  # library instead of the conda environment's own, which fails with an
  # incompatible-library error when the system and conda R installations use
  # different library trees -- see the Troubleshooting appendix entry below.
  # These variables are stripped from the child process's environment via
  # `env -u` so any nested Rscript call is forced to resolve its own R
  # installation from scratch rather than inheriting RStudio's.
  # ---------------------------------------------------------------------------
  r_env_vars_to_strip <- c(
    "R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE",
    "R_ENVIRON", "R_ENVIRON_USER", "R_PROFILE", "R_PROFILE_USER"
  )
  env_unset_flags <- as.vector(rbind("-u", r_env_vars_to_strip))

  start_time <- Sys.time()
  stderr_log_path <- paste0(log_path, ".stderr")
  if (file.exists(stderr_log_path)) unlink(stderr_log_path)
  exit_code <- system2(
    "/usr/bin/env",
    args   = c(env_unset_flags, shQuote(conda_executable), command_args),
    stdout = log_path,
    stderr = stderr_log_path
  )
  if (file.exists(stderr_log_path)) {
    stderr_lines <- readLines(stderr_log_path, warn = FALSE)
    if (length(stderr_lines)) write(stderr_lines, file = log_path, append = TRUE)
    unlink(stderr_log_path)
  }
  end_time <- Sys.time()
  duration_minutes <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)

  output_looks_complete <- file.exists(expected_output_path) &&
    isTRUE(file.info(expected_output_path)$size > 0)
  succeeded <- isTRUE(exit_code == 0) && output_looks_complete

  cat(
    "\n", step_label, "finished in", duration_minutes, "minutes ",
    "(reported exit code:", exit_code, "). Verifying output...\n"
  )

  if (succeeded) {
    cat(step_label, "verified — found expected output file:\n ", expected_output_path, "\n")
  } else {
    warning(
      step_label, " did not produce the expected output file (", expected_output_path,
      "). Reported exit code: ", exit_code, ". See the log excerpt printed below for the underlying error."
    )
    if (file.exists(log_path)) {
      log_lines  <- readLines(log_path, warn = FALSE)
      tail_lines <- utils::tail(log_lines, 40)
      cat(
        "\n--- Last", length(tail_lines), "line(s) of", basename(log_path), "---\n",
        paste(tail_lines, collapse = "\n"),
        "\n--- End of log excerpt (full log at ", log_path, ") ---\n\n"
      )
    } else {
      cat("\nNo log file found at:", log_path, "\n")
    }
  }

  list(exit_code = exit_code, duration_minutes = duration_minutes, succeeded = succeeded)
}
```

## Phylogenetic Placement

Places every ASV into PICRUSt2’s reference tree using
[EPA-NG](http://github.com/pierrebarbera/epa-ng).

``` r
placement_checkpoint <- if (!save_intermediate_files) {
  load_stage_checkpoint(
    "1_placement", c("placed_tree_lines", "place_log_lines", "placement_result")
  )
} else NULL
if (!is.null(placement_checkpoint)) {
  writeLines(placement_checkpoint$placed_tree_lines, placed_tree_path)
  writeLines(placement_checkpoint$place_log_lines, place_seqs_log_path)
  placement_result <- placement_checkpoint$placement_result
} else {
# Never let files from an interrupted placement become inputs to this run.
if (dir.exists(intermediate_folder)) unlink(intermediate_folder, recursive = TRUE)

# NOTE: two unrelated "-p" flags appear below:
#   - the first "-p" belongs to `conda run` and takes the environment's
#     directory (conda_env_prefix);
#   - the second "-p" belongs to place_seqs.py itself and takes the number
#     of CPU threads (num_threads).
place_seqs_args <- c(
  "run", "-p", shQuote(conda_env_prefix), "--no-capture-output",
  "place_seqs.py",
  "-s", shQuote(asv_fasta_path),
  "-o", shQuote(placed_tree_path),
  "-p", as.character(num_threads),
  "--intermediate", shQuote(intermediate_folder),
  "--verbose"
)

placement_result <- run_and_verify_step(
  step_label            = "Phylogenetic placement (place_seqs.py)",
  command_args          = place_seqs_args,
  log_path              = place_seqs_log_path,
  expected_output_path  = placed_tree_path
)

if (!placement_result$succeeded) {
  stop(
    "Phylogenetic placement failed — see the log excerpt above (full log at ",
    place_seqs_log_path, ") for the underlying error before proceeding."
  )
}

checkpoint_placement_path <- save_stage_checkpoint(
  "1_placement",
  list(
    placed_tree_lines = readLines(placed_tree_path, warn = FALSE),
    place_log_lines = readLines(place_seqs_log_path, warn = FALSE),
    placement_result = placement_result
  )
)
}
```

## Predict 16S rRNA Copy Number

Runs hidden-state prediction for the `16S` trait against the placement
tree produced above, using the
[castor](https://cran.r-project.org/package=castor) R package. The `-n`
flag additionally computes each ASV’s Nearest Sequenced Taxon Index
(NSTI) as a free byproduct.

``` r
hsp_checkpoint <- load_stage_checkpoint(
  "2_hsp", c("marker_prediction_lines", "hsp_log_lines", "hsp_result")
)
if (!is.null(hsp_checkpoint)) {
  writeLines(hsp_checkpoint$marker_prediction_lines, marker_prediction_path)
  writeLines(hsp_checkpoint$hsp_log_lines, hsp_log_path)
  hsp_result <- hsp_checkpoint$hsp_result
} else {
hsp_args <- c(
  "run", "-p", shQuote(conda_env_prefix), "--no-capture-output",
  "hsp.py",
  "-i", "16S",
  "-t", shQuote(placed_tree_path),
  "-o", shQuote(marker_prediction_path),
  "-p", as.character(num_threads),
  "-n"
)

hsp_result <- run_and_verify_step(
  step_label            = "16S copy number prediction (hsp.py)",
  command_args          = hsp_args,
  log_path              = hsp_log_path,
  expected_output_path  = marker_prediction_path
)
if (!hsp_result$succeeded) {
  stop(
    "16S copy number prediction failed — see the log excerpt above (full log at ",
    hsp_log_path, ") for the underlying error before proceeding."
  )
}
checkpoint_hsp_path <- save_stage_checkpoint(
  "2_hsp",
  list(
    marker_prediction_lines = readLines(marker_prediction_path, warn = FALSE),
    hsp_log_lines = readLines(hsp_log_path, warn = FALSE),
    hsp_result = hsp_result
  )
)
}
```

------------------------------------------------------------------------

# Data Import

## Import Predicted Copy Numbers

This chunk reads `hsp.py`’s raw copy-number prediction output and
validates it before use: every predicted copy number must be a finite,
positive number, since a genome cannot carry zero or a negative number
of 16S rRNA gene copies.

``` r
predicted_copy_number_raw <- read_tsv(marker_prediction_path, show_col_types = FALSE)

required_marker_columns <- c("sequence", "16S_rRNA_Count")
missing_marker_columns <- setdiff(required_marker_columns, colnames(predicted_copy_number_raw))
if (length(missing_marker_columns) > 0) {
  stop(
    "hsp.py's 16S copy number output is missing expected column(s): ",
    paste(missing_marker_columns, collapse = ", "),
    ". Found columns: ", paste(colnames(predicted_copy_number_raw), collapse = ", "), "."
  )
}

predicted_copy_number <- predicted_copy_number_raw %>%
  select(ASV_ID = sequence, Predicted_16S_Copy_Number = `16S_rRNA_Count`)

if (anyNA(predicted_copy_number$ASV_ID) || any(predicted_copy_number$ASV_ID == "") ||
    anyDuplicated(predicted_copy_number$ASV_ID)) {
  stop("hsp.py returned missing, empty, or duplicated ASV identifiers in: ",
       marker_prediction_path)
}

# Defensive check: every genome carries at least one copy of the 16S rRNA
# gene, so a non-positive or missing predicted copy number would make
# correction undefined (division by zero or a negative number).
invalid_copy_numbers <- predicted_copy_number %>%
  filter(!is.finite(Predicted_16S_Copy_Number) | Predicted_16S_Copy_Number <= 0)

if (nrow(invalid_copy_numbers) > 0) {
  stop(
    nrow(invalid_copy_numbers), " ASV(s) have a non-positive or missing predicted 16S copy ",
    "number, which would make copy-number correction undefined: ",
    paste(head(invalid_copy_numbers$ASV_ID, 10), collapse = ", "),
    if (nrow(invalid_copy_numbers) > 10) ", ... (truncated)" else "",
    ". Inspect ", marker_prediction_path, " directly for these ASVs."
  )
}

missing_prediction_ids <- setdiff(colnames(asv_count_matrix), predicted_copy_number$ASV_ID)
unexpected_prediction_ids <- setdiff(predicted_copy_number$ASV_ID, colnames(asv_count_matrix))
if (length(missing_prediction_ids) || length(unexpected_prediction_ids)) {
  stop(
    "PICRUSt2 predictions do not exactly match the current Step 5 ASV set. ",
    "Missing predictions: ",
    if (length(missing_prediction_ids)) paste(missing_prediction_ids, collapse = ", ") else "(none)",
    "; unexpected predictions: ",
    if (length(unexpected_prediction_ids)) paste(unexpected_prediction_ids, collapse = ", ") else "(none)",
    ". No corrected table will be published. Review placement logs and inputs."
  )
}

if (!"metadata_NSTI" %in% names(predicted_copy_number_raw)) {
  stop("hsp.py was run with -n but its output does not contain metadata_NSTI: ",
       marker_prediction_path)
}
if (any(!is.finite(predicted_copy_number_raw$metadata_NSTI)) ||
    any(predicted_copy_number_raw$metadata_NSTI < 0)) {
  stop("hsp.py returned missing, non-finite, or negative NSTI values.")
}
high_nsti_ids <- predicted_copy_number_raw$sequence[
  predicted_copy_number_raw$metadata_NSTI > max_nsti
]
if (length(high_nsti_ids)) {
  warning(length(high_nsti_ids), " ASV(s) exceed max_nsti = ", max_nsti,
          " and have predictions based on distant references: ",
          paste(head(high_nsti_ids, 10), collapse = ", "),
          if (length(high_nsti_ids) > 10) ", ..." else "", call. = FALSE)
}

cat(
  "Imported predicted 16S copy numbers for", nrow(predicted_copy_number), "ASVs.\n",
  "- Copy number range:", round(min(predicted_copy_number$Predicted_16S_Copy_Number), 2), "-",
  round(max(predicted_copy_number$Predicted_16S_Copy_Number), 2), "\n",
  "- Mean copy number:", round(mean(predicted_copy_number$Predicted_16S_Copy_Number), 2), "\n"
  , "- NSTI range:", round(min(predicted_copy_number_raw$metadata_NSTI), 3), "-",
  round(max(predicted_copy_number_raw$metadata_NSTI), 3), "\n",
  "- ASVs above max_nsti (", max_nsti, "):", length(high_nsti_ids), "\n"
)

run_provenance <- data.frame(
  PICRUSt2_Version = picrust2_version,
  PICRUSt2_Build = picrust2_build,
  Conda_Environment = normalizePath(conda_env_prefix, mustWork = TRUE),
  ASV_Sequences_MD5 = unname(tools::md5sum(input_asv_sequences_path)),
  ASV_Count_Table_MD5 = unname(tools::md5sum(input_asv_count_table_path)),
  Maximum_NSTI_Guidance = max_nsti,
  ASVs_Above_Maximum_NSTI = length(high_nsti_ids),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
```

## Import Step 5 ASV Count Table

This chunk reads Step 5’s raw ASV count table and reshapes it into a
numeric matrix (samples x ASVs) with `SampleID` moved into the row names
– the format the correction step below expects.

``` r
# The table was imported and cross-validated with asv_sequences.csv before
# PICRUSt2 was invoked; report the already validated matrix here.
cat(
  "Imported Step 5 ASV count table:\n",
  "- Samples:", nrow(asv_count_matrix), "\n",
  "- ASVs:", ncol(asv_count_matrix), "\n",
  "- Total raw reads:", sum(asv_count_matrix), "\n"
)
```

------------------------------------------------------------------------

# Apply Copy Number Correction

## Identify ASVs Excluded From Correction

Every Step 5 ASV must receive exactly one prediction. If placement
excludes an ASV, the notebook stops in the import section and does not
publish a partial corrected table. This preserves the exact Step 5 ASV
schema and keeps per-sample comparisons attributable solely to
copy-number correction.

``` r
asv_ids_in_counts        <- colnames(asv_count_matrix)
asv_ids_with_prediction  <- predicted_copy_number$ASV_ID

asvs_excluded_from_correction <- setdiff(asv_ids_in_counts, asv_ids_with_prediction)
asvs_used_for_correction <- asv_ids_in_counts
stopifnot(length(asvs_excluded_from_correction) == 0L)
cat("Every Step 5 ASV received exactly one predicted 16S copy number.\n")
```

## Compute the Corrected Abundance Table

This chunk divides each ASV’s raw per-sample counts by that ASV’s
predicted 16S copy number, producing the copy-number-corrected abundance
table – the core computation this notebook exists to perform.

``` r
# Build a copy-number lookup vector indexed by ASV_ID, then subset and
# reorder it to exactly match the COLUMN order of asv_count_matrix (ASVs are
# columns here, unlike PICRUSt2's own ASV-x-sample/BIOM convention) -- this
# guarantees the element-wise division below pairs each ASV with its own
# predicted copy number rather than risking a silent misalignment.
copy_number_lookup <- setNames(predicted_copy_number$Predicted_16S_Copy_Number, predicted_copy_number$ASV_ID)
copy_number_vector <- copy_number_lookup[asvs_used_for_correction]

raw_count_matrix_for_correction <- asv_count_matrix[, asvs_used_for_correction, drop = FALSE]

# Divide every ASV's per-sample counts by that ASV's predicted 16S copy
# number. sweep() with MARGIN = 2 applies copy_number_vector column-wise (one
# value per ASV), dividing every sample row in that column by the same
# number -- the classic copy-number-correction operation (Kembel et al.
# 2012; Langille et al. 2013), and numerically identical to what PICRUSt2's
# own metagenome_pipeline.py computes internally.
corrected_count_matrix <- sweep(raw_count_matrix_for_correction, MARGIN = 2, STATS = copy_number_vector, FUN = "/")

cat(
  "Copy number correction applied:\n",
  "- ASVs corrected:", ncol(corrected_count_matrix), "\n",
  "- ASVs excluded (no placement/prediction):", length(asvs_excluded_from_correction), "\n",
  "- Samples:", nrow(corrected_count_matrix), "\n",
  "- Total raw reads (corrected ASVs only):", round(sum(raw_count_matrix_for_correction), 1), "\n",
  "- Total corrected reads:", round(sum(corrected_count_matrix), 1), "\n"
)
```

## Export the Corrected Abundance Table (CSV)

This chunk writes the corrected ASV table to
[results/7_copy_number_correction/](../../results/7_copy_number_correction/)
as CSV, keeping the same sample x ASV layout as Step 5’s own output so
it can be used as a drop-in replacement anywhere downstream; the
exported file is linked immediately below via `render_output_links()`.

``` r
# Sample-rows / ASV-columns orientation, matching Step 5's asv_count_table.csv
# exactly, so this file can be read as a drop-in replacement anywhere Step 5's
# output is used downstream (including Steps 8 and 9).
corrected_asv_count_table_export <- data.frame(
  SampleID = rownames(corrected_count_matrix),
  corrected_count_matrix,
  check.names = FALSE,
  row.names = NULL
)

write_csv(corrected_asv_count_table_export, output_csv_path)

cat(
  "Copy-number-corrected ASV count table written to:\n", output_csv_path, "\n",
  "(", nrow(corrected_asv_count_table_export), "samples x", ncol(corrected_count_matrix), "ASVs )\n"
)
```

------------------------------------------------------------------------

# Results Summary: Before vs. After

## Per-Sample Comparison

Total read count per sample, before (raw) and after
(copy-number-corrected).

``` r
common_samples <- rownames(asv_count_matrix)

raw_total_per_sample       <- rowSums(asv_count_matrix)[common_samples]
corrected_total_per_sample <- rowSums(corrected_count_matrix)[common_samples]

# Guard against samples with zero total raw reads -- dividing by zero here
# would silently produce Inf/NaN. Reported as NA instead.
zero_count_samples <- names(raw_total_per_sample)[raw_total_per_sample == 0]
if (length(zero_count_samples) > 0) {
  warning(
    length(zero_count_samples), " sample(s) have zero total raw reads (",
    paste(zero_count_samples, collapse = ", "),
    "), so their Ratio (After / Before) and Percent Change are undefined and reported as NA."
  )
}

ratio_after_before <- ifelse(
  raw_total_per_sample == 0, NA_real_, corrected_total_per_sample / raw_total_per_sample
)
percent_change <- ifelse(is.na(ratio_after_before), NA_real_, round(100 * (ratio_after_before - 1), 2))

per_sample_comparison <- data.frame(
  `Sample ID`                     = common_samples,
  `Total Raw Reads (Before)`      = unname(raw_total_per_sample),
  `Total Corrected Reads (After)` = round(unname(corrected_total_per_sample), 2),
  `Ratio (After / Before)`        = round(unname(ratio_after_before), 4),
  `Percent Change`                = unname(percent_change),
  check.names = FALSE
)

# DT is namespace-qualified (not library()-loaded) here since this is its
# only use in this notebook; the fixed-height scroll box matches the
# interactive DT::datatable() tables used elsewhere in this pipeline.
head(per_sample_comparison, 10) %>%
  DT::datatable(options = list(pageLength = 10, scrollX = TRUE,
                                scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
                rownames = FALSE,
                caption = "First 10 samples (full table in the Excel export)")
```

## Per-ASV Comparison

Predicted copy number and relative-abundance shift per ASV, sorted by
the magnitude of that shift (largest movers first).

``` r
raw_total_per_asv       <- colSums(raw_count_matrix_for_correction)
corrected_total_per_asv <- colSums(corrected_count_matrix)

total_raw_reads_corrected_asvs <- sum(raw_count_matrix_for_correction)
total_corrected_reads          <- sum(corrected_count_matrix)

per_asv_comparison <- data.frame(
  `ASV_ID`                            = asvs_used_for_correction,
  `Predicted 16S Copy Number`         = round(unname(copy_number_vector[asvs_used_for_correction]), 3),
  `Total Raw Reads (Before)`          = unname(raw_total_per_asv),
  `Total Corrected Reads (After)`     = round(unname(corrected_total_per_asv), 2),
  `Relative Abundance Raw (%)`        = round(100 * unname(raw_total_per_asv) / total_raw_reads_corrected_asvs, 4),
  `Relative Abundance Corrected (%)`  = round(100 * unname(corrected_total_per_asv) / total_corrected_reads, 4),
  check.names = FALSE
)

per_asv_comparison$`Relative Abundance Shift (pp)` <- round(
  per_asv_comparison$`Relative Abundance Corrected (%)` - per_asv_comparison$`Relative Abundance Raw (%)`, 4
)

per_asv_comparison <- per_asv_comparison[order(-abs(per_asv_comparison$`Relative Abundance Shift (pp)`)), ]
rownames(per_asv_comparison) <- NULL

head(per_asv_comparison, 15) %>%
  DT::datatable(options = list(pageLength = 15, scrollX = TRUE,
                                scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
                rownames = FALSE,
                caption = "15 ASVs with the largest relative-abundance shift (full table in the Excel export)")
```

## Prediction Completeness

The pipeline requires complete prediction coverage. Reaching this
section confirms that no ASVs were excluded.

``` r
cat("Prediction coverage: 100% (", length(asvs_used_for_correction), " of ",
    length(asv_ids_in_counts), " ASVs).\n", sep = "")
```

------------------------------------------------------------------------

# Export to Excel

## Write Data and Comparison Sheets

The `Predicted_16S_Copy_Numbers` sheet contains `hsp.py`’s raw output
exactly as PICRUSt2 produced it. The original TSV is also retained as an
independently reusable machine-readable file.

## Document and Export Column Dictionary

Every column written to this workbook is documented in a trailing
`Column_Dictionary` sheet (Sheet / Column / Explanation).

## Preserve Raw Prediction File

The original `hsp.py` TSV is retained as a machine-readable provenance
artefact in addition to its workbook sheet.

``` r
cat("Retained raw PICRUSt2 prediction file:\n", marker_prediction_path, "\n")
render_output_links(marker_prediction_path, labels = "Raw PICRUSt2 16S prediction and NSTI table (TSV)")
```

------------------------------------------------------------------------

# Output File Summary

The tree below lists the current run’s outputs. Final downstream
deliverables are retired before external processing begins, command
outputs are recreated or restored only from signature-validated
checkpoints, and incompatible checkpoints are removed.

------------------------------------------------------------------------

# Recommended Next Step

The copy-number-corrected ASV table exported above enables two optional
downstream paths, and which one to take depends on whether you have
microbial-load data. If you have an independent, flow-cytometry-derived
microbial-load (cell-count) measurement for your samples, proceed to
[Step 8 — Microbial Load Correction](8_microbial_load_correction.md),
which specifically requires this notebook’s output as its input to
compute a Quantitative Microbiome Profile. If you do not have (or do not
need) microbial-load data, you can instead skip directly to [Step 9 —
Phyloseq Object Assembly](9_phyloseq_object.md), which automatically
detects and incorporates this notebook’s output
(`results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv`)
if present, building an additional copy-number-corrected phyloseq object
alongside its raw-counts object – with no error if this notebook was
skipped instead.

------------------------------------------------------------------------

# Session Information

------------------------------------------------------------------------

# References

## PICRUSt2

- [PICRUSt2 GitHub repository](https://github.com/picrust/picrust2)
- [PICRUSt2 wiki (user
  manual)](https://github.com/picrust/picrust2/wiki)
- Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for
  prediction of metagenome functions. *Nat Biotechnol* 38, 685-688.
  <https://doi.org/10.1038/s41587-020-0548-6>

## Copy Number Correction

- Langille MGI, Zaneveld J, Caporaso JG, et al. (2013). Predictive
  functional profiling of microbial communities using 16S rRNA marker
  gene sequences. *Nat Biotechnol* 31, 814-821.
  <https://doi.org/10.1038/nbt.2676>
- Kembel SW, Wu M, Eisen JA, Green JL (2012). Incorporating 16S Gene
  Copy Number Information Improves Estimates of Microbial Diversity and
  Abundance. *PLoS Comput Biol* 8(10):e1002743.
  <https://doi.org/10.1371/journal.pcbi.1002743>
- Louca S, Doebeli M (2018). Efficient comparative phylogenetics on
  large trees. *Bioinformatics* 34(6):1053-1055. (`castor`, hidden-state
  prediction) — <https://cran.r-project.org/package=castor>
- Barbera P, et al. (2019). EPA-ng: Massively Parallel Evolutionary
  Placement of Genetic Sequences. *Systematic Biology* 68(2):365-369. —
  <https://github.com/pierrebarbera/epa-ng>

## Related

- [Step 5 — DADA2 Pipeline](5_dada2_pipeline.md) — this notebook’s
  required input.
- [Step 8 — Microbial Load Correction](8_microbial_load_correction.md) —
  optional next step; requires this notebook’s copy-number-corrected
  output as its own input.
- [Step 9 — Phyloseq Object Assembly](9_phyloseq_object.md) —
  automatically incorporates this notebook’s output if present, whether
  or not Step 8 is also run.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### conda / `picrust2` Environment Issues

**Symptom**:
`Could not automatically determine which conda executable to use`,
`No conda environment named 'picrust2' was found`, or a conda
`EnvironmentLocationNotFound` error.

**Actions**:

- Confirm which conda installation and environment you expect to use by
  running `echo $CONDA_EXE` and `conda env list` directly in a terminal.
- If the environment does not exist yet, run
  [install_picrust2.sh](../../setup/install_picrust2.sh) from the
  repository root to create it.
- If `resolve_conda_executable()` is still resolving the wrong
  installation, override it manually in
  [Configuration](#adjust-run-parameters):
  `conda_executable <- "/full/path/to/your/conda"`.

### `hsp.py` Fails With an Rcpp / Architecture-Mismatch Error

**Symptom**: `place_seqs.py` succeeds (`placed_seqs.tre` is produced),
but the [Predict 16S rRNA Copy Number](#run-hsp) chunk stops with
`16S copy number prediction failed`, and
[results/7_copy_number_correction/logs/hsp_16S.log](../../results/7_copy_number_correction/logs/hsp_16S.log)
contains something like:

    WARNING: ignoring environment value of R_HOME
    Error: package or namespace load failed for 'Rcpp' in dyn.load(file, DLLpath = DLLpath, ...):
     unable to load shared object '.../R.framework/Versions/.../Rcpp/libs/Rcpp.so':
      ... (mach-o file, but is an incompatible architecture (have 'arm64', need 'x86_64'))

**Cause**: `hsp.py` internally runs `castor_nsti.R` inside the PICRUSt2
environment. R-specific variables inherited from another Linux R
installation can make that child process resolve packages from an
incompatible library tree. This notebook strips those variables from the
child command.

**Fix**: [Define Helper Function](#func-run-verify)
(`run_and_verify_step()`) runs every conda-wrapped command through
`/usr/bin/env -u R_HOME -u R_LIBS -u R_LIBS_USER -u R_LIBS_SITE -u R_ENVIRON -u R_ENVIRON_USER`
`-u R_PROFILE -u R_PROFILE_USER <conda executable> ...`, stripping these
variables before `conda run` ever sees them. If this error still occurs
after updating to this version of the notebook, re-run this notebook
from the top so the fix is re-sourced

### `hsp.py`’s Log File Is Empty

**Symptom**:
[results/7_copy_number_correction/logs/hsp_16S.log](../../results/7_copy_number_correction/logs/hsp_16S.log)
is empty (or a few bytes) after a run.

**Cause**: This is expected, not a bug — `hsp.py` has no documented
`--verbose` flag of its own. Confirm success via the ASV count reported
in [Import Predicted Copy Numbers](#import-copy-numbers) instead of via
this log’s length.

### Sequence Placement Failures / Very Few ASVs Retained

**Symptom**: `place_seqs.py` completes, but Step 7 stops because one or
more Step 5 ASVs have no copy-number prediction.

**Possible causes**:

- Input sequences are too short, too degenerate, or not actually 16S
  rRNA sequences.
- Too few PICRUSt2 reference sequences are sufficiently similar to your
  study sequences (see the [PICRUSt2 Key
  Limitations](https://github.com/picrust/picrust2/wiki/Key-Limitations)
  reference) — amplicons from unusual or under-studied environments are
  more likely to see high exclusion rates.

**Actions**:

- Check `place_seqs_log_path` (linked in [Output File
  Summary](#output-file-summary)) for the specific EPA-NG/gappa warning
  about excluded sequences.
- `save_intermediate_files` defaults to `FALSE` in [Adjust Run
  Parameters](#adjust-run-parameters), so the alignment/placement
  working files will not be available for further diagnosis unless you
  set it to `TRUE` before re-running.
