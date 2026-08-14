Step 8: Microbial Load Correction (Quantitative Microbiome Profiling)
================

- [Introduction](#introduction)
  - [Purpose](#purpose)
  - [Prerequisites](#prerequisites)
  - [Cell Count File Format](#cell-count-format)
  - [What This Notebook Does](#what-this-notebook-does)
  - [Expected Input](#expected-input)
  - [Expected Output](#expected-output)
- [Environment Setup](#environment-setup)
  - [Load Required Packages](#load-packages)
- [Configuration](#configuration)
  - [Adjust Run Parameters](#adjust-run-parameters)
  - [Define Path Parameters](#define-paths)
- [Data Import](#data-import)
  - [Import Copy-Number-Corrected ASV Table](#import-corrected-table)
  - [Import Cell Count File](#import-cell-counts)
  - [Validate Matching Sample Sets](#validate-samples)
- [Compute Quantitative Microbiome
  Profile](#compute-quantitative-microbiome-profile)
  - [Rarefy to a Common Sampling Depth Per Cell](#rarefy-per-cell)
  - [Rescale by Cell Count](#rescale-by-cell-count)
  - [Export the Microbial-Load-Corrected Table (CSV)](#export-csv)
- [Results Summary](#results-summary)
  - [Sampling Depth Diagnostics](#sampling-depth-summary)
- [Export to Excel](#export-to-excel)
  - [Write Summary Sheet](#export-summary-excel)
  - [Document and Export Column Dictionary](#column-dictionary)
- [Output File Summary](#output-file-summary)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [Methods](#methods)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [Sample Sets Do Not Match](#sample-sets-do-not-match)
    - [`package 'phyloseq' not found`](#package-phyloseq-not-found)
    - [Very Low Minimum Sampling Depth / Most Data
      Discarded](#very-low-minimum-sampling-depth--most-data-discarded)
    - [Negative or Zero Cell Counts](#negative-or-zero-cell-counts)

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

This notebook is **Step 8** of the 16S rRNA sequencing data processing
pipeline, and is **entirely optional**. It converts the
copy-number-corrected ASV abundance table from [Step
7](7_copy_number_correction.md) into an *absolute*, cell-count-scaled
**Quantitative Microbiome Profile (QMP)**, following the method of
[Vandeputte et al. (2017), *Quantitative microbiome profiling links gut
community variation to microbial load*, *Nature*
551:507-511](https://doi.org/10.1038/nature24460).

<div class="alert alert-info">

**This step requires data this pipeline cannot generate for you**: an
independent measurement of **microbial load** (cells per gram/mL of
sample material) for every sample, typically obtained by flow cytometry.
If you do not have this measurement, skip this notebook — everything
downstream of Step 5 still works on relative abundances without it.

</div>

Standard amplicon sequencing output is inherently **compositional**:
every sample’s read counts are constrained to sum to (approximately) the
same sequencing depth, so read counts only ever carry information about
a taxon’s abundance *relative to every other taxon in the same sample*.
If total microbial load genuinely differs between samples or groups
(which it does substantially in, for example, gut microbiome studies of
disease vs. health, or across diet interventions), naive
relative-abundance analysis can produce **spurious associations** — a
taxon can appear to “increase” in relative terms purely because other
taxa dropped in absolute number, without that taxon’s own absolute
abundance having changed at all. Vandeputte et al. (2017) showed this
effect is large enough to change the direction of associations in real
gut microbiome cohorts.

Quantitative Microbiome Profiling (QMP) addresses this by rescaling each
sample’s relative abundances by its own independently-measured
**microbial load** (total cells per unit sample), so the resulting table
reflects genuine absolute-abundance differences between samples rather
than compositional artifacts. The method rarefies every sample to a
common **sampling depth per cell** (not a common read depth) before
rescaling, so that samples with different sequencing depths and
different microbial loads become directly comparable.

## Prerequisites

Before running this notebook, ensure that:

1.  **Step 7 (Copy Number Correction)** has been completed successfully,
    producing
    [copy_number_corrected_asv_count_table.csv](../../results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv)
    in
    [results/7_copy_number_correction/](../../results/7_copy_number_correction/).
    QMP specifically requires a **copy-number-corrected** table as input
    (Vandeputte et al. 2017) — do not use Step 5’s raw
    [asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv)
    here!
2.  You have an independent **microbial load measurement** (cells per
    gram or per mL) for every sample in your study, prepared as a TSV
    file (see [Cell Count File Format](#cell-count-format) below).
3.  Required R packages are installed — run
    [setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
    once per R environment if you have not already done so (this
    notebook additionally requires
    [phyloseq](https://joey711.github.io/phyloseq/), already a
    dependency of this workflow for Step 9).

## Cell Count File Format

Prepare a two-column TSV (tab-separated) file with one row per sample:

| SampleID | Cell_Count |
|:---------|:-----------|
| Sample_A | 3200000000 |
| Sample_B | 1850000000 |
| Sample_C | 4100000000 |

- `SampleID` must exactly match the sample identifiers used in Step
  5/Step 7’s ASV tables (the same identifiers derived from your FASTQ
  file names).
- `Cell_Count` is the estimated absolute microbial load for that sample
  (e.g. cells per gram of stool from flow cytometry), as a positive
  number. The unit itself does not matter to the calculation below as
  long as it is used consistently across all samples — the output table
  will be expressed in that same unit.

The project cell-count file is
[data/cell_count/cell_count.tsv](../../data/cell_count/cell_count.tsv),
which is also this notebook’s default expected input path. Replace its
values for a different dataset, or edit `cell_counts_path` in
[Configuration](#adjust-run-parameters) if you keep the file elsewhere.

<div class="alert alert-info">

**Provenance**: the Quantitative Microbiome Profiling (QMP) algorithm
implemented in this notebook follows the original
[QMP.R](https://github.com/raeslab/QMP) reference script (original
concept: Gwen Falony; original script contributors: Doris Vandeputte,
Gunter Kathagen, Kevin d’Hoe, Joao Sabino, Mireia Valles-Colomer, Sara
Vieira-Silva — [raeslab/QMP](https://github.com/raeslab/QMP)),
reimplemented here to match this project’s notebook conventions (path
configuration, logging, Excel export, Column_Dictionary documentation).
The underlying rarefaction-to-even-sampling-depth algorithm is unchanged
from the original.

</div>

## What This Notebook Does

1.  **Imports** Step 7’s copy-number-corrected ASV count table and your
    cell-count file.
2.  **Validates** that every sample in the corrected abundance table has
    a cell count; extra cell-count rows are reported and ignored.
3.  **Rarefies every sample to a common sampling depth per cell**, using
    [phyloseq](https://joey711.github.io/phyloseq/)’s
    `rarefy_even_depth()`, seeded for reproducibility.
4.  **Rescales** the rarefied, normalized relative abundances by each
    sample’s own cell count, producing the final
    microbial-load-corrected abundance table.
5.  **Exports** the corrected table as CSV, plus an Excel summary
    (per-sample sampling depth diagnostics) with a trailing
    `Column_Dictionary` sheet.

## Expected Input

- [results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv](../../results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv)
  — copy-number-corrected, sample x ASV abundance table from Step 7.
- [data/cell_count/cell_count.tsv](../../data/cell_count/cell_count.tsv)
  (default path; user-provided) — one row per sample, `SampleID` +
  `Cell_Count` columns.

## Expected Output

- [results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv](../../results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv)
  — the final Quantitative Microbiome Profile: sample x ASV absolute
  abundances.
- [results/8_microbial_load_correction/microbial_load_correction_summary.xlsx](../../results/8_microbial_load_correction/microbial_load_correction_summary.xlsx)
  — `Run_Provenance` (inputs, checksums, seed, threshold, and software
  version), `Sampling_Depth_Summary` (per-sample diagnostics), and a
  trailing `Column_Dictionary` sheet.

------------------------------------------------------------------------

# Environment Setup

<div class="alert alert-info">

**Before running this notebook**: install the required R packages by
running
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R).
Complete [Step 7](7_copy_number_correction.md) and prepare your
cell-count TSV file (see [Cell Count File Format](#cell-count-format))
before running.

</div>

## Load Required Packages

The chunk below loads every R package this notebook depends on, plus
three helper functions shared across the whole pipeline: Excel writing,
column-dictionary documentation, and clickable output links (each
sourced from [R/functions/](../functions/)).

``` r
# here: Project-relative file paths
library(here)

# openxlsx: Excel file creation and manipulation
library(openxlsx)

# readr: Fast, consistent delimited-file import/export
library(readr)

# dplyr: Data manipulation grammar
library(dplyr)

# phyloseq: Provides rarefy_even_depth(), the core rarefaction routine QMP is
# built on (already a dependency of this workflow for Step 9)
library(phyloseq)

# Source the custom Excel utility function from the project's function library
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source the custom column-dictionary utility function from the project's
# function library.
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source the custom output-links utility function from the project's
# function library.
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

Set the path to your cell count file and the rarefaction random seed
here before running the rest of the notebook. `cell_counts_path` is the
one setting most users need to change (see [Cell Count File
Format](#cell-count-format)); the rarefaction seed only needs to change
if you deliberately want a different (but still reproducible) random
subsample.

``` r
# ------------------------------------------------------------------
# Path to your cell-count TSV file (see "Cell Count File Format" above).
# ------------------------------------------------------------------
cell_counts_path <- here("data", "cell_count", "cell_count.tsv")  # <-- EDIT THIS if your file lives elsewhere

# ------------------------------------------------------------------
# Random seed for phyloseq::rarefy_even_depth(), so rarefaction is
# reproducible across runs. 711 matches the seed used in the original QMP.R
# reference script (https://github.com/raeslab/QMP) this notebook is adapted from.
# ------------------------------------------------------------------
microbial_load_rarefaction_seed <- 711

if (length(microbial_load_rarefaction_seed) != 1L ||
    !is.numeric(microbial_load_rarefaction_seed) ||
    !is.finite(microbial_load_rarefaction_seed) ||
    microbial_load_rarefaction_seed < 0 ||
    microbial_load_rarefaction_seed != round(microbial_load_rarefaction_seed) ||
    microbial_load_rarefaction_seed > (.Machine$integer.max - 1000000L)) {
  stop("microbial_load_rarefaction_seed must be one non-negative whole number safely within R's integer range.")
}
microbial_load_rarefaction_seed <- as.integer(microbial_load_rarefaction_seed)

# Warn when rarefaction retains less than this fraction of a sample's
# integerized copy-number-corrected abundance.
minimum_rarefaction_retention <- 0.10

cat(
  "Configuration:\n",
  "- Cell count file:", cell_counts_path, "\n",
  "- Rarefaction seed:", microbial_load_rarefaction_seed, "\n"
)
```

## Define Path Parameters

This chunk resolves every input and output path relative to the [project
root](../../) via `here()`, creates this step’s own output folder if it
does not already exist, and fails fast with an informative error if
either Step 7’s output or the cell count file cannot be found – so a
missing prerequisite is caught immediately, before any computation
begins.

``` r
# Base results folder for all pipeline outputs
results_folder <- here("results")

# Folder containing Step 7's copy-number-corrected output (this step's input)
step7_output_folder <- here(results_folder, "7_copy_number_correction")

# Specific output folder for this step (Step 8: Microbial Load Correction)
output_folder <- here(results_folder, "8_microbial_load_correction")
dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
checkpoints_folder <- here(output_folder, "checkpoints")
dir.create(checkpoints_folder, recursive = TRUE, showWarnings = FALSE)

# Path to Step 7's copy-number-corrected ASV count table
input_corrected_table_path <- here(step7_output_folder, "copy_number_corrected_asv_count_table.csv")

# This notebook's primary deliverable: the microbial-load-corrected abundance table
output_csv_path <- here(output_folder, "microbial_load_corrected_abundance_table.csv")

# Summary workbook
output_excel_path <- here(output_folder, "microbial_load_correction_summary.xlsx")

# Retire prior final deliverables immediately. If current input validation or
# computation fails, Step 9 must not consume an older QMP table as current.
previous_final_outputs <- c(output_csv_path, output_excel_path)
unlink(previous_final_outputs[file.exists(previous_final_outputs)])

# Confirm inputs actually exist before going any further
if (!file.exists(input_corrected_table_path)) {
  stop(
    "Could not find Step 7's copy-number-corrected ASV count table: ", input_corrected_table_path, "\n",
    "Run Step 7 (7_copy_number_correction.md) first -- QMP specifically requires a ",
    "copy-number-corrected table, not Step 5's raw asv_count_table.csv."
  )
}
if (!file.exists(cell_counts_path)) {
  stop(
    "Could not find the cell-count file: ", cell_counts_path, "\n",
    "Prepare it first -- see the \"Cell Count File Format\" section in the Introduction above."
  )
}

step7_provenance_path <- here(step7_output_folder, "copy_number_correction_summary.xlsx")
step5_asv_sequences_path <- here(results_folder, "5_dada2_pipeline", "asv_sequences.csv")
step5_asv_counts_path <- here(results_folder, "5_dada2_pipeline", "asv_count_table.csv")
if (!file.exists(step7_provenance_path) ||
    !all(file.exists(c(step5_asv_sequences_path, step5_asv_counts_path)))) {
  stop("Step 7 provenance or the current Step 5 ASV inputs are missing. Rerun Steps 5 and 7 before Step 8.")
}

cat(
  "Path Configuration:\n",
  "- Input corrected ASV table:", input_corrected_table_path, "\n",
  "- Input cell count file:", cell_counts_path, "\n",
  "- Output folder:", output_folder, "\n"
)
```

------------------------------------------------------------------------

# Data Import

## Import Copy-Number-Corrected ASV Table

This chunk reads Step 7’s copy-number-corrected ASV table and reshapes
it into a numeric matrix (samples x ASVs) with `SampleID` moved into the
row names – the format `rarefy_even_depth()` and the rescaling step
below both expect.

``` r
corrected_table_raw <- read_csv(
  input_corrected_table_path,
  col_types = cols(SampleID = col_character(), .default = col_double()),
  name_repair = "minimal",
  show_col_types = FALSE
)

if (anyDuplicated(names(corrected_table_raw))) {
  stop("Step 7's corrected table contains duplicated column names: ",
       paste(unique(names(corrected_table_raw)[duplicated(names(corrected_table_raw))]),
             collapse = ", "))
}

if (!"SampleID" %in% colnames(corrected_table_raw)) {
  stop(
    "Step 7's copy_number_corrected_asv_count_table.csv is missing the expected 'SampleID' column. ",
    "Found columns: ", paste(head(colnames(corrected_table_raw), 10), collapse = ", "), "..."
  )
}
if (nrow(corrected_table_raw) == 0L || anyNA(corrected_table_raw$SampleID) ||
    any(corrected_table_raw$SampleID == "") || anyDuplicated(corrected_table_raw$SampleID)) {
  stop("Step 7's corrected ASV table must contain at least one non-missing, unique SampleID.")
}

corrected_abundance_matrix <- as.matrix(corrected_table_raw %>% select(-SampleID))
rownames(corrected_abundance_matrix) <- corrected_table_raw$SampleID

if (!is.numeric(corrected_abundance_matrix) || ncol(corrected_abundance_matrix) == 0L) {
  stop("Step 7's corrected ASV columns must form a non-empty numeric matrix.")
}
if (any(!is.finite(corrected_abundance_matrix)) || any(corrected_abundance_matrix < 0)) {
  stop("Step 7's corrected ASV table contains non-finite or negative abundances.")
}

if (sum(corrected_abundance_matrix) <= 0) {
  stop("Step 7's corrected ASV table contains no positive abundance.")
}

# Validate Step 7 against the current Step 5 inputs and its own provenance.
if (!("Run_Provenance" %in% openxlsx::getSheetNames(step7_provenance_path))) {
  stop("Step 7's workbook lacks Run_Provenance. Rerun the current Step 7 notebook.")
}
step7_provenance <- openxlsx::read.xlsx(step7_provenance_path, sheet = "Run_Provenance")
required_step7_provenance <- c("ASV_Sequences_MD5", "ASV_Count_Table_MD5")
if (nrow(step7_provenance) != 1L ||
    !all(required_step7_provenance %in% names(step7_provenance))) {
  stop("Step 7's Run_Provenance sheet is missing required input checksums.")
}
current_step5_checksums <- c(
  ASV_Sequences_MD5 = unname(tools::md5sum(step5_asv_sequences_path)),
  ASV_Count_Table_MD5 = unname(tools::md5sum(step5_asv_counts_path))
)
recorded_step5_checksums <- c(
  ASV_Sequences_MD5 = as.character(step7_provenance$ASV_Sequences_MD5[[1]]),
  ASV_Count_Table_MD5 = as.character(step7_provenance$ASV_Count_Table_MD5[[1]])
)
if (!identical(as.character(recorded_step5_checksums),
               as.character(current_step5_checksums))) {
  stop("Step 7 was not generated from the current Step 5 ASV files. Rerun Step 7.")
}

step5_asv_sequences <- read_csv(
  step5_asv_sequences_path,
  col_types = cols(ASV_ID = col_character(), Sequence = col_character()),
  name_repair = "minimal", show_col_types = FALSE
)
if (!all(c("ASV_ID", "Sequence") %in% names(step5_asv_sequences))) {
  stop("Step 5's ASV sequence file must contain ASV_ID and Sequence columns.")
}
if (nrow(step5_asv_sequences) == 0L || anyNA(step5_asv_sequences$ASV_ID) ||
    any(step5_asv_sequences$ASV_ID == "") || anyDuplicated(step5_asv_sequences$ASV_ID) ||
    !setequal(colnames(corrected_abundance_matrix), step5_asv_sequences$ASV_ID)) {
  stop("Step 7's corrected ASV columns do not exactly match the current Step 5 ASV set.")
}
corrected_abundance_matrix <- corrected_abundance_matrix[, step5_asv_sequences$ASV_ID, drop = FALSE]

cat(
  "Imported Step 7 copy-number-corrected ASV table:\n",
  "- Samples:", nrow(corrected_abundance_matrix), "\n",
  "- ASVs:", ncol(corrected_abundance_matrix), "\n"
)
```

## Import Cell Count File

This chunk reads your cell count TSV file and validates it before use:
every `Cell_Count` must be a finite, positive number (a zero, missing,
or negative value would make the sampling-depth-per-cell calculation
below undefined), and every `SampleID` must be unique.

``` r
cell_counts_raw <- read_tsv(
  cell_counts_path,
  col_types = cols(SampleID = col_character(), Cell_Count = col_double()),
  name_repair = "minimal",
  show_col_types = FALSE
)

if (anyDuplicated(names(cell_counts_raw))) {
  stop("The cell count file contains duplicated column names: ",
       paste(unique(names(cell_counts_raw)[duplicated(names(cell_counts_raw))]),
             collapse = ", "))
}

required_cell_count_columns <- c("SampleID", "Cell_Count")
missing_cell_count_columns <- setdiff(required_cell_count_columns, colnames(cell_counts_raw))
if (length(missing_cell_count_columns) > 0) {
  stop(
    "Cell count file is missing expected column(s): ",
    paste(missing_cell_count_columns, collapse = ", "),
    ". Found columns: ", paste(colnames(cell_counts_raw), collapse = ", "),
    ". See the \"Cell Count File Format\" section in the Introduction."
  )
}

# Defensive check: a non-positive or missing cell count would make the
# sampling-depth-per-cell calculation below undefined (division by zero or a
# negative number).
invalid_cell_counts <- cell_counts_raw %>%
  filter(!is.finite(Cell_Count) | Cell_Count <= 0)

if (nrow(invalid_cell_counts) > 0) {
  stop(
    nrow(invalid_cell_counts), " sample(s) in the cell count file have a non-positive or ",
    "missing Cell_Count, which would make microbial load correction undefined: ",
    paste(invalid_cell_counts$SampleID, collapse = ", ")
  )
}

if (any(duplicated(cell_counts_raw$SampleID))) {
  stop(
    "The cell count file has duplicate SampleID entries: ",
    paste(unique(cell_counts_raw$SampleID[duplicated(cell_counts_raw$SampleID)]), collapse = ", ")
  )
}
if (anyNA(cell_counts_raw$SampleID) || any(cell_counts_raw$SampleID == "")) {
  stop("The cell count file contains a missing or empty SampleID.")
}

cat("Imported cell counts for", nrow(cell_counts_raw), "samples.\n")
```

## Validate Matching Sample Sets

QMP requires a cell count for every sample in the corrected abundance
table so every sample’s rescaling factor is well-defined. A cell-count
file may also contain additional samples that were excluded upstream;
those rows are reported and ignored.

``` r
samples_in_abundance_table <- rownames(corrected_abundance_matrix)
samples_in_cell_counts     <- cell_counts_raw$SampleID

missing_cell_counts   <- setdiff(samples_in_abundance_table, samples_in_cell_counts)
missing_from_abundance <- setdiff(samples_in_cell_counts, samples_in_abundance_table)

if (length(missing_cell_counts) > 0) {
  stop(
    "The corrected ASV table contains sample(s) without a cell count: ",
    paste(missing_cell_counts, collapse = ", "), ".\n",
    "Add those samples to the cell count file (or confirm sample identifiers match exactly) before proceeding."
  )
}

if (length(missing_from_abundance) > 0) {
  warning(
    "Ignoring ", length(missing_from_abundance),
    " cell-count sample(s) not present in the corrected ASV table: ",
    paste(missing_from_abundance, collapse = ", "),
    call. = FALSE
  )
}

# Retain and order only the measurements used by the corrected abundance table.
cell_counts_raw <- cell_counts_raw[
  match(samples_in_abundance_table, cell_counts_raw$SampleID), , drop = FALSE
]

cat(
  "Cell counts matched for all", length(samples_in_abundance_table),
  "abundance-table samples.\n"
)

pipeline_signature <- list(
  inputs = c(
    step7_corrected_table = unname(tools::md5sum(input_corrected_table_path)),
    step7_provenance = unname(tools::md5sum(step7_provenance_path)),
    cell_counts = unname(tools::md5sum(cell_counts_path))
  ),
  parameters = list(
    rarefaction_seed = microbial_load_rarefaction_seed,
    minimum_rarefaction_retention = minimum_rarefaction_retention
  ),
  phyloseq_version = as.character(packageVersion("phyloseq"))
)

checkpoint_path <- here(checkpoints_folder, "checkpoint_qmp.rds")
load_qmp_checkpoint <- function() {
  if (!file.exists(checkpoint_path)) return(NULL)
  checkpoint <- tryCatch(readRDS(checkpoint_path), error = function(e) NULL)
  required <- c("corrected_abundance_matrix_int", "cell_count_vector",
                "sample_read_totals", "sampling_depths",
                "minimum_sampling_depth", "rarefy_to_reads", "rarefaction_retention",
                "rarefied_matrix", "microbial_load_corrected_abundance_matrix")
  if (is.null(checkpoint) || !identical(checkpoint$signature, pipeline_signature) ||
      !all(required %in% names(checkpoint$objects))) {
    message("Ignoring incompatible QMP checkpoint.")
    unlink(checkpoint_path)
    return(NULL)
  }
  message("Resuming from validated QMP checkpoint: ", checkpoint_path)
  checkpoint$objects
}

save_qmp_checkpoint <- function(objects) {
  temporary_path <- paste0(checkpoint_path, ".tmp")
  saveRDS(list(signature = pipeline_signature, objects = objects),
          temporary_path, compress = TRUE)
  if (file.exists(checkpoint_path)) unlink(checkpoint_path)
  if (!file.rename(temporary_path, checkpoint_path)) {
    stop("Could not finalize QMP checkpoint: ", checkpoint_path)
  }
  checkpoint_path
}
```

------------------------------------------------------------------------

# Compute Quantitative Microbiome Profile

## Rarefy to a Common Sampling Depth Per Cell

Following Vandeputte et al. (2017): every sample’s copy-number-corrected
abundance total (rounded up to whole counts for rarefaction) is divided
by its own cell count to get a **sampling depth per cell**; every sample
is then rarefied to the minimum sampling depth observed across samples,
so samples differing in both sequencing depth and microbial load become
comparable before the final rescaling step.

``` r
qmp_checkpoint_objects <- load_qmp_checkpoint()
qmp_checkpoint_loaded <- !is.null(qmp_checkpoint_objects)
if (qmp_checkpoint_loaded) {
  list2env(qmp_checkpoint_objects, envir = environment())
} else {
# Round counts up to whole reads (rarefy_even_depth() requires integer
# counts); ceiling (rather than floor/round) matches the original QMP.R
# reference implementation exactly.
corrected_abundance_matrix_int <- ceiling(corrected_abundance_matrix)

# Re-order the cell count vector to exactly match the abundance table's
# sample (row) order -- guards against a silent sample mix-up even though
# Validate Matching Sample Sets above already confirmed both sets are equal.
cell_count_vector <- setNames(cell_counts_raw$Cell_Count, cell_counts_raw$SampleID)
cell_count_vector <- cell_count_vector[rownames(corrected_abundance_matrix_int)]

# Sampling depth per cell = total reads for this sample / this sample's cell count
sample_read_totals   <- rowSums(corrected_abundance_matrix_int)
sampling_depths       <- sample_read_totals / cell_count_vector
minimum_sampling_depth <- min(sampling_depths)

# Target read count to rarefy each sample down to, so every sample ends up
# at the same sampling-depth-per-cell (the minimum observed)
rarefy_to_reads <- floor(cell_count_vector * minimum_sampling_depth)

if (any(sample_read_totals <= 0) || any(rarefy_to_reads < 1)) {
  problem_samples <- names(sample_read_totals)[sample_read_totals <= 0 | rarefy_to_reads < 1]
  stop("Every sample must contain at least one corrected read and have a rarefaction target of at least one read. Problem samples: ",
       paste(problem_samples, collapse = ", "))
}

rarefaction_retention <- rarefy_to_reads / sample_read_totals
low_retention_samples <- names(rarefaction_retention)[
  rarefaction_retention < minimum_rarefaction_retention
]
if (length(low_retention_samples)) {
  warning(
    length(low_retention_samples), " sample(s) retain less than ",
    100 * minimum_rarefaction_retention, "% of integerized abundance during rarefaction: ",
    paste(low_retention_samples, collapse = ", "), call. = FALSE
  )
}

cat(
  "Sampling depth per cell (reads / cell count):\n",
  "- Range across samples:", signif(min(sampling_depths), 4), "-", signif(max(sampling_depths), 4), "\n",
  "- Minimum (target for all samples):", signif(minimum_sampling_depth, 4), "\n",
  "- Corresponding rarefaction target range (reads):", round(min(rarefy_to_reads)), "-", round(max(rarefy_to_reads)), "\n"
)

# Build a phyloseq otu_table (samples as rows, matching taxa_are_rows = FALSE)
# and rarefy one sample at a time, each to its own target read count
# (rarefy_to_reads[i]), exactly as in the original QMP.R reference script.
# The configured base seed plus each sample's alphabetical rank makes every
# sample's random subsample reproducible even if input rows are reordered.
corrected_abundance_phyloseq <- otu_table(corrected_abundance_matrix_int, taxa_are_rows = FALSE)

rarefied_matrix <- matrix(
  nrow = nrow(corrected_abundance_matrix_int),
  ncol = ncol(corrected_abundance_matrix_int),
  dimnames = dimnames(corrected_abundance_matrix_int)
)

for (sample_index in seq_len(nrow(corrected_abundance_matrix_int))) {
  # Seed assignment is based on alphabetic sample rank, not input row order,
  # so rearranging otherwise identical input rows does not change a sample's
  # random rarefaction realization.
  sample_seed_offset <- match(
    rownames(corrected_abundance_matrix_int)[sample_index],
    sort(rownames(corrected_abundance_matrix_int))
  ) - 1L
  rarefied_matrix[sample_index, ] <- rarefy_even_depth(
    corrected_abundance_phyloseq[sample_index, ],
    sample.size = rarefy_to_reads[sample_index],
    rngseed     = microbial_load_rarefaction_seed + sample_seed_offset,
    replace     = FALSE,
    trimOTUs    = FALSE,
    verbose     = FALSE
  )
}

if (!identical(dimnames(rarefied_matrix), dimnames(corrected_abundance_matrix_int)) ||
    any(!is.finite(rarefied_matrix)) || any(rarefied_matrix < 0) ||
    !identical(as.numeric(rowSums(rarefied_matrix)), as.numeric(rarefy_to_reads))) {
  stop("Rarefaction output failed sample/ASV alignment, value, or requested-depth validation.")
}

cat("Rarefaction complete for", nrow(rarefied_matrix), "samples.\n")
}
```

## Rescale by Cell Count

This step converts each sample’s rarefied counts to relative abundance,
then multiplies by that sample’s own cell count – turning “relative
abundance at a common sampling-depth-per-cell” into the final, absolute
Quantitative Microbiome Profile.

``` r
if (!qmp_checkpoint_loaded) {
# Normalize each sample's rarefied counts to relative abundance, then rescale
# by that sample's own cell count -- converting from "relative abundance at a
# common sampling-depth-per-cell" to an absolute, quantitative abundance.
normalised_rarefied_matrix <- rarefied_matrix / rowSums(rarefied_matrix)
microbial_load_corrected_abundance_matrix <- sweep(normalised_rarefied_matrix, MARGIN = 1, STATS = cell_count_vector, FUN = "*")
if (any(!is.finite(microbial_load_corrected_abundance_matrix)) ||
    any(microbial_load_corrected_abundance_matrix < 0) ||
    !isTRUE(all.equal(unname(rowSums(microbial_load_corrected_abundance_matrix)),
                     unname(cell_count_vector), tolerance = 1e-8))) {
  stop("The final QMP matrix failed finite-value or cell-count scaling validation.")
}

checkpoint_qmp <- save_qmp_checkpoint(list(
  corrected_abundance_matrix_int = corrected_abundance_matrix_int,
  cell_count_vector = cell_count_vector,
  sample_read_totals = sample_read_totals,
  sampling_depths = sampling_depths,
  minimum_sampling_depth = minimum_sampling_depth,
  rarefy_to_reads = rarefy_to_reads,
  rarefaction_retention = rarefaction_retention,
  rarefied_matrix = rarefied_matrix,
  microbial_load_corrected_abundance_matrix = microbial_load_corrected_abundance_matrix
))
}

# Validate loaded checkpoints as rigorously as newly computed results. The
# signature proves compatibility; these checks also reject damaged RDS data.
if (!identical(dimnames(rarefied_matrix), dimnames(corrected_abundance_matrix_int)) ||
    any(!is.finite(rarefied_matrix)) || any(rarefied_matrix < 0) ||
    !identical(as.numeric(rowSums(rarefied_matrix)), as.numeric(rarefy_to_reads))) {
  stop("Rarefaction checkpoint failed sample/ASV alignment, value, or requested-depth validation.")
}
if (!identical(dimnames(microbial_load_corrected_abundance_matrix),
               dimnames(corrected_abundance_matrix_int)) ||
    any(!is.finite(microbial_load_corrected_abundance_matrix)) ||
    any(microbial_load_corrected_abundance_matrix < 0) ||
    !isTRUE(all.equal(unname(rowSums(microbial_load_corrected_abundance_matrix)),
                     unname(cell_count_vector), tolerance = 1e-8))) {
  stop("The final QMP checkpoint failed alignment, finite-value, or cell-count scaling validation.")
}

cat(
  "Quantitative Microbiome Profile computed:\n",
  "- Samples:", nrow(microbial_load_corrected_abundance_matrix), "\n",
  "- ASVs:", ncol(microbial_load_corrected_abundance_matrix), "\n"
)
```

## Export the Microbial-Load-Corrected Table (CSV)

This chunk writes the final QMP abundance table to
[results/8_microbial_load_correction/](../../results/8_microbial_load_correction/)
as CSV, using the same sample order and ASV columns established above;
the exported file is linked immediately below via
`render_output_links()`.

``` r
microbial_load_corrected_abundance_table_export <- data.frame(
  SampleID = rownames(microbial_load_corrected_abundance_matrix),
  microbial_load_corrected_abundance_matrix,
  check.names = FALSE,
  row.names = NULL
)

write_csv(microbial_load_corrected_abundance_table_export, output_csv_path)

cat(
  "Quantitative Microbiome Profile written to:\n", output_csv_path, "\n",
  "(", nrow(microbial_load_corrected_abundance_table_export), "samples x", ncol(microbial_load_corrected_abundance_matrix), "ASVs )\n"
)
```

------------------------------------------------------------------------

# Results Summary

## Sampling Depth Diagnostics

The table below previews, per sample, the diagnostic values computed
during rarefaction above (read totals, cell counts, sampling depth, and
rarefaction target); only the first 10 rows are shown here, the full
table is written to the Excel workbook next.

``` r
sampling_depth_summary <- data.frame(
  `Sample ID`                    = rownames(corrected_abundance_matrix_int),
  `Total Corrected Reads`        = unname(sample_read_totals),
  `Cell Count`                   = unname(cell_count_vector),
  `Sampling Depth (reads/cell)`  = signif(unname(sampling_depths), 6),
  `Rarefied To (reads)`          = round(unname(rarefy_to_reads)),
  `Rarefaction Retention (%)`    = round(100 * unname(rarefaction_retention), 2),
  check.names = FALSE
)

# DT is namespace-qualified (not library()-loaded) here since this is its
# only use in this notebook; the fixed-height scroll box matches the
# interactive DT::datatable() tables used elsewhere in this pipeline.
head(sampling_depth_summary, 10) %>%
  DT::datatable(options = list(pageLength = 10, scrollX = TRUE,
                                scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
                rownames = FALSE,
                caption = "First 10 samples (full table in the Excel export)")
```

------------------------------------------------------------------------

# Export to Excel

## Write Summary Sheet

This chunk writes the sampling depth summary above to the output
workbook using the project’s shared Excel-writing helper,
[add_sheet_to_excel_function.R](../functions/add_sheet_to_excel_function.R).

## Document and Export Column Dictionary

This chunk documents every column of the provenance and sampling-depth
sheets in a trailing `Column_Dictionary` sheet, using the project’s
shared
[build_column_dictionary_function.R](../functions/build_column_dictionary_function.R)
helper – this project’s standard convention of shipping a data
dictionary alongside every exported workbook. The chunk’s own code and
output are hidden from the rendered report (`include=FALSE`), since only
the resulting workbook matters to the reader.

------------------------------------------------------------------------

# Output File Summary

The tree below lists outputs from the validated current run. Previous
final deliverables are retired when a run starts, and the stochastic
result is restored only from a checkpoint whose input checksums,
parameters, and phyloseq version match exactly.

------------------------------------------------------------------------

# Recommended Next Step

The Quantitative Microbiome Profile exported above is ready to be
incorporated into the rest of the pipeline. Proceed to [Step 9 —
Phyloseq Object Assembly](9_phyloseq_object.md), which automatically
detects and imports this notebook’s output
([results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv](../../results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv))
if present – alongside the outputs of Steps 5-7 – to build the final,
annotated `phyloseq` object(s) used for downstream statistical analysis
and visualization. If you skipped this notebook, Step 9 simply falls
back to the relative-abundance (or copy-number-corrected, if Step 7 was
run) ASV table instead, with no error.

------------------------------------------------------------------------

# Session Information

------------------------------------------------------------------------

# References

## Methods

- Vandeputte D, Kathagen G, D’hoe K, et al. (2017). Quantitative
  microbiome profiling links gut community variation to microbial load.
  *Nature* 551, 507-511. <https://doi.org/10.1038/nature24460>
- [raeslab/QMP GitHub repository](https://github.com/raeslab/QMP) —
  original reference implementation (`QMP.R`) this notebook’s algorithm
  follows.
- McMurdie PJ, Holmes S (2013). phyloseq: An R Package for Reproducible
  Interactive Analysis and Graphics of Microbiome Census Data. *PLoS
  ONE*, 8(4):e61217. (provides `rarefy_even_depth()`)
  <https://doi.org/10.1371/journal.pone.0061217>

## Related

- [Step 7 — 16S rRNA Gene Copy Number
  Correction](7_copy_number_correction.md) — this notebook’s required
  input.
- [Step 9 — Phyloseq Object Assembly](9_phyloseq_object.md) — this
  notebook’s recommended next step.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### Sample Sets Do Not Match

**Error**:
`The corrected ASV table contains sample(s) without a cell count.`

**Cause**: Every `SampleID` in the corrected abundance table must have
an exactly matching row in the cell count TSV. A common cause is a
spreadsheet program silently altering identifiers (e.g. converting a
purely numeric sample ID into a number, dropping a leading zero) when
the cell count TSV was created or edited. Extra rows in the cell count
file are allowed, reported as a warning, and ignored.

**Actions**:

- Compare the identifiers listed in the error message directly against
  [results/5_dada2_pipeline/asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv)’s
  `SampleID` column.
- Re-save the cell count TSV with `SampleID` explicitly formatted as
  text if your spreadsheet program is altering numeric-looking
  identifiers.

### `package 'phyloseq' not found`

**Cause**: `phyloseq` is a Bioconductor package installed by
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R),
shared with Step 9.

**Actions**:

- Run
  [setup/install_R_dependencies.R](../../setup/install_R_dependencies.R)
  if you have not already done so.

### Very Low Minimum Sampling Depth / Most Data Discarded

**Symptom**: The [Rarefy to a Common Sampling Depth Per
Cell](#rarefy-per-cell) chunk prints a warning about a very low minimum
sampling depth, or the resulting QMP table looks implausibly sparse.

**Cause**: Every sample is rarefied down to match the sample with the
*lowest* reads-per-cell ratio — one poorly-sequenced sample (low read
depth) or one sample with an implausibly high cell count can therefore
force heavy data loss across the *entire* dataset, not just that one
sample.

**Actions**:

- Identify the limiting sample from the `Sampling Depth (reads/cell)`
  column in the `Sampling_Depth_Summary` Excel sheet (the row with the
  smallest value).
- Consider whether that sample’s sequencing depth or cell count
  measurement is reliable; if not, exclude it from both the ASV table
  and the cell count file and re-run this notebook.
- If the low value is genuine (a real sample with unusually low absolute
  abundance relative to its sequencing depth), this data loss is an
  expected consequence of the QMP method itself (Vandeputte et al. 2017)
  rather than a bug.

### Negative or Zero Cell Counts

**Error**:
`sample(s) in the cell count file have a non-positive or missing Cell_Count`

**Cause**: Flow cytometry cell counts must be positive; a zero or
negative value makes the sampling-depth-per-cell ratio undefined
(division by zero) or nonsensical.

**Actions**:

- Re-check the source flow cytometry data for the listed sample(s); a
  zero value usually indicates a failed or missing measurement rather
  than a genuine zero microbial load.
