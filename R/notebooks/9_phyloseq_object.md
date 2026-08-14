Step 9: Phyloseq Object Construction
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
  - [Define Data Selection Parameters](#define-data-params)
- [Import Data](#import-data)
  - [Load Abundance Tables](#load-abundance)
  - [Load Taxonomy Table(s)](#load-taxonomy)
  - [Load Sample Metadata (Optional)](#load-metadata)
  - [Load Phylogenetic Tree (Optional)](#load-tree)
- [Construct Phyloseq Object(s)](#construct-phyloseq-objects)
  - [Build Phyloseq Objects](#build-phyloseq)
  - [Explore Phyloseq Objects](#explore-phyloseq)
  - [Save Phyloseq Objects](#save-phyloseq)
- [Taxonomic Aggregation](#taxonomic-aggregation)
  - [Aggregate at Genus Level](#aggregate-genus)
  - [Identify Top Genera](#top-genera)
- [Visualization](#visualization)
  - [Define Barplot Function](#define-barplot)
  - [Assign Consistent Genus Colors](#genus-colors)
  - [Generate Genus-Level Barplot](#genus-barplot)
- [Sample Statistics](#sample-statistics)
  - [Per-Sample Abundance](#sample-stats)
  - [Taxonomy Summary](#taxonomy-summary)
  - [Record Run Provenance](#run-provenance)
  - [Document and Export Column Dictionary](#column-dictionary)
- [Output File Summary](#output-file-summary)
- [Interactive Exploration with
  Shiny-Phyloseq](#interactive-exploration-with-shiny-phyloseq)
- [Recommended Next Step](#recommended-next-step)
- [Session Information](#session-information)
- [References](#references)
  - [Methods](#methods)
  - [Databases](#databases)
  - [Related](#related)
- [Appendix: Troubleshooting Guide](#appendix-troubleshooting-guide)
  - [Common Issues and Solutions](#common-issues-and-solutions)
    - [Fewer Phyloseq Objects Than
      Expected](#fewer-phyloseq-objects-than-expected)
    - [ASV Names Don’t Match](#asv-names-dont-match)
    - [Phylogenetic Tree Tips Don’t
      Match](#phylogenetic-tree-tips-dont-match)
    - [Memory Issues with Large
      Datasets](#memory-issues-with-large-datasets)
    - [Missing Taxonomy Ranks](#missing-taxonomy-ranks)

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
</style>

# Introduction

## Purpose

This notebook is **Step 9**, the final step, of the 16S rRNA sequencing
data processing pipeline. It constructs one or more
[phyloseq](https://joey711.github.io/phyloseq/) objects from the ASV
count table(s) and taxonomy assignments generated earlier in the
pipeline. The phyloseq object is a powerful data structure that
integrates sample data, taxonomy, and optionally phylogenetic trees for
comprehensive microbiome analysis.

Because this notebook runs last, it automatically detects and picks up
the output of every optional upstream step that was run: for each
taxonomy database detected on disk ([SILVA](https://www.arb-silva.de/)
and/or [GTDB](https://gtdb.ecogenomic.org/) – whichever one or both were
produced by Step 5), it builds **up to three phyloseq objects** — one
from **raw counts** (Step 5, always available), one from
**copy-number-corrected counts** (Step 7, if that optional notebook was
run), and one from **microbial-load-corrected counts** (Step 8, if that
optional notebook was run). Any combination whose input file is not
found is skipped automatically, with a clear message explaining why.
Every database’s results — phyloseq objects, barplots, and summary
workbook — are kept in their own subfolder, so SILVA- and GTDB-based
outputs never mix.

## Prerequisites

Before running this notebook, ensure that:

1.  **[Step 5 (DADA2 Pipeline)](5_dada2_pipeline.md)** has been
    completed successfully.
2.  The ASV count table
    ([asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv))
    is present in
    [results/5_dada2_pipeline/](../../results/5_dada2_pipeline/).
3.  At least one taxonomy table is present
    ([silva_taxonomy_table.csv](../../results/5_dada2_pipeline/silva_taxonomy_table.csv)
    and/or
    [gtdb_taxonomy_table.csv](../../results/5_dada2_pipeline/gtdb_taxonomy_table.csv)).
4.  **Optional**: [Step 7 (Copy Number
    Correction)](7_copy_number_correction.md) has been run, to
    additionally build a copy-number-corrected phyloseq object per
    taxonomy database.
5.  **Optional**: [Step 8 (Microbial Load
    Correction)](8_microbial_load_correction.md) has been run, to
    additionally build a microbial-load-corrected phyloseq object per
    taxonomy database.
6.  **Optional**: Sample metadata file in TSV (or CSV) format – see
    [data/README.md](../../data/README.md) for the required layout. The
    ready-to-edit [data/metadata.tsv](../../data/metadata.tsv) contains
    ten generic placeholder samples (`S01`–`S10`); replace its rows and
    values with your own before drawing conclusions from metadata-based
    groupings.
7.  **Optional**: Phylogenetic trees from [Step
    6](6_phylogenetic_tree.md).

## What This Notebook Does

The workflow accomplishes the following tasks:

1.  **Import Data**: Loads whichever abundance table(s) are available
    (raw, copy-number-corrected, microbial-load-corrected), taxonomy
    assignments, and optional metadata/phylogenetic tree.
2.  **Construct Phyloseq Objects**: Creates one integrated data
    structure per (taxonomy database x available abundance source)
    combination.
3.  **Data Exploration**: Summarizes the composition and taxonomy of
    every phyloseq object built.
4.  **Taxonomic Aggregation**: Aggregates data at the Genus level for
    every phyloseq object.
5.  **Relative Abundance**: Calculates relative abundance
    transformations.
6.  **Visualization**: Generates an interactive top-taxa Genus-level
    barplot for each phyloseq object that contains at least one
    classified genus, using a shared genus-to-color mapping.
7.  **Export**: Saves every phyloseq object for downstream analysis.

## Expected Input

- **Always required**, from
  [results/5_dada2_pipeline/](../../results/5_dada2_pipeline/) (Step 5):
  - [asv_count_table.csv](../../results/5_dada2_pipeline/asv_count_table.csv)
    — raw ASV abundance matrix.
  - [silva_taxonomy_table.csv](../../results/5_dada2_pipeline/silva_taxonomy_table.csv)
    and/or
    [gtdb_taxonomy_table.csv](../../results/5_dada2_pipeline/gtdb_taxonomy_table.csv)
    — taxonomy assignments.
- **Optional**, automatically included if present:
  - [results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv](../../results/7_copy_number_correction/copy_number_corrected_asv_count_table.csv)
    (Step 7).
  - [results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv](../../results/8_microbial_load_correction/microbial_load_corrected_abundance_table.csv)
    (Step 8).
  - [results/6_phylogenetic_tree/phylogenetic_tree.nwk](../../results/6_phylogenetic_tree/phylogenetic_tree.nwk)
    (Step 6).
  - Sample metadata TSV (or CSV) file – see
    [data/README.md](../../data/README.md).

## Expected Output

Every output below is split into one self-contained subfolder per
taxonomy database under
[results/9_phyloseq_object/](../../results/9_phyloseq_object/) –
e.g. [results/9_phyloseq_object/SILVA/](../../results/9_phyloseq_object/SILVA/)
and/or
[results/9_phyloseq_object/GTDB/](../../results/9_phyloseq_object/GTDB/)
– rather than a single flat folder mixing both databases together, since
a phyloseq object (and everything derived from one) is always specific
to exactly one taxonomy database. Within each database’s own subfolder,
the `.RData` objects and the interactive barplots are further split into
their own
[phyloseq_objects/](../../results/9_phyloseq_object/SILVA/phyloseq_objects/)
and [barplots/](../../results/9_phyloseq_object/SILVA/barplots/)
subfolders.

- Up to three phyloseq objects per taxonomy database, saved as `.RData`
  files (one per available abundance source),
  e.g. [SILVA/phyloseq_objects/phyloseq_object_silva_raw_counts.RData](../../results/9_phyloseq_object/SILVA/phyloseq_objects/phyloseq_object_silva_raw_counts.RData).
- One interactive HTML barplot per phyloseq object containing classified
  genera – the top `n_top_taxa` genera by relative abundance, excluding
  unclassified/NA genera. Objects with no classified genus are retained,
  but their top-genus table and plot are skipped with a warning.
- A summary statistics Excel workbook per taxonomy database,
  e.g. [SILVA/phyloseq_summary_silva.xlsx](../../results/9_phyloseq_object/SILVA/phyloseq_summary_silva.xlsx),
  containing `Tax_Summary`, per-source `Top_Genera_*`/`Sample_Stats_*`,
  `Run_Provenance`, and a trailing `Column_Dictionary` sheet.
- Genus-level aggregated objects are intermediate in-memory objects used
  to generate the summary tables and barplots; they are not separate
  disk outputs.

------------------------------------------------------------------------

# Environment Setup

<div class="alert alert-info">

**Before running this notebook**: install the required R packages by
running
[setup/install_R_dependencies.R](../../setup/install_R_dependencies.R).
Complete [Step 5](5_dada2_pipeline.md) at minimum; [Step
6](6_phylogenetic_tree.md), [Step 7](7_copy_number_correction.md), and
[Step 8](8_microbial_load_correction.md) are optional and are picked up
automatically if present.

</div>

## Load Required Packages

The chunk below loads every R package this notebook depends on, plus
three helper functions shared across the whole pipeline: Excel writing,
column-dictionary documentation, and clickable output links (each
sourced from [R/functions/](../functions)).

``` r
# phyloseq: Core package for microbiome data analysis
# Provides data structures and functions for microbiome analysis
library(phyloseq)

# ape: Analyses of Phylogenetics and Evolution
# Required for phylogenetic tree handling in phyloseq
library(ape)

# data.table: High-performance data manipulation
# Provides fast operations for handling large data tables
library(data.table)

# DT: Interactive tables in R Markdown
# Creates searchable, sortable HTML tables
library(DT)

# dplyr: Data manipulation grammar
# Provides select() and the %>% pipe used throughout this notebook
library(dplyr)

# fs: Cross-platform filesystem operations
# Consistent functions for file and directory manipulation
library(fs)

# ggplot2: Grammar of graphics plotting
# Used for creating publication-quality visualizations
library(ggplot2)

# here: Project-relative file paths
# Enables reproducible path construction regardless of working directory
library(here)

# htmlwidgets: Framework for creating HTML widgets
# Used for exporting interactive plots
library(htmlwidgets)

# openxlsx: Excel file creation and manipulation
# Read and write Excel files without Java dependencies
library(openxlsx)

# plotly: Interactive plotting
# Converts ggplot objects to interactive HTML visualizations
library(plotly)

# RColorBrewer: Color palettes for visualization
# Provides aesthetically pleasing color schemes
library(RColorBrewer)

# readr: Fast, consistent delimited-file import
# Used to read every abundance table in this notebook's standard "SampleID +
# one column per ASV" layout
library(readr)

# stringr: Consistent string manipulation
# Intuitive string processing functions
library(stringr)

# Source our custom Excel utility function from the project's function library
# This function handles creating or appending sheets to Excel workbooks
source(here("R", "functions", "add_sheet_to_excel_function.R"))

# Source our custom column dictionary builder from the project's function library
# This function documents every column of a sheet before it is exported to Excel
source(here("R", "functions", "build_column_dictionary_function.R"))

# Source our custom output-links helper from the project's function library
# This function renders clickable Markdown links to output files/folders
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
# Input folder containing DADA2 pipeline results (Step 5, always required)
dada2_results_folder <- here("results", "5_dada2_pipeline")

# Input folder containing the phylogenetic tree (optional, from Step 6)
tree_results_folder <- here("results", "6_phylogenetic_tree")

# Input folder containing copy-number-corrected counts (optional, from Step 7)
copy_number_results_folder <- here("results", "7_copy_number_correction")

# Input folder containing microbial-load-corrected counts (optional, from Step 8)
microbial_load_results_folder <- here("results", "8_microbial_load_correction")

# Base results folder for all pipeline outputs
results_folder <- here("results")

# Specific output folder for this step (Step 9: Phyloseq Object)
final_output_folder <- here(results_folder, "9_phyloseq_object")
output_folder <- here(results_folder, ".9_phyloseq_object_staging")

# Create the output directory if it does not already exist
if (dir_exists(output_folder)) dir_delete(output_folder)
dir_create(output_folder, recurse = TRUE)
cat("Created isolated staging directory for the current Step 9 run.\n")

# Every output this notebook writes -- phyloseq .RData objects, barplots, and
# the summary workbook -- is split into one self-contained subfolder per
# taxonomy database (e.g. results/9_phyloseq_object/SILVA/,
# results/9_phyloseq_object/GTDB/) rather than a single flat folder mixing
# both databases' files together, since a phyloseq object (and everything
# derived from one) is always specific to exactly one taxonomy database.
# Within each database's own subfolder, the .RData phyloseq objects and the
# interactive HTML barplots are further split into their own "phyloseq_objects/"
# and "barplots/" subfolders, so the two file types never mix even within a
# single database. The summary Excel workbook is neither a raw object nor a
# plot, so it stays directly in the database's own top-level subfolder rather
# than nested in either.
# All helpers below are functions of db_name (not fixed path variables) --
# called with the current combo's meta$db_name everywhere a
# database-specific path is needed -- rather than fixed variables, because
# which database(s) are actually processed is only known once Load Taxonomy
# Table(s) below has run. Each folder helper creates its folder on first use,
# so every chunk that writes a database-specific file can call it directly
# without a separate setup step.

# Returns (and creates, if necessary) the top-level per-database output
# subfolder for db_name (e.g. "SILVA" or "GTDB").
db_output_folder <- function(db_name) {
    folder <- here(output_folder, db_name)
    dir_create(folder)
    folder
}

# Returns (and creates) db_name's "phyloseq_objects/" subfolder, for its
# saved .RData phyloseq objects.
db_phyloseq_folder <- function(db_name) {
    folder <- here(db_output_folder(db_name), "phyloseq_objects")
    dir_create(folder)
    folder
}

# Returns (and creates) db_name's "barplots/" subfolder, for its interactive
# HTML barplots.
db_barplot_folder <- function(db_name) {
    folder <- here(db_output_folder(db_name), "barplots")
    dir_create(folder)
    folder
}

# Returns the path to db_name's own summary Excel workbook, directly inside
# its top-level per-database output subfolder (not nested in either
# phyloseq_objects/ or barplots/).
db_excel_path <- function(db_name) {
    here(db_output_folder(db_name), paste0("phyloseq_summary_", tolower(db_name), ".xlsx"))
}
```

## Define Data Selection Parameters

Configure the optional components to include. Which **taxonomy
database(s)** (SILVA / GTDB) and **abundance sources** (raw /
copy-number-corrected / microbial-load-corrected) get included is not
configured here — both are detected automatically, in [Load Taxonomy
Table(s)](#load-taxonomy) and [Load Abundance Tables](#load-abundance)
below respectively, based on which output files are actually present on
disk.

``` r
# Path to sample metadata file (optional)
# Set to NULL if no metadata is available. Accepts either a tab-separated
# (.tsv/.txt) or comma-separated (.csv) file, detected automatically from the
# file extension -- sample identifiers must be in the first column, matching
# the SampleID convention used throughout this workflow (the text before the
# first underscore in each FASTQ filename; see data/README.md). See
# data/README.md for the full required format.
#
# Defaults to the ready-to-edit template shipped at data/metadata.tsv, which
# contains generic S01-S10 sample IDs and placeholder values. Replace its rows
# and values with your own before drawing conclusions from metadata-based
# groupings. Set to NULL to build phyloseq
# objects without any metadata.
metadata_path <- here("data", "metadata.tsv")

# Whether to include the phylogenetic tree (if available from Step 6)
include_tree <- TRUE

# Number of top taxa to display in visualizations
n_top_taxa <- 10

if (length(include_tree) != 1L || !is.logical(include_tree) || is.na(include_tree)) {
    stop("include_tree must be exactly TRUE or FALSE.")
}
if (length(n_top_taxa) != 1L || !is.numeric(n_top_taxa) || is.logical(n_top_taxa) ||
    !is.finite(n_top_taxa) || n_top_taxa < 1 ||
    n_top_taxa != as.integer(n_top_taxa)) {
    stop("n_top_taxa must be one positive integer.")
}

# Display configuration
cat("Data Selection Configuration:\n",
    "- Taxonomy database(s): auto-detected below, in Load Taxonomy Table(s)\n",
    "- Metadata file:", ifelse(is.null(metadata_path), "Not provided", metadata_path), "\n",
    "- Include phylogenetic tree:", include_tree, "\n",
    "- Top taxa for visualization:", n_top_taxa, "\n")
```

------------------------------------------------------------------------

# Import Data

## Load Abundance Tables

This notebook can automatically build phyloseq objects from up to three
abundance tables. **Raw counts** (Step 5) are always required;
**copy-number-corrected** (Step 7) and **microbial-load-corrected**
(Step 8) counts are picked up automatically if their output files exist,
and skipped – with a message, not an error – otherwise.

``` r
# Registry of every abundance-table source this notebook knows how to build a
# phyloseq object from. `required = TRUE` means the notebook stops if the
# file is missing; `required = FALSE` means it is silently skipped (with a
# message) if the corresponding optional upstream step was not run.
# `label` is a short, Excel-sheet-name-safe tag (sheet names cannot exceed 31
# characters) used only when constructing per-sheet names later in this
# notebook -- file names, list/variable names, and console messages all use
# the full `display` name instead.
abundance_source_definitions <- list(
    raw_counts = list(
        display  = "Raw Counts",
        label    = "Raw",
        path     = here(dada2_results_folder, "asv_count_table.csv"),
        required = TRUE,
        provenance = NULL
    ),
    copy_number_corrected = list(
        display  = "Copy-Number Corrected",
        label    = "CN_Corrected",
        path     = here(copy_number_results_folder, "copy_number_corrected_asv_count_table.csv"),
        required = FALSE,
        provenance = here(copy_number_results_folder, "copy_number_correction_summary.xlsx")
    ),
    microbial_load_corrected = list(
        display  = "Microbial-Load Corrected",
        label    = "ML_Corrected",
        path     = here(microbial_load_results_folder, "microbial_load_corrected_abundance_table.csv"),
        required = FALSE,
        provenance = here(microbial_load_results_folder, "microbial_load_correction_summary.xlsx")
    )
)
```

``` r
# Read a Sample x ASV Abundance Table
#
# Reads a CSV abundance table written in this workflow's standard layout (a
# leading `SampleID` column followed by one column per ASV) and returns it as
# a numeric matrix with sample identifiers as row names.
#
# Arguments:
#   path - Path to the CSV file.
# Returns:
#   A numeric matrix (samples x ASVs).
read_abundance_table <- function(path) {
    table_raw <- read_csv(
        path,
        col_types = cols(SampleID = col_character(), .default = col_double()),
        name_repair = "minimal", show_col_types = FALSE
    )

    if (anyNA(names(table_raw)) || any(trimws(names(table_raw)) == "") ||
        anyDuplicated(names(table_raw))) {
        stop("Abundance table contains missing, empty, or duplicated column names: ", path)
    }

    if (!"SampleID" %in% colnames(table_raw)) {
        stop("Abundance table is missing the expected 'SampleID' column: ", path)
    }
    if (nrow(table_raw) == 0L || anyNA(table_raw$SampleID) ||
        any(table_raw$SampleID == "") || anyDuplicated(table_raw$SampleID)) {
        stop("Abundance table must contain at least one non-missing, unique SampleID: ", path)
    }

    abundance_matrix <- as.matrix(table_raw %>% select(-SampleID))
    rownames(abundance_matrix) <- table_raw$SampleID

    # Step 5 already writes canonical sample identifiers. Preserve them
    # verbatim here because underscores can be legitimate identifier content.
    if (!is.numeric(abundance_matrix) || ncol(abundance_matrix) == 0L) {
        stop("Abundance table must contain at least one numeric ASV column: ", path)
    }
    if (any(!is.finite(abundance_matrix)) || any(abundance_matrix < 0)) {
        stop("Abundance table contains non-finite or negative values: ", path)
    }
    if (any(rowSums(abundance_matrix) <= 0)) {
        stop("Every sample must have positive total abundance in: ", path)
    }

    abundance_matrix
}

# Load every abundance source that is available, skipping optional ones that
# are not, and stopping only if a required one (raw counts) is missing.
abundance_tables <- list()

step5_asv_sequences_path <- here(dada2_results_folder, "asv_sequences.csv")
if (!all(file_exists(c(step5_asv_sequences_path,
                       abundance_source_definitions$raw_counts$path)))) {
    stop("Required Step 5 ASV sequence/count inputs are missing. Rerun Step 5.")
}
current_step5_checksums <- c(
    ASV_Sequences_MD5 = unname(tools::md5sum(step5_asv_sequences_path)),
    ASV_Count_Table_MD5 = unname(tools::md5sum(abundance_source_definitions$raw_counts$path))
)

validate_optional_provenance <- function(source_key, source_def) {
    if (is.null(source_def$provenance) || !file_exists(source_def$provenance)) {
        stop("The ", source_key, " table exists but its provenance workbook is missing: ",
             source_def$provenance, ". Rerun its upstream notebook.")
    }
    sheets <- openxlsx::getSheetNames(source_def$provenance)
    if (!"Run_Provenance" %in% sheets) {
        stop("The ", source_key, " provenance workbook lacks Run_Provenance. Rerun its upstream notebook.")
    }
    provenance <- openxlsx::read.xlsx(source_def$provenance, sheet = "Run_Provenance")
    if (source_key == "copy_number_corrected") {
        required <- names(current_step5_checksums)
        if (nrow(provenance) != 1L || !all(required %in% names(provenance)) ||
            !identical(as.character(unlist(provenance[1, required], use.names = FALSE)),
                       as.character(unname(current_step5_checksums)))) {
            stop("Step 7 outputs were not generated from the current Step 5 inputs. Rerun Step 7.")
        }
    }
    if (source_key == "microbial_load_corrected") {
        required <- c("Step7_Corrected_Table_MD5", "Step7_Provenance_MD5")
        step7_table <- abundance_source_definitions$copy_number_corrected$path
        step7_book <- abundance_source_definitions$copy_number_corrected$provenance
        if (nrow(provenance) != 1L || !all(required %in% names(provenance)) ||
            !all(file_exists(c(step7_table, step7_book))) ||
            !identical(as.character(provenance$Step7_Corrected_Table_MD5[[1]]),
                       as.character(unname(tools::md5sum(step7_table)))) ||
            !identical(as.character(provenance$Step7_Provenance_MD5[[1]]),
                       as.character(unname(tools::md5sum(step7_book))))) {
            stop("Step 8 outputs are stale relative to the current Step 7 outputs. Rerun Step 8.")
        }
    }
}

for (source_key in names(abundance_source_definitions)) {

    source_def <- abundance_source_definitions[[source_key]]

    if (file_exists(source_def$path)) {

        if (!source_def$required) validate_optional_provenance(source_key, source_def)

        abundance_tables[[source_key]] <- read_abundance_table(source_def$path)

        cat(source_def$display, "abundance table found:\n",
            "   ", source_def$path, "\n",
            "  Samples:", nrow(abundance_tables[[source_key]]), "\n",
            "  ASVs:", ncol(abundance_tables[[source_key]]), "\n")

    } else if (source_def$required) {

        stop(
            source_def$display, " abundance table not found: ", source_def$path,
            "\nPlease ensure Step 5 (DADA2 Pipeline) has been completed."
        )

    } else {

        cat("-", source_def$display, "abundance table not found -- skipping.\n",
            "   Expected at:", source_def$path, "\n",
            "   Run the corresponding optional notebook first if you want this phyloseq object.\n")

    }
}
```

## Load Taxonomy Table(s)

This chunk automatically detects which taxonomy database(s) Step 5
produced – SILVA, GTDB, or both – and loads whichever are actually
present, exactly mirroring the auto-detection already used for the
abundance sources above ([Load Abundance Tables](#load-abundance)). A
database whose file is missing is skipped with a message, not an error;
the chunk only stops if *neither* database’s file can be found, since a
phyloseq object cannot be built without at least one taxonomy
assignment.

``` r
# Registry of every taxonomy database this notebook knows how to look for.
# Mirrors the abundance_source_definitions registry above: each entry is
# skipped, with a message, if its file is not found -- there is no
# "required" database, since any one of SILVA or GTDB is sufficient to
# build a phyloseq object.
taxonomy_database_paths <- list(
    SILVA = here(dada2_results_folder, "silva_taxonomy_table.csv"),
    GTDB  = here(dada2_results_folder, "gtdb_taxonomy_table.csv")
)

taxonomy_files <- list()

for (db_name in names(taxonomy_database_paths)) {

    db_path <- taxonomy_database_paths[[db_name]]

    if (file_exists(db_path)) {
        taxonomy_files[[db_name]] <- db_path
        cat(db_name, "taxonomy file found:", db_path, "\n")
    } else {
        cat("-", db_name, "taxonomy file not found -- skipping.\n",
            "   Expected at:", db_path, "\n")
    }
}
```

``` r
# Check that at least one taxonomy file was found
if (length(taxonomy_files) == 0) {
    stop("No taxonomy files found! Expected at least one of:\n",
         paste0("  - ", unlist(taxonomy_database_paths), collapse = "\n"),
         "\nPlease ensure Step 5 (DADA2 Pipeline) has been completed with at least one taxonomy database.")
}

cat("\nTaxonomy database(s) available for this run:", paste(names(taxonomy_files), collapse = ", "), "\n")
```

## Load Sample Metadata (Optional)

If `metadata_path` was set above, this chunk loads it (auto-detecting a
tab- or comma-separated file from its extension – see
[data/README.md](../../data/README.md) for the required format) and
checks that its sample identifiers actually match the abundance table’s
samples – a phyloseq object is still built without metadata, but any
downstream grouping/faceting by a metadata variable requires this to
succeed.

``` r
sample_metadata <- NULL

if (!is.null(metadata_path) && !file_exists(metadata_path)) {
    stop("Configured metadata file does not exist: ", metadata_path,
         ". Set metadata_path <- NULL to run without metadata.")
}

if (!is.null(metadata_path) && file_exists(metadata_path)) {

    # Delimiter is inferred from the file extension, so this notebook accepts
    # either a tab-separated (.tsv/.txt, this project's convention -- see
    # data/README.md) or comma-separated (.csv) metadata file without any
    # extra configuration.
    if (!grepl("\\.(tsv|txt|csv)$", metadata_path, ignore.case = TRUE)) {
        stop("Unsupported metadata extension. Use .tsv, .txt, or .csv: ", metadata_path)
    }
    metadata_delimiter <- if (grepl("\\.(tsv|txt)$", metadata_path, ignore.case = TRUE)) "\t" else ","

    metadata_raw <- read_delim(
        metadata_path, delim = metadata_delimiter,
        col_types = cols(.default = col_character()),
        name_repair = "minimal", show_col_types = FALSE, trim_ws = TRUE
    )
    if (ncol(metadata_raw) < 2L || names(metadata_raw)[[1]] != "SampleID" ||
        anyDuplicated(names(metadata_raw))) {
        stop("Metadata must begin with a unique SampleID column and contain no duplicated headers.")
    }
    metadata_ids <- metadata_raw$SampleID
    if (!length(metadata_ids) || anyNA(metadata_ids) ||
        any(trimws(metadata_ids) == "") || anyDuplicated(metadata_ids)) {
        stop("Metadata must contain at least one non-empty, unique SampleID.")
    }
    sample_metadata <- as.data.frame(metadata_raw[-1], stringsAsFactors = FALSE,
                                     check.names = FALSE)
    rownames(sample_metadata) <- metadata_ids

    if (identical(basename(metadata_path), "metadata.tsv") &&
        identical(rownames(sample_metadata), sprintf("S%02d", 1:10)) &&
        "SubjectID" %in% names(sample_metadata) &&
        identical(as.character(sample_metadata$SubjectID), sprintf("Subject%02d", 1:10))) {
        stop("data/metadata.tsv is still the generic S01-S10 template. Replace it with experimental metadata or set metadata_path <- NULL.")
    }

    cat("\nSample metadata loaded:\n",
        "  Samples:", nrow(sample_metadata), "\n",
        "  Variables:", ncol(sample_metadata), "\n",
        "  Columns:", paste(colnames(sample_metadata), collapse = ", "), "\n")

    # Align metadata to the complete abundance-table sample set. Missing
    # metadata rows are filled with NA rather than allowing phyloseq() to
    # prune those biological samples to the metadata/abundance intersection.
    abundance_sample_ids <- rownames(abundance_tables[["raw_counts"]])
    matching_samples <- intersect(abundance_sample_ids, rownames(sample_metadata))
    missing_metadata_samples <- setdiff(abundance_sample_ids, rownames(sample_metadata))
    extra_metadata_samples <- setdiff(rownames(sample_metadata), abundance_sample_ids)

    if (length(matching_samples) == 0L) {
        stop("Metadata has no SampleID values matching the current abundance table. ",
             "Replace the template rows with your samples or set metadata_path <- NULL.")
    }

    if (length(missing_metadata_samples) > 0L) {
        warning(
            "Metadata is missing for ", length(missing_metadata_samples), " abundance-table sample(s): ",
            paste(missing_metadata_samples, collapse = ", "),
            ". Their metadata fields will be filled with NA; the samples will remain in every phyloseq object."
        )
    }
    if (length(extra_metadata_samples) > 0L) {
        warning(
            "Ignoring ", length(extra_metadata_samples), " metadata row(s) not present in the abundance table: ",
            paste(extra_metadata_samples, collapse = ", ")
        )
    }

    sample_metadata <- sample_metadata[abundance_sample_ids, , drop = FALSE]
    rownames(sample_metadata) <- abundance_sample_ids

    datatable(
        sample_metadata,
        options = list(pageLength = 5, scrollX = TRUE,
                        scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
        caption = "Sample metadata preview"
    )
} else {
    cat("No sample metadata provided. Phyloseq objects will be created without sample data.\n")
}
```

## Load Phylogenetic Tree (Optional)

If `include_tree` is `TRUE` and Step 6’s output is present, this chunk
loads the Newick tree so it can be incorporated into the phyloseq
object(s) below (enabling phylogeny-aware analyses like UniFrac
downstream); it is skipped gracefully, with an explanatory message, if
Step 6 was not run.

``` r
phylo_tree <- NULL

if (include_tree) {
    tree_path <- here(tree_results_folder, "phylogenetic_tree.nwk")

    if (file_exists(tree_path)) {
        phylo_tree <- read.tree(tree_path)

        if (!inherits(phylo_tree, "phylo") || length(phylo_tree$tip.label) < 2L ||
            anyNA(phylo_tree$tip.label) || any(trimws(phylo_tree$tip.label) == "") ||
            anyDuplicated(phylo_tree$tip.label) ||
            (!is.null(phylo_tree$edge.length) &&
             any(!is.finite(phylo_tree$edge.length) | phylo_tree$edge.length < 0))) {
            stop("Step 6 tree is structurally invalid: tip labels must be unique/non-empty and branch lengths finite/non-negative.")
        }

        cat("\nPhylogenetic tree loaded:\n",
            "  Tips:", length(phylo_tree$tip.label), "\n",
            "  Internal nodes:", phylo_tree$Nnode, "\n",
            "  Is rooted:", is.rooted(phylo_tree), "\n")

        # Check for tip label matching against the raw-counts table (a
        # representative abundance source)
        matching_tips <- intersect(colnames(abundance_tables[["raw_counts"]]), phylo_tree$tip.label)
        if (length(matching_tips) < ncol(abundance_tables[["raw_counts"]])) {
            warning("Only ", length(matching_tips), "/",
                    ncol(abundance_tables[["raw_counts"]]),
                    " raw ASVs have matching tree tips. Objects that include the tree will discard unmatched ASVs.")
        }
    } else {
        cat("Phylogenetic tree not found at:", tree_path, "\n")
        cat("Phyloseq objects will be created without a tree.\n")
        cat("Run Step 6 (Phylogenetic Tree) to generate the tree.\n")
    }
}
```

------------------------------------------------------------------------

# Construct Phyloseq Object(s)

## Build Phyloseq Objects

Every combination of taxonomy database x available abundance source gets
its own phyloseq object, keyed by a `combo_id` of the form
`<db_name>_<abundance_source>` (e.g. `SILVA_raw_counts`,
`SILVA_copy_number_corrected`, `GTDB_microbial_load_corrected`). A
parallel `phyloseq_combo_meta` list stores each combo’s display metadata
(taxonomy database, full abundance-source name, and its Excel-safe short
label) so every downstream section can look up how to label its output
without recomputing anything.

``` r
# Initialize lists to store phyloseq objects and their associated metadata,
# both keyed by combo_id ("<db_name>_<abundance_source>")
phyloseq_objects <- list()
phyloseq_combo_meta <- list()

step9_signature <- list(
    abundance = vapply(abundance_source_definitions[names(abundance_tables)],
                       function(x) unname(tools::md5sum(x$path)), character(1)),
    taxonomy = vapply(taxonomy_files, function(x) unname(tools::md5sum(x)), character(1)),
    metadata = if (is.null(metadata_path)) NA_character_ else unname(tools::md5sum(metadata_path)),
    tree = if (is.null(phylo_tree)) NA_character_ else unname(tools::md5sum(tree_path)),
    include_tree = include_tree,
    phyloseq_version = as.character(packageVersion("phyloseq"))
)
step9_checkpoint_path <- here(results_folder, ".9_phyloseq_object_checkpoint.rds")
checkpoint <- if (file_exists(step9_checkpoint_path)) {
    tryCatch(readRDS(step9_checkpoint_path), error = function(e) NULL)
} else NULL
checkpoint_loaded <- !is.null(checkpoint) && identical(checkpoint$signature, step9_signature) &&
    is.list(checkpoint$phyloseq_objects) && is.list(checkpoint$phyloseq_combo_meta)

if (checkpoint_loaded) {
    phyloseq_objects <- checkpoint$phyloseq_objects
    phyloseq_combo_meta <- checkpoint$phyloseq_combo_meta
    message("Resuming Step 9 object construction from a compatible checkpoint.")
} else {
    if (file_exists(step9_checkpoint_path)) file_delete(step9_checkpoint_path)

for (db_name in names(taxonomy_tables)) {

    # Get taxonomy table for this database
    tax_df <- taxonomy_tables[[db_name]]

    # Remove non-taxonomic columns if present (e.g., Unique_Tax)
    tax_columns <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    available_tax_columns <- intersect(tax_columns, colnames(tax_df))
    tax_matrix <- as.matrix(tax_df[, available_tax_columns, drop = FALSE])

    for (source_key in names(abundance_tables)) {

        combo_id <- paste0(db_name, "_", source_key)
        source_def <- abundance_source_definitions[[source_key]]

        cat("\nBuilding phyloseq object:", db_name, "taxonomy x", source_def$display, "\n\n")

        otu_matrix <- abundance_tables[[source_key]]

        # Ensure ASV order matches between OTU and taxonomy tables
        common_asvs <- intersect(colnames(otu_matrix), rownames(tax_matrix))

        if (length(common_asvs) == 0) {
            warning("No matching ASVs between ", source_def$display, " counts and ", db_name, " taxonomy! Skipping ", combo_id, ".")
            next
        }

        cat("  Matching ASVs:", length(common_asvs), "\n")
        taxonomy_loss <- setdiff(colnames(otu_matrix), common_asvs)
        if (length(taxonomy_loss)) {
            warning(combo_id, " will discard ", length(taxonomy_loss), "/",
                    ncol(otu_matrix), " ASVs without matching taxonomy (",
                    round(100 * length(taxonomy_loss) / ncol(otu_matrix), 1), "%).")
        }

        # Subset to common ASVs
        otu_subset <- otu_matrix[, common_asvs, drop = FALSE]
        tax_subset <- tax_matrix[common_asvs, , drop = FALSE]

        # Build phyloseq components
        OTU <- otu_table(otu_subset, taxa_are_rows = FALSE)
        TAX <- tax_table(tax_subset)

        # Add sample data if available. It was aligned to the complete raw
        # abundance sample set above, so this subset preserves every sample in
        # the current abundance source, including NA-filled metadata rows.
        if (!is.null(sample_metadata)) {
            SAM <- sample_data(sample_metadata[rownames(otu_subset), , drop = FALSE])
        } else {
            SAM <- NULL
        }

        # Add phylogenetic tree if available
        if (!is.null(phylo_tree)) {
            # Prune tree to matching tips
            matching_tips <- intersect(common_asvs, phylo_tree$tip.label)
            if (length(matching_tips) > 0) {
                tree_loss <- setdiff(common_asvs, matching_tips)
                if (length(tree_loss)) {
                    warning(combo_id, " will discard ", length(tree_loss), "/",
                            length(common_asvs), " taxonomy-matched ASVs without tree tips (",
                            round(100 * length(tree_loss) / length(common_asvs), 1), "%).")
                }
                TREE <- prune_taxa(matching_tips, phylo_tree)
                # Also subset OTU and TAX to matching tips
                OTU <- otu_table(otu_subset[, matching_tips, drop = FALSE], taxa_are_rows = FALSE)
                TAX <- tax_table(tax_subset[matching_tips, , drop = FALSE])
                cat("  Tree tips included:", length(matching_tips), "\n")
            } else {
                TREE <- NULL
                cat("  No matching tree tips found\n")
            }
        } else {
            TREE <- NULL
        }

        final_otu_matrix <- as(OTU, "matrix")
        if (taxa_are_rows(OTU)) final_otu_matrix <- t(final_otu_matrix)
        zero_samples <- rownames(final_otu_matrix)[rowSums(final_otu_matrix) <= 0]
        if (length(zero_samples)) {
            stop(combo_id, " has samples with zero abundance after taxonomy/tree pruning: ",
                 paste(zero_samples, collapse = ", "))
        }

        # Create the phyloseq object
        if (!is.null(SAM) && !is.null(TREE)) {
            ps <- phyloseq(OTU, TAX, SAM, TREE)
        } else if (!is.null(SAM)) {
            ps <- phyloseq(OTU, TAX, SAM)
        } else if (!is.null(TREE)) {
            ps <- phyloseq(OTU, TAX, TREE)
        } else {
            ps <- phyloseq(OTU, TAX)
        }

        # Store the phyloseq object and its display metadata
        phyloseq_objects[[combo_id]] <- ps
        phyloseq_combo_meta[[combo_id]] <- list(
            db_name           = db_name,
            abundance_source  = source_key,
            source_display    = source_def$display,
            source_label      = source_def$label
        )
        phyloseq_combo_meta[[combo_id]]$input_asvs <- ncol(otu_matrix)
        phyloseq_combo_meta[[combo_id]]$taxonomy_matched_asvs <- length(common_asvs)
        phyloseq_combo_meta[[combo_id]]$final_asvs <- ntaxa(ps)

        cat("\n  Phyloseq object created successfully! (", combo_id, ")\n")
    }
}

    temporary_checkpoint <- paste0(step9_checkpoint_path, ".tmp")
    saveRDS(list(signature = step9_signature,
                 phyloseq_objects = phyloseq_objects,
                 phyloseq_combo_meta = phyloseq_combo_meta), temporary_checkpoint)
    if (file_exists(step9_checkpoint_path)) file_delete(step9_checkpoint_path)
    file_move(temporary_checkpoint, step9_checkpoint_path)
}
```

## Explore Phyloseq Objects

For every phyloseq object built above, this chunk prints basic
composition statistics (taxa/sample counts, total reads) and, per
taxonomic rank, how many ASVs could not be classified – a quick sanity
check before moving on to aggregation and visualization.

``` r
for (combo_id in names(phyloseq_objects)) {

    ps <- phyloseq_objects[[combo_id]]
    meta <- phyloseq_combo_meta[[combo_id]]

    cat("\n", meta$db_name, "-", meta$source_display, "phyloseq summary:\n")

    # Basic statistics
    cat("  Taxonomic ranks:", paste(rank_names(ps), collapse = ", "), "\n")
    cat("  Number of taxa (ASVs):", ntaxa(ps), "\n")
    cat("  Number of samples:", nsamples(ps), "\n")
    cat("  Sample names:", paste(head(sample_names(ps), 5), collapse = ", "),
        ifelse(nsamples(ps) > 5, "...", ""), "\n")
    cat("  Total abundance:", format(sum(sample_sums(ps)), big.mark = ","), "\n")
    cat("  Has sample data:", !is.null(sample_data(ps, errorIfNULL = FALSE)), "\n")
    cat("  Has phylogenetic tree:", !is.null(phy_tree(ps, errorIfNULL = FALSE)), "\n")

    # Get a sense of how many NA values you have in each taxonomic rank:
    tax_na_counts <- apply(tax_table(ps), 2, function(x) sum(is.na(x)))
    cat("\n  NA counts by taxonomic rank:\n")
    for (rank in names(tax_na_counts)) {
        cat("    ", rank, ":", tax_na_counts[rank], "/", ntaxa(ps), "\n")
    }
}
```

## Save Phyloseq Objects

This chunk writes each phyloseq object built above to its own `.RData`
file, inside that combination’s database-specific `phyloseq_objects/`
subfolder
(e.g. `SILVA/phyloseq_objects/phyloseq_object_silva_raw_counts.RData`,
`GTDB/phyloseq_objects/phyloseq_object_gtdb_raw_counts.RData`), so each
one can be reloaded independently in downstream analysis scripts without
re-running this notebook.

``` r
# results='asis' + cat() (rather than print()) is required for
# render_output_links()'s Markdown output to actually render as a clickable
# link -- render_output_links() returns a knitr::asis_output() object (a
# character string tagged for raw-Markdown insertion), and that tagging is
# only honoured when knitr auto-prints a chunk's own top-level return value.
# An explicit print() call inside a for loop bypasses that auto-print path
# entirely and falls back to base R's print.default(), which dumps the
# object's raw character content AND its class/attribute metadata as
# literal visible text (e.g. `[1] "- [text](url)" attr(,"class") ...`)
# instead of rendering the link. cat() has no such special-casing: it just
# writes the object's character content straight to the chunk's output
# stream, which results='asis' then passes through as raw Markdown, exactly
# as intended.
for (combo_id in names(phyloseq_objects)) {

    ps <- phyloseq_objects[[combo_id]]
    meta <- phyloseq_combo_meta[[combo_id]]

    # Create filename based on database and abundance source, and save into
    # that database's own phyloseq_objects/ subfolder
    phyloseq_filename <- paste0("phyloseq_object_", tolower(meta$db_name), "_", meta$abundance_source, ".RData")
    phyloseq_path <- here(db_phyloseq_folder(meta$db_name), phyloseq_filename)

    # Save atomically within staging, then reopen and validate the object.
    temporary_phyloseq_path <- paste0(phyloseq_path, ".tmp")
    save(ps, file = temporary_phyloseq_path)
    verification_environment <- new.env(parent = emptyenv())
    loaded_names <- load(temporary_phyloseq_path, envir = verification_environment)
    if (!identical(loaded_names, "ps") ||
        !inherits(verification_environment$ps, "phyloseq") ||
        !identical(sample_names(verification_environment$ps), sample_names(ps)) ||
        !identical(taxa_names(verification_environment$ps), taxa_names(ps))) {
        stop("Saved phyloseq object failed verification: ", phyloseq_path)
    }
    if (file_exists(phyloseq_path)) file_delete(phyloseq_path)
    file_move(temporary_phyloseq_path, phyloseq_path)

    cat("\n**", meta$db_name, "-", meta$source_display, "**\n\n", sep = "")
    cat("Object staged and verified; final links are published after the complete run is validated.\n\n")
}
```

------------------------------------------------------------------------

# Taxonomic Aggregation

## Aggregate at Genus Level

For every phyloseq object built above, this chunk collapses ASVs sharing
the same Genus assignment into a single row using phyloseq’s
[tax_glom()](https://joey711.github.io/phyloseq/) (`NArm = FALSE`, so
unclassified ASVs are kept as their own group rather than dropped), then
also computes a genus-level relative-abundance version for later
visualization.

``` r
# Store aggregated objects, keyed by combo_id
genus_objects <- list()
genus_relative <- list()

for (combo_id in names(phyloseq_objects)) {

    ps <- phyloseq_objects[[combo_id]]
    meta <- phyloseq_combo_meta[[combo_id]]

    cat("\nAggregating", combo_id, "at Genus level...\n")

    # Aggregate the OTU table at the genus level
    ps_genus <- tax_glom(ps, taxrank = "Genus", NArm = FALSE)
    genus_objects[[combo_id]] <- ps_genus

    cat("  Original ASVs:", ntaxa(ps), "\n")
    cat("  Unique genera:", ntaxa(ps_genus), "\n")

    # Calculate relative abundances at the genus level
    if (any(sample_sums(ps_genus) <= 0)) {
        stop(combo_id, " has a zero-total sample after genus aggregation.")
    }
    ps_genus_rel <- transform_sample_counts(ps_genus, function(x) x / sum(x))
    genus_relative[[combo_id]] <- ps_genus_rel

    cat("  Relative abundance calculated\n")
}
```

## Identify Top Genera

This chunk ranks each combination’s genera by total abundance in that
source table and keeps the top `n_top_taxa` (set in
[Configuration](#define-data-params)), previewing them as an interactive
table and exporting the full ranking to a per-combination Excel sheet.
**Ranking basis**: “total abundance” here means each genus’s abundance
**summed across every sample** within that taxonomy database x abundance
source combination – one shared, dataset-wide top `n_top_taxa`, not a
separate top `n_top_taxa` chosen independently for each individual
sample. A genus that is highly abundant in only one compositionally
distinct sample (e.g. a control or an outlier) can therefore still fall
outside this shared ranking, even though it dominates that one sample’s
own composition – such a genus is simply absent from that combination’s
barplot and Excel sheet, rather than shown with a near-zero value.
Genera with an unclassified (`NA`) Genus assignment are excluded from
the ranking/selection itself – `tax_glom(..., NArm = FALSE)` in
[Aggregate at Genus Level](#aggregate-genus) deliberately collapses
every unclassified ASV into one row rather than dropping it, and that
row is often abundant enough to otherwise occupy a top-N slot that a
named genus should have. `Relative_Abundance` is still computed against
the *total* community (including unclassified reads), so the percentages
shown remain a true share of the whole sample rather than being inflated
by excluding the unclassified fraction from the denominator too.

``` r
top_genera_list <- list()
top_genera_df_list <- list()

# Interactive tables are accumulated here (one per combo_id) rather than
# printed inside the loop -- explicit print() on an htmlwidget inside a for
# loop bypasses knitr's dependency-aware knit_print() method, so the
# widget's JS/CSS dependencies never get registered and it renders blank.
# Auto-printing one combined tagList() as the chunk's own final top-level
# expression (below, after the loop) is what makes knit_print() fire
# correctly for every table.
report_tables <- list()

# Iterate over the genus-level objects built in the preceding section. The
# output lists are empty at this point and are populated inside this loop, so
# iterating over names(top_genera_list) here would execute zero times and leave
# both the top-genus selections and shared color palette empty.
for (combo_id in names(genus_objects)) {

    ps_genus <- genus_objects[[combo_id]]
    meta <- phyloseq_combo_meta[[combo_id]]

    # Sum abundance for each genus (raw reads, copy-number-corrected values,
    # or microbial-load values, depending on this object's source).
    genus_sums <- taxa_sums(ps_genus)

    # Total community abundance, including any unclassified (NA Genus) reads
    # -- kept as the Relative_Abundance denominator below so percentages
    # reflect a true share of the whole sample, even after the unclassified
    # genus itself is excluded from ranking/selection just below.
    total_abundance <- sum(genus_sums)

    # Exclude the unclassified (NA Genus) row, if present, before ranking --
    # see this section's intro for why. Sorting by abundance only among
    # classified genera means "top N" always means the top N *named* genera.
    classified_asv_ids <- rownames(tax_table(ps_genus))[!is.na(tax_table(ps_genus)[, "Genus"])]
    genus_sums_classified <- genus_sums[intersect(names(genus_sums), classified_asv_ids)]

    # Sort the classified genera by abundance
    sorted_genus <- sort(genus_sums_classified, decreasing = TRUE)

    if (!length(sorted_genus)) {
        warning(combo_id, " contains no classified genera; its top-genus table and barplot will be skipped.")
        next
    }

    # Select top wanted genera
    top_genus <- names(sorted_genus)[1:min(n_top_taxa, length(sorted_genus))]
    top_genera_list[[combo_id]] <- top_genus

    # Copy-number-corrected and microbial-load-corrected values are
    # fractional (unlike raw counts, which are already whole numbers of
    # reads), so Total_Abundance is rounded to 1 decimal place for those two
    # sources only -- Relative_Abundance below is still computed from the
    # unrounded values, so its percentages are unaffected by this rounding.
    genus_total_abundance <- as.numeric(sorted_genus[top_genus])
    if (meta$abundance_source != "raw_counts") {
        genus_total_abundance <- round(genus_total_abundance, 1)
    }

    # Create summary table
    top_genera_df <- data.frame(
        Rank = seq_along(top_genus),
        ASV_ID = top_genus,
        Genus = as.character(tax_table(ps_genus)[top_genus, "Genus"]),
        Total_Abundance = genus_total_abundance,
        Relative_Abundance = round(as.numeric(sorted_genus[top_genus]) / total_abundance * 100, 2),
        stringsAsFactors = FALSE
    )
    top_genera_df_list[[combo_id]] <- top_genera_df

    cat("\nTop", n_top_taxa, "genera (", meta$db_name, "-", meta$source_display, "):\n")

    report_tables[[combo_id]] <- datatable(
        top_genera_df,
        options = list(pageLength = 10, scrollX = TRUE,
                        scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
        caption = paste("Top", n_top_taxa, "genera by abundance -", meta$db_name, "-", meta$source_display)
    )

    # Export to Excel, into this combination's own database-specific workbook
    # (short source_label keeps the sheet name within Excel's 31-character
    # limit)
    add_sheet_to_excel(
        workbook_path = db_excel_path(meta$db_name),
        sheet_name = paste0("Top_Genera_", meta$source_label),
        data = top_genera_df,
        rownames = FALSE,
        overwrite = TRUE
    )
}

# Auto-print every accumulated table as one combined tagList() -- this is the
# chunk's own top-level return value, so knitr's evaluate::evaluate() invokes
# knit_print() on it correctly (see comment above the report_tables list).
htmltools::tagList(report_tables)
```

------------------------------------------------------------------------

# Visualization

## Define Barplot Function

This chunk defines a reusable helper, `create_taxonomy_barplot()`,
called once per taxonomy database x abundance source combination in the
section below – centralizing the plot styling and the interactive HTML
export in one place instead of repeating it at the call site.

``` r
# Create Interactive Taxonomy Barplot
#
# Generates an interactive barplot of taxonomic composition using ggplot2 and
# plotly. Color-to-taxon mapping is supplied externally via `color_palette`
# (built once, in Assign Consistent Genus Colors, from the union of every
# taxon appearing across ALL combinations processed in this run) rather than
# computed locally from this one call's own taxa -- that is what makes a
# given genus draw in the same color across every plot this notebook
# generates. Previously, each call derived its own palette from only its own
# subset of taxa, so the same genus could land on a different color in every
# combination's plot (and, coincidentally, the same color could be reused
# for two different genera across plots), making cross-plot color
# comparisons meaningless.
#
# Arguments:
#   ps            - A phyloseq object
#   tax_rank      - Taxonomic rank to display (default: "Genus")
#   title         - Plot title
#   y_title       - Y-axis title
#   legend_title  - Legend (color key) title
#   color_palette - Named character vector mapping taxon name -> hex color,
#                   shared across every call to this function in this run
#                   (see Assign Consistent Genus Colors). scale_fill_manual()
#                   matches palette entries to the data by NAME, not by
#                   position, so this stays correct even though each call's
#                   ps object only contains a SUBSET of the palette's full
#                   taxon universe -- unused names in the palette are
#                   silently ignored by ggplot2.
#   output_path   - Path to save HTML file (optional)
# Returns:
#   A plotly object
create_taxonomy_barplot <- function(ps, tax_rank = "Genus", title, y_title, legend_title, color_palette, output_path = NULL) {

    # Generate barplot. The externally-supplied, run-wide color_palette (see
    # Assign Consistent Genus Colors) -- not a palette derived from this
    # call's own taxa -- is what keeps a given genus' color identical across
    # every combination's plot.
    plot <- plot_bar(ps, fill = tax_rank) +
        scale_fill_manual(values = color_palette, name = legend_title) +
        ggtitle(title) +
        labs(x = "", y = y_title) +
        theme_minimal() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right",
            legend.text = element_text(size = 8)
        )

    # Convert ggplot object to plotly for interactivity. ggplotly() does not
    # reliably carry over a ggplot2 legend title set via scale_*(name = ...)
    # -- a well-documented ggplotly() conversion quirk -- so the legend title
    # is set a second time here, directly on the plotly object via layout(),
    # to guarantee it actually appears as intended in the rendered widget.
    interactive_plot <- ggplotly(plot) %>%
        layout(legend = list(title = list(text = legend_title)))

    # Export the interactive plot to an HTML file if path provided
    if (!is.null(output_path)) {
        saveWidget(interactive_plot, output_path, selfcontained = TRUE)
        cat("  Plot saved to:", output_path, "\n")
    }

    return(interactive_plot)
}
```

## Assign Consistent Genus Colors

So the same genus is always drawn in the same color in every barplot
below – rather than each combination’s plot independently re-deriving
its own color assignment from only its own subset of taxa – this chunk
builds ONE color palette shared across every call to
`create_taxonomy_barplot()`, keyed by genus name. It is built from the
union of every genus appearing in ANY combination’s top `n_top_taxa`
selection ([Identify Top Genera](#top-genera)), sorted alphabetically so
the color assigned to each genus is reproducible and does not depend on
which combination happened to be processed first.

``` r
# Union of every genus name appearing in any combination's top-N selection,
# across every taxonomy database x abundance source combination processed in
# this run -- the full set of genera that could ever need a color below.
# Sorting alphabetically (rather than leaving insertion/discovery order,
# which would depend on combo processing order) makes the color assigned to
# each genus reproducible from one run to the next.
all_top_genera <- sort(unique(unlist(lapply(top_genera_df_list, function(df) df$Genus))))
n_unique_genera_total <- length(all_top_genera)

cat("Building a shared color palette for", n_unique_genera_total,
    "unique genera across", length(top_genera_df_list), "combination(s)...\n")

# RColorBrewer's "Paired" palette is a well-distinguished 12-color
# qualitative palette; colorRampPalette() interpolates it up to however many
# unique genera actually need a color across the whole run, so every genus
# receives its own color regardless of how many total genera are involved.
# The resulting NAMED vector (genus -> hex color) is what
# create_taxonomy_barplot() below matches to each plot's data by name.
genus_color_palette <- setNames(
    if (n_unique_genera_total) colorRampPalette(brewer.pal(12, "Paired"))(n_unique_genera_total) else character(),
    all_top_genera
)
```

## Generate Genus-Level Barplot

For every combination, this chunk generates a single interactive barplot
– relative abundance restricted to the top `n_top_taxa` genera
(excluding any unclassified/NA genus – see [Identify Top
Genera](#top-genera)) – exporting it as a self-contained HTML widget and
displaying it inline. Every plot uses the same genus-to-color mapping
built in [Assign Consistent Genus Colors](#genus-colors), so a given
genus is always the same color no matter which combination’s plot it
appears in, making a direct visual, genus-by-genus comparison across
plots meaningful. Absolute-abundance, full (non-restricted)
relative-abundance, and any species-level barplot are intentionally not
generated: absolute counts mostly reflect sequencing depth rather than
biology in compositional 16S data, an unrestricted genus- or
species-level barplot becomes visually unreadable once the categorical
color palette runs out of distinguishable hues, and species-level
assignments from short-read 16S amplicons are typically too fragmented
to be worth plotting routinely (the underlying Genus- and Species-level
data are still fully available – for Genus, in the [Top
Genera](#top-genera) Excel sheet and `genus_objects`/`genus_relative`;
for Species, by running `tax_glom(ps, taxrank = "Species")` yourself on
the saved [phyloseq object](#save-phyloseq)). Barplots are written into
that combination’s database-specific `barplots/` subfolder
(e.g. `SILVA/barplots/`, `GTDB/barplots/`) under this notebook’s output
folder, so each database’s set of interactive files can be browsed or
shared independently.

``` r
# This chunk deliberately does NOT use the results='asis' + print(htmltools::
# tagList(widget)) pattern used elsewhere in this project (e.g.
# 5_dada2_pipeline.md's display-taxonomy-tables chunk) for widgets emitted
# from inside a loop. That pattern renders the widget's container <div> and
# data correctly, but does NOT actually register the widget's own JS/CSS
# library with the document -- verified empirically (rendered to real HTML
# and inspected in an actual browser): the plot area comes out as blank
# space, because base R's print() generic dispatches an htmlwidget-wrapped
# htmltools::tagList() to htmltools' own plain print.shiny.tag.list()
# /print.shiny.tag(), NOT to knitr's dependency-aware knit_print.shiny.tag.list().
# knit_print() -- the method that actually registers a widget's dependencies
# with the document via knitr::knit_meta_add() -- is only invoked
# automatically by knitr's own auto-print machinery, and only for a chunk's
# own top-level, non-nested return value; an explicit print() call buried
# inside a for loop never goes through that path, no matter the chunk's
# results option.
#
# The fix verified to actually work: accumulate every combination's heading,
# plot, and download link as htmltools tags in a single list across the
# whole loop, then auto-print ONE combined htmltools::tagList() as this
# chunk's own last top-level expression (see below, after the loop). That is
# a genuine top-level return value, so knitr's auto-print correctly routes
# it through knit_print.shiny.tag.list() and registers every widget's
# dependencies -- confirmed by rendering this exact pattern and checking
# window.Plotly was defined and the chart's <svg> actually drew in a real
# browser.
report_tags <- list()

# Plot only combinations for which the top-genus selection above succeeded.
# A combination with no classified genera is deliberately skipped there and
# therefore has no entry in top_genera_list.
for (combo_id in names(top_genera_list)) {

    meta <- phyloseq_combo_meta[[combo_id]]

    # Per-database barplots/ subfolder (created on first use by
    # db_barplot_folder(), defined in Define Path Parameters). Since the
    # database name is now encoded by the enclosing folders, it is no longer
    # repeated in the filename itself below -- only the abundance source is.
    barplot_db_folder <- db_barplot_folder(meta$db_name)

    combo_label <- meta$abundance_source
    barplot_path <- here(barplot_db_folder, paste0("barplot_genus_top", n_top_taxa, "_", combo_label, ".html"))

    # Barplot for top N genera (relative abundance). The plot title
    # intentionally omits "Top N Genera" and "Relative Abundance" -- that
    # information now lives in the legend title and y-axis label
    # respectively (see create_taxonomy_barplot()'s legend_title argument),
    # so the title itself only needs to identify this specific combination.
    ps_top_genus <- prune_taxa(top_genera_list[[combo_id]], genus_relative[[combo_id]])
    plot_top <- create_taxonomy_barplot(
        ps = ps_top_genus,
        tax_rank = "Genus",
        title = paste(meta$db_name, "-", meta$source_display),
        y_title = "Relative Abundance",
        legend_title = paste("Top", n_top_taxa, "Genera"),
        color_palette = genus_color_palette,
        output_path = barplot_path
    )
    if (!file_exists(barplot_path) || file_info(barplot_path)$size <= 0) {
        stop("Interactive barplot was not written successfully: ", barplot_path)
    }

    # Portable link text/href, matching render_output_links()'s own
    # convention (display path relative to the project root; href relative
    # to R/notebooks/, where this notebook's own rendered .html lives) --
    # built natively as htmltools tags here, rather than by calling
    # render_output_links() itself, since that function returns Markdown
    # text meant for results='asis' + cat(), not for embedding inside an
    # htmltools tag tree.
    final_barplot_path <- here(final_output_folder, meta$db_name, "barplots",
                               path_file(barplot_path))
    barplot_display_path <- as.character(path_rel(final_barplot_path, start = here()))
    barplot_href <- as.character(path_rel(final_barplot_path, start = here("R", "notebooks")))

    report_tags[[length(report_tags) + 1]] <- htmltools::h4(paste(meta$db_name, "-", meta$source_display))
    report_tags[[length(report_tags) + 1]] <- plot_top
    report_tags[[length(report_tags) + 1]] <- if (file_exists(barplot_path)) {
        htmltools::tags$ul(htmltools::tags$li(htmltools::tags$a(
            href = barplot_href,
            paste0("Genus barplot (top ", n_top_taxa, "): ", barplot_display_path)
        )))
    } else {
        htmltools::tags$p(htmltools::strong("Not found: "), barplot_display_path,
            " (expected but was not created -- check the chunk above for warnings or errors)")
    }
}
```

------------------------------------------------------------------------

# Sample Statistics

## Per-Sample Abundance

For every combination, this chunk tabulates each sample’s total
abundance in the units of its source table and its observed ASV
richness, previews the table sorted by total abundance, and exports the
full table to Excel. For raw counts the total is reads; for corrected
sources it is a fractional copy-number-corrected total or a
microbial-load-scaled total, rounded to 1 decimal place.

``` r
sample_stats_list <- list()

# Interactive tables are accumulated here (one per combo_id) rather than
# printed inside the loop -- explicit print() on an htmlwidget inside a for
# loop bypasses knitr's dependency-aware knit_print() method, so the
# widget's JS/CSS dependencies never get registered and it renders blank.
# Auto-printing one combined tagList() as the chunk's own final top-level
# expression (below, after the loop) is what makes knit_print() fire
# correctly for every table.
report_tables <- list()

for (combo_id in names(phyloseq_objects)) {

    ps <- phyloseq_objects[[combo_id]]
    meta <- phyloseq_combo_meta[[combo_id]]

    # Copy-number-corrected and microbial-load-corrected values are
    # fractional (unlike raw counts, which are already whole numbers of
    # reads), so Total_Abundance is rounded to 1 decimal place for those two
    # sources only when it goes into the table below. Percentage, sort
    # order, and the console summary further down are all computed from
    # these unrounded per-sample totals, so this display-only rounding does
    # not affect them.
    sample_total_abundance <- sample_sums(ps)
    if (meta$abundance_source != "raw_counts") {
        rounded_total_abundance <- round(sample_total_abundance, 1)
    } else {
        rounded_total_abundance <- sample_total_abundance
    }

    # Calculate per-sample statistics
    sample_stats <- data.frame(
        SampleID = sample_names(ps),
        Total_Abundance = rounded_total_abundance,
        Observed_ASVs = apply(otu_table(ps), 1, function(x) sum(x > 0)),
        stringsAsFactors = FALSE
    )

    # Add percentage of total
    sample_stats$Pct_Total_Abundance <- round(sample_total_abundance / sum(sample_total_abundance) * 100, 2)

    # Sort by total abundance
    sample_stats <- sample_stats[order(sample_total_abundance, decreasing = TRUE), ]
    sample_stats_list[[combo_id]] <- sample_stats

    cat("\nSample statistics (", meta$db_name, "-", meta$source_display, "):\n",
        "  Total samples:", nrow(sample_stats), "\n",
        "  Mean abundance/sample:", format(round(mean(sample_total_abundance), 2), big.mark = ","), "\n",
        "  Median abundance/sample:", format(round(median(sample_total_abundance), 2), big.mark = ","), "\n",
        "  Min abundance:", format(round(min(sample_total_abundance), 2), big.mark = ","), "\n",
        "  Max abundance:", format(round(max(sample_total_abundance), 2), big.mark = ","), "\n",
        "  Mean ASVs/sample:", round(mean(sample_stats$Observed_ASVs), 1), "\n")

    report_tables[[combo_id]] <- datatable(
        sample_stats,
        options = list(pageLength = 10, scrollX = TRUE,
                        scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
        caption = paste("Per-sample statistics -", meta$db_name, "-", meta$source_display)
    )

    # Export to Excel, into this combination's own database-specific workbook
    add_sheet_to_excel(
        workbook_path = db_excel_path(meta$db_name),
        sheet_name = paste0("Sample_Stats_", meta$source_label),
        data = sample_stats,
        rownames = FALSE,
        overwrite = TRUE
    )
}

# Auto-print every accumulated table as one combined tagList() -- this is the
# chunk's own top-level return value, so knitr's evaluate::evaluate() invokes
# knit_print() on it correctly (see comment above the report_tables list).
htmltools::tagList(report_tables)
```

## Taxonomy Summary

For every taxonomy database, this chunk counts how many distinct taxa
are observed at each taxonomic rank and what fraction of ASVs were
successfully classified at that rank, previews the result, and exports
it to Excel as a single `Tax_Summary` sheet. The summary describes the
ASVs retained in that database’s raw-count phyloseq object after any
taxonomy and tree matching, rather than claiming to describe ASVs that
were pruned during construction.

``` r
tax_summary_list <- list()

# Interactive tables are accumulated here (one per db_name) rather than
# printed inside the loop -- explicit print() on an htmlwidget inside a for
# loop bypasses knitr's dependency-aware knit_print() method, so the
# widget's JS/CSS dependencies never get registered and it renders blank.
# The combined tagList() is auto-printed at the very end of this chunk
# (after the worksheet-reordering loop below), which is what makes
# knit_print() fire correctly for every table -- see that loop's own
# comment for why it has to stay between the two.
report_tables <- list()

for (db_name in names(taxonomy_tables)) {

    # Raw counts is the representative object for this database. The summary
    # explicitly describes its retained ASVs after taxonomy/tree matching.
    representative_combo_id <- paste0(db_name, "_raw_counts")

    if (!representative_combo_id %in% names(phyloseq_objects)) {
        warning("No raw-counts phyloseq object found for ", db_name,
                " -- skipping its Tax_Summary sheet.")
        next
    }

    ps <- phyloseq_objects[[representative_combo_id]]

    # Count unique taxa at each rank
    tax_summary <- data.frame(
        Rank = rank_names(ps),
        stringsAsFactors = FALSE
    )

    tax_summary$Unique_Taxa <- sapply(rank_names(ps), function(rank) {
        length(unique(na.omit(tax_table(ps)[, rank])))
    })

    tax_summary$NA_Count <- sapply(rank_names(ps), function(rank) {
        sum(is.na(tax_table(ps)[, rank]))
    })

    tax_summary$Classified_Pct <- round(
        (ntaxa(ps) - tax_summary$NA_Count) / ntaxa(ps) * 100, 1
    )
    tax_summary_list[[db_name]] <- tax_summary

    cat("\nTaxonomy summary (", db_name, ", from its raw-counts phyloseq object):\n")

    report_tables[[db_name]] <- datatable(
        tax_summary,
        options = list(pageLength = 10, dom = 't',
                        scrollY = "400px", scrollCollapse = TRUE, paging = FALSE),
        caption = paste("Taxonomy classification summary -", db_name)
    )

    # Export to Excel, into this database's own summary workbook -- ONE
    # "Tax_Summary" sheet per database (not one per abundance source).
    add_sheet_to_excel(
        workbook_path = db_excel_path(db_name),
        sheet_name = "Tax_Summary",
        data = tax_summary,
        rownames = FALSE,
        overwrite = TRUE
    )
}

# Move "Tax_Summary" to be the first sheet (leftmost tab) in each database's
# workbook. It is physically written after every Top_Genera_* and
# Sample_Stats_* sheet above, but worksheetOrder() controls Excel's TAB
# order independently of the physical order sheets were added in, so
# reordering here does not disturb any sheet's contents. This ordering also
# survives every later loadWorkbook()/saveWorkbook() cycle in this notebook
# (e.g. when the header-cell popup comments and the trailing
# Column_Dictionary sheet are added below), so sheets written after this
# point still land after Tax_Summary rather than before it.
#
# This loop -- not report_tables -- has to be the second-to-last statement
# in the chunk: only a chunk's own final top-level expression gets
# auto-printed, so the combined tagList() below has to come after all of
# the Excel bookkeeping, not before it.
for (db_name in names(tax_summary_list)) {

    workbook_path <- db_excel_path(db_name)
    sheet_names_in_workbook <- getSheetNames(workbook_path)
    tax_summary_position <- which(sheet_names_in_workbook == "Tax_Summary")

    if (length(tax_summary_position) == 1L) {
        wb <- loadWorkbook(workbook_path)
        worksheetOrder(wb) <- c(
            tax_summary_position,
            setdiff(seq_along(sheet_names_in_workbook), tax_summary_position)
        )
        saveWorkbook(wb, workbook_path, overwrite = TRUE)
    }
}

# Auto-print every accumulated table as one combined tagList() -- this is the
# chunk's own top-level return value, so knitr's evaluate::evaluate() invokes
# knit_print() on it correctly (see comment above the report_tables list).
htmltools::tagList(report_tables)
```

## Record Run Provenance

Each database workbook records the exact inputs, checksums,
configuration, software versions, and ASV retention for every abundance
object produced in this run.

## Document and Export Column Dictionary

Every `Top_Genera_<source_label>` and `Sample_Stats_<source_label>`
sheet (one per abundance source processed) shares an identical column
layout within a database’s own workbook — only the underlying values
differ; `Tax_Summary` is written once per database and needs no such
per-source repetition. Build one `descriptions` lookup per sheet family
and reuse it for every abundance source’s copy of that sheet, then
combine everything into a trailing `Column_Dictionary` sheet, appended
separately to *each* taxonomy database’s own summary workbook, so every
workbook documents itself independently.

------------------------------------------------------------------------

# Output File Summary

The tree below lists every file this notebook has written to its own
output folder,
[results/9_phyloseq_object/](../../results/9_phyloseq_object/), as a
clickable, portable link (relative to this notebook’s own location) –
built live from what is actually on disk at knit time via the project’s
shared
[render_output_tree_function.R](../functions/render_output_tree_function.R)
helper. The exact set of `.RData` phyloseq objects, interactive `.html`
barplots, and per-database Excel workbooks depends on which taxonomy
database(s) and optional upstream steps (6-8) were used in this run, so
most filenames here are dynamically generated rather than fixed – only
the database subfolders (top-level, `phyloseq_objects/`, and
`barplots/`) and each database’s summary workbook get an inline
description; every other file’s name is already self-descriptive
(e.g. `phyloseq_object_silva_raw_counts.RData`,
`barplot_genus_top10_raw_counts.html`).

------------------------------------------------------------------------

# Interactive Exploration with Shiny-Phyloseq

<div class="alert alert-info">

**Optional**: For interactive exploration of a phyloseq object, you can
use the Shiny-Phyloseq web application. Uncomment and run the code below
to launch it.

</div>

``` r
# Launch Shiny-Phyloseq for interactive data exploration
# This will open a web browser with the interactive application
shiny::runGitHub("shiny-phyloseq", "joey711")
```

------------------------------------------------------------------------

# Recommended Next Step

This is the last notebook of the DADA2 16S pipeline itself, so there is
no further DADA2 notebook to run – the phyloseq object(s) exported above
are analysis-ready for whatever community-ecology workflow you use next
(alpha/beta diversity, ordination, differential abundance testing,
etc.). If your next question is about the community’s *functional*
potential rather than its taxonomic composition, this pipeline’s natural
companion project is the separate
**PICRUSt2_16S_Functional_Inference_Workflow** repository: it takes this
pipeline’s ASV table and representative sequences as input and uses
[PICRUSt2](https://github.com/picrust/picrust2/wiki) to predict each
sample’s functional gene content (e.g. KEGG pathways, EC numbers) via
phylogenetic placement and hidden-state prediction. It is a separate,
standalone repository rather than part of this one, so it is named here
rather than linked.

------------------------------------------------------------------------

# Session Information

Record the R environment for reproducibility.

------------------------------------------------------------------------

# References

## Methods

- McMurdie PJ, Holmes S (2013). phyloseq: An R Package for Reproducible
  Interactive Analysis and Graphics of Microbiome Census Data. *PLoS
  ONE*, 8(4):e61217. <https://doi.org/10.1371/journal.pone.0061217>
- [phyloseq: Explore microbiome profiles using
  R](https://joey711.github.io/phyloseq/index.html) — official
  documentation for the package this entire notebook is built on.
- [phyloseq
  tutorials](https://joey711.github.io/phyloseq/tutorials-index.html)
- [Shiny-phyloseq: An interactive web
  application](https://joey711.github.io/shiny-phyloseq/) — the optional
  interactive explorer launched in [Interactive Exploration with
  Shiny-Phyloseq](#interactive-exploration-with-shiny-phyloseq) above.

## Databases

- Quast C, Pruesse E, Yilmaz P, et al. (2013). The SILVA ribosomal RNA
  gene database project: improved data processing and web-based tools.
  *Nucleic Acids Research*, 41(D1):D590-D596.
  <https://doi.org/10.1093/nar/gks1219> —
  [SILVA](https://www.arb-silva.de/) is one of the two taxonomy
  databases this notebook can build a phyloseq object from.
- Parks DH, Chuvochina M, Rinke C, et al. (2022). GTDB: an ongoing
  census of bacterial and archaeal diversity through a phylogenetically
  consistent, rank normalized and complete genome-based taxonomy.
  *Nucleic Acids Research*, 50(D1):D199-D207.
  <https://doi.org/10.1093/nar/gkab776> —
  [GTDB](https://gtdb.ecogenomic.org/) is this notebook’s other
  supported taxonomy database.

## Related

- [Step 5 — DADA2 Pipeline](5_dada2_pipeline.md) — this notebook’s
  required input (ASV count table and taxonomy assignments).
- [Step 6 — Phylogenetic Tree](6_phylogenetic_tree.md) — optional input,
  incorporated into the phyloseq object(s) if present.
- [Step 7 — 16S rRNA Gene Copy Number
  Correction](7_copy_number_correction.md) — optional input, adds a
  copy-number-corrected phyloseq object per taxonomy database.
- [Step 8 — Microbial Load Correction](8_microbial_load_correction.md) —
  optional input, adds a microbial-load-corrected phyloseq object per
  taxonomy database.

------------------------------------------------------------------------

# Appendix: Troubleshooting Guide

## Common Issues and Solutions

### Fewer Phyloseq Objects Than Expected

**Symptom**: The [Build Phyloseq Objects](#build-phyloseq) section
reports fewer combinations than `number of taxonomy databases x 3`.

**Cause**: `copy_number_corrected` and/or `microbial_load_corrected`
abundance tables were not found, because [Step
7](7_copy_number_correction.md) and/or [Step
8](8_microbial_load_correction.md) have not been run. This is expected –
both are optional – and is reported as a skip message, not an error, in
[Load Abundance Tables](#load-abundance).

**Actions**:

- Run [Step 7](7_copy_number_correction.md) and/or [Step
  8](8_microbial_load_correction.md) first if you want the additional
  phyloseq object(s), then re-run this notebook.
- If you only need raw-counts phyloseq objects, no action is needed.

### ASV Names Don’t Match

**Error**: No matching ASVs between OTU table and taxonomy

**Cause**:

- ASV count table uses sample names as row names instead of column
  names.
- Taxonomy table has different ASV IDs.

**Actions**:

- Verify the format of your input files.
- Check that ASV IDs are consistent between files.
- Ensure ASV count table has samples as rows and ASVs as columns.

### Phylogenetic Tree Tips Don’t Match

**Warning**: Not all ASVs have matching tree tips

**Cause**:

- Tree was built with different ASV sequences.
- ASV IDs changed between steps.

**Actions**:

- Rebuild the tree using the same ASV sequences.
- Verify ASV ID consistency across all input files.

### Memory Issues with Large Datasets

**Error**: Cannot allocate memory

**Actions**:

- Process one taxonomy database at a time: temporarily rename or move
  the other database’s taxonomy CSV
  ([silva_taxonomy_table.csv](../../results/5_dada2_pipeline/silva_taxonomy_table.csv)
  or
  [gtdb_taxonomy_table.csv](../../results/5_dada2_pipeline/gtdb_taxonomy_table.csv),
  in [results/5_dada2_pipeline/](../../results/5_dada2_pipeline/)) out
  of that folder before running this notebook, so [Load Taxonomy
  Table(s)](#load-taxonomy) only detects the one you want – restore the
  file afterward to include it on a later run.
- Filter low-abundance ASVs before creating phyloseq object.
- Increase system memory.

### Missing Taxonomy Ranks

**Warning**: NA values in taxonomy table

**Cause**: This is common and expected. Many ASVs cannot be classified
to all taxonomic levels.

**Actions**:

- Use `NArm = FALSE` in `tax_glom()` to retain unclassified taxa.
- Consider using less stringent classification thresholds in [Step
  5](5_dada2_pipeline.md).
