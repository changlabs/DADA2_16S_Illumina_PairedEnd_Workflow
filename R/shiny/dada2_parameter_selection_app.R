# =============================================================================
# DADA2 PARAMETER EXPLORER
# Step 4 of the paired-end Illumina 16S rRNA sequencing workflow
# =============================================================================
#
# Run after primer trimming (Step 3) and before the DADA2 pipeline (Step 5).
# The app follows one four-tab workflow:
#
#   Visualize > Select > Validate > Export
#
# Visualize documents the amplified target, primer coordinates, and expected
# primer-trimmed length range for the assay that was actually used. It is not an
# experimental-design, primer-recommendation, or specificity-prediction tool.
#
# Select loads primer-trimmed paired FASTQs and helps the user choose forward
# and reverse truncLen/maxEE values from observed quality, estimated retention,
# and overlap at the maximum expected target length. truncQ is fixed at the
# DADA2 default of 2.
#
# Validate optionally runs the real DADA2 workflow on representative samples
# and compares predicted merged-read retention with observed merging.
#
# Export reviews the current assay, filtering parameters, estimates, and run
# metadata, then writes:
#   results/4_dada2_parameter_selection/dada2_filter_parameters.xlsx
# Step 5 imports the Parameters sheet from this workbook automatically.
#
# The current user guide is:
#   R/notebooks/4_dada2_parameter_selection.Rmd
#
# AUTHOR: Amro Abbas
# FILENAME: dada2_parameter_selection_app.R
# =============================================================================


# =============================================================================
# SECTION 1: LIBRARY IMPORTS
# =============================================================================
# Each library serves a specific purpose. We load all dependencies upfront
# to ensure they're available when needed and to fail fast if missing.

library(shiny)          # Core web framework - provides reactive programming model
                        # for building interactive web applications in R

library(shinyFiles)     # Cross-platform file/directory browser dialogs
                        # Allows users to select local directories for FASTQ files

library(shinyjs)        # JavaScript utilities for Shiny applications
                        # Used for showing/hiding UI elements, custom JS handlers

library(ShortRead)      # Bioconductor package for FASTQ file manipulation
                        # Provides FastqStreamer for memory-efficient reading

library(plotly)         # Interactive JavaScript plotting via R
                        # Enables zoom, pan, hover tooltips on all plots

library(dplyr)          # Data manipulation verbs (filter, mutate, summarize)
                        # The core of our data processing pipelines

library(DT)             # Interactive DataTables for R
                        # Provides sortable, filterable tables with column formatting

library(bslib)          # Bootstrap 5 theming for Shiny
                        # Modern, responsive UI components and layout

library(future)         # Unified parallel processing framework
                        # Manages worker processes for parallel file processing

library(htmltools)      # HTML generation utilities
                        # Used for building custom UI elements

library(openxlsx)       # Excel workbook creation/writing
                        # Writes the workbook saved from Export

# Shared workbook helpers keep the Export workbook consistent with the other
# Excel outputs in this pipeline. Paths are relative to R/shiny/, the app's
# working directory when launched with shiny::runApp().
source(file.path("..", "functions", "add_sheet_to_excel_function.R"), local = TRUE)
source(file.path("..", "functions", "build_column_dictionary_function.R"), local = TRUE)

# Paired FASTQ sampling, paired-retention calculations, representative-sample
# selection, and empirical DADA2 validation used by Select and Validate.
source(file.path("functions", "paired_read_retention_engine_function.R"), local = TRUE)


# =============================================================================
# SECTION 2: PARALLEL PROCESSING CONFIGURATION
# =============================================================================
# Processing FASTQ files is both I/O-intensive (reading compressed files)
# and CPU-intensive (quality calculations). Parallel processing dramatically
# reduces total time for datasets with many samples.

# Detect available CPU cores and reserve one for the UI thread
# This prevents the application from becoming unresponsive during processing
detected_cores <- suppressWarnings(parallel::detectCores())
if (length(detected_cores) != 1L || is.na(detected_cores) ||
    !is.finite(detected_cores) || detected_cores < 1) {
    detected_cores <- 2L
}
n_cores <- min(4L, max(1L, as.integer(detected_cores) - 1L))
# - detectCores(): Returns number of logical cores (may include hyperthreads)
# - We subtract 1 to leave resources for Shiny UI and operating system
# - min(4, ...) avoids spawning an excessive number of persistent R processes
# - max(1, ...) ensures we always have at least one worker

# The plan is activated only while the application is running (Section 9), so
# sourcing this script does not alter the caller's global future configuration.


# =============================================================================
# SECTION 3: CONFIGURATION CONSTANTS
# =============================================================================
# These constants define file discovery, default locations, plotting, overlap
# classification, and reproducible paired-read sampling.

# -----------------------------------------------------------------------------
# 3.1 File Pattern Choices
# -----------------------------------------------------------------------------
# Paired-end sequencing produces two files per sample (forward and reverse).
# Different sequencing facilities use different naming conventions.
# Users select the pattern matching their files.

PATTERN_CHOICES <- list(
    # Pattern 1: Full Illumina naming with lane info
    # Example: Sample1_S1_L001_R1_001.fastq.gz
    "_L001_R1_001 / _L001_R2_001" = list(fwd = "_L001_R1_001", rev = "_L001_R2_001"),
    
    # Pattern 2: Simple R1/R2 suffix (very common)
    # Example: Sample1_R1.fastq.gz, Sample1_R2.fastq.gz
    "_R1 / _R2" = list(fwd = "_R1", rev = "_R2"),
    
    # Pattern 3: Numeric suffix (SRA/ENA convention)
    # Example: Sample1_1.fastq.gz, Sample1_2.fastq.gz
    "_1 / _2" = list(fwd = "_1", rev = "_2")
)

# -----------------------------------------------------------------------------
# 3.2 Filtering Constants
# -----------------------------------------------------------------------------
# truncQ is fixed at the DADA2 default (2) and is not user-tunable. DADA2's
# filterAndTrim() applies truncQ by truncating each read at the first base
# with Q <= truncQ), so the parameter cannot be genuinely eliminated from the
# filtering math; instead it is pinned to the engine's
# RETENTION_DEFAULT_TRUNCQ constant everywhere it is needed.

# -----------------------------------------------------------------------------
# Direction Identity Colors
# -----------------------------------------------------------------------------
# Reused in card headers, plots, summaries, and table bars.
FWD_COLOR        <- "#1B9E77"              # teal-green   (Forward)
REV_COLOR        <- "#D95F02"              # orange       (Reverse)
PAIRED_COLOR     <- "#7570B3"              # muted purple (Paired / combined)
# Default retention-vs-maxEE curve color.
MAXEE_LINE_COLOR <- "#6f42c1"

# -----------------------------------------------------------------------------
# Overlap-Quality Color Scale
# -----------------------------------------------------------------------------
# Shared by the overlap bar and status message. Adjacent tiers use related hues
# while GOOD and FAIL use darker variants than MODERATE and LOW YIELD.
OVERLAP_STATUS_COLORS <- c(
    "overlap-good"     = "#2980b9",   # darker blue -- "Belize Hole"
    "overlap-moderate" = "#3498db",   # theme info (blue) -- "Peter River"
    "overlap-critical" = "#f39c12",   # theme warning (orange)
    "overlap-lowyield" = "#e74c3c",   # theme danger (red) -- "alizarin"
    "overlap-fail"     = "#c0392b"    # darker red -- "pomegranate"
)

# Darker variant of each color above, for anything that must be read as TEXT
# or a thin stroke on a light/white background: the overlap box's border,
# its "Overlap X bp" in-plot annotation font, and the "Overlap check"
# message's body text. The raw theme colors above (success teal, warning
# orange especially) fall well short of WCAG AA's 4.5:1 text-contrast
# minimum against white (~2.2-2.4:1 as measured); these variants are each
# darkened just enough to clear 4.5:1 while staying visibly the same hue, so
# the box/caption stay legible instead of washing out. "Fail" needs no
# separate darkening -- "#c0392b" already measures ~5.4:1 on white. "Good"'s
# own base color ("#2980b9") only measures ~4.3:1, just short of 4.5:1, so it
# still needs a slightly darker text variant like every other tier.
OVERLAP_STATUS_TEXT_COLORS <- c(
    "overlap-good"     = "#2471a3",
    "overlap-moderate" = "#2a7aaf",
    "overlap-critical" = "#a0670c",
    "overlap-lowyield" = "#d04436",
    "overlap-fail"     = "#c0392b"
)

# Hex color for a classify_overlap_status() class string (falls back to the
# danger red for any unrecognized class, defensively) -- used for translucent
# fills/tints only; see overlap_status_text_color() for stroke/font/text.
overlap_status_color <- function(status_class) {
    col <- OVERLAP_STATUS_COLORS[[status_class]]
    if (is.null(col)) "#e74c3c" else col
}

# WCAG-AA-legible darker variant, for strokes/annotation fonts/body text.
overlap_status_text_color <- function(status_class) {
    col <- OVERLAP_STATUS_TEXT_COLORS[[status_class]]
    if (is.null(col)) "#c0392b" else col
}

# overlap_status_color()'s result as an "r, g, b" string (via
# grDevices::col2rgb(), so this can never drift out of sync with
# OVERLAP_STATUS_COLORS above), for building rgba() fills/tints at whatever
# alpha the caller needs.
overlap_status_rgb <- function(status_class) {
    paste(grDevices::col2rgb(overlap_status_color(status_class))[, 1], collapse = ", ")
}

# -----------------------------------------------------------------------------
# Quality-Plot / Slider Alignment Margins
# -----------------------------------------------------------------------------
# Single source of truth for the plotly plot margins used by the read-quality
# plots AND the horizontal padding of the truncLen slider that visually acts as
# their x-axis. The SAME numbers are interpolated into both the plotly
# layout(margin = ...) call (make_quality_position_plot) and the .axis-slider-
# chart CSS variables (--plot-left / --plot-right, tags$style below), so the
# slider track and the plot's data region begin/end at the same pixel and can
# never drift apart -- across window resizes, sidebar toggles, or narrow screens.
# The slider sits BELOW each plot as a row: [label][-][ track ][+]. The label +
# minus button occupy a fixed-width LEFT block, the plus button a fixed-width
# RIGHT block, and the slider track fills the middle. For the track to line up
# with the plot's data region, the plotly left/right margins are set to those
# same block widths PLUS the ionRangeSlider handle inset (the draggable track sits
# a few px inside its own container on each end). Deriving the plot margins from
# the control-block widths + inset -- in one place -- keeps the slider handle and
# the plot's dashed line at the same pixel across resizes. QP_IRS_INSET is the
# single number to nudge if the line still sits slightly off the handle.
QP_CTRL_LEFT  <- 84   # px: label + minus-button block (left of the slider track)
QP_CTRL_RIGHT <- 30   # px: plus-button block (right of the slider track)
QP_IRS_INSET  <- 9    # px: ionRangeSlider handle half-width inset (tune if off)
QP_MARGIN_L   <- QP_CTRL_LEFT  + QP_IRS_INSET   # plotly left margin  (= track start px)
QP_MARGIN_R   <- QP_CTRL_RIGHT + QP_IRS_INSET   # plotly right margin (= track end px)
QP_MARGIN_T   <- 10   # px: top margin
QP_MARGIN_B   <- 2    # px: bottom margin (tight -- no x tick labels, and the slider
                      #     is pulled up close underneath via .axis-slider-row)
TRUNC_LINE_COLOR <- "#e74c3c"   # dashed truncLen line drawn through the quality plot

# -----------------------------------------------------------------------------
# 3.3 Paired-Read Sampling Reproducibility
# -----------------------------------------------------------------------------
# The paired FASTQ sampler uses this base seed plus the sample index. The value
# is recorded in Export so the quality-profile sampling is reproducible.
SAMPLING_BASE_SEED <- 20260806L

# -----------------------------------------------------------------------------
# 3.4 Default FASTQ Directory
# -----------------------------------------------------------------------------
# Pre-populates the directory chooser with the pipeline's standard
# primer-trimmed-reads location so the app is immediately usable on load,
# without requiring a manual Browse click. The Browse button remains
# available to point the app at a different location if needed (e.g. a
# different sequencing run, or a project laid out differently from the
# repository template).
#
# shiny::runApp() sets the working directory to the folder containing this
# script (R/shiny/) for the duration of the app session, so the pipeline's
# results/ folder is reached by going up two levels to the repository root.
# If the app is launched from a different working directory, or the
# expected folder does not exist yet (e.g. before Step 3 of the pipeline
# has been run), DEFAULT_FASTQ_DIR safely resolves to NULL and the app
# falls back to its original behaviour of requiring the user to Browse
# manually.

CANDIDATE_DEFAULT_FASTQ_DIR <- file.path(
    "..", "..", "results", "3_cutadapt_primer_trimming", "primer_trimmed_reads"
)

DEFAULT_FASTQ_DIR <- if (dir.exists(CANDIDATE_DEFAULT_FASTQ_DIR)) {
    normalizePath(CANDIDATE_DEFAULT_FASTQ_DIR, mustWork = FALSE)
} else {
    NULL   # Folder not found yet -- fall back to requiring manual Browse
}

# -----------------------------------------------------------------------------
# 3.5 Project Root Directory + Path Relativizer
# -----------------------------------------------------------------------------
# Resolve the project root once so Export can record project-relative input and
# output locations whenever possible.
PROJECT_ROOT_DIR <- normalizePath(file.path("..", ".."), mustWork = FALSE)

# Strips PROJECT_ROOT_DIR off the front of an absolute path, returning a
# forward-slash, project-root-relative path (e.g.
# "results/3_cutadapt_primer_trimming/primer_trimmed_reads"). Falls back
# to returning `path` unchanged if it isn't actually under the project
# root -- e.g. the user Browsed to a FASTQ directory living somewhere
# else entirely -- since there is no meaningful project-relative form to
# show in that case, and silently truncating an unrelated path would be
# misleading rather than clarifying.
to_project_relative_path <- function(path) {
    if (is.null(path) || !nzchar(path)) return(path)

    normalized <- normalizePath(path, mustWork = FALSE)
    project_prefix <- paste0(PROJECT_ROOT_DIR, .Platform$file.sep)
    if (identical(normalized, PROJECT_ROOT_DIR) || startsWith(normalized, project_prefix)) {
        # substring()/nchar() rather than sub() with PROJECT_ROOT_DIR
        # interpolated into a regex pattern -- a real filesystem path can
        # contain regex metacharacters (e.g. a folder literally named
        # "results (v2)"), which would silently corrupt a regex-based
        # strip; plain character-position slicing has no such risk.
        relative <- substring(normalized, nchar(PROJECT_ROOT_DIR) + 1)
        relative <- sub("^[/\\\\]+", "", relative)   # drop the leading slash left behind
        gsub("\\\\", "/", relative)                  # Windows backslashes -> forward slashes
    } else {
        path
    }
}

# -----------------------------------------------------------------------------
# 3.6 Sequencing Platform Read Lengths
# -----------------------------------------------------------------------------
# Nominal maximum per-read cycle count for common Illumina paired-end kits.
# input$vis_platform is selected in Select and caps the truncLen controls at the
# stated platform length. Defensive downstream fallbacks handle a missing or
# unmatched selection.

PLATFORM_READ_LENGTHS <- c(
    "MiSeq v2 (2×250)" = 250,
    "MiSeq v3 (2×300)" = 300,
    "NovaSeq (2×150)"  = 150,
    "iSeq (2×150)"     = 150
)

# -----------------------------------------------------------------------------
# 3.7 Export Workbook Location
# -----------------------------------------------------------------------------
# Export writes the workbook where Step 5 expects it. The workbook contains an
# Info sheet, a numeric Parameters sheet consumed by Step 5, and a generated
# Column_Dictionary sheet. The Info sheet mixes text and numeric values, so its
# numeric cells are restored after the shared table writer has saved the file;
# see fix_excel_numeric_typed_cells().

REPORT_EXCEL_OUTPUT_DIR <- file.path("..", "..", "results", "4_dada2_parameter_selection")
REPORT_EXCEL_FILENAME <- "dada2_filter_parameters.xlsx"
REPORT_EXCEL_PARAMETERS_SHEET <- "Parameters"
REPORT_EXCEL_INFO_SHEET <- "Info"

# Project-relative label shown in the Export table. File operations use
# REPORT_EXCEL_OUTPUT_DIR above.
REPORT_EXCEL_SAVE_LOCATION_LABEL <- file.path("results", "4_dada2_parameter_selection", REPORT_EXCEL_FILENAME)

# Column descriptions shared by the Parameters and Info sheets and used to
# build the trailing Column_Dictionary sheet.
REPORT_EXCEL_COLUMN_DESCRIPTIONS <- c(
  Parameter   = "For the Parameters sheet: the exact variable name R/notebooks/5_dada2_pipeline.Rmd (Step 5) assigns this value to. For the Info sheet: a descriptive key for reference only, not read by Step 5.",
  Value       = "The value for this parameter or reference field.",
  Description = "A short human-readable explanation of this row."
)

# Restores real Excel numeric typing to specific cells in an already-written
# sheet's Value column -- the workaround for openxlsx's inability to write a
# single column with genuinely mixed per-cell types (see the long comment
# above REPORT_EXCEL_OUTPUT_DIR for the full investigation). Reloads the
# saved workbook and issues one single-cell writeData() call per numeric
# row, then re-saves; every other cell already written by
# add_sheet_to_excel() (including the still-text cells in the same Value
# column) is left untouched.
#
# workbook_path: path to the already-saved .xlsx -- add_sheet_to_excel()
#   must have already written `sheet_name` into it.
# sheet_name: the sheet to patch (e.g. REPORT_EXCEL_INFO_SHEET).
# ordered_table: the EXACT data.frame already written to that sheet (i.e.
#   after apply_report_row_order() and any column drops) -- used only to
#   look up which Excel row each Parameter ended up on. Row 1 of the sheet
#   is the header, so data row i is Excel row i + 1.
# numeric_values: a named numeric vector (e.g. report_info_numeric_values())
#   whose names match values in ordered_table$Parameter, giving the real
#   numeric value that Parameter's Value cell should hold. A name not
#   present in ordered_table$Parameter is silently skipped (keeps this
#   function reusable even if a future numeric_values entry doesn't apply
#   to every sheet it might be used on); an NA value is also skipped,
#   leaving that cell exactly as add_sheet_to_excel() already wrote it
#   (blank -- writeDataTable()'s default keepNA = FALSE already omits NA
#   values rather than writing the literal text "NA").
fix_excel_numeric_typed_cells <- function(workbook_path, sheet_name, ordered_table, numeric_values) {
    value_col <- which(names(ordered_table) == "Value")
    if (length(value_col) != 1) {
        stop("fix_excel_numeric_typed_cells(): expected exactly one 'Value' column in ordered_table.")
    }

    wb <- loadWorkbook(workbook_path)
    for (parameter_name in names(numeric_values)) {
        raw_value <- numeric_values[[parameter_name]]
        if (is.na(raw_value)) next
        data_row <- which(ordered_table$Parameter == parameter_name)
        if (length(data_row) != 1) next
        writeData(wb, sheet_name, x = raw_value,
                  startRow = data_row + 1, startCol = value_col, colNames = FALSE)
    }
    destination_directory <- dirname(workbook_path)
    temporary_path <- tempfile(".workbook-", tmpdir = destination_directory, fileext = ".xlsx")
    on.exit(unlink(temporary_path), add = TRUE)
    saveWorkbook(wb, temporary_path, overwrite = TRUE)
    if (!file.rename(temporary_path, workbook_path)) {
        stop("fix_excel_numeric_typed_cells(): could not atomically replace workbook: ", workbook_path)
    }
    invisible(workbook_path)
}


# =============================================================================
# SECTION 4: 16S rRNA REFERENCE DATA
# =============================================================================
# The 16S ribosomal RNA gene is the gold standard for bacterial identification.
# Its structure of conserved and variable regions enables universal PCR with
# species-specific discrimination.

# -----------------------------------------------------------------------------
# 4.1 Variable Region Coordinates
# -----------------------------------------------------------------------------
# These coordinates are based on the E. coli 16S rRNA gene (GenBank J01859),
# which is the standard reference for 16S position numbering.
#
# IMPORTANT: Actual positions vary between bacterial species due to indels.
# These coordinates are approximate guides, not absolute positions.

variable_regions <- data.frame(
    region = c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9"),
    
    # Start positions (5' boundary of each variable region)
    start = c(69, 137, 433, 576, 822, 986, 1117, 1243, 1435),
    
    # End positions (3' boundary of each variable region)
    end = c(99, 242, 497, 682, 879, 1043, 1173, 1294, 1465),
    
    stringsAsFactors = FALSE
) %>%
    mutate(
        # Calculate derived properties
        length = end - start + 1,       # Length of each region in bp
        midpoint = (start + end) / 2,   # Center point for label placement
        
        # Assign visually distinct colors for the visualization
        # Using a carefully chosen palette for accessibility
        color = c("#E63946",   # V1: Red
                  "#F4A261",   # V2: Orange
                  "#E9C46A",   # V3: Yellow
                  "#2A9D8F",   # V4: Teal (most common target)
                  "#264653",   # V5: Dark blue
                  "#8338EC",   # V6: Purple
                  "#FF006E",   # V7: Pink
                  "#3A86FF",   # V8: Blue
                  "#06D6A0")   # V9: Green
    )

# -----------------------------------------------------------------------------
# 4.2 Conserved Region Coordinates
# -----------------------------------------------------------------------------
# Conserved regions flank the variable regions and serve as primer binding sites.
# Their high sequence similarity across bacteria enables "universal" primers.

conserved_regions <- data.frame(
    name = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9"),
    
    # Conserved regions fill the gaps between variable regions
    start = c(1, 100, 243, 498, 683, 880, 1044, 1174, 1295),
    end = c(68, 136, 432, 575, 821, 985, 1116, 1242, 1434)
)

# Total length of the 16S gene (E. coli reference)
FULL_16S_LENGTH <- 1542

# -----------------------------------------------------------------------------
# 4.3 Primer Database
# -----------------------------------------------------------------------------
# Common primer pairs used in 16S studies, with literature references.
# This helps users choose appropriate primers for their research goals.

primer_database <- data.frame(
    # Human-readable name for the primer pair. Rows 1 and 2's forward
    # primer was renamed from the more commonly seen "27F" to "8F" -- see
    # the forward_start comment below and the Visualize tab's
    # collapsible Help pane for the full explanation (both names refer to the
    # same E. coli positions 8-27 footprint; "8F" is used here because it
    # matches forward_start's 5'-start convention, avoiding the apparent
    # "27F but forward_start=8" mismatch that prompted this rename).
    #
    # Full-length V1-V9 primer pairs are outside this paired-end Illumina
    # short-read workflow and are therefore not offered.
    name = c(
        "8F / 338R",       # Classic V1-V2 primers (commonly cited as "27F/338R")
        "8F / 519R",       # Extended V1-V3 coverage (commonly cited as "27F/519R")
        "341F / 785R",     # V3-V4 (Herlemann)
        "341F / 806R",     # V3-V4 (widely used)
        "515F / 806R",     # V4 only (Earth Microbiome Project)
        "515F / 926R",     # V4-V5 (improved Archaea)
        "515F / 944R",     # V4-V5 (marine)
        "784F / 1061R",    # V5-V6
        "926F / 1392R",    # V6-V8
        "967F / 1391R"     # V6-V8 (deep sequencing)
    ),

    # Individual primer names
    forward_name = c("8F", "8F", "341F", "341F", "515F", "515F", "515F",
                     "784F", "926F", "967F"),
    reverse_name = c("338R", "519R", "785R", "806R", "806R", "926R", "944R",
                     "1061R", "1392R", "1391R"),

    # Binding positions on 16S gene (5'-most base for forward_start, 3'-most
    # base for reverse_end, both in standard E. coli J01859 16S rRNA gene
    # numbering).
    #
    # Note on rows 1 and 2 (8F/338R, 8F/519R) -- both share forward_start =
    # 8 and the forward_name "8F". This primer is far more commonly cited
    # in the literature as "27F" (its older, Lane 1991-era name, which
    # numbers the primer by its 3'-most base instead of its 5'-most base).
    # "27F" and "8F" are the exact same 20-nt primer footprint -- E. coli
    # positions 8-27 -- confirmed via Frank et al. 2008 (Appl Environ
    # Microbiol 74(8):2461-2470, doi:10.1128/aem.02272-07), which
    # explicitly states "27f (spanning positions 8 to 27 in Escherichia
    # coli rRNA coordinates)". This table uses "8F" for both rows instead,
    # deliberately, so forward_name always matches forward_start directly
    # (as it already does for every other forward primer here -- 341F=341,
    # 515F=515, 784F=784, 926F=926, 967F=967) -- see Visualize Help for the
    # full explanation aimed at end
    # users who may only know the "27F" name.
    forward_start = c(8, 8, 341, 341, 515, 515, 515, 784, 926, 967),
    reverse_end = c(338, 519, 785, 806, 806, 926, 944, 1061, 1392, 1391),

    # Which variable regions are targeted
    target_regions = c("V1-V2", "V1-V3", "V3-V4", "V3-V4", "V4", "V4-V5",
                       "V4-V5", "V5-V6", "V6-V8", "V6-V8"),

    # Forward/reverse primer nucleotide sequences (5' -> 3', standard IUPAC
    # degenerate-base codes preserved as published -- M=A/C, W=A/T, R=A/G,
    # Y=C/T, N=any, V=A/C/G, H=A/C/T, K=G/T). Verified individually against
    # either the primer's originating publication or a well-established
    # protocol/database (not from memory alone). Two entries (944R,
    # 784F/1061R) rest only on secondary sources -- flagged via
    # sequence_note below, surfaced in the Primer Database table as a
    # footnote marker with its explanation as a footnote beneath the
    # table, and called out in that section's caption -- since the
    # primary publications were inaccessible (CAPTCHA/paywall) at
    # verification time and should be independently confirmed before
    # ordering physical oligos from either.
    forward_seq = c(
        "AGAGTTTGATCMTGGCTCAG",   # 8F, degenerate-M variant (8F/338R; commonly cited as 27F/338R)
        "AGAGTTTGATCMTGGCTCAG",   # 8F, degenerate-M variant (8F/519R; commonly cited as 27F/519R)
        "CCTACGGGNGGCWGCAG",      # 341F (341F/785R)
        "CCTACGGGNGGCWGCAG",      # 341F (341F/806R)
        "GTGCCAGCMGCCGCGGTAA",    # 515F, original/Caporaso form (515F/806R)
        "GTGYCAGCMGCCGCGGTAA",    # 515F, Parada-modified form (515F/926R)
        "GTGCCAGCMGCCGCGGTAA",    # 515F, original/Caporaso form (515F/944R)
        "AGGATTAGATACCCTGGTA",    # 784F
        "AAACTYAAAKGAATTGACGG",   # 926F
        "CAACGCGAAGAACCTTACC"     # 967F (single-oligo representation; see note)
    ),
    reverse_seq = c(
        "TGCTGCCTCCCGTAGGAGT",    # 338R
        "GWATTACCGCGGCKGCTG",     # 519R
        "GACTACHVGGGTATCTAATCC",  # 785R
        "GGACTACHVGGGTWTCTAAT",   # 806R
        "GGACTACHVGGGTWTCTAAT",   # 806R (identical sequence, confirmed consistent)
        "CCGYCAATTYMTTTRAGTTT",   # 926R
        "GAATTAAACCACATGCTC",     # 944R (low confidence -- see sequence_note)
        "CRRCACGAGCTGACGAC",      # 1061R
        "ACGGGCGGTGTGTRC",        # 1392R
        "GACGGGCGGTGWGTRCA"       # 1391R
    ),

    # Short, table-facing caveat shown only for the lower-confidence
    # sequences above (empty string for every well-corroborated entry).
    # Note: text below must avoid literal double-quote characters -- it is
    # interpolated into an HTML title='...' attribute (see vis_primer_table
    # below) delimited by double quotes, so an embedded " would truncate
    # the tooltip early. Use single quotes for any quoting within the text.
    sequence_note = c(
        "Labeled '8F' here to match forward_start=8, but far more commonly cited in the literature as '27F' (Lane 1991-era name, numbered by its 3' end instead). See this tab's own Help section (Primer Database) below for the full naming history.",
        "Labeled '8F' here to match forward_start=8, but far more commonly cited in the literature as '27F' (Lane 1991-era name, numbered by its 3' end instead). See this tab's own Help section (Primer Database) below for the full naming history.",
        "", "", "", "",
        "No verified primary publication defines this primer pair (see Reference column) -- sequence traces only to a compiled review table; confirm before use.",
        "Sequence verified only against a secondary source citing Andersson et al. 2008 -- confirm against the original Table 3 before use.",
        "",
        "Shown as a single representative oligo; the Sogin lab's own protocol deploys 967F as a 4-variant degenerate mix -- confirm which form your workflow needs."
    ),
    
    # Literature citations for each primer pair. 338R is sourced to Amann et
    # al. 1990 and Daims et al. 1999; 785R to Herlemann et al. 2011 and
    # Klindworth et al. 2013. No verified primary source was found for 515F/944R,
    # so the table explicitly identifies its compiled-review provenance.
    reference = c(
        "Lane 1991; Weisburg et al. 1991 (8F, commonly cited as 27F); Amann et al. 1990; Daims et al. 1999 (338R)",
        "Lane 1991; Turner et al. 1999",
        "Herlemann et al. 2011; Klindworth et al. 2013",
        "Caporaso et al. 2011; Klindworth et al. 2013",
        "Caporaso et al. 2012 (Earth Microbiome Project)",
        "Parada et al. 2016; Walters et al. 2016",
        "No verified primary source (see notes)",
        "Andersson et al. 2008",
        "Engelbrektson et al. 2010; Haas et al. 2011",
        "Sogin et al. 2006"
    ),

    # HTML version of `reference` above, with each citation hyperlinked to
    # its DOI (or, where a citation has no DOI -- Lane 1991 is a book
    # chapter -- to a stable lookup page instead). Used only by the Primer
    # Database table (vis_primer_table_help in Visualize Help, rendered with
    # sanitize.text.function = identity so this markup
    # is not HTML-escaped); the plain-text `reference` column above is left
    # untouched, since nothing else in the app uses it. Verified against
    # the publisher/PubMed record for each paper individually rather than
    # guessed.
    reference_html = c(
        # 8F / 338R, commonly cited as 27F/338R (8F: Lane 1991 + Weisburg et
        # al. 1991; 338R: Amann et al. 1990's original EUB338 probe + Daims
        # et al. 1999's reassessment)
        paste0('<a href="https://scholar.google.com/scholar?q=Lane+1991+16S%2F23S+rRNA+sequencing" target="_blank">Lane 1991</a>; ',
               '<a href="https://doi.org/10.1128/jb.173.2.697-703.1991" target="_blank">Weisburg et al. 1991</a> (8F, commonly cited as 27F); ',
               '<a href="https://doi.org/10.1128/aem.56.6.1919-1925.1990" target="_blank">Amann et al. 1990</a>; ',
               '<a href="https://doi.org/10.1016/S0723-2020(99)80053-8" target="_blank">Daims et al. 1999</a> (338R)'),
        # 8F / 519R, commonly cited as 27F/519R
        paste0('<a href="https://scholar.google.com/scholar?q=Lane+1991+16S%2F23S+rRNA+sequencing" target="_blank">Lane 1991</a>; ',
               '<a href="https://doi.org/10.1111/j.1550-7408.1999.tb04612.x" target="_blank">Turner et al. 1999</a>'),
        # 341F / 785R: 785R originates with Herlemann et al. 2011 and was
        # standardized as S-D-Bact-0785-a-A-21 by Klindworth et al. 2013.
        paste0('<a href="https://doi.org/10.1038/ismej.2011.41" target="_blank">Herlemann et al. 2011</a>; ',
               '<a href="https://doi.org/10.1093/nar/gks808" target="_blank">Klindworth et al. 2013</a>'),
        # 341F / 806R
        paste0('<a href="https://doi.org/10.1073/pnas.1000080107" target="_blank">Caporaso et al. 2011</a>; ',
               '<a href="https://doi.org/10.1093/nar/gks808" target="_blank">Klindworth et al. 2013</a>'),
        # 515F / 806R
        '<a href="https://doi.org/10.1038/ismej.2012.8" target="_blank">Caporaso et al. 2012</a> (Earth Microbiome Project)',
        # 515F / 926R
        paste0('<a href="https://doi.org/10.1111/1462-2920.13023" target="_blank">Parada et al. 2016</a>; ',
               '<a href="https://doi.org/10.1128/msystems.00009-15" target="_blank">Walters et al. 2016</a>'),
        # 515F / 944R -- no verifiable primary source found (see notes
        # column and the citation-correction comment above); the sequence
        # traces only to a compiled review-table entry (Fuks et al. 2018,
        # Microbiome) rather than an original Fuhrman-lab/Needham publication
        paste0('<span class="text-danger">No verified primary source</span> -- compiled in ',
               '<a href="https://doi.org/10.1186/s40168-017-0396-x" target="_blank">Fuks et al. 2018</a>'),
        # 784F / 1061R
        '<a href="https://doi.org/10.1371/journal.pone.0002836" target="_blank">Andersson et al. 2008</a>',
        # 926F / 1392R
        paste0('<a href="https://doi.org/10.1038/ismej.2009.153" target="_blank">Engelbrektson et al. 2010</a>; ',
               '<a href="https://genome.cshlp.org/content/21/3/494" target="_blank">Haas et al. 2011</a>'),
        # 967F / 1391R
        '<a href="https://doi.org/10.1073/pnas.0605127103" target="_blank">Sogin et al. 2006</a>'
    ),

    # Practical notes about each primer pair
    notes = c(
        # Naming history for 8F (commonly cited as 27F) lives in the
        # forward_start comment above, the sequence_note column (rendered
        # in-app as a footnote marker), and the Visualize tab's
        # collapsible Help pane -- not duplicated here since this notes
        # column is omitted from output$vis_primer_table_help.
        "Classic universal primer pair, good for V1-V2",
        "Extended V1-V3 coverage",
        "Popular for Illumina, misses some taxa",
        "Widely used V3-V4, good taxonomic resolution",
        "Earth Microbiome Project standard, V4 only",
        "Modified 515F with better Archaea coverage",
        # Updated to disclose the citation-verification finding above: no
        # Fuhrman-lab or Needham publication actually defines this primer
        # pair, so its use should be treated with more caution than the
        # other entries in this table.
        "Marine microbiome studies -- primer sequence/pairing not traceable to a verified primary publication",
        "Good for specific environments",
        "Broad coverage V6-V8",
        "Deep sequencing studies"
    ),
    
    stringsAsFactors = FALSE
) %>%
    mutate(
        # Calculate expected amplicon size
        amplicon_size = reverse_end - forward_start + 1,

        # Primer-trimmed insert length bounds (bp) -- the sequence length
        # DADA2 actually sees in colnames(seqtab) after cutadapt (Step 3)
        # has removed both primers, NOT the raw primer-to-primer span
        # above (amplicon_size). These are the real, per-primer-pair
        # values that drive Visualize > Target Length and Export's
        # amplicon_min_length/amplicon_max_length filter
        # bounds (Section 8.14).
        #
        # Derivation: for each primer pair, in silico PCR was run against four
        # reference databases --
        # SILVA 138.2 SSURef NR99 (~200k bacterial sequences), a
        # phylogenetically stratified NCBI sample (500 sequences x 10
        # major phyla = 5,000 total), GTDB r220 bacterial species
        # representatives (66k full-length sequences), and NCBI RefSeq
        # Targeted Loci 16S bacteria (26k sequences). For each database,
        # the primer-trimmed amplicon lengths were extracted, the 1st and
        # 99th percentiles computed, and the filter bounds set as p1-2
        # (min) and p99+2 (max). The final bounds below are the UNION
        # across all four sources -- the lowest minimum and the highest
        # maximum from any single database -- so a legitimate amplicon
        # captured by any one database is not discarded. See the Help
        # tab's / Visualizer help pane's Primer Database section for the
        # same explanation surfaced in Visualize Help.
        amplicon_min_length = c(273, 430, 399, 401, 249, 362, 401, 242, 440, 387),
        amplicon_max_length = c(345, 525, 436, 438, 260, 384, 425, 269, 483, 419)
    )

# -----------------------------------------------------------------------------
# 4.4 Primer Database Footnote Marker Assignment
# -----------------------------------------------------------------------------
# Maps each distinct non-empty sequence_note to a superscript letter in first-
# appearance order. Rows sharing the same note share one marker.
# Shared by build_primer_database_table() (which appends each row's marker
# to its Fwd/Rev sequence) and build_primer_database_footnote_html() (which
# lists each DISTINCT note once, keyed by this same marker) so a marker in
# the table always points to exactly one matching footnote entry below it.
primer_database_footnote_markers <- function() {
    distinct_notes <- unique(primer_database$sequence_note[primer_database$sequence_note != ""])
    setNames(letters[seq_along(distinct_notes)], distinct_notes)
}

# -----------------------------------------------------------------------------
# 4.5 Primer Database Display-Table Builder
# -----------------------------------------------------------------------------
# Builds the formatted data frame rendered in Visualize Help.
build_primer_database_table <- function() {
    footnote_markers <- primer_database_footnote_markers()
    # as.integer() calls must live in mutate(), not in select().
    # select() only accepts column selection expressions; passing function
    # calls like as.integer(col) causes an "object not found" error in
    # current dplyr because the expression is interpreted as a column name,
    # not a transformation.
    primer_database %>%
        mutate(
            forward_start = as.integer(forward_start),
            reverse_end   = as.integer(reverse_end),
            amplicon_size = as.integer(amplicon_size),
            # Fwd/Rev primer length (bp), each primer's own end
            # position, and the reverse primer's start position --
            # computed from the actual verified sequences (forward_seq/
            # reverse_seq, Section 4.3) rather than stored separately,
            # so they can never drift out of sync with the sequences
            # shown in the Fwd/Rev Primer columns below. forward_start
            # is already the forward primer's 5' start and reverse_end
            # is already the reverse primer's 3' end (both standard
            # E. coli J01859 numbering, see the forward_start comment
            # in Section 4.3) -- forward_end/reverse_start fill in the
            # other two corners of each primer's own footprint.
            forward_length = nchar(forward_seq),
            reverse_length = nchar(reverse_seq),
            forward_end    = forward_start + forward_length - 1L,
            reverse_start  = reverse_end - reverse_length + 1L,
            amplicon_min_length = as.integer(amplicon_min_length),
            amplicon_max_length = as.integer(amplicon_max_length),
            # Per-row footnote marker letter (Section 4.4) -- "" for rows
            # with no caveat (sequence_note == ""). unname() avoids the
            # named-character-vector names leaking into this column.
            footnote_marker = unname(ifelse(
                sequence_note == "", "",
                footnote_markers[sequence_note]
            )),
            # Sequences rendered in a monospace <code> tag (readable as a
            # literal oligo string), with a superscript footnote marker
            # appended for entries whose sequence_note flags a caveat (a
            # secondary-source-only verification, a naming-convention
            # mismatch, or a variant/mixture caveat -- see Section 4).
            # Distinct caveats receive distinct superscript letters; the full
            # note is always visible below the table rather than on hover.
            fwd_seq_html = ifelse(
                sequence_note == "",
                paste0("<code>", forward_seq, "</code>"),
                paste0("<code>", forward_seq, "</code> <sup style='color:#c77700;'>", footnote_marker, "</sup>")
            ),
            rev_seq_html = ifelse(
                sequence_note == "",
                paste0("<code>", reverse_seq, "</code>"),
                paste0("<code>", reverse_seq, "</code> <sup style='color:#c77700;'>", footnote_marker, "</sup>")
            )
        ) %>%
        select(
            `Target`      = target_regions,
            `Primer Pair` = name,
            `Fwd Primer (5'->3')` = fwd_seq_html,
            # Exact footprint of each primer in E. coli reference coordinates.
            `Fwd Start (bp)` = forward_start,
            `Fwd End (bp)`   = forward_end,
            `Fwd Len (bp)`   = forward_length,
            `Rev Primer (5'->3')` = rev_seq_html,
            `Rev Start (bp)` = reverse_start,
            `Rev End (bp)`   = reverse_end,
            `Rev Len (bp)`   = reverse_length,
            # Amplicon (bp): the full primer-to-primer span (Fwd Start
            # to Rev End inclusive) as conventionally reported in the
            # primer's own literature -- i.e. it INCLUDES both primers'
            # lengths above, unlike Visualize's "Amplicon (bp,
            # w/o primers)" metric/Gene Map, which excludes them.
            `Amplicon (bp)` = amplicon_size,
            # Primer-trimmed insert bounds used to populate Visualize >
            # Target Length; see the derivation in Section 4.3.
            `Target Min (bp)` = amplicon_min_length,
            `Target Max (bp)` = amplicon_max_length,
            # Hyperlinked citations (reference_html, Section 4) rather
            # than the plain-text `reference` column, so each citation
            # is clickable. sanitize.text.function below is required for
            # the <a href> markup to render instead of being escaped.
            `Reference`   = reference_html
        )
}

# -----------------------------------------------------------------------------
# 4.6 Primer Database Footnote Builder
# -----------------------------------------------------------------------------
# Returns one always-visible footnote per distinct caveat, listing every primer
# pair to which the note applies.
build_primer_database_footnote_html <- function() {
    flagged <- primer_database[primer_database$sequence_note != "", , drop = FALSE]
    if (nrow(flagged) == 0) return(NULL)

    marker_map <- primer_database_footnote_markers()
    distinct_notes <- names(marker_map)

    tagList(
        p(class = "text-muted small mb-1 mt-2", tags$strong("Footnotes:")),
        tagList(lapply(seq_along(distinct_notes), function(i) {
            note_text <- distinct_notes[i]
            marker <- marker_map[[note_text]]
            # Every primer-pair name sharing this exact note text (usually
            # one, sometimes more -- e.g. 8F/338R and 8F/519R) is listed
            # together, so the footnote explicitly states which row(s) it
            # applies to.
            matching_names <- flagged$name[flagged$sequence_note == note_text]
            p(class = "text-muted small mb-1",
              tags$sup(style = "color:#c77700;", marker), " ",
              tags$strong(paste(matching_names, collapse = ", ")), ": ", note_text)
        }))
    )
}

# -----------------------------------------------------------------------------
# 4.7 Target Region <-> Literature Primer Helpers
# -----------------------------------------------------------------------------
# Translate a selected target (for example V3-V4) into available primer pairs
# and the individual region codes highlighted on the Gene Map.

# Distinct Target values actually present in primer_database, in the
# table's row order. This exposes only combinations represented in the primer
# database and avoids redundant target choices.
VIS_TARGET_REGION_CHOICES <- unique(primer_database$target_regions)

# Named startup defaults keep the two initial controls consistent. The primer
# is the first V3-V4 database entry, matching the auto-selection rule.
DEFAULT_TARGET_REGION <- "V3-V4"
DEFAULT_PRIMER_PAIR   <- "341F / 785R"

# Literature Primers dropdown choices for a given Target Region Selection
# value, preserving primer_database row order. The caller auto-selects the first
# match after a target-region change.
vis_primer_choices_for_target <- function(target) {
    matches <- primer_database[primer_database$target_regions == target, , drop = FALSE]
    setNames(matches$name, matches$name)
}

# Individual variable-region codes (e.g. c("V3", "V4")) spanned by a given
# Target Region Selection value -- used only to translate the single
# selected target string into the set of variable_regions$region codes the
# Gene Map (output$vis_main_plot, Section 8.6) draws with a highlighted
# selection border.
target_to_region_codes <- function(target) {
    switch(target,
           "V1-V2" = c("V1", "V2"),
           "V1-V3" = c("V1", "V2", "V3"),
           "V3-V4" = c("V3", "V4"),
           "V4"    = c("V4"),
           "V4-V5" = c("V4", "V5"),
           "V5-V6" = c("V5", "V6"),
           "V6-V8" = c("V6", "V7", "V8"),
           c("V4"))
}


# =============================================================================
# SECTION 5: HELPER FUNCTIONS - FILE OPERATIONS
# =============================================================================

# -----------------------------------------------------------------------------
# 5.1 Paired File Detection
# -----------------------------------------------------------------------------
# Exact forward/reverse pairing and sample-count helpers are defined in
# functions/paired_read_retention_engine_function.R. They match complete stems before
# deriving unique sample IDs and are shared with the regression tests.


# =============================================================================
# SECTION 6: QUALITY PROFILE SUMMARIZATION AND VALIDATION HELPERS
# =============================================================================

# -----------------------------------------------------------------------------
# Derive the fields used by the quality plots and retention calculations from
# one direction of the row-aligned paired quality sample.
summarize_quality_matrix <- function(qual_matrix) {
    if (is.null(qual_matrix) || nrow(qual_matrix) == 0) return(NULL)

    n_sampled <- nrow(qual_matrix)
    max_len   <- ncol(qual_matrix)

    # Per-position median quality (robust to outliers).
    position_stats <- data.frame(
        Position = 1:max_len,
        Median   = apply(qual_matrix, 2, median, na.rm = TRUE)
    )

    # Per-position expected error and its cumulative sum. NA-padded positions
    # beyond a read's own length remain NA and are handled downstream.
    error_matrix   <- 10^(-qual_matrix / 10)
    cum_ee_matrix  <- t(apply(error_matrix, 1, cumsum))
    if (n_sampled == 1) cum_ee_matrix <- matrix(cumsum(error_matrix[1, ]), nrow = 1)

    list(
        n_sampled      = n_sampled,
        position_stats = position_stats,
        cum_ee_matrix  = cum_ee_matrix,
        qual_matrix    = qual_matrix
    )
}


# -----------------------------------------------------------------------------
# Adaptive default for Validate's sample-count field, mirroring
# recommended_pairs_per_sample() from the retention engine but tuned for a much
# more expensive per-unit cost: each validation sample runs the REAL DADA2 pipeline
# (filterAndTrim -> learnErrors -> dada -> mergePairs -> removeBimeraDenovo),
# minutes rather than milliseconds. Small cohorts are validated in full since
# the real cost stays low; larger cohorts are capped well below their total
# size, because select_validation_samples()'s retention-stratified spread
# already gives good coverage of the cohort's quality range without
# validating every single sample.
recommended_validation_samples <- function(n_total_samples) {
    if (!is.finite(n_total_samples) || n_total_samples <= 0) return(3L)
    if (n_total_samples <= 5)   return(as.integer(n_total_samples))  # small cohort: validate all of it
    if (n_total_samples <= 20)  return(5L)
    if (n_total_samples <= 50)  return(6L)
    if (n_total_samples <= 100) return(8L)
    10L  # very large cohorts: cap runtime, stratified coverage is already good
}

# -----------------------------------------------------------------------------
# select_validation_samples: choose which samples the empirical DADA2 stage
# runs on, using a retention-stratified spread to avoid bias from file order.
#
# The samples are ranked by their per-sample TRUE paired retention (paired_pass
# = fwd_pass & rev_pass) at a fixed reference parameter set (the top-ranked
# candidate), then the selection walks evenly spaced quantile positions across
# that ranking so the chosen subset spans the cohort's quality range -- always
# anchoring the lowest-retention (p0, worst-case) and highest-retention (p1)
# samples when n >= 2, with the interior filled at even quantiles (so n = 3
# gives lowest / median / highest). This tests the quality range instead of one
# corner of the distribution. For n = 1 the most representative sample (median)
# is used rather than the worst, so a lone-sample sanity check is not alarmingly
# pessimistic. Samples whose retention cannot be computed (NA) sink to the low
# end so they are stress-tested rather than silently dropped.
#
# retention_by_sample : named numeric vector, per-sample paired retention (%),
#                       aligned to sample order; NA allowed.
# n_val               : desired number of validation samples (clamped to cohort).
# Returns a data.frame (one row per chosen sample, ordered lowest -> highest
# retention) with: index (position in the original sample order), name,
# retention (%), rank (1 = lowest retention in the cohort), and role
# ("lowest retention" / "median" / "highest retention" / "pNN retention tier").
select_validation_samples <- function(retention_by_sample, n_val) {

    if (!is.numeric(retention_by_sample)) {
        stop("select_validation_samples(): `retention_by_sample` must be numeric.")
    }
    if (!is.numeric(n_val) || length(n_val) != 1L || is.na(n_val) ||
        !is.finite(n_val) || n_val < 1) {
        stop("select_validation_samples(): `n_val` must be one positive number.")
    }

    # ---- Clamp the request to what the cohort actually offers -----------------
    n_total <- length(retention_by_sample)                 # samples available
    if (n_total == 0L) {
        return(data.frame(index = integer(0), name = character(0),
                          retention = numeric(0), rank = integer(0),
                          role = character(0), stringsAsFactors = FALSE))
    }
    n_val <- max(1L, min(as.integer(n_val), n_total))      # 1 .. n_total

    # ---- Rank samples ascending by finite retention -----------------------------
    # Missing estimates are eligible only if the requested subset is larger than
    # the finite cohort; otherwise they cannot displace an informative sample.
    finite_indices <- which(is.finite(retention_by_sample))
    missing_indices <- which(!is.finite(retention_by_sample))
    finite_order <- finite_indices[order(retention_by_sample[finite_indices], decreasing = FALSE)]
    ord <- if (length(finite_order) >= n_val) finite_order else c(finite_order, missing_indices)
    selection_total <- length(ord)

    # ---- Target evenly spaced quantile positions + their intended roles --------
    # n_val == 1 -> the median (most representative single sample); n_val >= 2 ->
    # seq(0, 1) which always includes both extremes (worst and best) plus evenly
    # spaced interior tiers. Each target carries the role it is meant to fill, so
    # the label reflects design intent (e.g. the central pick is "median") rather
    # than being re-derived from a rank that even-sized cohorts never land on.
    target_q <- if (n_val == 1L) 0.5 else seq(0, 1, length.out = n_val)

    # Role that a given target quantile is meant to represent.
    label_for_q <- function(q) {
        if (isTRUE(all.equal(q, 0)))   return("lowest retention")
        if (isTRUE(all.equal(q, 1)))   return("highest retention")
        if (isTRUE(all.equal(q, 0.5))) return("median")
        sprintf("p%02d retention tier", round(q * 100))
    }
    target_role <- vapply(target_q, label_for_q, character(1))

    # Map each target quantile to a 1-based position in the sorted-by-retention
    # order (round to the nearest rank).
    target_pos <- round(target_q * (selection_total - 1)) + 1L

    # ---- Resolve collisions so exactly n_val DISTINCT samples are returned -----
    # With few samples (or quantile rounding) two targets can land on the same
    # sorted position; search outward (up first, then down) for the nearest
    # still-unused position so the returned set always has n_val members. The
    # role of each target travels with the position it ends up claiming.
    chosen_pos  <- integer(0)
    chosen_role <- character(0)
    for (k in seq_along(target_pos)) {
        p      <- target_pos[k]                            # preferred position
        cand_p <- p
        offset <- 0L
        while (cand_p %in% chosen_pos || cand_p < 1L || cand_p > selection_total) {
            offset <- offset + 1L
            up   <- p + offset                             # search above first
            down <- p - offset                             # then below
            cand_p <- if (up   <= selection_total && !(up   %in% chosen_pos)) up
                      else if (down >= 1L    && !(down %in% chosen_pos)) down
                      else NA_integer_                      # nothing left to take
            if (is.na(cand_p)) break
        }
        if (!is.na(cand_p)) {
            chosen_pos  <- c(chosen_pos, cand_p)
            chosen_role <- c(chosen_role, target_role[k])
        }
    }

    # ---- Order lowest -> highest retention, carrying roles ---------------------
    ord_pos     <- order(chosen_pos)
    chosen_pos  <- chosen_pos[ord_pos]
    chosen_role <- chosen_role[ord_pos]

    # ---- Truthfulness override for the extremes --------------------------------
    # After collision resolution, force whichever sample actually sits at rank 1
    # / rank n_total to read "lowest" / "highest retention" so those labels are
    # never wrong; a lone-cohort sample is simply "only sample".
    if (selection_total == 1L) {
        chosen_role[] <- "only sample"
    } else {
        chosen_role[chosen_pos == 1L]      <- "lowest retention"
        chosen_role[chosen_pos == selection_total] <- "highest retention"
    }

    # ---- Translate sorted positions back to original-order sample indices ------
    idx  <- ord[chosen_pos]                                # position in input order
    rank <- chosen_pos                                     # 1 = lowest retention

    # ---- Assemble the result, ordered lowest -> highest retention --------------
    nm  <- names(retention_by_sample)
    result <- data.frame(
        index     = idx,
        name      = if (is.null(nm)) as.character(idx) else nm[idx],
        retention = unname(retention_by_sample[idx]),
        rank      = rank,
        role      = chosen_role,
        stringsAsFactors = FALSE
    )
    result$role[!is.finite(result$retention)] <- "retention unavailable"
    result
}

# validate_candidates_with_dada2: empirical validation stage (brief section 16).
#
# Runs the REAL DADA2 pipeline (filterAndTrim -> learnErrors -> dada ->
# mergePairs -> makeSequenceTable -> removeBimeraDenovo) on a handful of
# representative samples for each of the top candidates, and returns a per-
# candidate track table (reads.in / filtered / denoisedF / denoisedR / merged /
# nonchim) plus a "Predicted merged %" (the surrogate: paired filter retention x
# overlap-mergeability) next to "Real merged %", and the gap between them. The
# comparison is deliberately against MERGED, not nonchim: the surrogate models
# filtering + overlap, both of which map to the merged step, whereas chimera
# removal is not controlled by truncLen/maxEE, so gapping against nonchim would
# unfairly penalize the surrogate for something it cannot predict. The nonchim
# column is retained purely as a reference. A large gap flags candidates where
# overlap-region mismatches (not modeled by the surrogate) cost real merges.
#
# Also returns a per-sample track table used by the Validate tab. reads.in and
# filtered come directly from
# filterAndTrim()'s own per-file rows; denoisedF/denoisedR/merged/nonchim are
# matched back to sample name (not position) via named dada2 objects, so a
# sample dropped during filtering can never silently misalign with its
# neighbor. This app passes one current parameter set, so $per_sample is one
# data frame.
#
# IMPORTANT: this function is heavy (minutes to tens of minutes on real data)
# and depends on the Bioconductor 'dada2' package, so it is NOT exercised by
# the sandbox unit tests -- only parse-checked and run on the user's machine.
# It is written to be safe if dada2 is absent (the caller guards first).
#
# candidates   : data frame with trunc_len_fwd/rev, max_ee_fwd/rev,
#                overall_paired_retention (the surrogate prediction).
# fnFs / fnRs  : full forward/reverse FASTQ paths (rv_exp$files$fnFs/fnRs).
# sample_names : sample labels aligned to fnFs/fnRs.
# n_samples    : how many samples to validate on (used only if sample_indices
#                is NULL, in which case it falls back to the leading n samples).
# sample_indices : explicit integer positions (into fnFs/fnRs/sample_names) of
#                the samples to validate on. Supplied by the caller from
#                select_validation_samples() so the empirical stage runs on a
#                retention-stratified spread rather than the first n files. When
#                NULL, the leading n samples are used as a standalone fallback.
# progress     : optional function(frac, msg) for a Shiny progress bar.
# Returns a list: $summary (the per-candidate aggregated table) and
# $per_sample (a data.frame, one row per validated sample, or NULL on error).
validate_candidates_with_dada2 <- function(candidates, fnFs, fnRs, sample_names,
                                           n_samples = 3, trunc_q = 2,
                                           sample_indices = NULL,
                                           progress = NULL) {
    if (!requireNamespace("dada2", quietly = TRUE)) {
        return(list(summary = data.frame(Note = "dada2 is not installed.", stringsAsFactors = FALSE),
                    per_sample = NULL))
    }
    note <- function(frac, msg) if (!is.null(progress)) progress(frac, msg)

    # Prefer the caller-supplied stratified indices; fall back to the leading n
    # samples only when none were provided. Both paths are clamped to valid,
    # in-range, de-duplicated positions so a bad index can never subscript past
    # the file list.
    if (!is.null(sample_indices)) {
        sel <- unique(as.integer(sample_indices))
        sel <- sel[is.finite(sel) & sel >= 1 & sel <= length(fnFs)]
        if (length(sel) == 0) sel <- seq_len(max(1, min(n_samples, length(fnFs))))
    } else {
        sel <- seq_len(max(1, min(n_samples, length(fnFs))))
    }
    subF   <- fnFs[sel]; subR <- fnRs[sel]; subN <- sample_names[sel]

    # A unique scratch directory prevents cross-session collisions; on.exit()
    # removes intermediate filtered FASTQs on success or error.
    work_dir <- tempfile(pattern = "dada2_validate_")
    dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
    on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

    n_cand <- nrow(candidates)
    results <- vector("list", n_cand)
    per_sample_track <- NULL

    for (ci in seq_len(n_cand)) {
        cand <- candidates[ci, ]
        note((ci - 1) / n_cand, sprintf("Candidate %d/%d: filtering...", ci, n_cand))

        filtF <- file.path(work_dir, paste0(subN, "_C", ci, "_F_filt.fastq.gz"))
        filtR <- file.path(work_dir, paste0(subN, "_C", ci, "_R_filt.fastq.gz"))
        # Name the filtered-file paths by sample (standard DADA2 tutorial
        # convention: names(filtFs) <- sample.names) so dada2::dada() and
        # dada2::mergePairs() return per-sample results keyed by sample name
        # whenever there is more than one sample -- this is what the
        # per-sample matching-by-name below relies on.
        names(filtF) <- subN
        names(filtR) <- subN

        one <- tryCatch({
            ft <- dada2::filterAndTrim(
                subF, filtF, subR, filtR,
                truncLen = c(cand$trunc_len_fwd, cand$trunc_len_rev),
                maxEE    = c(cand$max_ee_fwd, cand$max_ee_rev),
                truncQ   = trunc_q, maxN = 0, rm.phix = TRUE,
                matchIDs = TRUE,
                compress = TRUE, multithread = FALSE, verbose = FALSE)

            keep <- file.exists(filtF) & file.exists(filtR)
            if (!any(keep)) stop("no reads survived filtering")
            kept_names <- subN[keep]

            note((ci - 0.7) / n_cand, sprintf("Candidate %d/%d: learning errors...", ci, n_cand))
            errF <- dada2::learnErrors(filtF[keep], multithread = FALSE, verbose = 0)
            errR <- dada2::learnErrors(filtR[keep], multithread = FALSE, verbose = 0)

            note((ci - 0.4) / n_cand, sprintf("Candidate %d/%d: denoising + merging...", ci, n_cand))
            ddF <- dada2::dada(filtF[keep], err = errF, multithread = FALSE, verbose = 0)
            ddR <- dada2::dada(filtR[keep], err = errR, multithread = FALSE, verbose = 0)
            mrg <- dada2::mergePairs(ddF, filtF[keep], ddR, filtR[keep], verbose = FALSE)

            getN <- function(x) sum(dada2::getUniques(x))

            # Normalize ddF/ddR/mrg to NAMED lists (one element per kept
            # sample), regardless of how many samples were validated.
            # dada2::dada() returns a single "dada"-class object directly
            # (not wrapped in a list) when given exactly one sample -- and a
            # dada-class object is ITSELF internally list-shaped, so a plain
            # is.list() check cannot tell "one dada object" apart from "a
            # list of several dada objects" (is.list() is TRUE for both).
            # inherits(x, "dada") disambiguates correctly via the object's
            # actual class instead. mergePairs() has no such ambiguity: its
            # single-sample result is a plain data.frame, a type multi-
            # sample results (a list of data.frames) can never be confused
            # with, so is.data.frame() alone is sufficient there.
            ddF_list <- if (inherits(ddF, "dada")) setNames(list(ddF), kept_names[1]) else ddF
            ddR_list <- if (inherits(ddR, "dada")) setNames(list(ddR), kept_names[1]) else ddR
            mrg_list <- if (is.data.frame(mrg)) setNames(list(mrg), kept_names[1]) else mrg

            # seqtab/seqtab2 are built from mrg_list (always a named list,
            # never a bare single data.frame) specifically so
            # makeSequenceTable()'s resulting matrix rownames are always the
            # sample names, for every sample count -- avoiding any ambiguity
            # in how it names a single-sample result.
            seqtab  <- dada2::makeSequenceTable(mrg_list)
            seqtab2 <- dada2::removeBimeraDenovo(seqtab, method = "consensus",
                                                 multithread = FALSE, verbose = FALSE)

            denoisedF_by_name <- vapply(ddF_list, getN, numeric(1))
            denoisedR_by_name <- vapply(ddR_list, getN, numeric(1))
            merged_by_name    <- vapply(mrg_list, getN, numeric(1))
            nonchim_by_name   <- setNames(rowSums(seqtab2), rownames(seqtab2))

            reads_in  <- sum(ft[, 1])
            filtered  <- sum(ft[, 2])
            denoisedF <- sum(denoisedF_by_name)
            denoisedR <- sum(denoisedR_by_name)
            merged    <- sum(merged_by_name)
            nonchim   <- sum(nonchim_by_name)

            # Per-sample breakdown -- reads.in/filtered are taken directly
            # from filterAndTrim()'s own per-file rows (always available,
            # aligned to subN, regardless of which samples survived
            # filtering); denoisedF/denoisedR/merged/nonchim are looked up
            # BY NAME so a dropped sample (kept == FALSE, rare given the
            # FASTQs are already primer-trimmed) gets NA there instead of
            # silently misaligning with a neighboring sample.
            per_sample_track <- data.frame(
                sample    = subN,
                reads_in  = as.numeric(ft[, 1]),
                filtered  = as.numeric(ft[, 2]),
                denoisedF = unname(denoisedF_by_name[subN]),
                denoisedR = unname(denoisedR_by_name[subN]),
                merged    = unname(merged_by_name[subN]),
                nonchim   = unname(nonchim_by_name[subN]),
                stringsAsFactors = FALSE
            )

            # Predicted merged fraction = predicted filter survival (paired
            # retention) x predicted overlap-mergeability (mergeable fraction).
            # This is the fair surrogate to compare against REAL merged: both
            # cover "survives filtering AND has enough overlap to merge". The
            # surrogate cannot model chimera removal, so nonchim is NOT the right
            # comparison (kept below only as a reference count/percentage).
            merge_frac <- if (is.finite(cand$mergeable_pair_fraction)) cand$mergeable_pair_fraction else 0
            predicted_merged_pct <- round(cand$overall_paired_retention * merge_frac * 100, 1)

            data.frame(
                Candidate       = paste0("truncLen ", cand$trunc_len_fwd, "/", cand$trunc_len_rev,
                                         ", maxEE ", cand$max_ee_fwd, "/", cand$max_ee_rev),
                reads.in        = reads_in,
                filtered        = filtered,
                denoisedF       = denoisedF,
                denoisedR       = denoisedR,
                merged          = merged,
                nonchim         = nonchim,
                `Predicted merged %` = predicted_merged_pct,
                `Real merged %`      = round(100 * merged / max(reads_in, 1), 1),
                `Real nonchim %`     = round(100 * nonchim / max(reads_in, 1), 1),
                check.names = FALSE, stringsAsFactors = FALSE)
        }, error = function(e) {
            merge_frac <- if (is.finite(cand$mergeable_pair_fraction)) cand$mergeable_pair_fraction else 0
            data.frame(
                Candidate = paste0("truncLen ", cand$trunc_len_fwd, "/", cand$trunc_len_rev,
                                   ", maxEE ", cand$max_ee_fwd, "/", cand$max_ee_rev),
                reads.in = NA, filtered = NA, denoisedF = NA, denoisedR = NA,
                merged = NA, nonchim = NA,
                `Predicted merged %` = round(cand$overall_paired_retention * merge_frac * 100, 1),
                `Real merged %` = NA, `Real nonchim %` = NA,
                Error = conditionMessage(e), check.names = FALSE, stringsAsFactors = FALSE)
        })
        results[[ci]] <- one
    }

    note(1, "Validation complete.")
    # Align columns across candidates (some rows may carry an Error column).
    all_cols <- unique(unlist(lapply(results, names)))
    results  <- lapply(results, function(df) {
        miss <- setdiff(all_cols, names(df))
        for (m in miss) df[[m]] <- NA
        df[, all_cols, drop = FALSE]
    })
    out <- do.call(rbind, results)
    # Gap vs MERGED (not nonchim): the surrogate predicts filtering + overlap
    # mergeability, which map to the real merged step; chimera removal is not
    # something these parameters control, so comparing to nonchim would be
    # unfair. A large gap flags where the surrogate over-predicts merging
    # (usually overlap-region mismatches, which the surrogate does not model).
    if (all(c("Real merged %", "Predicted merged %") %in% names(out))) {
        out$`Gap vs merged (pp)` <- ifelse(is.finite(out$`Real merged %`),
                                           round(out$`Predicted merged %` - out$`Real merged %`, 1), NA)
    }
    rownames(out) <- NULL
    list(summary = out, per_sample = per_sample_track)
}


# Vectorized retention calculation for per-sample table
#
# Calculates retention for one sample at the selected parameter values.

calculate_retention_vectorized <- function(qual_matrix, cum_ee_matrix, 
                                           truncLen, maxEE, truncQ) {
    if (is.null(qual_matrix) || is.null(cum_ee_matrix)) return(NA)
    if (truncLen > ncol(cum_ee_matrix) || truncLen < 1) return(NA)
    
    n_reads <- nrow(qual_matrix)
    if (n_reads == 0) return(NA)
    
    qual_subset <- qual_matrix[, 1:truncLen, drop = FALSE]
    qual_subset[is.na(qual_subset)] <- 40
    
    low_qual <- qual_subset <= truncQ
    
    first_low_qual <- apply(low_qual, 1, function(row) {
        pos <- which(row)
        if (length(pos) == 0) return(truncLen + 1)
        return(pos[1])
    })
    
    effective_len <- pmin(first_low_qual - 1, truncLen)
    fails_immediately <- effective_len < 1
    
    row_idx <- 1:n_reads
    col_idx <- pmax(effective_len, 1)
    ee_vals <- cum_ee_matrix[cbind(row_idx, col_idx)]
    
    passes <- !fails_immediately & !is.na(ee_vals) & ee_vals <= maxEE

    return(round(100 * sum(passes) / n_reads, 1))
}

# read_ee_at_truncLen: return EACH read's total expected errors at a given
# truncLen (with truncQ applied), rather than a single retention scalar. Reads
# that fail for length/truncQ reasons are set to Inf so they never pass any
# finite maxEE. This is the per-read basis for the retention-vs-maxEE curve
# retention(maxEE) = mean(ee <= maxEE) over the pooled reads,
# which is exactly the empirical CDF of read EE -- the correct object for
# choosing maxEE (the median cumulative-EE curve only shows the 50th-percentile
# read and cannot say what fraction a maxEE threshold retains). The pass logic
# mirrors calculate_retention_vectorized() exactly, so the aggregate retention
# derived here matches the per-sample table's Forward/Reverse % values.
read_ee_at_truncLen <- function(qual_matrix, cum_ee_matrix, truncLen, truncQ) {
    if (is.null(qual_matrix) || is.null(cum_ee_matrix)) return(numeric(0))
    if (truncLen > ncol(cum_ee_matrix) || truncLen < 1)  return(numeric(0))
    n_reads <- nrow(qual_matrix)
    if (n_reads == 0) return(numeric(0))

    qual_subset <- qual_matrix[, 1:truncLen, drop = FALSE]
    qual_subset[is.na(qual_subset)] <- 40
    low_qual <- qual_subset <= truncQ
    first_low_qual <- apply(low_qual, 1, function(row) {
        pos <- which(row); if (length(pos) == 0) truncLen + 1 else pos[1]
    })
    effective_len <- pmin(first_low_qual - 1, truncLen)
    fails         <- effective_len < 1
    col_idx       <- pmax(effective_len, 1)
    ee_vals       <- cum_ee_matrix[cbind(seq_len(n_reads), col_idx)]
    ee_vals[fails | is.na(ee_vals)] <- Inf   # never pass any finite maxEE
    ee_vals
}

# =============================================================================
# SECTION 7: USER INTERFACE DEFINITION
# =============================================================================
# The UI is built using bslib (Bootstrap 5) components for a modern, responsive
# design. The top-level navbar presents the four workflow steps.

ui <- page_navbar(
    id = "main_nav",
    
    # =========================================================================
    # 7.1 Application Title and Theme
    # =========================================================================
    
    title = tags$span(
        icon("dna"), " ",  # DNA helix icon for branding
        "DADA2 Parameter Explorer"
    ),
    
    # Bootstrap 5 theme with custom colors
    # Flatly is a clean, professional theme suitable for scientific applications
    theme = bs_theme(
        version = 5,                    # Bootstrap 5 for modern components
        bootswatch = "flatly",          # Clean, professional look
        primary = "#2c3e50",            # Dark blue for primary actions
        secondary = "#95a5a6",          # Gray for secondary elements
        success = "#18bc9c",            # Teal for success states
        info = "#3498db",               # Blue for information
        warning = "#f39c12",            # Orange for warnings
        danger = "#e74c3c",             # Red for errors/danger
        font_scale = 0.9                # Slightly smaller text for data-dense UI
    ),

    # =========================================================================
    # 7.1b Processing Console (relocated per-tab, see Section 7.4)
    # =========================================================================
    # The Processing Console is authored once, inside the Visualize tab's own
    # sidebar (the default active tab -- see Section 7.4 below), exactly
    # where it originally lived at the bottom of the Select tab's sidebar.
    # To keep it visible regardless of which of the four tabs is active, an
    # empty `console_slot_<value>` anchor div sits at the bottom of every
    # OTHER tab's own sidebar, and a server-side `observeEvent` on
    # `input$main_nav` (see Section 8.2b) tells the client, via a
    # `move_console` custom message handled in `header` below, to relocate
    # the single #exp_console_accordion DOM node to sit directly before the
    # newly active tab's anchor. It is inserted as a direct sidebar child in
    # every case (never appended inside the anchor div itself), because
    # bslib only gives a sidebar's own direct-child accordions its
    # edge-to-edge "flush" negative margin -- an intermediate wrapper would
    # leave the console visibly narrower than the sidebar. The console
    # therefore always appears "in the sidebar, at the bottom, full width"
    # of whichever tab is showing -- never as a separate page-level sidebar
    # -- and its message history survives the relocation because the same
    # DOM node (and #console_container within it) simply moves, rather than
    # being duplicated.

    # Cross-cutting setup that must appear on every page regardless of which
    # nav panel is active (shinyjs initialisation, custom CSS, custom JS
    # message handlers) is passed via `header` rather than as bare
    # positional arguments. page_navbar()'s ... expects only
    # nav_panel()/nav_menu() (plus nav_spacer()/nav_item()); passing
    # anything else positionally makes bslib log a "Navigation containers
    # expect ..." warning once per unrecognised element -- exactly what
    # useShinyjs() and the two tags$head() blocks below were causing (3
    # warnings, one per element, every app launch).
    header = tagList(

        # Enable shinyjs for JavaScript operations
        useShinyjs(),

        # =========================================================================
        # 7.2 Custom CSS Styles
        # =========================================================================
        # Additional styling for custom components not covered by Bootstrap

        tags$head(
        tags$style(HTML("
            /* ============================================================
               GLOBAL STYLES
               ============================================================ */
            
            /* App header styling with gradient background */
            .app-header {
                background: linear-gradient(135deg, #264653 0%, #2A9D8F 100%);
                color: white;
                padding: 1.5rem;
                margin: -1rem -1rem 1.5rem -1rem;
                border-radius: 0 0 0.75rem 0.75rem;
                box-shadow: 0 4px 15px rgba(0,0,0,0.1);
            }
            
            .app-header h2 {
                font-weight: 700;
                margin-bottom: 0.25rem;
            }
            
            .app-header p {
                opacity: 0.9;
                margin-bottom: 0;
            }

            /* ============================================================
               QUALITY PROFILES SLIDER ROWS (truncLen / maxEE)
               ============================================================ */
            /* Each row is: bold label | - button | ionRangeSlider | + button.
               align-items:center lines the label and buttons up with the
               slider TRACK. The ion slider by default prints its min/max as
               grey boxes ABOVE the track (.irs-min / .irs-max) AND again as the
               first/last numbers of the tick scale BELOW it (.irs-grid) --
               duplicate 0 / 250 labels. Hiding the top boxes keeps only the
               bottom tick scale. */
            .qp-slider { align-items: center; }
            .qp-slider .irs-min,
            .qp-slider .irs-max { display: none !important; }
            /* Raise the flanking label + buttons onto the slider TRACK. The
               ionRangeSlider box includes the tick grid below the track, so pure
               vertical centering lands them below the line; position:relative +
               a negative top pulls them back up level with the track without
               reflowing the row. (Estimated offset -- tweak if still off.) */
            .qp-slider > strong,
            .qp-slider > .btn { position: relative; top: -14px; }

            /* Processing notifications at the BOTTOM-LEFT (default is top-right). */
            #shiny-notification-panel {
                top: auto !important;
                bottom: 12px !important;
                left: 12px !important;
                right: auto !important;
            }

            /* ============================================================
               16S VISUALIZER SPECIFIC STYLES
               ============================================================ */
            
            /* Metric display boxes */
            .metric-box {
                background: linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%);
                border-radius: 10px;
                padding: 1rem;
                text-align: center;
                border: 1px solid #dee2e6;
                transition: transform 0.2s ease;
            }
            
            .metric-box:hover {
                transform: translateY(-2px);
            }
            
            .metric-value {
                font-size: 1.75rem;
                font-weight: 700;
                color: #2A9D8F;
                line-height: 1.2;
            }
            
            .metric-label {
                font-size: 0.8rem;
                color: #6c757d;
                text-transform: uppercase;
                letter-spacing: 0.5px;
            }
            
            /* Explanation sections */
            .explanation-section {
                background: #f8f9fa;
                border-radius: 10px;
                padding: 1.25rem;
                margin-top: 1rem;
            }
            
            .explanation-section h5 {
                color: #264653;
                font-weight: 700;
                margin-bottom: 0.75rem;
            }
            /* Help is intentionally long-form. The outer accordion is the only
               pane; headings and whitespace provide hierarchy inside it. */
            #vis_help_accordion .explanation-section,
            #exp_help_accordion .explanation-section,
            #val_help_accordion .explanation-section,
            #report_help_accordion .explanation-section {
                background: transparent;
                border-radius: 0;
                padding: 0;
                margin-top: 1.4rem;
            }
            #vis_help_accordion .explanation-section:first-child,
            #exp_help_accordion .explanation-section:first-child,
            #val_help_accordion .explanation-section:first-child,
            #report_help_accordion .explanation-section:first-child { margin-top: 0; }
            .help-longform > h5 {
                color: #264653;
                border-bottom: 1px solid #dce3e6;
                padding-bottom: 0.4rem;
                margin-top: 1.5rem;
            }
            #vis_help_accordion .alert,
            #exp_help_accordion .alert,
            #val_help_accordion .alert,
            #report_help_accordion .alert {
                background: transparent;
                color: var(--bs-body-color) !important;
                border-width: 0 0 0 3px;
                border-radius: 0;
                padding: 0.25rem 0 0.25rem 0.8rem !important;
            }
            #report_help_accordion .reference-item {
                background: transparent;
                box-shadow: none;
                border-left: 0;
                border-bottom: 1px solid #edf0f1;
                border-radius: 0;
                padding: 0.45rem 0;
            }
            /* Info cards */
            .info-card {
                background: white;
                border-radius: 10px;
                padding: 1rem;
                box-shadow: 0 2px 8px rgba(0,0,0,0.06);
                border-left: 3px solid #2A9D8F;
                margin-bottom: 1rem;
            }
            
            /* Reference citations */
            .reference-item {
                padding: 0.5rem 0.75rem;
                background: white;
                border-radius: 6px;
                margin-bottom: 0.5rem;
                border-left: 3px solid #3A86FF;
            }

            .page-intro { margin-bottom: 0.85rem; }
            .page-intro h4 { color: #264653; margin-bottom: 0.15rem; }
            .page-intro p { color: #66747b; margin-bottom: 0; }

            .empty-state {
                min-height: 360px;
                display: flex;
                flex-direction: column;
                align-items: center;
                justify-content: center;
                text-align: center;
                color: #6c757d;
                background: repeating-linear-gradient(135deg, #fff, #fff 12px,
                            #fafbfb 12px, #fafbfb 24px);
                border: 1px dashed #cbd4d8;
                border-radius: 0.5rem;
                padding: 2rem;
            }
            .empty-state .fa, .empty-state .fas { font-size: 2rem; color: #95a5a6; }

            .status-strip {
                display: flex;
                flex-wrap: wrap;
                gap: 0.5rem;
                margin-bottom: 0.75rem;
            }
            .status-chip {
                background: #eef2f3;
                border: 1px solid #d9e0e3;
                border-radius: 999px;
                padding: 0.3rem 0.65rem;
                font-size: 0.78rem;
            }
            
            /* ============================================================
               DADA2 PARAMETER EXPLORER SPECIFIC STYLES
               ============================================================ */
            
            /* Console styling for real-time feedback. Hosted inside the
               page-level global console sidebar (see the `sidebar`
               argument of page_navbar() below), so it keeps the original
               tall/narrow sidebar sizing. */
            .console-container {
                background: #1e1e1e;
                color: #d4d4d4;
                font-family: 'Consolas', 'Monaco', monospace;
                font-size: 11px;
                padding: 0.5rem;
                border-radius: 4px;
                box-sizing: border-box;
                width: 100%;
                max-width: 100%;
                height: clamp(320px, calc(100vh - 500px), 480px);
                overflow-y: auto;
                overflow-x: hidden;
                white-space: normal;
                overflow-wrap: anywhere;
                word-break: break-word;
            }

            /* Form controls within cards */
            .card .form-control, .card .form-select {
                font-size: 0.9rem;
            }
            /* Long project paths must wrap inside the 300px parameter margin
               instead of creating a distracting horizontal scrollbar. */
            .bslib-sidebar-layout > .sidebar { overflow-x: hidden; }
            #exp_selected_dir, #exp_selected_dir pre {
                white-space: normal !important;
                overflow-wrap: anywhere;
                word-break: break-word;
                max-width: 100%;
            }
            
            /* Compact accordion panels */
            .accordion-button {
                padding: 0.75rem 1rem;
            }

            /* Keep the Export preview within the main-panel width while allowing
               long paths and descriptions to wrap. */
            #report_table_preview table.dataTable {
                table-layout: fixed;
                width: 100% !important;
            }
            #report_table_preview table.dataTable td,
            #report_table_preview table.dataTable th {
                white-space: normal;
                word-wrap: break-word;
                overflow-wrap: break-word;
            }

            /* Retained-reads table: DataTables' automatic column sizing can
               make the generated table a few pixels wider than its output
               container, which is enough for bslib's card body to show a
               redundant horizontal scrollbar. Keep this four-column table
               pinned to the available width; its existing in-cell bars and
               column sizing remain intact. */
            #exp_sample_table {
                width: 100% !important;
                max-width: 100%;
            }
            #exp_sample_table table.dataTable {
                table-layout: fixed;
                width: 100% !important;
            }

            /* Keep the Validate results table readable beside its sidebar. */
            #optim_validation_table table.dataTable {
                table-layout: fixed;
                width: 100% !important;
            }
            #optim_validation_table table.dataTable td,
            #optim_validation_table table.dataTable th {
                white-space: normal;
                word-wrap: break-word;
                overflow-wrap: break-word;
                font-size: 0.78rem;
            }

            /* Join the Select sidebar panels into one compact accordion. */
            #exp_controls_accordion.accordion > .accordion-item {
                margin-bottom: 0;
            }

            /* Give collapsed Help accordions a consistent visual identity. */
            #vis_help_accordion .accordion-button.collapsed,
            #exp_help_accordion .accordion-button.collapsed,
            #val_help_accordion .accordion-button.collapsed,
            #report_help_accordion .accordion-button.collapsed {
                background-color: #eaf4fb;
                color: #0b5583;
            }
            #vis_help_accordion .accordion-button.collapsed:hover,
            #exp_help_accordion .accordion-button.collapsed:hover,
            #val_help_accordion .accordion-button.collapsed:hover,
            #report_help_accordion .accordion-button.collapsed:hover {
                background-color: #dcedf7;
            }

            /* ============================================================
               SHARED NAVIGATION STYLES
               ============================================================ */
            
            /* Active nav tab highlighting */
            .nav-pills .nav-link.active {
                background: linear-gradient(135deg, #264653 0%, #2A9D8F 100%);
            }
            
            /* Checkbox styling */
            .form-check-input:checked {
                background-color: #2A9D8F;
                border-color: #2A9D8F;
            }

            /* ============================================================
               TOP-LEVEL NAVBAR: always-visible tabs
               ============================================================ */
            /* Keep compact desktop tabs while preserving Bootstrap's
               responsive navbar collapse at narrower widths. */
            .navbar-nav .nav-link {
                font-size: 0.8rem;
                padding-left: 0.6rem !important;
                padding-right: 0.6rem !important;
                white-space: nowrap;
            }
            #main_nav > .nav-item {
                display: flex !important;
                align-items: center;
            }
            #main_nav .workflow-separator {
                color: rgba(255,255,255,0.55);
                padding: 0 0.15rem;
                font-weight: 700;
            }
            /* Keep the workflow navigation reachable on long pages. */
            .navbar {
                position: sticky;
                top: 0;
                z-index: 1030;
            }

            @media (max-width: 991.98px) {
                .navbar-collapse { max-height: calc(100vh - 4rem); overflow-y: auto; }
                .navbar-nav .nav-link { padding: 0.55rem 0.75rem !important; }
            }

            @media (max-width: 575.98px) {
                .metric-value { font-size: 1.35rem; }
                .metric-label { font-size: 0.68rem; }
            }

        ")),
        # Per-direction Quality Profiles slider colors.
        # are tinted with the Forward / Reverse identity colors -- the filled bar
        # and the value bubble. Built with the FWD_COLOR / REV_COLOR constants.
        tags$style(HTML(paste0(
            ".qp-fwd .irs-bar,.qp-fwd .irs-single,.qp-fwd .irs-from,.qp-fwd .irs-to",
            "{background-color:", FWD_COLOR, " !important;border-color:", FWD_COLOR, " !important;}",
            ".qp-rev .irs-bar,.qp-rev .irs-single,.qp-rev .irs-from,.qp-rev .irs-to",
            "{background-color:", REV_COLOR, " !important;border-color:", REV_COLOR, " !important;}",
            # Slider-as-x-axis layout (slider below its plot). Each
            # control is a flex row [label][-][ track ][+]. The left/right blocks
            # are fixed-width (--ctrl-left / --ctrl-right, from the shared R
            # constants) and equal the plotly margins MINUS the handle inset, so the
            # slider track spans exactly the plot's data region. Driven by the same
            # numbers as the plot margins, they cannot drift apart on resize.
            ".axis-slider-chart{--ctrl-left:", QP_CTRL_LEFT, "px;--ctrl-right:", QP_CTRL_RIGHT, "px;}",
            # margin-top negative pulls the whole slider up so it sits as close as
            # possible under the plot's x-axis.
            ".axis-slider-row{display:flex;align-items:center;margin-top:-6px;margin-bottom:2px;}",
            # padding-right leaves a couple of px between the minus button and the
            # track start; the track (asr-mid) still begins at --ctrl-left, so the
            # alignment with the plot's data region is preserved.
            ".axis-slider-row .asr-left{flex:0 0 var(--ctrl-left);width:var(--ctrl-left);",
            "display:flex;align-items:center;justify-content:flex-end;gap:4px;padding-right:4px;}",
            ".axis-slider-row .asr-mid{flex:1 1 auto;min-width:0;}",
            # padding-left leaves the matching gap between the track end and the plus.
            ".axis-slider-row .asr-right{flex:0 0 var(--ctrl-right);width:var(--ctrl-right);",
            "display:flex;align-items:center;justify-content:flex-start;padding-left:4px;}",
            ".axis-slider-row .asr-mid .form-group{margin-bottom:0;}",
            ".axis-slider-row .asr-mid .shiny-input-container{width:100% !important;}",
            ".axis-slider-row .asr-label{font-size:0.78rem;font-weight:600;white-space:nowrap;}",
            ".axis-slider-row .btn{padding:0 0.4rem;line-height:1.25;}",
            # Hide the redundant min/max end labels above the track (the 0 / 250
            # grey boxes); the numeric grid BELOW the track is the visible scale,
            # and the value bubble stays above the handle.
            ".axis-slider-row .irs-min,.axis-slider-row .irs-max{display:none !important;}"
        )))
    ),

    # =========================================================================
    # 7.3 JavaScript Handlers for Real-Time Console
    # =========================================================================
    # Custom message handlers allow the server to push updates to the UI
    # without waiting for the reactive cycle to complete
    
    tags$head(
        tags$script(HTML("
            // Handler for appending console messages in real-time
            // Called from server via session$sendCustomMessage('console_msg', ...)
            Shiny.addCustomMessageHandler('console_msg', function(msg) {
                var container = document.getElementById('console_container');
                if (container) {
                    // Create new message element
                    var div = document.createElement('div');
                    div.style.color = msg.color;
                    div.style.marginBottom = '2px';
                    // Format: [timestamp] message
                    var stamp = document.createElement('span');
                    stamp.style.color = '#6c757d';
                    stamp.textContent = '[' + msg.time + '] ';
                    div.appendChild(stamp);
                    div.appendChild(document.createTextNode(String(msg.text)));
                    // Append and scroll to bottom
                    container.appendChild(div);
                    container.scrollTop = container.scrollHeight;
                }
            });
            
            // Handler for clearing console
            Shiny.addCustomMessageHandler('console_clear', function(msg) {
                var container = document.getElementById('console_container');
                if (container) {
                    container.innerHTML = '<div style=\"color:#6c757d;font-style:italic;\">Select a FASTQ folder, then load quality profiles.</div>';
                }
            });

            // Handler that relocates the single Processing Console accordion
            // node next to the sidebar anchor belonging to whichever tab
            // just became active (sent by the `main_nav` observeEvent on
            // the server -- see Section 9). Uses insertBefore against the
            // console_slot_<value> anchor -- never appendChild into it --
            // so the accordion always ends up a DIRECT child of that tab's
            // sidebar-content, exactly like the sidebar's own controls
            // accordion above it. That direct-child relationship is what
            // gives it bslib's edge-to-edge flush negative margin (see
            // the comment above the Visualize tab's Processing Console);
            // appending it inside the anchor div instead would leave it
            // one level too deep and visibly narrower than the sidebar.
            Shiny.addCustomMessageHandler('move_console', function(msg) {
                var panel = document.getElementById('exp_console_accordion');
                var anchor = document.getElementById(msg.slot);
                if (panel && anchor && anchor.parentElement) {
                    anchor.parentElement.insertBefore(panel, anchor);
                }
            });
        "))
        )
    ),

    # =========================================================================
    # 7.4 Visualize Tab
    # =========================================================================
    # Documents the amplified target and primer coordinates.
    
    nav_panel(
        title = tags$span(icon("dna"), " Visualize"),
        value = "visualizer",

        # The navbar label identifies the tab; no duplicate header is needed.

        # Main layout with sidebar for controls
        layout_sidebar(
            fillable = FALSE,
            
            # -----------------------------------------------------------------
            # Sidebar: amplified-target controls
            # -----------------------------------------------------------------
            sidebar = sidebar(
                width = 340,
                # The sidebar stays visible; its individual sections remain
                # independently collapsible.
                open = "always",

                # Collapsible accordion for organized controls
                accordion(
                    id = "vis_controls_accordion",
                    open = c("vis_region_panel", "vis_primer_panel", "vis_target_length_panel"),

                    # .............................................................
                    # Target Region Selection Panel
                    # .............................................................
                    accordion_panel(
                        title = "Amplified Region",
                        value = "vis_region_panel",
                        icon = icon("crosshairs"),

                        # Choices are the distinct targets represented in the
                        # primer database and filter the primer-pair choices.
                        p(class = "text-muted small",
                          "Select your target region(s)."),

                        # Single-select (was a 9-checkbox checkboxGroupInput
                        # letting users pick any individual V1-V9
                        # combination) -- a "target" as the Primer Database
                        # defines it is always exactly one of these 7
                        # values, never an arbitrary combination, so a
                        # single choice is the correct control here.
                        selectInput(
                            "vis_target_region_select",
                            label = NULL,
                            choices = VIS_TARGET_REGION_CHOICES,
                            selected = DEFAULT_TARGET_REGION
                        )
                    ),
                    
                    # .............................................................
                    # Primer Pair Panel
                    # .............................................................
                    accordion_panel(
                        title = "Primer Pair",
                        value = "vis_primer_panel",
                        icon = icon("arrows-left-right-to-line"),

                        # Published pairs are the primary path; manual fields
                        # remain directly below for custom coordinates.
                        p(class = "text-muted small",
                          "Select a primer pair or enter coordinates below."),

                        h6(icon("flask"), " Literature Primers:", class = "mb-2"),
                        # Choices restricted to primer pairs matching the
                        # selected Target Region Selection value
                        # (Section 4.7's vis_primer_choices_for_target()). The
                        # server-side observer for
                        # input$vis_target_region_select (Section 8.5)
                        # keeps these choices in sync whenever the target
                        # region changes; initial choices here just need to
                        # match DEFAULT_TARGET_REGION so the two panels
                        # agree on app startup.
                        selectInput(
                            "vis_primer_preset_select",
                            label = NULL,
                            choices = vis_primer_choices_for_target(DEFAULT_TARGET_REGION),
                            # Startup selection; its database row fills all
                            # coordinate, primer-length, and target-length fields.
                            selected = DEFAULT_PRIMER_PAIR
                        ),
                        # Manual coordinates are normal, always-visible
                        # parameters again. Selecting a literature pair fills
                        # them automatically, while they remain editable for
                        # assays with verified custom coordinates.
                        hr(class = "my-3"),
                        h6("Primer Coordinates", class = "mb-1"),
                        p(class = "text-muted small mb-2",
                          "Fields auto-fill when primers are selected, but values can be entered manually!"),

                        # Numeric inputs for primer positions
                        layout_columns(
                            col_widths = c(6, 6),
                            numericInput(
                                "vis_forward_primer_start",
                                label = tags$span(icon("arrow-right", style = sprintf("color:%s;", FWD_COLOR)), " Fwd start"),
                                value = 341, min = 1, max = 1500, step = 1
                            ),
                            numericInput(
                                "vis_reverse_primer_end",
                                label = tags$span(icon("arrow-left", style = sprintf("color:%s;", REV_COLOR)), " Rev end"),
                                value = 806, min = 1, max = 1542, step = 1
                            )
                        ),

                        # Primer LENGTH inputs (bp), added alongside the
                        # position inputs above so the Gene Map can subtract
                        # each primer's own footprint from the primer-to-
                        # primer span and show the true amplicon length left
                        # once cutadapt (Step 3) trims both primers off --
                        # not just the raw genomic distance between them.
                        # These four numericInputs' own R-level defaults
                        # (341/806/17/20 bp) are only ever visible for an
                        # instant: the Literature Primers dropdown above
                        # defaults to "341F / 785R" (V3-V4), and selecting
                        # any preset -- including the default one, on
                        # startup -- immediately overwrites all four via
                        # vis_update_primers_and_regions() with that pair's
                        # real, verified start/end/lengths. Left as generic
                        # round numbers here rather than duplicating
                        # 785R's real values, to avoid two places that
                        # could drift out of sync with primer_database.
                        layout_columns(
                            col_widths = c(6, 6),
                            numericInput(
                                "vis_forward_primer_length",
                                label = tags$span(icon("ruler-horizontal", style = sprintf("color:%s;", FWD_COLOR)), " Fwd length"),
                                value = 17, min = 1, max = 50, step = 1
                            ),
                            numericInput(
                                "vis_reverse_primer_length",
                                label = tags$span(icon("ruler-horizontal", style = sprintf("color:%s;", REV_COLOR)), " Rev length"),
                                value = 20, min = 1, max = 50, step = 1
                            )
                        )

                    ),

                    # .............................................................
                    # Target Length Panel
                    # .............................................................
                    # Primer-trimmed post-merge length bounds are distinct from
                    # reference coordinates and have their own section.
                    accordion_panel(
                        title = "Target Length",
                        value = "vis_target_length_panel",
                        icon = icon("ruler-horizontal"),

                        # These inputs hold the primer-trimmed insert
                        # length filter bounds (amplicon_min_length/
                        # amplicon_max_length, Section 4.3) that
                        # 5_dada2_pipeline.Rmd applies AFTER mergePairs() to
                        # discard chimeric/off-target merges. Same
                        # auto-fill-but-editable pattern as the Fwd/Rev
                        # start/end/length inputs in the Primer Pair
                        # panel above: selecting a literature primer pair
                        # overwrites both fields via
                        # vis_update_primers_and_regions() (Section 8.5)
                        # with that pair's real, in-silico-PCR-derived
                        # bounds (see the Primer Database's Target Min
                        # (bp)/Target Max (bp) columns, or the methodology
                        # write-up below the table on this tab's own Help
                        # pane), but the values stay fully editable
                        # afterward for a manually-entered primer pair, or
                        # to override a literature pair's bounds with the
                        # user's own knowledge of their sample.
                        # amplicon_length_range() reads both inputs directly.
                        p(class = "text-muted small mb-1",
                          "Select the expected primer-trimmed length range."),
                        p(class = "text-muted small mb-2",
                          "Fields auto-fill when primers are selected, but values can be entered manually!"),
                        layout_columns(
                            col_widths = c(6, 6),
                            numericInput(
                                "vis_target_length_min",
                                label = tags$span(icon("ruler-horizontal", style = sprintf("color:%s;", FWD_COLOR)), " Min (bp)"),
                                # Generic V3-V4-ish placeholder, deliberately
                                # not 341F/785R's real value (399) -- same
                                # rationale as the Fwd/Rev length defaults
                                # in the Primer Pair panel above: this
                                # value is only ever visible for an instant
                                # before the Literature Primers dropdown's
                                # startup selection overwrites it.
                                value = 380, min = 1, max = 2000, step = 1
                            ),
                            numericInput(
                                "vis_target_length_max",
                                label = tags$span(icon("ruler-horizontal", style = sprintf("color:%s;", REV_COLOR)), " Max (bp)"),
                                value = 460, min = 1, max = 2000, step = 1
                            )
                        )
                    )
                ),

                # -------------------------------------------------------------
                # Processing Console
                # -------------------------------------------------------------
                # Authored here (rather than in Select's sidebar, where it
                # used to live) because Visualize is the default active tab
                # on load, so it needs to render correctly here with no JS
                # relocation required. Deliberately a DIRECT child of this
                # sidebar (not wrapped in a placeholder div) because bslib
                # only gives a sidebar's own accordions its edge-to-edge
                # "flush" negative margin via the CSS rule
                # `.sidebar-content > .accordion { margin: 0 calc(-1 *
                # var(--_padding)) }` -- a wrapper div in between breaks that
                # direct-child match and leaves the accordion visibly
                # narrower than its sibling controls accordion. When another
                # tab becomes active, the `main_nav` observeEvent in the
                # server tells the client (via the `move_console` message
                # handled in `header` below) to move this entire
                # #exp_console_accordion node to sit directly before that
                # tab's own console_slot_<value> anchor -- keeping it a
                # direct sidebar child (and therefore correctly flush) in
                # every tab, while the same DOM node (and #console_container
                # inside it) simply relocates, so console history is
                # preserved across tab switches and no ID is ever
                # duplicated.
                accordion(
                    id = "exp_console_accordion",
                    open = "exp_console_panel",
                    accordion_panel(
                        title = "Processing Details",
                        value = "exp_console_panel",
                        icon = icon("terminal"),
                        div(class = "d-flex justify-content-end mb-1",
                            actionButton("exp_clear_console", "Clear",
                                         class = "btn-outline-secondary btn-sm py-0")),
                        div(id = "console_container", class = "console-container",
                            div(style = "color:#6c757d;font-style:italic;",
                                "Select a FASTQ folder, then load quality profiles."))
                    )
                ),

                # Zero-size anchor marking where the console accordion above
                # belongs when this tab is active again. The `move_console`
                # handler in `header` re-inserts #exp_console_accordion
                # immediately before this marker (never appends into it), so
                # the marker itself never becomes the accordion's parent --
                # it stays a direct sidebar child throughout.
                tags$div(id = "console_slot_visualizer")
            ),

            # -----------------------------------------------------------------
            # Main panel: amplified-target visualization
            # -----------------------------------------------------------------
            tagList(
                div(class = "page-intro",
                    h4("Visualize the amplified 16S target"),
                    p("View the amplified region, primer coordinates, and expected primer-trimmed target length on the 16S rRNA gene.")),
                # Main 16S gene map -- shown first so the map itself
                # (the primary content of this tab) is immediately
                # visible, with the summary metrics for the current
                # selection reinforcing it right below.
                card(
                    card_header(tags$span(icon("dna"), " 16S rRNA Gene Map"), class = "bg-light"),
                    card_body(plotlyOutput("vis_main_plot", height = "300px"))
                ),

                br(),

                # Assay coordinates and primer-trimmed target-length summary.
                layout_columns(
                    col_widths = c(4, 4, 4),

                    div(class = "metric-box",
                        div(class = "metric-value", textOutput("vis_amplicon_size", inline = TRUE)),
                        # Label clarifies this is the primer-free length
                        # (excludes both primers' own footprints) -- see
                        # output$vis_amplicon_size above -- since the Gene
                        # Map below shows the primer regions as separate
                        # colored segments, outside this number.
                        div(class = "metric-label", "Amplicon (bp, w/o primers)")
                    ),
                    div(class = "metric-box",
                        div(class = "metric-value", textOutput("vis_regions_covered", inline = TRUE)),
                        div(class = "metric-label", "Regions")
                    ),
                    div(class = "metric-box",
                        div(class = "metric-value", textOutput("vis_coverage_pct", inline = TRUE)),
                        div(class = "metric-label", "16S Coverage %")
                    )
                ),

                br(),

                # Supporting biology and primer-database reference. Collapsed by
                # default because the Gene Map is the tab's primary content and
                # the app documents an existing assay rather than designing one.
                accordion(
                    id = "vis_help_accordion",
                    open = FALSE,
                    accordion_panel(
                        title = "How to interpret the visualization",
                        icon = icon("question-circle"),

                        div(class = "help-longform mt-0",
                            h5("What this tab does"),

                        div(class = "explanation-section",
                            h5("Purpose"),
                            p("The Visualize tab places the known primers and resulting amplified target on the ",
                              tags$em("Escherichia coli"), " 16S rRNA reference sequence (GenBank J01859). It shows which variable regions are covered, the primer-to-primer span, and the expected primer-trimmed insert length. These coordinates provide the amplicon context used by the overlap assessment in Select. The tab does not design primers, recommend an assay, or estimate taxonomic specificity.")
                        ),

                        div(class = "explanation-section",
                            h5("Quick guide"),
                            tags$ol(class = "mb-2",
                                tags$li("Select the amplified variable region used by the assay."),
                                tags$li("Select its primer pair, or enter the known coordinates directly."),
                                tags$li("Confirm the expected primer-trimmed target-length range."),
                                tags$li("Read the primer positions, amplified span, and variable-region coverage from the gene map.")),
                            p(class = "small text-muted mb-0",
                              "This tab documents and visualizes an existing assay. It does not design or recommend an experiment.")),

                        h5("Scientific background"),

                        div(class = "explanation-section",
                            h5("16S rRNA variable regions"),
                            p("The 16S ribosomal RNA (rRNA) gene is a highly conserved component of
                              the bacterial ribosome, essential for protein synthesis. This ~1,500 bp
                              gene contains ", tags$strong("nine variable regions (V1-V9)"), "
                              interspersed with conserved regions. This mosaic structure makes it
                              useful for bacterial identification and phylogenetic analysis."),
                            # Coordinates use one E. coli reference frame; biological
                            # insertions/deletions explain the target-length ranges.
                            p(class = "text-muted small fst-italic",
                              "Note: all region and primer-binding coordinates shown in this app",
                              " (Gene Map, Primer Database) are given in standard ",
                              tags$em("Escherichia coli"), " 16S rRNA gene reference numbering",
                              " (GenBank J01859). Actual amplicon length still varies somewhat",
                              " between different bacterial species/strains, due to small",
                              " insertions and deletions within the variable region(s) -- see",
                              " \"How were Target Min (bp) / Target Max (bp) derived?\" below the",
                              " Primer Database table for how this variation is quantified."),
                            tagList(
                                div(
                                    h6("Variable regions (V1-V9)"),
                                    tags$ul(
                                        tags$li("Evolve rapidly, accumulating mutations"),
                                        tags$li("Can provide genus- or, for some taxa and regions, species-level discrimination"),
                                        tags$li("Different regions offer varying taxonomic resolution"),
                                        tags$li("V4 and V3-V4 are widely used targets in short-read surveys")
                                    )
                                ),
                                div(
                                    h6("Conserved regions"),
                                    tags$ul(
                                        tags$li("Relatively conserved across broad bacterial groups"),
                                        tags$li("Essential for ribosome function"),
                                        tags$li("Serve as universal primer binding sites"),
                                        tags$li("Enable PCR amplification across diverse taxa")
                                    )
                                )
                            )
                        ),

                        div(class = "explanation-section",
                            h5("How amplicon sequencing works"),
                            tags$ol(
                                tags$li(tags$strong("Primer Design: "), "Primers are designed to bind
                                        conserved regions flanking the variable region(s) of
                                        interest."),
                                tags$li(tags$strong("PCR Amplification: "), "The target region is
                                        amplified from all bacteria in the sample, creating millions
                                        of copies."),
                                tags$li(tags$strong("Sequencing: "), "Amplicons are sequenced using
                                        paired-end reads (e.g., 2", tags$span(HTML("&times;")), "250 bp
                                        on Illumina MiSeq)."),
                                tags$li(tags$strong("Read Merging: "), "Forward and reverse reads are
                                        merged using the overlapping region to create full-length
                                        amplicon sequences -- the basis for the ",
                                        tags$strong("Select"), " tab's overlap assessment."),
                                tags$li(tags$strong("Taxonomy Assignment: "), "Sequences are compared
                                        to reference databases (SILVA, Greengenes, RDP) for taxonomic
                                        classification.")
                            )
                        ),

                        div(class = "explanation-section",
                            h5("Primer database"),
                            p(class = "text-muted small",
                              "Published primer pairs available in the ",
                              tags$strong("Literature Primers"), " dropdown above, with their
                              binding positions, primer sequences, target region(s), amplicon
                              size, primer-trimmed target length filter bounds, and source
                              citation."),
                            div(class = "table-responsive", tableOutput("vis_primer_table_help")),
                            uiOutput("vis_primer_table_footnote_help"),
                            p(class = "text-muted small mt-2 mb-0",
                              "Sequences use standard IUPAC degenerate-base codes (e.g. M=A/C, W=A/T,",
                              " R=A/G, Y=C/T, N=any, V=A/C/G, H=A/C/T, K=G/T) exactly as published,",
                              " and were individually verified against each primer pair's originating",
                              " publication or an established protocol/database. A superscript",
                              " letter (", tags$span(style = "color:#c77700;", tags$sup("a"), ", ",
                              tags$sup("b"), ", ..."), ") next to a sequence means it has a caveat --",
                              " each letter is a distinct caveat, listed in its own footnote below",
                              " the table."),

                            # These bounds auto-populate the Target Length panel
                            # when a literature primer pair is selected.
                            div(class = "alert alert-info py-2 mt-2 mb-0",
                                tags$strong("How were Target Min (bp) / Target Max (bp) derived?"),
                                p(class = "small mb-1 mt-1",
                                  "These are ", tags$strong("primer-trimmed insert lengths"),
                                  " -- the sequence length DADA2 sees in ",
                                  tags$code("colnames(seqtab)"), " after cutadapt has removed the",
                                  " primers -- not the raw primer-to-primer span (",
                                  tags$code("Amplicon (bp)"), " column above)."),
                                p(class = "small mb-1",
                                  "For each primer pair, ", tags$em("in silico"), " PCR was run",
                                  " against four reference databases:"),
                                tags$ul(class = "small mb-1",
                                    tags$li("SILVA 138.2 SSURef NR99 (~200k bacterial sequences)"),
                                    tags$li("A phylogenetically stratified NCBI sample (500 sequences",
                                            " x 10 major phyla = 5,000 total)"),
                                    tags$li("GTDB r220 bacterial species representatives (66k",
                                            " full-length sequences)"),
                                    tags$li("NCBI RefSeq Targeted Loci 16S bacteria (26k sequences)")
                                ),
                                p(class = "small mb-0",
                                  "For each database, the primer-trimmed amplicon lengths were",
                                  " extracted, the 1st and 99th percentiles computed, and the filter",
                                  " bounds set as p1", tags$code("-2"), " (min) and p99",
                                  tags$code("+2"), " (max). The final bounds are the ",
                                  tags$strong("union"), " across all four sources -- the lowest",
                                  " minimum and the highest maximum from any single database. This",
                                  " means a legitimate amplicon captured by any one database will",
                                  " not be discarded.")
                            ),

                            # Explain the 8F naming next to the affected rows.
                            div(class = "alert alert-info py-2 mt-2 mb-0",
                                tags$strong("Why \"8F\" instead of the more familiar \"27F\"?"),
                                " The classic universal forward primer for the 16S V1 region has two names in",
                                " circulation for the exact same 20-nt binding site -- E. coli J01859",
                                " positions 8-27: ", tags$strong("\"27F\""), " (Lane 1991-era name, numbering the",
                                " primer by its 3'-most base) and ", tags$strong("\"8F\""),
                                " (Turner et al. 1999-era name, numbering it by its 5'-most base instead). \"27F\"",
                                " is by far the more commonly seen name in the literature. This table labels the",
                                " primer \"8F\" in the two rows above, deliberately, so the name always matches",
                                " the numeric ", tags$code("Fwd Start (bp)"), " column next to it -- exactly as it",
                                " already does for every other forward primer in this table (341F starts at 341,",
                                " 515F at 515, 784F at 784, 926F at 926, 967F at 967). Using \"27F\" here would have",
                                " made it look, incorrectly, like the table had an off-by-19 error. If your own",
                                " protocol or prior work refers to this primer as \"27F\", know that it is the",
                                " same oligo pairing (8F/338R and 8F/519R above correspond to what is usually",
                                " written as 27F/338R and 27F/519R). Source: Frank et al. 2008 (",
                                tags$em("Appl Environ Microbiol"), " 74(8):2461-2470), which documents",
                                " \"27f\" as spanning E. coli positions 8 to 27.")
                        )
                        )
                    )
                )
            )
        )
    ),
    nav_item(tags$span(class = "workflow-separator", ">")),

    # =========================================================================
    # 7.5 Select Tab
    # =========================================================================
    # FASTQ discovery, quality profiles, filtering controls, and retention.

    nav_panel(
        title = tags$span(icon("sliders-h"), " Select"),
        value = "explorer",

        # Natural-height layout prevents nested page/console scrollbars.
        layout_sidebar(
            fillable = FALSE,

            # -----------------------------------------------------------------
            # Sidebar: data and processing controls
            # -----------------------------------------------------------------
            sidebar = sidebar(
                width = 320,
                padding = c(10, 10),

                # Keep the sidebar visible and use only the section chevrons.
                open = "always",

                # All Select controls start open. Processing Details/console
                # now lives in the global console bar (see the `header`
                # argument of page_navbar() above), not in this accordion.
                accordion(
                    id = "exp_controls_accordion",
                    open = c("exp_platform_panel", "exp_data_panel"),

                # .............................................................
                # Sequencing Platform Panel
                # .............................................................
                # The sequencing kit is selected once and caps both truncLen
                # controls at its nominal read length.
                accordion_panel(
                    title = "Read Length",
                    value = "exp_platform_panel",
                    icon = icon("satellite-dish"),

                    selectInput("vis_platform", label = "Sequencing kit",
                                choices = names(PLATFORM_READ_LENGTHS),
                                selected = "MiSeq v3 (2×300)", width = "100%")
                ),

                # .............................................................
                # FASTQ Data Panel
                # .............................................................
                accordion_panel(
                    title = "FASTQ Data",
                    value = "exp_data_panel",
                    icon = icon("folder-open"),

                    tags$label(class = "form-label fw-semibold", "Primer-trimmed FASTQ folder"),
                    p(class = "small text-muted mb-2",
                      "The default project folder is selected automatically; choose another folder if needed."),
                    shinyDirButton("exp_dir_select", "Choose folder...",
                                   "Select FASTQ Directory",
                                   class = "btn-outline-primary btn-sm w-100"),
                    tags$div(class = "form-text mt-2", "Selected folder"),
                    verbatimTextOutput("exp_selected_dir", placeholder = TRUE),

                    selectInput("exp_pattern_choice", "Forward / reverse filename suffix",
                                choices = names(PATTERN_CHOICES),
                                selected = names(PATTERN_CHOICES)[1]),

                    # Load button — bottom of panel
                    actionButton("exp_load_data", "Load quality profiles",
                                 class = "btn-primary btn-sm w-100 mt-3",
                                 icon = icon("chart-line"))
                )

                ),

                # Zero-size anchor marking where the Processing Console
                # accordion belongs when this tab is active. The panel
                # itself is authored once, in the Visualize tab's sidebar
                # (the default active tab), and relocated to sit directly
                # before this anchor by JS whenever this tab becomes active
                # -- see the `main_nav` observeEvent in the server and the
                # `move_console` message handler in `header` below. The
                # accordion is inserted as a direct sidebar child (never
                # appended inside this anchor div), which is what keeps it
                # visually flush with the sidebar -- see the comment above
                # the Visualize tab's Processing Console for why that
                # matters. Keeps this exactly where "Processing Details"
                # used to sit (bottom of this sidebar) while staying visible
                # on every tab, without duplicating #console_container.
                tags$div(id = "console_slot_explorer")
            ),
            
            # -----------------------------------------------------------------
            # Main panel: quality analysis and parameter controls
            # -----------------------------------------------------------------
            tagList(
                div(class = "page-intro",
                    h4("Choose filtering parameters"),
                    p("Load paired FASTQs, trim poor-quality read tails, and preserve enough overlap for reliable merging.")),
            div(id = "exp_empty_state", class = "empty-state",
                icon("chart-line"),
                h5(class = "mt-3", "No quality data loaded"),
                p(class = "mb-1", "Choose the primer-trimmed FASTQ folder and click Load quality profiles."),
                p(class = "small mb-0", "Quality plots and parameter controls will appear when processing completes.")),
            div(id = "exp_analysis_content", style = "display:none;",
            uiOutput("exp_status_summary"),
            navset_card_pill(
                id = "exp_main_tabs",

                # Quality profiles and live filtering controls
                nav_panel(
                    title = tags$span(icon("chart-area"), " Quality Profiles"),

                    # Forward and reverse quality/retention panes flank a compact
                    # three-bar retained-reads summary.
                    layout_columns(
                        col_widths = c(5, 2, 5),

                        card(
                            card_header(class = "text-white py-2",
                                        style = sprintf("background-color:%s;", FWD_COLOR),
                                        "Forward Reads"),
                            card_body(
                                class = "p-2",
                                # Quality plot with truncLen as its aligned x-axis:
                                # the row is
                                # [label][-][ track ][+]; the fixed-width left/right
                                # blocks are sized (--ctrl-left/right) to the plot's
                                # margins so the track spans exactly the data region.
                                div(class = "axis-slider-chart qp-fwd",
                                    plotlyOutput("exp_plot_fwd_quality", height = "240px"),
                                    div(class = "axis-slider-row",
                                        div(class = "asr-left",
                                            tags$span("truncLen", class = "asr-label"),
                                            actionButton("exp_truncLenF_dec", "−",
                                                         class = "btn btn-outline-secondary btn-sm")),
                                        div(class = "asr-mid",
                                            sliderInput("exp_truncLenF", NULL, min = 0, max = 300,
                                                        value = 200, step = 1, ticks = TRUE, width = "100%")),
                                        div(class = "asr-right",
                                            actionButton("exp_truncLenF_inc", "+",
                                                         class = "btn btn-outline-secondary btn-sm")))
                                ),
                                # Retention plot with its maxEE slider below it, same
                                # aligned [label][-][ track ][+] layout.
                                div(class = "axis-slider-chart qp-fwd",
                                    plotlyOutput("exp_plot_fwd_retmaxee", height = "200px"),
                                    div(class = "axis-slider-row",
                                        div(class = "asr-left",
                                            tags$span("maxEE", class = "asr-label"),
                                            actionButton("exp_maxEEF_dec", "−",
                                                         class = "btn btn-outline-secondary btn-sm")),
                                        div(class = "asr-mid",
                                            sliderInput("exp_maxEEF", NULL, min = 0, max = 10,
                                                        value = 2, step = 0.1, ticks = TRUE, width = "100%")),
                                        div(class = "asr-right",
                                            actionButton("exp_maxEEF_inc", "+",
                                                         class = "btn btn-outline-secondary btn-sm")))
                                )
                            )
                        ),

                        card(
                            card_header(class = "bg-secondary text-white py-2", "Retained reads"),
                            card_body(
                                class = "p-2 d-flex align-items-center",
                                plotlyOutput("exp_retbar_plot", height = "470px")
                            )
                        ),

                        card(
                            card_header(class = "text-white py-2",
                                        style = sprintf("background-color:%s;", REV_COLOR),
                                        "Reverse Reads"),
                            card_body(
                                class = "p-2",
                                # Mirror of the Forward card: reverse quality plot with
                                # its truncLen slider below, reverse retention plot with
                                # its maxEE slider below -- both aligned x-axes.
                                div(class = "axis-slider-chart qp-rev",
                                    plotlyOutput("exp_plot_rev_quality", height = "240px"),
                                    div(class = "axis-slider-row",
                                        div(class = "asr-left",
                                            tags$span("truncLen", class = "asr-label"),
                                            actionButton("exp_truncLenR_dec", "−",
                                                         class = "btn btn-outline-secondary btn-sm")),
                                        div(class = "asr-mid",
                                            sliderInput("exp_truncLenR", NULL, min = 0, max = 300,
                                                        value = 200, step = 1, ticks = TRUE, width = "100%")),
                                        div(class = "asr-right",
                                            actionButton("exp_truncLenR_inc", "+",
                                                         class = "btn btn-outline-secondary btn-sm")))
                                ),
                                div(class = "axis-slider-chart qp-rev",
                                    plotlyOutput("exp_plot_rev_retmaxee", height = "200px"),
                                    div(class = "axis-slider-row",
                                        div(class = "asr-left",
                                            tags$span("maxEE", class = "asr-label"),
                                            actionButton("exp_maxEER_dec", "−",
                                                         class = "btn btn-outline-secondary btn-sm")),
                                        div(class = "asr-mid",
                                            sliderInput("exp_maxEER", NULL, min = 0, max = 10,
                                                        value = 2, step = 0.1, ticks = TRUE, width = "100%")),
                                        div(class = "asr-right",
                                            actionButton("exp_maxEER_inc", "+",
                                                         class = "btn btn-outline-secondary btn-sm")))
                                )
                            )
                        )
                    ),
                    
                    # Read coverage and overlap at the current truncLen values.
                    card(
                        card_header(class = "bg-secondary text-white py-2",
                                    tags$span(icon("arrows-left-right"),
                                              " Paired-Read Overlap at Current truncLen")),
                        card_body(
                            class = "p-2",
                            # Fixed-width status column beside the flexible plot.
                            div(class = "d-flex", style = "gap:10px; align-items:center;",
                                div(style = "flex:0 0 168px;",
                                    uiOutput("exp_overlap_msg_left")),
                                div(style = "flex:1 1 auto; min-width:0;",
                                    plotlyOutput("exp_overlap_plot", height = "260px")))
                        )
                    )
                ),

                # Per-sample retained reads
                nav_panel(
                    title = tags$span(icon("table"), " Retained reads"),

                    card(
                        card_header(
                            class = "bg-secondary text-white py-2",
                            div(class = "d-flex justify-content-between align-items-center",
                                # me-2 guarantees a visible gap before the
                                # sample-count badge even if the flex
                                # justify-content-between spacing reads as
                                # too tight against the badge in practice.
                                span("Retained reads", class = "me-2"),
                                span(class = "badge bg-light text-dark",
                                     textOutput("exp_n_samples", inline = TRUE))
                            )
                        ),
                        card_body(
                            class = "p-2",
                            # Fill available viewport space without capping tall
                            # tables; the page itself scrolls when needed.
                            style = "min-height: calc(100vh - 300px);",

                            # Aggregate retention summary above the per-sample table.
                            fluidRow(class = "mb-2",
                                column(4, div(class = "text-center",
                                    strong("Forward", class = "small d-block"),
                                    span(textOutput("exp_retention_F", inline = TRUE),
                                         class = "fs-5 fw-bold",
                                         style = sprintf("color:%s;", FWD_COLOR)))),
                                column(4, div(class = "text-center",
                                    strong("Reverse", class = "small d-block"),
                                    span(textOutput("exp_retention_R", inline = TRUE),
                                         class = "fs-5 fw-bold",
                                         style = sprintf("color:%s;", REV_COLOR)))),
                                column(4, div(class = "text-center",
                                    strong("Paired", class = "small d-block"),
                                    span(textOutput("exp_retention_C", inline = TRUE),
                                         class = "fs-5 fw-bold",
                                         style = sprintf("color:%s;", PAIRED_COLOR))))
                            ),
                            hr(class = "my-2"),

                            # Natural height avoids an internal horizontal or
                            # vertical scroll container.
                            DTOutput("exp_sample_table", fill = FALSE)
                        )
                    )
                )
            )),

            br(),

            # Parameter definitions, overlap thresholds, and practical guidance.
            accordion(
                id = "exp_help_accordion",
                open = FALSE,
                accordion_panel(
                    title = "How to choose parameters",
                    icon = icon("question-circle"),

                    div(class = "help-longform mt-0",
                        h5("What this tab does"),

                    div(class = "explanation-section",
                        h5("Purpose"),
                        p("The Select tab uses observed forward- and reverse-read quality profiles to choose DADA2 filtering and truncation parameters. It estimates how many primer-trimmed reads satisfy the selected truncLen and maxEE thresholds, and checks whether the retained read lengths can overlap across the expected amplicon-length range. The estimate supports parameter selection; it does not reproduce denoising, merging, or chimera removal.")),

                    div(class = "explanation-section",
                        h5("Quick guide"),
                        tags$ol(class = "mb-2",
                            tags$li("Load the primer-trimmed paired FASTQ files."),
                            tags$li("Set each truncLen near the point where its quality profile begins to deteriorate, while preserving sufficient overlap."),
                            tags$li("Use maxEE to remove reads with excessive accumulated error; 2 is a common starting value, not a universal optimum."),
                            tags$li("Keep the live overlap assessment at MODERATE or GOOD whenever possible."),
                            tags$li("Inspect per-sample retention before optional validation.")),
                        p(class = "small text-muted mb-0",
                          "A useful setting retains enough high-quality sequence to merge the pair; no single truncLen or maxEE is correct for every run.")),

                    h5("Scientific background"),

                    div(class = "explanation-section",
                        h5("Parameter definitions and interpretation"),

                        h6("truncLen", class = "text-primary fw-bold"),
                        p("truncLen specifies the retained length of the forward and reverse reads. Bases beyond the selected position are removed, and reads shorter than the requested length are discarded. Choose values that remove unreliable read tails without sacrificing the sequence required for paired-read overlap. The sum of the retained lengths must exceed the expected amplicon length, and a practical overlap margin should remain."),
                        div(class = "alert alert-warning py-2",
                            tags$strong("Important:"), " Ensure overlap for merging! DADA2's mergePairs() enforces
                            a hard floor of ", tags$code("minOverlap = 12"), " bp, but the recommended minimum is
                            amplicon-specific and usually higher -- see ",
                            tags$strong("\"Why Read Overlap Matters\""), " below for the full thresholds. The
                            current overlap status for your selection is shown live below the plot on the
                            Quality Profiles tab."),

                        h6("amplicon_min_length / amplicon_max_length", class = "text-primary fw-bold"),
                        p("Lower and upper bounds for the post-merge amplicon length filter, applied after
                        mergePairs() to discard chimeric or off-target merges whose length falls outside the
                        expected range for the selected primer pair. These come from the Primer Database's own ",
                        tags$code("Target Min (bp)"), " / ", tags$code("Target Max (bp)"),
                        " columns (see Visualize Help for the full derivation
                        methodology -- in short, the union of 1st/99th-percentile primer-trimmed amplicon
                        lengths observed via ", tags$em("in silico"), " PCR against four reference
                        databases), auto-filled into the Visualize sidebar's ",
                        tags$strong("Target Length"), " panel whenever a
                        literature primer pair is selected there -- ", tags$strong("but still manually
                        editable"), " afterward, for a manually-entered primer pair or to override a
                        literature pair's bounds with your own knowledge of the sample."),
                        p(class = "small text-muted mb-0",
                          "The resulting min/max values are shown in the Export table
                          (parameters ", tags$code("amplicon_min_length"), " and ",
                          tags$code("amplicon_max_length"), "). The upper bound (",
                          tags$code("Target Max (bp)"), ", \"amp_p99\") is what the Quality Profiles tab's
                          overlap status note below uses as the amplicon length -- see \"Why Read Overlap
                          Matters\" below for the full formula."),

                        h6("maxEE", class = "text-primary fw-bold"),
                        p("maxEE stands for maximum expected allowed errors. The expected error of a sequence is
                        calculated based on the quality scores of each base in that sequence -- essentially, a sum
                        of the probabilities of each base being an error, as inferred from the quality scores (Q)
                        provided by the sequencing platform: maxEE = sum(P). Q is related to the probability of an
                        error (P) by the formula P = 10^(-Q/10), so a higher Q-score indicates a lower probability
                        of error: Q = -10log10(P). Sequences with a total expected error rate exceeding the maxEE
                        threshold are discarded."),
                        tags$ul(
                            tags$li(tags$span(class = "text-success", tags$strong("Lower maxEE:")), " stricter filtering; fewer reads pass."),
                            tags$li(tags$span(class = "text-warning", tags$strong("maxEE around 1-2:")), " a common starting range for Illumina data."),
                            tags$li(tags$span(class = "text-danger", tags$strong("Higher maxEE:")), " more reads pass, including reads with greater expected error.")
                        )
                    ),

                    div(class = "explanation-section",
                        h5("Why read overlap matters"),
                        p("For paired-end Illumina sequencing, reads must overlap sufficiently to be
                        merged. Rather than judging overlap against a single average amplicon length,
                        this app classifies it against ", tags$strong("overlap at p99"), " -- the read
                        overlap computed at the ", tags$em("longest"), " amplicon length this primer
                        pair is expected to produce (", tags$code("amp_p99"), ", the Primer Database's
                        Target Max (bp)): ",
                        tags$code("overlap_at_p99 = truncLen_fwd + truncLen_rev - amp_p99"),
                        ". This worst-case framing means the classification already accounts for
                        biological length variation, rather than needing a separate margin added on top."),
                        tags$table(class = "table table-sm table-bordered small mb-2",
                            tags$thead(class = "table-light",
                                tags$tr(tags$th("Status"), tags$th("Overlap at p99"), tags$th("Basis"))
                            ),
                            tags$tbody(
                                tags$tr(
                                    tags$td(tags$span(class = "badge bg-primary", "GOOD")),
                                    tags$td(HTML("&ge; 50 bp")),
                                    tags$td("Comfortable; truncLen can be guided by quality scores alone.")
                                ),
                                tags$tr(
                                    tags$td(tags$span(class = "badge bg-info text-dark", "MODERATE")),
                                    tags$td("20-49 bp"),
                                    tags$td("Meets the DADA2 tutorial rule: ",
                                            tags$code("truncLen_fwd + truncLen_rev - amp_p99 >= 20"),
                                            " -- equivalent to \"20 + biological.length.variation\" at
                                            the shortest amplicon (p1).")
                                ),
                                tags$tr(
                                    tags$td(tags$span(class = "badge bg-warning text-dark", "CRITICAL")),
                                    tags$td("12-19 bp"),
                                    tags$td("Above the ", tags$code("minOverlap = 12"), " software default
                                            so merging proceeds, but below the tutorial recommendation;
                                            reduced accuracy due to short overlap near read ends.")
                                ),
                                tags$tr(
                                    tags$td(tags$span(class = "badge bg-danger", "LOW YIELD")),
                                    tags$td("0-11 bp"),
                                    tags$td("Reads physically overlap but fall below ",
                                            tags$code("minOverlap = 12"), "; ", tags$code("mergePairs()"),
                                            " attempts alignment, marks pairs as ",
                                            tags$code("accept = FALSE"), ", and silently drops them from
                                            output by default -- merging yield will be severely reduced.
                                            Lowering ", tags$code("minOverlap"), " is not recommended.")
                                ),
                                tags$tr(
                                    tags$td(tags$span(class = "badge bg-dark", "FAIL")),
                                    tags$td("< 0 bp"),
                                    tags$td("Reads do not reach each other; no overlap region exists;
                                            merging is impossible regardless of ",
                                            tags$code("minOverlap"), " setting.")
                                )
                            )
                        ),
                        p(class = "text-muted fst-italic",
                          "Note: Actual per-read overlap also depends on the sample's real insert size
                        distribution, not just the primer pair's expected range. Use quality trimming
                        to remove low-quality bases from read ends before merging. The live status for
                        your selected primer pair and truncation lengths is shown beside the overlap plot in Quality Profiles.")
                    ),

                    div(class = "explanation-section",
                        h5("Practical guidance"),
                        tags$ol(
                            tags$li("Choose truncation positions from the shape of the observed quality profiles; Q30 is a useful reference, not a mandatory cutoff."),
                            tags$li("Assess forward and reverse reads independently because their quality profiles can differ."),
                            tags$li("Check per-sample retention so pooled summaries do not hide problematic samples."),
                            tags$li("Confirm that both truncLen values preserve adequate overlap at the expected maximum amplicon length."),
                            tags$li("Use Validate on representative samples, then inspect read tracking in the full pipeline.")
                        )
                    ),

                    )
                )
            )
            )
        )
    ),
    nav_item(tags$span(class = "workflow-separator", ">")),

    # =========================================================================
    # 7.6 Validate Tab
    # =========================================================================
    nav_panel(
        title = tags$span(icon("flask"), " Validate"),
        value = "validation",

        # Validation inputs remain in the sidebar; results use the main panel.
        layout_sidebar(
            fillable = FALSE,

            sidebar = sidebar(
                width = 300,
                open = "always",

                accordion(
                    id = "val_controls_accordion",
                    open = c("val_settings_panel"),
                    accordion_panel(
                        title = "Validation Settings",
                        value = "val_settings_panel",
                        icon = icon("sliders-h"),

                        uiOutput("val_current_params"),

                        tags$p(class = "text-muted small mb-1",
                               "Choose the number of representative samples."),
                        tags$p(class = "text-muted small mb-2",
                               "The recommended value is automatically filled based on the dataset size!"),

                        numericInput("optim_validate_n_samples", "Samples to validate",
                                     value = 3, min = 1, max = 50, step = 1, width = "100%"),

                        # Select-all disables the numeric field and validates
                        # the complete loaded cohort.
                        checkboxInput("optim_validate_use_all",
                                      "Validate on all loaded samples", value = FALSE),

                        actionButton("optim_validate", "Validate with real DADA2",
                                     class = "btn-success btn-sm w-100 mt-2",
                                     icon = icon("flask"))
                    )
                ),

                # Zero-size anchor marking where the Processing Console
                # accordion belongs when this tab is active -- see Section
                # 7.1b and the longer comment on the Visualize tab's
                # Processing Console. The console panel itself is authored
                # once in the Visualize tab's sidebar and relocated to sit
                # directly before this anchor, as a direct sidebar child, by
                # JS whenever this tab becomes active.
                tags$div(id = "console_slot_validation")
            ),

            tagList(
                div(class = "page-intro",
                    h4("Validate with the real DADA2 pipeline"),
                    p("Optionally compare the fast retention estimate with observed through-pipeline retention on representative samples.")),
                tags$p(class = "text-muted mb-2",
                       style = "font-size:0.82rem;",
                       "Run the complete DADA2 workflow on a representative subset and compare observed through-pipeline retention with the fast quality-based estimate from Select. Runtime depends on sample size and the number selected."),

                uiOutput("optim_validation_status"),

                # One row per validated sample; percentage columns use in-cell
                # bars for quick comparison.
                card(
                    card_header(class = "bg-success text-white py-2",
                                tags$span(icon("table"),
                                          " Predicted vs. Real Retention per Sample")),
                    card_body(
                        class = "p-2",
                        uiOutput("optim_validation_table_ui")
                    )
                ),

                br(),

                # Validation interpretation and scientific context.
                accordion(
                    id = "val_help_accordion",
                    open = FALSE,
                    accordion_panel(
                        title = "How to interpret validation",
                        icon = icon("question-circle"),
                        div(class = "help-longform mt-0",
                            h5("What this tab does"),
                            div(class = "explanation-section",
                                h5("Purpose"),
                                p("The Validate tab runs the selected settings through the actual DADA2 workflow on representative samples. It reports retention after filtering, denoising, paired-read merging, and chimera removal. The fair estimate comparison is Predicted Merged % versus Real Merged %; Real Non-chim % is reported separately for final-yield context.")),

                            div(class = "explanation-section",
                                h5("Quick guide"),
                                tags$ol(
                                  tags$li("Choose a small representative subset; validating every sample is usually unnecessary."),
                                  tags$li("Run validation after the Quality Profiles and per-sample retention look acceptable."),
                                  tags$li("Compare predicted retention with filtered, denoised, merged, and non-chimeric retention.")),
                                div(class = "alert alert-info py-2",
                                  tags$strong("Interpretation: "),
                                  "small differences are expected because the fast estimate does not model denoising or overlap-region mismatches. Because Difference (pp) is Predicted Merged % minus Real Merged %, large or consistently positive differences mean the estimate was optimistic and the parameters or overlap should be reviewed."),
                                p(class = "small text-muted mb-0",
                                  "Validation is optional and may take several minutes to tens of minutes. Treat the comparison as a diagnostic check rather than a formal pass/fail test. Parameter definitions and overlap thresholds are documented in Select Help.")),

                            h5("Scientific background"),
                            div(class = "explanation-section",
                                h5("Why observed and predicted retention differ"),
                                p("The Select estimate evaluates read-level quality filters. DADA2 subsequently learns run-specific error rates, infers exact amplicon sequence variants, merges pairs only when their overlap is consistent, and removes sequences identified as chimeric. Losses at these later stages are biological and algorithmic outcomes that a quality-only estimate cannot predict."),
                                p(class = "mb-0", "Interpret differences across several representative samples. A small discrepancy is expected; a large systematic shortfall can indicate insufficient overlap, poor reverse-read quality, atypical error profiles, non-specific amplification, or substantial chimera formation. Validation supports judgment but does not establish a universal acceptance threshold."))
                        )
                    )
                )
            )
        )
    ),
    nav_item(tags$span(class = "workflow-separator", ">")),

    # =========================================================================
    # 7.7 Export Tab
    # =========================================================================
                nav_panel(
                  title = tags$span(icon("file-export"), " Export"),
                  value = "report",

                  # The sidebar contains only the save action and its status.
                  layout_sidebar(
                    fillable = FALSE,

                    sidebar = sidebar(
                      width = 300,
                      open = "always",
                      p(class = "small mb-2",
                        "Review the table and save the selected parameters for subsequent pipeline steps."),
                      actionButton("report_save_project", "Save parameters",
                                   class = "btn-success btn-sm w-100",
                                   icon = icon("floppy-disk")),
                      uiOutput("report_save_status"),

                      # Zero-size anchor marking where the Processing
                      # Console accordion belongs when this tab is active --
                      # see Section 7.1b and the longer comment on the
                      # Visualize tab's Processing Console. The console
                      # panel itself is authored once in the Visualize tab's
                      # sidebar and relocated to sit directly before this
                      # anchor, as a direct sidebar child, by JS whenever
                      # this tab becomes active.
                      tags$div(id = "console_slot_report")
                    ),

                    tagList(
                      card(
                        card_body(
                          class = "p-2",
                          # Natural height avoids a nested table scrollbar.
                          DTOutput("report_table_preview", fill = FALSE)
                        )
                      ),

                      br(),

                      # Export purpose and concise instructions.
                      accordion(
                        id = "report_help_accordion",
                        open = FALSE,
                        accordion_panel(
                          title = "How export works",
                          icon = icon("question-circle"),

                          div(class = "help-longform mt-0",
                              h5("What this tab does"),

                          div(class = "explanation-section",
                              h5("Purpose"),
                              p("The Export tab provides a final, read-only review of the assay definition, DADA2 parameters, overlap assessment, and available retention results. It writes these values to the workbook consumed by pipeline Step 5. Exporting records the current choices; it does not validate them or run the analysis.")),

                          div(class = "explanation-section",
                              h5("Quick guide"),
                              tags$ol(class = "mb-2",
                                tags$li("Review the grouped parameter table."),
                                tags$li("Return to Visualize or Select if any value needs adjustment."),
                                tags$li("Click Save parameters; Step 5 detects the saved workbook automatically.")),
                              p(class = "small text-muted mb-0",
                                "Load paired FASTQ quality profiles before export; saving remains disabled until data-derived retention values are available."))
                      )
                    )
                  )
                    )
                  )
                )

)
# =============================================================================
# SECTION 8: SERVER LOGIC
# =============================================================================
# The server function contains all reactive logic, event handlers, and output
# renderers. It is organized into sections for each tab, with shared utilities
# at the top.

server <- function(input, output, session) {
    
    # =========================================================================
    # 8.1 SHARED: Reactive Values Storage
    # =========================================================================
    # ReactiveValues stores loaded FASTQ summaries and validation results.
    
    # Select/Validate data
    rv_exp <- reactiveValues(
        files = NULL,                # Paired file information

        # Use the standard primer-trimmed directory when it exists.
        selected_path = DEFAULT_FASTQ_DIR,

        # Updated when the selected directory changes.
        n_detected_samples = if (!is.null(DEFAULT_FASTQ_DIR)) {
            default_patterns <- PATTERN_CHOICES[[1L]]
            estimate_paired_sample_count(
                DEFAULT_FASTQ_DIR,
                default_patterns$fwd,
                default_patterns$rev
            )
        } else {
            NULL
        },

        quality_data_fwd = NULL,     # Per-sample forward quality data
        quality_data_rev = NULL,     # Per-sample reverse quality data
        aggregated_fwd = NULL,       # Aggregated forward statistics
        aggregated_rev = NULL,       # Aggregated reverse statistics
        # Row-aligned paired quality matrices keyed by sample name.
        paired_quality_data = NULL,
        # Sampling metadata recorded in Export.
        profile_sampling = NULL,
        # Empirical DADA2 validation result, if run.
        validation_results = NULL,
        # Monotonic tokens prevent completed background work from being attached
        # after the user has changed the inputs that launched it.
        quality_request_id = 0L,
        validation_request_id = 0L
    )
    
    # =========================================================================
    # 8.2 SHARED: Console Message Helper
    # =========================================================================
    # Helper function to send real-time messages to the console
    # Uses custom JavaScript handler for immediate display
    
    add_console_msg <- function(msg, type = "info") {
        timestamp <- format(Sys.time(), "%H:%M:%S")
        
        # Color coding by message type
        color <- switch(type,
                        "info" = "#18bc9c",      # Teal for general info
                        "success" = "#2ecc71",   # Green for success
                        "warning" = "#f39c12",   # Orange for warnings
                        "error" = "#e74c3c",     # Red for errors
                        "progress" = "#3498db",  # Blue for progress
                        "validation" = "#9b59b6",# Purple for validation
                        "#95a5a6")               # Gray default
        
        # Send to JavaScript handler
        session$sendCustomMessage("console_msg",
                                  list(text = msg, time = timestamp, color = color))
    }

    # =========================================================================
    # 8.2b SHARED: Processing Console Relocation
    # =========================================================================
    # Whenever the active workflow tab changes, tell the client (via the
    # `move_console` custom message, handled in Section 7.3) to relocate the
    # single #exp_console_accordion node -- authored once in the Visualize
    # tab's sidebar (Section 7.4) -- into the newly active tab's own
    # `console_slot_<value>` placeholder (Sections 7.4/7.5/7.6/7.7). This is
    # what keeps the console visible, at the bottom of the sidebar, on every
    # tab, without duplicating its DOM node or losing message history.

    observeEvent(input$main_nav, {
        session$sendCustomMessage("move_console",
                                  list(slot = paste0("console_slot_", input$main_nav)))
    }, ignoreInit = FALSE)

    # =========================================================================
    # 8.3 SHARED: File Browser Setup
    # =========================================================================
    # Configure available volumes for the directory browser
    # This allows users to navigate their filesystem
    
    volumes <- c(
        Home = path.expand("~"),     # User's home directory
        Working = getwd(),            # Current working directory
        Root = "/"                    # Filesystem root
    )
    
    # Add common data directories if they exist
    if (dir.exists("/data")) volumes <- c(volumes, Data = "/data")
    if (dir.exists("/mnt")) volumes <- c(volumes, Mount = "/mnt")
    
    # Register the directory chooser
    shinyDirChoose(input, "exp_dir_select", roots = volumes, session = session)

    # Announce the pre-populated default directory (if found, Section 3.4)
    # so the console history mirrors what a manual Browse selection would
    # have logged, and it's obvious to the user why data appears loaded
    # without having clicked Browse.
    #
    # Read DEFAULT_FASTQ_DIR (a plain constant) rather than
    # rv_exp$selected_path here: reactiveValues fields can only be read
    # from within a reactive consumer (reactive()/observe()/render*/
    # isolate()), and this code runs once as plain top-level statements in
    # the server function body, before any reactive context exists. At this
    # exact point rv_exp$selected_path == DEFAULT_FASTQ_DIR regardless,
    # since no Browse event can have fired yet.
    if (!is.null(DEFAULT_FASTQ_DIR)) {
        add_console_msg(paste("Default directory pre-selected:",
                               DEFAULT_FASTQ_DIR), "info")
    }

    
    # =========================================================================
    # VISUALIZE SERVER LOGIC
    # =========================================================================
    
    # -------------------------------------------------------------------------
    # 8.4 VIS: Reactive Expressions
    # -------------------------------------------------------------------------
    
    # Calculate amplicon coordinates based on primer positions
    vis_amplicon_coords <- reactive({
        fwd_start <- input$vis_forward_primer_start
        rev_end <- input$vis_reverse_primer_end
        fwd_len <- input$vis_forward_primer_length
        rev_len <- input$vis_reverse_primer_length

        # Treat invalid primer lengths as zero so coordinate-based outputs remain
        # available while a field is being edited.
        if (is.null(fwd_len) || is.na(fwd_len) || fwd_len < 0) fwd_len <- 0
        if (is.null(rev_len) || is.na(rev_len) || rev_len < 0) rev_len <- 0

        # Validate inputs
        if (is.null(fwd_start) || is.null(rev_end) ||
            is.na(fwd_start) || is.na(rev_end) ||
            fwd_start >= rev_end) {
            return(list(start = 0, end = 0, size = 0,
                        forward_length = 0, reverse_length = 0,
                        insert_start = 0, insert_end = 0, insert_size = 0))
        }

        full_span <- rev_end - fwd_start + 1

        # Clamp so an oversized manual length entry can never claim more of
        # the primer-to-primer span than actually exists -- insert_size is
        # floored at 0 (via the clamp below) rather than going negative.
        fwd_len <- min(fwd_len, full_span)
        rev_len <- min(rev_len, full_span - fwd_len)

        list(
            start = fwd_start,
            end = rev_end,
            # Full primer-to-primer genomic span.
            size = full_span,
            # Primer lengths and the retained insert shown on the gene map.
            forward_length = fwd_len,
            reverse_length = rev_len,
            insert_start = fwd_start + fwd_len,
            insert_end = rev_end - rev_len,
            insert_size = max(0, full_span - fwd_len - rev_len)
        )
    })
    
    # Determine which variable regions are covered by the amplicon
    vis_covered_regions <- reactive({
        amp <- vis_amplicon_coords()
        if (amp$size == 0) return(character(0))

        # A region is covered if it overlaps with the amplicon
        variable_regions %>%
            filter(start <= amp$end & end >= amp$start) %>%
            pull(region)
    })

    # Expected primer-trimmed length range. Visualize auto-fills it from the
    # selected database row and still permits manual edits. Export saves both
    # bounds; exp_overlap_info() classifies overlap against the maximum.
    amplicon_length_range <- reactive({
        min_val <- input$vis_target_length_min
        max_val <- input$vis_target_length_max

        if (is.null(min_val) || is.null(max_val) || is.na(min_val) || is.na(max_val)) {
            # Defensive fallback while a field is temporarily cleared.
            return(list(min = NA, max = NA))
        }

        list(min = min_val, max = max_val)
    })

    # Five-tier status driven by "overlap at p99" --
    # the read overlap computed against the LONGEST amplicon length this
    # primer pair is expected to produce (amp_p99, the Primer Database's
    # own Target Max (bp) / amplicon_max_length, Section 4.3), i.e. the
    # worst-case (smallest-overlap) scenario within the expected length
    # distribution:
    #     overlap_at_p99 = truncLen_fwd + truncLen_rev - amp_p99
    # This is mathematically equivalent to a 20 bp rule evaluated at the
    # SHORTEST amplicon (p1) -- overlap_p1 >= 20 + (amp_p99 - amp_p1)
    # simplifies to (truncLen_fwd + truncLen_rev) - amp_p99 >= 20, i.e.
    # overlap_at_p99 >= 20 -- but is a single, directly-interpretable
    # number instead of two separate terms (a point-estimate overlap plus
    # a separate "biological length variation" margin).
    #
    #   Status     Overlap at p99   Basis
    #   GOOD       >= 50 bp         Comfortable; truncLen can be guided by
    #                               quality scores alone.
    #   MODERATE   20-49 bp         Meets the DADA2 tutorial rule
    #                               (truncLen_fwd + truncLen_rev - amp_p99
    #                               >= 20), equivalent to "20 +
    #                               biological.length.variation" evaluated
    #                               at the shortest amplicon (p1).
    #   CRITICAL   12-19 bp         Above mergePairs()'s minOverlap = 12 bp
    #                               software default so merging proceeds,
    #                               but below the tutorial recommendation;
    #                               reduced accuracy from short overlap
    #                               near read ends.
    #   LOW YIELD  0-11 bp          Reads physically overlap but fall
    #                               below minOverlap = 12; mergePairs()
    #                               attempts alignment, marks these pairs
    #                               accept = FALSE, and silently drops
    #                               them from output by default -- merging
    #                               yield will be severely reduced.
    #                               Lowering minOverlap is not recommended.
    #   FAIL       < 0 bp           Reads do not reach each other; no
    #                               overlap region exists; merging is
    #                               impossible regardless of minOverlap.
    classify_overlap_status <- function(overlap_at_p99) {
        if (overlap_at_p99 >= 50) {
            list(class = "overlap-good", label = "GOOD", icon = "check-circle",
                 text = paste0("GOOD (", overlap_at_p99, " bp overlap at p99). Comfortable -- ",
                               "truncLen can be guided by quality scores alone."))
        } else if (overlap_at_p99 >= 20) {
            list(class = "overlap-moderate", label = "MODERATE", icon = "check",
                 text = paste0("MODERATE (", overlap_at_p99, " bp overlap at p99). Meets the ",
                               "DADA2 tutorial rule (truncLen_fwd + truncLen_rev - amp_p99 >= 20 bp)."))
        } else if (overlap_at_p99 >= 12) {
            list(class = "overlap-critical", label = "CRITICAL", icon = "exclamation-triangle",
                 text = paste0("CRITICAL (", overlap_at_p99, " bp overlap at p99). Above ",
                               "mergePairs()'s minOverlap = 12 bp default so merging proceeds, but ",
                               "below the tutorial recommendation -- reduced accuracy from short ",
                               "overlap near read ends."))
        } else if (overlap_at_p99 >= 0) {
            list(class = "overlap-lowyield", label = "LOW YIELD", icon = "exclamation-circle",
                 text = paste0("LOW YIELD (", overlap_at_p99, " bp overlap at p99). Reads physically ",
                               "overlap but fall below minOverlap = 12 -- mergePairs() attempts ",
                               "alignment, marks these pairs accept = FALSE, and silently drops them ",
                               "from output by default; merging yield will be severely reduced. ",
                               "Lowering minOverlap is not recommended."))
        } else {
            list(class = "overlap-fail", label = "FAIL", icon = "times-circle",
                 text = paste0("FAIL (", overlap_at_p99, " bp overlap at p99). Reads do not reach ",
                               "each other -- no overlap region exists; merging is impossible ",
                               "regardless of minOverlap. Use longer reads or a shorter amplicon."))
        }
    }

    # Calculate overlap for the truncLen values selected in Select.
    exp_overlap_info <- reactive({
        amp <- vis_amplicon_coords()
        # Debounced (Section 8.8) -- this reactive drives a full plotly
        # re-render (output$exp_overlap_plot below) on every change, so it
        # should only recompute once the user pauses on a truncLen value,
        # not on every intermediate tick while dragging.
        fwd <- exp_truncLenF_db()
        rev <- exp_truncLenR_db()
        req(!is.null(fwd), !is.null(rev), !is.na(fwd), !is.na(rev))

        # Total read coverage
        total_read_length <- fwd + rev

        # Point-estimate overlap = total reads - amplicon size. Uses
        # amp$insert_size (the primer-free length of the CURRENTLY entered
        # primer start/end/lengths), NOT amp$size (the full primer-to-
        # primer span): truncLenF/R (fwd/rev above) are measured from the
        # start of the loaded FASTQ reads, and those reads already had both
        # primers removed by cutadapt (Step 3) before DADA2 ever saw them.
        # Comparing them against the full primer-inclusive span would
        # subtract off the primer lengths a second time, understating the
        # true overlap by roughly amp$forward_length + amp$reverse_length
        # bp. Still used for the plot's shaded overlap region and the
        # summary strip above the plot (a concrete, single-number
        # "estimated overlap" for THIS primer pair) -- positive = overlap,
        # negative = gap. Positive = overlap, Negative = gap.
        overlap <- total_read_length - amp$insert_size

        # Overlap at p99 -- the read overlap against the LONGEST amplicon
        # length this primer pair is expected to produce (amp_p99, the
        # Target Length: min max sidebar's upper bound). This, not the
        # point-estimate `overlap` above, drives the qualitative status
        # classification below -- see classify_overlap_status()'s comment
        # for the full rationale. Falls back to amp$insert_size if the
        # Target Length inputs are momentarily NA (defensive only, see
        # amplicon_length_range()'s own comment).
        len_range <- amplicon_length_range()
        amp_p99 <- if (!is.na(len_range$max)) len_range$max else amp$insert_size
        overlap_at_p99 <- total_read_length - amp_p99

        status <- classify_overlap_status(overlap_at_p99)

        list(
            overlap = max(0, overlap),
            gap = max(0, -overlap),
            overlap_at_p99 = overlap_at_p99,
            status = status
        )
    })
    
    # -------------------------------------------------------------------------
    # 8.5 VIS: Target Region / Literature Primers Dropdown Handlers
    # -------------------------------------------------------------------------
    # Data flow: target region -> matching primer-pair choices -> coordinates,
    # primer lengths, and target-length bounds. A target change selects the
    # first matching database row; manual coordinate edits remain possible.

    # Helper function to update primers (position, length, and target
    # length bounds) for a given literature primer pair name.
    vis_update_primers_and_regions <- function(primer_name) {
        primer_info <- primer_database %>% filter(name == primer_name)
        if (nrow(primer_info) == 1) {
            updateNumericInput(session, "vis_forward_primer_start",
                               value = primer_info$forward_start)
            updateNumericInput(session, "vis_reverse_primer_end",
                               value = primer_info$reverse_end)
            # Auto-fill verified primer lengths and target-length bounds.
            updateNumericInput(session, "vis_forward_primer_length",
                               value = nchar(primer_info$forward_seq))
            updateNumericInput(session, "vis_reverse_primer_length",
                               value = nchar(primer_info$reverse_seq))
            updateNumericInput(session, "vis_target_length_min",
                               value = primer_info$amplicon_min_length)
            updateNumericInput(session, "vis_target_length_max",
                               value = primer_info$amplicon_max_length)
        }
    }

    # Primer preset dropdown
    observeEvent(input$vis_primer_preset_select, {
        req(input$vis_primer_preset_select != "")
        vis_update_primers_and_regions(input$vis_primer_preset_select)
    })

    # Target Region Selection dropdown -- re-filters the Literature
    # Primers dropdown's choices to the primer pair(s) targeting the newly
    # selected region. ignoreInit = TRUE: the initial choices/selection
    # for both dropdowns are already set consistently in the UI
    # (DEFAULT_TARGET_REGION/DEFAULT_PRIMER_PAIR, Sections 4.7 and 7.4), so
    # this observer should only act on a genuine user change, not
    # re-fire (and potentially clobber the deliberate startup pairing)
    # during Shiny's initial reactive flush.
    observeEvent(input$vis_target_region_select, {
        new_choices <- vis_primer_choices_for_target(input$vis_target_region_select)
        matches <- primer_database %>% filter(target_regions == input$vis_target_region_select)

        # Always auto-select the FIRST matching primer pair, in
        # primer_database's own row order -- same behavior whether this
        # target has one published pair or several (e.g. "V3-V4" ->
        # 341F/785R, listed before 341F/806R). Cascades into
        # vis_update_primers_and_regions() above via the
        # input$vis_primer_preset_select observer, since updateSelectInput
        # changing the value fires that observer same as a manual pick.
        new_selected <- matches$name[1]

        updateSelectInput(session, "vis_primer_preset_select",
                          choices = new_choices, selected = new_selected)
    }, ignoreInit = TRUE)

    # -------------------------------------------------------------------------
    # 8.6 VIS: Output Renderers
    # -------------------------------------------------------------------------
    
    # Metric outputs
    # Shows the primer-free amplicon length (insert_size), not the raw
    # primer-to-primer genomic span (size) -- i.e. the length the sequence
    # actually has once cutadapt (Step 3) trims both primers off, using
    # the Fwd/Rev length inputs above (or the lengths auto-filled from a
    # selected literature primer pair). This is the one place in the app
    # where the primer-length-aware amplicon length is surfaced as a
    # number; the Gene Map (output$vis_main_plot, below) visualizes the
    # same breakdown directly.
    output$vis_amplicon_size <- renderText({
        vis_amplicon_coords()$insert_size
    })
    
    output$vis_regions_covered <- renderText({
        covered <- vis_covered_regions()
        if (length(covered) == 0) "None" else paste(covered, collapse = "–")
    })
    
    output$vis_coverage_pct <- renderText({
        amp <- vis_amplicon_coords()
        round(amp$size / FULL_16S_LENGTH * 100, 1)
    })
    
    # Main 16S plot
    output$vis_main_plot <- renderPlotly({
        amp <- vis_amplicon_coords()
        covered <- vis_covered_regions()
        
        # Create base plot

        # Keep every variable-region start/end coordinate on one baseline. A
        # compact horizontal font reduces crowding while preserving exact labels.
        region_boundary_ticks <- sort(unique(c(variable_regions$start,
                                                variable_regions$end)))

        p <- plot_ly() %>%
            layout(
                # Boundary labels use a compact font so all 18 coordinates
                # can remain horizontal on one baseline.
                xaxis = list(title = "Position (bp)", range = c(0, FULL_16S_LENGTH + 50),
                             showgrid = TRUE, gridcolor = "#e9ecef", zeroline = FALSE,
                             tickmode = "array", tickvals = region_boundary_ticks,
                             tickangle = 0, tickfont = list(size = 10),
                             titlefont = list(size = 15)),
                yaxis = list(title = "", showticklabels = FALSE, range = c(-0.6, 2.2),
                             showgrid = FALSE, zeroline = FALSE),
                showlegend = FALSE,
                # Horizontal labels need less bottom space than the former
                # staggered/two-line presentation.
                margin = list(l = 50, r = 50, t = 20, b = 45),
                paper_bgcolor = "white", plot_bgcolor = "white",
                hovermode = "closest"
            )

        # Add conserved regions (background)
        #
        # Note on mode = "lines": every filled shape below (fill = "toself")
        # also carries a `line` style for its border. plotly.js requires
        # "lines" in `mode` for that border to actually render -- omitting
        # it (e.g. mode = "none") makes plotly.js silently rewrite the mode
        # itself at draw time and log a console warning each time
        # ("A line object has been specified, but lines is not in the
        # mode..."). Setting mode = "lines" explicitly renders identically
        # but avoids the console noise, and documents the intent.
        #
        # Conserved regions are structural context and are always drawn.
        # hoverinfo = "skip" disables hover popups entirely for these region
        # rectangles -- both conserved and variable regions are already
        # labeled directly on the map (region name text for variable
        # regions; the grey conserved blocks are visually self-evident
        # between them), so a hover popup added no information a user
        # cannot already see and could expose raw styling values.
        for (i in 1:nrow(conserved_regions)) {
            p <- p %>% add_trace(
                type = "scatter", mode = "lines",
                # The 5th x/y pair repeats the 1st (start, bottom) point,
                # explicitly closing the rectangle's path. Without it, the
                # `mode = "lines"` stroke only connects the 4 listed points
                # (bottom edge -> right edge -> top edge), leaving the left
                # edge unstroked -- `fill = "toself"` still closes that edge
                # for FILLING purposes, but plotly.js does not also draw a
                # stroke along a fill-only closing segment.
                x = c(conserved_regions$start[i], conserved_regions$end[i],
                      conserved_regions$end[i], conserved_regions$start[i],
                      conserved_regions$start[i]),
                y = c(0.7, 0.7, 1.3, 1.3, 0.7),
                fill = "toself", fillcolor = "rgba(200, 200, 200, 0.3)",
                line = list(width = 0),
                hoverinfo = "skip"
            )
        }

        # Add variable regions
        # Selected-region codes derived from the single Target Region
        # Selection value (e.g. "V3-V4" -> c("V3", "V4")) via
        # target_to_region_codes() (Section 4.7).
        selected_region_codes <- target_to_region_codes(input$vis_target_region_select)
        for (i in 1:nrow(variable_regions)) {
            vr <- variable_regions[i, ]
            is_selected <- vr$region %in% selected_region_codes
            is_covered <- vr$region %in% covered

            opacity <- if (is_selected && is_covered) 0.9 else if (is_covered) 0.6 else 0.25
            # A thick dark-grey outline makes the selection unambiguous.
            border_width <- if (is_selected) 5 else 1
            border_color <- if (is_selected) "#4D4D4D" else vr$color

            p <- p %>% add_trace(
                type = "scatter", mode = "lines",
                # 5th point repeats the 1st to close the rectangle's stroke
                # (see the matching comment in the conserved-regions loop
                # above) -- without it the left edge of every region box
                # renders unstroked.
                x = c(vr$start, vr$end, vr$end, vr$start, vr$start),
                y = c(0.7, 0.7, 1.3, 1.3, 0.7),
                fill = "toself",
                fillcolor = paste0(vr$color, sprintf("%02X", round(opacity * 255))),
                line = list(color = border_color, width = border_width),
                # Labels are always visible, so region hover is unnecessary.
                hoverinfo = "skip"
            )
            
            # Region labels are always shown; the selected label is larger.
            p <- p %>% add_annotations(
                x = vr$midpoint, y = 1,
                text = vr$region, showarrow = FALSE,
                font = list(size = if (is_selected) 15 else 12,
                            color = if (is_selected && is_covered) "white"
                            else if (is_covered) "#333" else "#888")
            )
        }
        
        # Add amplicon: drawn as up to three adjoining segments rather than
        # one solid bar, so the forward/reverse primers' own footprints
        # (amp$forward_length/reverse_length bp each, from the Primer
        # Positions panel's length inputs or a selected literature pair)
        # are visually distinct from the actual sequenced insert between
        # them -- a teal segment (forward primer) and an orange segment
        # (reverse primer) flank a grey middle segment (amp$insert_start
        # to amp$insert_end -- the true amplicon once cutadapt, Step 3,
        # trims both primers off). Colors match the app's direction palette.
        # Either flank segment is skipped entirely if its primer length is
        # 0 (no length entered) -- the grey segment then simply extends to
        # that side's outer edge instead.
        if (amp$size > 0) {
            # Forward primer segment (blue), amp$start to amp$insert_start.
            if (amp$forward_length > 0) {
                p <- p %>% add_trace(
                    type = "scatter", mode = "lines",
                    # 5th point repeats the 1st to close the rectangle's
                    # stroke (see the matching comment in the conserved/
                    # variable region loops above) -- without it the left
                    # edge of this box renders unstroked.
                    x = c(amp$start, amp$insert_start, amp$insert_start, amp$start, amp$start),
                    y = c(0.35, 0.35, 0.55, 0.55, 0.35),
                    fill = "toself", fillcolor = "rgba(27, 158, 119, 0.45)",
                    line = list(color = "#1B9E77", width = 2),
                    hoverinfo = "skip"
                )
            }

            # True amplicon segment (grey), amp$insert_start to
            # amp$insert_end -- the sequence actually left after both
            # primers are trimmed off. Only drawn if any of it remains
            # (insert_size could be 0 if the entered primer lengths
            # consume the entire primer-to-primer span).
            if (amp$insert_size > 0) {
                p <- p %>% add_trace(
                    type = "scatter", mode = "lines",
                    x = c(amp$insert_start, amp$insert_end, amp$insert_end, amp$insert_start, amp$insert_start),
                    y = c(0.35, 0.35, 0.55, 0.55, 0.35),
                    fill = "toself", fillcolor = "rgba(38, 70, 83, 0.3)",
                    line = list(color = "#264653", width = 2),
                    # hoverinfo = "skip" (rather than a hovertemplate +
                    # <extra></extra>, this plot's original fix) -- the
                    # amplicon length breakdown is already shown directly
                    # on the map via the add_annotations() label just
                    # below, so no hover popup is needed on any element of
                    # this static diagram.
                    hoverinfo = "skip"
                )
            }

            # Reverse primer segment (red), amp$insert_end to amp$end.
            if (amp$reverse_length > 0) {
                p <- p %>% add_trace(
                    type = "scatter", mode = "lines",
                    x = c(amp$insert_end, amp$end, amp$end, amp$insert_end, amp$insert_end),
                    y = c(0.35, 0.35, 0.55, 0.55, 0.35),
                    fill = "toself", fillcolor = "rgba(217, 95, 2, 0.45)",
                    line = list(color = "#D95F02", width = 2),
                    hoverinfo = "skip"
                )
            }

            # Primer arrows retain hover text because their complete footprints
            # are not otherwise printed on the map. <extra></extra> suppresses
            # Plotly's secondary hover box.
            p <- p %>% add_trace(
                type = "scatter", mode = "markers",
                x = c(amp$start), y = c(0.1),
                marker = list(symbol = "triangle-right", size = 12, color = "#1B9E77"),
                hovertemplate = paste0("Forward primer: ", amp$start, "-",
                                       max(amp$start, amp$insert_start - 1),
                                       " bp (", amp$forward_length, " bp)<extra></extra>")
            )

            p <- p %>% add_trace(
                type = "scatter", mode = "markers",
                x = c(amp$end), y = c(0.1),
                marker = list(symbol = "triangle-left", size = 12, color = "#D95F02"),
                hovertemplate = paste0("Reverse primer: ", min(amp$end, amp$insert_end + 1), "-",
                                       amp$end, " bp (", amp$reverse_length, " bp)<extra></extra>")
            )

            # Amplicon length label -- leads with the true, primer-free
            # amplicon length (amp$insert_size, matching the "Amplicon (bp,
            # w/o primers)" metric above the map), with the primer-inclusive
            # total and each primer's own length spelled out alongside it
            # so the full breakdown is visible without needing to hover.
            p <- p %>% add_annotations(
                x = (amp$start + amp$end) / 2, y = -0.1,
                text = paste0("<b>", amp$insert_size, " bp</b> amplicon w/o primers  |  ",
                              amp$size, " bp total (Fwd ", amp$forward_length,
                              " + Rev ", amp$reverse_length, " bp)"),
                showarrow = FALSE, font = list(size = 13, color = "#264653")
            )
        }

        # Add 16S backbone
        p <- p %>% add_trace(
            type = "scatter", mode = "lines",
            x = c(1, FULL_16S_LENGTH), y = c(1.5, 1.5),
            line = list(color = "#264653", width = 3),
            hoverinfo = "skip"
        )
        
        # Add 16S backbone end markers.
        p <- p %>% add_annotations(
            x = c(1, FULL_16S_LENGTH), y = c(1.7, 1.7),
            text = c("5'", "3'"), showarrow = FALSE,
            font = list(size = 14, color = "#264653")
        )
        
        p
    })
    
    # Target Amplicon visualization driven by the Select truncLen controls.
    output$exp_overlap_plot <- renderPlotly({
        amp <- vis_amplicon_coords()
        ov <- exp_overlap_info()
        # Debounced (Section 8.8), and matched to ov's own debounced
        # truncLen reads above -- this is a full plotly re-render (many
        # add_trace()/add_annotations() calls), and using the raw, live
        # input$ here while ov above used the debounced value would let the
        # plot's read bars and ov's own overlap-bp figure disagree
        # mid-drag.
        fwd <- exp_truncLenF_db()
        rev <- exp_truncLenR_db()
        # Regions covered, shown as an in-plot label.
        covered <- vis_covered_regions()
        region_label <- if (length(covered) > 0) paste(covered, collapse = ", ") else "none selected"

        if (amp$size == 0 || amp$insert_size <= 0) {
            return(plot_ly(type = "scatter", mode = "markers") %>% layout(
                annotations = list(list(x = 0.5, y = 0.5,
                                        text = "Set valid primer positions",
                                        showarrow = FALSE, xref = "paper", yref = "paper"))
            ))
        }

        # Compute read extents before plotting so the x-axis ticks can include
        # them: fwd_end is where the
        # forward read's own truncLen ends (amp$insert_start + fwd_len),
        # rev_start is where the reverse read's own truncLen ends in the
        # 5' direction it's sequenced from (amp$insert_end - rev_len).
        # Both clamp against amp$insert_size (not amp$size) since these
        # reads only ever cover the primer-trimmed insert.
        fwd_len <- min(fwd, amp$insert_size)
        fwd_end <- amp$insert_start + fwd_len
        rev_len <- min(rev, amp$insert_size)
        rev_start <- amp$insert_end - rev_len

        # All x-coordinates below are absolute 16S gene positions, spanning
        # the full primer-to-primer footprint (amp$start to amp$end). The
        # forward/reverse primer segments are drawn alongside the amplicon
        # bar so the reader can distinguish primers from the true amplicon.
        # The forward/
        # reverse READ bars and the overlap region still anchor on
        # amp$insert_start/amp$insert_end (not amp$start/amp$end), since
        # the loaded FASTQ files already had both primers removed by
        # cutadapt (Step 3) before DADA2 ever saw them -- a read's first
        # base sits at amp$insert_start, not amp$start.
        # x-axis tick marks are pinned to all six boundary positions --
        # primer start, amplicon start, forward read end, reverse read
        # end, amplicon end, primer end -- via tickmode = "array" -- same
        # approach used for the 16S rRNA Gene Map's own x-axis (Section
        # 11.6) -- so every coordinate (including where each read's own
        # truncLen actually ends) is readable directly off the axis
        # without hovering. tickangle is left at 0 (horizontal, not
        # tilted) per explicit preference; labels may crowd when primer
        # lengths (or truncLen values close to the insert's own edges)
        # are short relative to the plot width, the same tradeoff already
        # used on the Gene Map.
        axis_tick_positions <- sort(unique(c(amp$start, amp$insert_start,
                                              fwd_end, rev_start,
                                              amp$insert_end, amp$end)))

        # Dashed vertical reference lines marking the true amplicon's
        # start (amp$insert_start, immediately after where the forward
        # primer was trimmed off) and end (amp$insert_end, immediately
        # before where the reverse primer was trimmed off) -- drawn as
        # layout "shapes" (rather than as their own scatter/hover trace)
        # so they read purely as static reference guides and never
        # themselves contribute an entry to the hover tooltip. Their bp
        # values are surfaced as small axis-adjacent annotations directly
        # below the plot area instead, so the exact coordinate is always
        # visible without requiring a hover at all.
        # Drawn only across the annotated bars, not the empty margins.
        primer_line_shapes <- list(
            list(type = "line", xref = "x", yref = "y",
                 x0 = amp$insert_start, x1 = amp$insert_start, y0 = 0.28, y1 = 1.38,
                 line = list(color = "#1B9E77", width = 1.5, dash = "dash")),
            list(type = "line", xref = "x", yref = "y",
                 x0 = amp$insert_end, x1 = amp$insert_end, y0 = 0.28, y1 = 1.38,
                 line = list(color = "#D95F02", width = 1.5, dash = "dash"))
        )

        p <- plot_ly() %>%
            layout(
                xaxis = list(title = "16S Gene Position (bp)",
                             range = c(amp$start - 10, amp$end + 10),
                             tickmode = "array", tickvals = axis_tick_positions,
                             tickangle = 0, tickfont = list(size = 13),
                             titlefont = list(size = 15),
                             showgrid = TRUE, gridcolor = "#e9ecef", zeroline = FALSE),
                # Tight range keeps the diagram prominent without empty bands.
                yaxis = list(title = "", showticklabels = FALSE, range = c(-1.05, 1.5),
                             showgrid = FALSE, zeroline = FALSE),
                showlegend = TRUE,
                # Vertical legend outside the plotting region.
                legend = list(orientation = "v", x = 1.02, xanchor = "left",
                              y = 0.5, yanchor = "middle", font = list(size = 13)),
                margin = list(l = 12, r = 150, t = 30, b = 34),
                paper_bgcolor = "white", plot_bgcolor = "white",
                shapes = primer_line_shapes,
                # Show only the nearest trace rather than stacked x-coordinate
                # tooltips at shared boundaries.
                hovermode = "closest"
            ) %>%
            add_annotations(
                # No bp value appended here -- amp$insert_start is already
                # readable directly off the x-axis (tickvals pinned to it
                # just above), so repeating it in this label would be
                # redundant.
                x = amp$insert_start, y = -0.9, xref = "x", yref = "y",
                text = "Amplicon start",
                showarrow = FALSE, font = list(size = 10, color = "#1B9E77"),
                xanchor = "left"
            ) %>%
            add_annotations(
                x = amp$insert_end, y = -0.9, xref = "x", yref = "y",
                text = "Amplicon end",
                showarrow = FALSE, font = list(size = 10, color = "#D95F02"),
                xanchor = "right"
            )

        # Amplicon bar -- drawn across the PRIMER-TRIMMED insert
        # (amp$insert_start to amp$insert_end), not the full primer-to-
        # primer span, since that insert is what the loaded FASTQ files
        # (and therefore DADA2) actually contain. Added first so it anchors the legend
        # and reads as the plot's central reference bar, with the primer
        # segments and reads ordered around it rather than wrapped
        # immediately on either side of it.
        # Note on mode = "lines": see output$vis_main_plot above -- these
        # fill = "toself" bars need "lines" in mode for their border to
        # render without plotly.js logging a console warning. hovertemplate
        # (rather than hoverinfo = "text" + text =) is used throughout this
        # plot, same as vis_main_plot above, to suppress plotly's default
        # "extra" secondary hover box via a trailing "<extra></extra>". Each
        # box's 5th x/y pair repeats its 1st point to close the rectangle's
        # stroke -- without it, `mode = "lines"` only strokes the 3 listed
        # segments and the left edge renders unstroked (fill = "toself"
        # closes the shape for FILLING purposes only, not for the stroke).
        #
        # hoveron = "fills" avoids duplicate hover text at the repeated closing
        # point and makes the whole rectangle interactive.
        p <- p %>% add_trace(
            type = "scatter", mode = "lines",
            x = c(amp$insert_start, amp$insert_end, amp$insert_end, amp$insert_start, amp$insert_start),
            y = c(0.3, 0.3, 0.6, 0.6, 0.3),
            fill = "toself", fillcolor = "rgba(200, 200, 200, 0.5)",
            line = list(width = 1, color = "#999"),
            # This bar spans amp$insert_start to amp$insert_end, i.e. the exact
            # primer-free target length DADA2's mergePairs() would produce
            # from a successfully merged forward/reverse read pair, so
            # "Merged Amplicon" describes what the bar represents (the
            # target/expected merged-sequence length), not an actual
            # merge that has been performed on real data at this point.
            name = "Merged Amplicon",
            hoveron = "fills",
            hovertemplate = paste0("Merged Amplicon: ", amp$insert_size, " bp<extra></extra>")
        )

        # Amplicon-size label on the amplicon bar.
        p <- p %>% add_annotations(
            x = (amp$insert_start + amp$insert_end) / 2, y = 0.45,
            xref = "x", yref = "y",
            text = paste0("Amplicon ", amp$insert_size, " bp"),
            showarrow = FALSE, font = list(size = 11, color = "#333333")
        )

        # Forward primer segment -- the primer itself (amp$start to
        # amp$insert_start), drawn on the SAME row as the forward/reverse
        # READ bars below (y = 0.75-1.05), not the grey amplicon bar's row
        # (y = 0.3-0.6) -- so it sits immediately to the left of, and
        # visually contiguous with, the Forward trimmed reads bar it precedes on
        # the gene, rather than appearing to belong with the amplicon bar
        # underneath it. Uses a lighter fill and a dotted (rather than
        # solid) border so it still reads as "this was here on the gene
        # but is no longer present in the loaded reads" instead of being
        # mistaken for sequence DADA2 actually processes. Colors match the
        # 16S rRNA Gene Map's own forward-primer segment (Section 8.6).
        # The legend groups trimmed primer sequence with read coverage. Only
        # drawn when a forward primer length is available.
        if (amp$forward_length > 0) {
            p <- p %>% add_trace(
                type = "scatter", mode = "lines",
                x = c(amp$start, amp$insert_start, amp$insert_start, amp$start, amp$start),
                y = c(0.75, 0.75, 1.05, 1.05, 0.75),
                fill = "toself", fillcolor = "rgba(27, 158, 119, 0.3)",
                line = list(width = 1.5, color = "#1B9E77", dash = "dot"),
                name = paste0("Fwd primer (", amp$forward_length, " bp)"),
                hoveron = "fills",
                hovertemplate = paste0("Forward primer (trimmed off before DADA2): ",
                                       amp$forward_length, " bp<extra></extra>")
            )
        }

        # Reverse primer segment -- the primer itself (amp$insert_end to
        # amp$end), same styling/rationale/placement as the forward primer
        # segment above -- also on the read bars' row (y = 0.75-1.05), so
        # it sits immediately to the right of, and visually contiguous
        # with, the Reverse trimmed reads bar it follows. Only drawn if a reverse
        # primer length is set.
        if (amp$reverse_length > 0) {
            p <- p %>% add_trace(
                type = "scatter", mode = "lines",
                x = c(amp$insert_end, amp$end, amp$end, amp$insert_end, amp$insert_end),
                y = c(0.75, 0.75, 1.05, 1.05, 0.75),
                fill = "toself", fillcolor = "rgba(217, 95, 2, 0.3)",
                line = list(width = 1.5, color = "#D95F02", dash = "dot"),
                name = paste0("Rev primer (", amp$reverse_length, " bp)"),
                hoveron = "fills",
                hovertemplate = paste0("Reverse primer (trimmed off before DADA2): ",
                                       amp$reverse_length, " bp<extra></extra>")
            )
        }

        # Forward trimmed reads -- starts at amp$insert_start (immediately after
        # where the forward primer was trimmed off by cutadapt), not
        # amp$start (the forward primer's own start), since the loaded
        # FASTQ reads no longer contain the primer. fwd_len/fwd_end were
        # already computed above (needed there for the x-axis tick
        # positions), so they're just reused here.
        p <- p %>% add_trace(
            type = "scatter", mode = "lines",
            x = c(amp$insert_start, fwd_end, fwd_end, amp$insert_start, amp$insert_start),
            y = c(0.75, 0.75, 1.05, 1.05, 0.75),
            fill = "toself", fillcolor = "rgba(27, 158, 119, 0.7)",
            line = list(width = 2, color = "#1B9E77"),
            # "Forward trimmed reads" distinguishes this bar from the primer
            # segment in the same legend and makes
            # explicit that this bar represents the already primer-trimmed
            # reads DADA2 actually loads, not the raw (primer-inclusive)
            # forward read.
            name = paste0("Forward trimmed reads (", fwd, " bp)"),
            hoveron = "fills",
            hovertemplate = paste0("Forward trimmed reads: ", fwd, " bp<extra></extra>")
        )

        # Reverse trimmed reads -- ends at amp$insert_end (immediately before where
        # the reverse primer was trimmed off), not amp$end (the reverse
        # primer's own start), for the same already-primer-trimmed reason.
        # rev_len/rev_start were already computed above (needed there for
        # the x-axis tick positions), so they're just reused here.
        p <- p %>% add_trace(
            type = "scatter", mode = "lines",
            x = c(rev_start, amp$insert_end, amp$insert_end, rev_start, rev_start),
            y = c(0.75, 0.75, 1.05, 1.05, 0.75),
            fill = "toself", fillcolor = "rgba(217, 95, 2, 0.7)",
            line = list(width = 2, color = "#D95F02"),
            # "Reverse trimmed reads" -- same reasoning as "Forward trimmed
            # reads" above.
            name = paste0("Reverse trimmed reads (", rev, " bp)"),
            hoveron = "fills",
            hovertemplate = paste0("Reverse trimmed reads: ", rev, " bp<extra></extra>")
        )

        # Overlap region
        if (ov$overlap > 0) {
            overlap_start <- amp$insert_end - rev_len
            overlap_end <- amp$insert_start + fwd_len

            if (overlap_start < overlap_end) {
                # Color the overlap bar from the shared status tier.
                # ov$status$class is the same
                # classify_overlap_status() tier (Section 8.4) used
                # everywhere else in the app for this metric, so this box and
                # the "Overlap check" caption (output$exp_overlap_msg_left,
                # just below) both agree on one color per tier -- see
                # OVERLAP_STATUS_COLORS' definition above.
                # box_rgb (the raw theme color) tints the translucent fill;
                # box_text_color (its darker, WCAG-AA-legible variant) draws
                # the border stroke and the "Overlap X bp" annotation font,
                # both of which need to stay readable against the plot's
                # white background.
                box_rgb        <- overlap_status_rgb(ov$status$class)
                box_text_color <- overlap_status_text_color(ov$status$class)
                p <- p %>% add_trace(
                    type = "scatter", mode = "lines",
                    x = c(overlap_start, overlap_end, overlap_end, overlap_start, overlap_start),
                    y = c(1.15, 1.15, 1.35, 1.35, 1.15),
                    fill = "toself", fillcolor = sprintf("rgba(%s, 0.7)", box_rgb),
                    line = list(width = 2, color = box_text_color),
                    name = paste0("Overlap (", ov$overlap, " bp)"),
                    hoveron = "fills",
                    hovertemplate = paste0("Overlap: ", ov$overlap, " bp<extra></extra>")
                )
                # Place the overlap label above the bar so it remains legible.
                p <- p %>% add_annotations(
                    x = (overlap_start + overlap_end) / 2, y = 1.46,
                    xref = "x", yref = "y",
                    text = paste0("Overlap ", ov$overlap, " bp"),
                    showarrow = FALSE, font = list(size = 11, color = box_text_color))
            }
        }

        # Direction arrows
        p <- p %>% add_annotations(x = (amp$insert_start + fwd_end) / 2, y = 0.9, text = "→",
                                   showarrow = FALSE, font = list(size = 14, color = "white"))
        p <- p %>% add_annotations(x = rev_start + (amp$insert_end - rev_start) / 2, y = 0.9,
                                   text = "←", showarrow = FALSE,
                                   font = list(size = 14, color = "white"))

        # Region label in the plot; overlap status is rendered beside it.
        p <- p %>% add_annotations(
            x = 0.0, xref = "paper", xanchor = "left",
            y = 1.0, yref = "paper", yanchor = "bottom",
            text = paste0("Region(s): ", region_label),
            showarrow = FALSE, align = "left",
            font = list(size = 11, color = "#333333")
        )

        p
    })

    # Overlap-quality message rendered beside the plot. Tinted lightly so it reads
    # as a caption rather than a solid alert) by overlap QUALITY, not a fixed
    # color -- ov$status$class drives the same color used for the overlap
    # box on the plot above; see OVERLAP_STATUS_COLORS.
    output$exp_overlap_msg_left <- renderUI({
        ov <- exp_overlap_info()
        # Text uses the darker, WCAG-AA-legible variant; the
        # background/border tints use the raw theme color, same as the
        # overlap box on the plot above.
        status_text_color <- overlap_status_text_color(ov$status$class)
        status_rgb         <- overlap_status_rgb(ov$status$class)
        div(style = paste0("padding:8px 10px;border-radius:6px;font-size:0.8rem;",
                           "line-height:1.35;color:", status_text_color, ";",
                           "background:rgba(", status_rgb, ",0.12);",
                           "border:1px solid rgba(", status_rgb, ",0.4);"),
            tags$strong("Overlap check"), tags$br(),
            ov$status$text)
    })

    output$exp_status_summary <- renderUI({
        ov <- exp_overlap_info()
        fmt_retention <- function(x) {
            value <- tryCatch(x(), error = function(e) NA_real_)
            if (is.finite(value)) paste0(round(value, 1), "%") else "--"
        }
        div(class = "status-strip", `aria-label` = "Current analysis status",
            span(class = "status-chip", icon("vials"), " ",
                 if (is.null(rv_exp$files)) "0 samples" else paste(length(rv_exp$files$fnFs), "samples")),
            span(class = "status-chip", "Fwd retained: ", fmt_retention(exp_retention_fwd)),
            span(class = "status-chip", "Paired retained: ", fmt_retention(exp_retention_combined)),
            span(class = "status-chip", "Rev retained: ", fmt_retention(exp_retention_rev)),
            span(class = "status-chip",
                 style = paste0("border-color:", overlap_status_color(ov$status$class),
                                ";color:", overlap_status_text_color(ov$status$class), ";"),
                 icon("arrows-left-right"), " Overlap: ", tags$strong(ov$status$label)))
    })

    # Primer database and its caveat footnotes in Visualize Help.
    output$vis_primer_table_help <- renderTable({
        build_primer_database_table()
    }, striped = TRUE, hover = TRUE, bordered = TRUE,
       sanitize.text.function = function(x) x)

    # Footnotes are UI paragraphs rather than a data frame.
    output$vis_primer_table_footnote_help <- renderUI({
        build_primer_database_footnote_html()
    })

    
    # =========================================================================
    # SELECT SERVER LOGIC
    # =========================================================================
    
    # -------------------------------------------------------------------------
    # 8.7 EXP: Directory Selection Handler
    # -------------------------------------------------------------------------
    
    observeEvent(input$exp_dir_select, {
        if (!is.integer(input$exp_dir_select)) {
            path <- parseDirPath(volumes, input$exp_dir_select)
            if (length(path) > 0 && nchar(path) > 0) {
                rv_exp$selected_path <- path
                add_console_msg(paste("Selected:", path), "info")

                # Pre-load sample estimate for Export metadata and defaults.
                selected_patterns <- PATTERN_CHOICES[[input$exp_pattern_choice]]
                rv_exp$n_detected_samples <- estimate_paired_sample_count(
                    path,
                    selected_patterns$fwd,
                    selected_patterns$rev
                )
            }
        }
    })

    # Keep the pre-load estimate aligned with the currently selected filename
    # convention. Invalid or incomplete pair sets correctly report zero.
    observeEvent(input$exp_pattern_choice, {
        req(!is.null(rv_exp$selected_path))
        selected_patterns <- PATTERN_CHOICES[[input$exp_pattern_choice]]
        rv_exp$n_detected_samples <- estimate_paired_sample_count(
            rv_exp$selected_path,
            selected_patterns$fwd,
            selected_patterns$rev
        )
    }, ignoreInit = TRUE)
    
    output$exp_selected_dir <- renderText({
        if (is.null(rv_exp$selected_path)) return("No directory selected")
        basename(rv_exp$selected_path)
    })
    
    # =========================================================================
    # 8.8 EXP: Debounced Slider Reactives
    # =========================================================================
    # Expensive retention/overlap calculations wait until the user pauses. The
    # following interactions intentionally use raw slider values:
    #   - The dashed truncLen line's live position on the quality plot (the
    #     observeEvent(input$exp_truncLenF/R, ...) plotlyProxy relayout pair
    #     just below output$exp_plot_fwd/rev_quality) -- this only moves an
    #     existing line via a lightweight plotlyProxyInvoke("relayout", ...)
    #     call, it does not recompute or re-render the underlying curve, so
    #     there is no wasted work to avoid; debouncing it would just make the
    #     line visibly lag behind the slider handle while dragging.
    #   - The -/+ stepper buttons (register_stepper(), Section 8.12) and the
    #     platform-change slider-bound clamps (updateSliderInput() calls
    #     above) -- both need the slider's true CURRENT value at the instant
    #     they fire, not a stale pre-debounce one.
    #   - The empirical DADA2 validation run (observeEvent(input$optim_validate,
    #     ...), Section 8.13) and output$val_current_params -- triggered by an
    #     explicit button click / display of state at that instant, not by the
    #     slider itself, so there is nothing to debounce.
    DEBOUNCE_MILLIS <- 800

    exp_truncLenF_db <- debounce(reactive(input$exp_truncLenF), DEBOUNCE_MILLIS)
    exp_truncLenR_db <- debounce(reactive(input$exp_truncLenR), DEBOUNCE_MILLIS)
    exp_maxEEF_db    <- debounce(reactive(as.numeric(input$exp_maxEEF)), DEBOUNCE_MILLIS)
    exp_maxEER_db    <- debounce(reactive(as.numeric(input$exp_maxEER)), DEBOUNCE_MILLIS)

    # -------------------------------------------------------------------------
    # 8.9 EXP: Console Clear Handler
    # -------------------------------------------------------------------------

    observeEvent(input$exp_clear_console, {
        session$sendCustomMessage("console_clear", list())
    })
    
    # -------------------------------------------------------------------------
    # 8.10 EXP: Load Quality Profiles
    # -------------------------------------------------------------------------

    quality_profile_task <- ExtendedTask$new(function(files, target_reads, base_seed,
                                                       request_id, selected_path, pattern_choice) {
        promises::future_promise({
            sample_indices <- seq_along(files$fnFs)
            paired_results <- lapply(sample_indices, function(i) {
                prof <- extract_paired_quality_profile(
                    files$fnFs[i], files$fnRs[i], n_pairs = target_reads,
                    seed = base_seed + i
                )
                if (is.null(prof)) return(NULL)
                list(
                    paired = prof,
                    fwd = summarize_quality_matrix(prof$quality_matrix_fwd),
                    rev = summarize_quality_matrix(prof$quality_matrix_rev)
                )
            })
            names(paired_results) <- files$sample_names
            list(files = files, target_reads = target_reads,
                 paired_results = paired_results, request_id = request_id,
                 selected_path = selected_path, pattern_choice = pattern_choice)
        })
    })
    bind_task_button(quality_profile_task, "exp_load_data")
    
    observeEvent(input$exp_load_data, {
        req(rv_exp$selected_path)

        rv_exp$quality_request_id <- rv_exp$quality_request_id + 1L
        request_id <- rv_exp$quality_request_id

        # A new load is a transaction: discard every result tied to the prior
        # directory before inspecting or processing the new selection. This
        # prevents a failed reload from pairing a new file list with stale
        # quality matrices, plots, sampling metadata, or validation results.
        rv_exp$files <- NULL
        rv_exp$quality_data_fwd <- NULL
        rv_exp$quality_data_rev <- NULL
        rv_exp$aggregated_fwd <- NULL
        rv_exp$aggregated_rev <- NULL
        rv_exp$paired_quality_data <- NULL
        rv_exp$profile_sampling <- NULL
        rv_exp$validation_results <- NULL
        shinyjs::hide("exp_analysis_content")

        shinyjs::html("exp_empty_state",
                      html = paste0(
                          "<div class='spinner-border text-secondary' role='status' aria-label='Processing'></div>",
                          "<h5 class='mt-3'>Processing paired FASTQ files</h5>",
                          "<p class='mb-0'>Sampling reads and building quality profiles. This may take a few minutes.</p>"
                      ))
        patterns <- PATTERN_CHOICES[[input$exp_pattern_choice]]

        add_console_msg(paste("Loading from:", rv_exp$selected_path), "info")
        add_console_msg(paste("Pattern:", patterns$fwd, "/", patterns$rev), "info")
        add_console_msg("", "info")

        # Detect paired files
        files <- detect_paired_files(rv_exp$selected_path, patterns$fwd, patterns$rev)

        if (!is.null(files$error)) {
            add_console_msg(paste("ERROR:", files$error), "error")
            showNotification(files$error, type = "error")
            return()
        }

        # Sampling depth adapts to the detected sample count.
        target_reads <- recommended_pairs_per_sample(length(files$fnFs))
        add_console_msg(paste("Found", length(files$fnFs), "paired samples"), "success")
        add_console_msg(paste0("Auto-selected sampling depth: ",
                               format(target_reads, big.mark = ","), " reads/file"), "info")
        add_console_msg("", "info")

        # Fill Validate's recommended subset and cap it at the loaded cohort.
        n_total_loaded <- length(files$fnFs)
        recommended_val_n <- recommended_validation_samples(n_total_loaded)
        updateNumericInput(session, "optim_validate_n_samples",
                            value = recommended_val_n, max = n_total_loaded)
        add_console_msg(paste0("Recommended validation sample count: ", recommended_val_n,
                               " of ", n_total_loaded,
                               " (Validate tab -- adjust or select all there)."), "info")
        
        add_console_msg("Sampling paired reads in a background worker...", "progress")
        quality_profile_task$invoke(
            files, target_reads, SAMPLING_BASE_SEED, request_id,
            rv_exp$selected_path, input$exp_pattern_choice
        )
    })

    observe({
        result <- quality_profile_task$result()
        if (!identical(result$request_id, rv_exp$quality_request_id) ||
            !identical(normalizePath(result$selected_path, mustWork = FALSE),
                       normalizePath(rv_exp$selected_path, mustWork = FALSE)) ||
            !identical(result$pattern_choice, input$exp_pattern_choice)) {
            return()
        }
        files <- result$files
        paired_results <- result$paired_results
        quality_fwd <- lapply(paired_results, function(x) if (is.null(x)) NULL else x$fwd)
        quality_rev <- lapply(paired_results, function(x) if (is.null(x)) NULL else x$rev)
        paired_data <- lapply(paired_results, function(x) if (is.null(x)) NULL else x$paired)

        aggregate_profiles <- function(quality_list) {
            valid <- quality_list[!vapply(quality_list, is.null, logical(1))]
            if (length(valid) == 0L) return(NULL)
            rows <- do.call(rbind, lapply(valid, function(profile) data.frame(
                Position = profile$position_stats$Position,
                Median = profile$position_stats$Median,
                Weight = profile$n_sampled
            )))
            rows %>% group_by(Position) %>%
                summarise(Median = weighted.mean(Median, Weight, na.rm = TRUE), .groups = "drop") %>%
                as.data.frame() %>% list(position_stats = .)
        }

        rv_exp$quality_data_fwd <- quality_fwd
        rv_exp$quality_data_rev <- quality_rev
        rv_exp$paired_quality_data <- paired_data
        rv_exp$aggregated_fwd <- aggregate_profiles(quality_fwd)
        rv_exp$aggregated_rev <- aggregate_profiles(quality_rev)
        rv_exp$files <- files
        rv_exp$profile_sampling <- list(
            mode = "Automatic", requested = result$target_reads,
            actual = vapply(paired_data, function(p) if (is.null(p)) NA_real_ else as.numeric(p$n_sampled), numeric(1)),
            base_seed = SAMPLING_BASE_SEED
        )
        rv_exp$validation_results <- NULL
        shinyjs::hide("exp_empty_state")
        shinyjs::show("exp_analysis_content")
        add_console_msg(paste("Processed", length(files$fnFs), "samples."), "success")
        showNotification(paste("Processed", length(files$fnFs), "samples."), type = "message")
    })

    # A directory or naming-pattern change invalidates all loaded summaries and
    # any in-flight task launched for the previous selection.
    observeEvent(list(rv_exp$selected_path, input$exp_pattern_choice), {
        rv_exp$quality_request_id <- rv_exp$quality_request_id + 1L
        rv_exp$files <- NULL
        rv_exp$quality_data_fwd <- NULL
        rv_exp$quality_data_rev <- NULL
        rv_exp$aggregated_fwd <- NULL
        rv_exp$aggregated_rev <- NULL
        rv_exp$paired_quality_data <- NULL
        rv_exp$profile_sampling <- NULL
        rv_exp$validation_results <- NULL
        rv_exp$validation_request_id <- rv_exp$validation_request_id + 1L
        shinyjs::hide("exp_analysis_content")
    }, ignoreInit = TRUE)
    
    # -------------------------------------------------------------------------
    # 8.11 EXP: Sequencing Platform Live Feedback
    # -------------------------------------------------------------------------
    # Re-cap both truncLen controls when the sequencing kit changes.

    # ignoreInit = TRUE: don't reclamp truncLenF/R the instant the app
    # loads (before the user has touched this dropdown at all) -- only
    # react to an actual platform change.
    observeEvent(input$vis_platform, {
        platform_max <- PLATFORM_READ_LENGTHS[[input$vis_platform]]
        req(is.finite(platform_max))
        updateSliderInput(session, "exp_truncLenF", max = platform_max,
                          value = min(input$exp_truncLenF, platform_max))
        updateSliderInput(session, "exp_truncLenR", max = platform_max,
                          value = min(input$exp_truncLenR, platform_max))
    }, ignoreInit = TRUE)

    # -------------------------------------------------------------------------
    # 8.12 EXP: Retention Calculations
    # -------------------------------------------------------------------------
    
    # Pooled per-read total EE at the current truncLen (all reads across all
    # samples), computed ONCE per truncLen change and reused for both the
    # Forward/Reverse % readouts and the retention-vs-maxEE curve. The empirical
    # CDF is exact at every 0.1-step maxEE value and keeps all readouts consistent.
    # Debounced (Section 8.8): this loops read_ee_at_truncLen() over every
    # loaded sample, so it should fire once after the user pauses on a
    # truncLen value, not once per intermediate tick while dragging.
    exp_read_ee_fwd <- reactive({
        req(rv_exp$quality_data_fwd)
        unlist(lapply(rv_exp$quality_data_fwd, function(qd)
            read_ee_at_truncLen(qd$qual_matrix, qd$cum_ee_matrix,
                                exp_truncLenF_db(), RETENTION_DEFAULT_TRUNCQ)),
            use.names = FALSE)
    })
    exp_read_ee_rev <- reactive({
        req(rv_exp$quality_data_rev)
        unlist(lapply(rv_exp$quality_data_rev, function(qd)
            read_ee_at_truncLen(qd$qual_matrix, qd$cum_ee_matrix,
                                exp_truncLenR_db(), RETENTION_DEFAULT_TRUNCQ)),
            use.names = FALSE)
    })

    # Debounced maxEE reads (Section 8.8) -- same reasoning: these are
    # cheap on their own, but they feed the retention-vs-maxEE plots and the
    # Retained-reads barplot, so keeping them in step with the debounced
    # exp_read_ee_fwd/rev() above avoids re-rendering those on every tick too.
    exp_retention_fwd <- reactive({
        ee <- exp_read_ee_fwd()
        if (length(ee) == 0) return(NA)
        round(100 * mean(ee <= exp_maxEEF_db()), 1)
    })

    exp_retention_rev <- reactive({
        ee <- exp_read_ee_rev()
        if (length(ee) == 0) return(NA)
        round(100 * mean(ee <= exp_maxEER_db()), 1)
    })

    # -/+ buttons nudge truncLen by 1 bp and maxEE by 0.1,
    # clamped to range. truncLen's upper bound tracks the selected platform's
    # read length.
    exp_platform_max <- reactive({
        m <- PLATFORM_READ_LENGTHS[[input$vis_platform]]
        if (is.null(m) || !is.finite(m)) 300 else m
    })
    local({
        register_stepper <- function(id, step, lo, hi_fun) {
            observeEvent(input[[paste0(id, "_dec")]], {
                v <- input[[id]]; if (is.null(v)) return()
                updateSliderInput(session, id, value = max(lo, round(v - step, 3)))
            })
            observeEvent(input[[paste0(id, "_inc")]], {
                v <- input[[id]]; if (is.null(v)) return()
                updateSliderInput(session, id, value = min(hi_fun(), round(v + step, 3)))
            })
        }
        # Fixed-width button blocks preserve alignment with the plot region.
        register_stepper("exp_truncLenF", 1,   0, function() exp_platform_max())
        register_stepper("exp_truncLenR", 1,   0, function() exp_platform_max())
        register_stepper("exp_maxEEF",    0.1, 0, function() 10)
        register_stepper("exp_maxEER",    0.1, 0, function() 10)
    })
    
    # True paired retention: evaluate pass/fail
    # per read PAIR (paired_pass = fwd_pass & rev_pass) on the row-aligned paired
    # matrices and read-weight across samples, instead of multiplying the two
    # marginal retention rates (which wrongly assumes forward/reverse failures
    # are independent). A defensive marginal-product fallback handles missing
    # paired matrices.
    exp_retention_combined <- reactive({
        paired <- rv_exp$paired_quality_data
        if (is.null(paired)) {
            fwd <- tryCatch(exp_retention_fwd(), error = function(e) NA)
            rev <- tryCatch(exp_retention_rev(), error = function(e) NA)
            if (is.na(fwd) || is.na(rev)) return(NA)
            return(round((fwd / 100) * (rev / 100) * 100, 1))
        }
        valid <- paired[!vapply(paired, is.null, logical(1))]
        if (length(valid) == 0) return(NA)
        # Debounced truncLen/maxEE reads (Section 8.8): this simulates
        # paired retention for every loaded sample, one of the heaviest
        # reactives driven by these sliders, so it must not re-run on every
        # intermediate value while the user is still dragging.
        per_sample <- lapply(valid, function(s) {
            simulate_paired_retention(
                s$quality_matrix_fwd, s$quality_matrix_rev,
                exp_truncLenF_db(), exp_truncLenR_db(),
                exp_maxEEF_db(), exp_maxEER_db(),
                trunc_q = RETENTION_DEFAULT_TRUNCQ
            )
        })
        n_sampled <- vapply(valid, function(s) as.numeric(s$n_sampled), numeric(1))
        summ <- summarize_paired_retention(per_sample, n_sampled)
        if (!is.finite(summ$read_weighted_retention)) return(NA)
        round(summ$read_weighted_retention * 100, 1)
    })

    # -------------------------------------------------------------------------
    # 8.13 VAL: Empirical DADA2 Validation
    # -------------------------------------------------------------------------
    # Runs the real DADA2 pipeline at the current Select settings on a
    # retention-stratified sample subset. The shared helper supplies the
    # surrogate prediction and representative-sample ranking.

    # Read-only fields showing the parameters the next validation run will use.
    # These mirror the live Select sliders but cannot be edited here.
    output$val_current_params <- renderUI({
        readonly_field <- function(label, value) {
            div(class = "mb-2",
                tags$label(class = "form-label small fw-semibold mb-1", label),
                tags$input(type = "text", class = "form-control form-control-sm",
                           value = as.character(value), readonly = "readonly",
                           tabindex = "-1"))
        }

        tagList(
            h6(class = "mb-2", "Parameters to validate"),
            layout_columns(
                col_widths = c(6, 6),
                readonly_field("truncLen Fwd", input$exp_truncLenF),
                readonly_field("truncLen Rev", input$exp_truncLenR),
                readonly_field("maxEE Fwd", input$exp_maxEEF),
                readonly_field("maxEE Rev", input$exp_maxEER)
            )
        )
    })

    # Selecting all samples disables the numeric field and uses the loaded
    # cohort size. Clearing it restores manual sample-count entry.
    observeEvent(input$optim_validate_use_all, {
        if (isTRUE(input$optim_validate_use_all)) {
            shinyjs::disable("optim_validate_n_samples")
            if (!is.null(rv_exp$files) && !is.null(rv_exp$files$fnFs)) {
                updateNumericInput(session, "optim_validate_n_samples",
                                    value = length(rv_exp$files$fnFs))
            }
        } else {
            shinyjs::enable("optim_validate_n_samples")
        }
    })

    validation_task <- ExtendedTask$new(function(cand, fnFs, fnRs, sample_names,
                                                  n_val, sample_indices, trunc_q,
                                                  paired_qd, amp_p99, tF, tR, mF, mR,
                                                  request_id) {
        promises::future_promise({
            res <- validate_candidates_with_dada2(
                candidates = cand, fnFs = fnFs, fnRs = fnRs,
                sample_names = sample_names, n_samples = n_val,
                trunc_q = trunc_q, sample_indices = sample_indices,
                progress = NULL
            )
            if (!is.null(res$per_sample)) {
                predicted <- vapply(res$per_sample$sample, function(sn) {
                    p <- if (!is.null(paired_qd)) paired_qd[[sn]] else NULL
                    if (is.null(p)) return(NA_real_)
                    r <- simulate_paired_retention(
                        p$quality_matrix_fwd, p$quality_matrix_rev, tF, tR, mF, mR,
                        trunc_q = trunc_q, amplicon_p99 = amp_p99,
                        min_overlap = DADA2_ABSOLUTE_MIN_OVERLAP
                    )
                    if (!is.finite(r$paired_retention) || !is.finite(r$mergeable_fraction)) return(NA_real_)
                    round(r$paired_retention * r$mergeable_fraction * 100, 1)
                }, numeric(1))
                res$per_sample$predicted_merged_pct <- predicted
                res$per_sample$real_merged_pct <- round(100 * res$per_sample$merged / pmax(res$per_sample$reads_in, 1), 1)
                res$per_sample$real_nonchim_pct <- round(100 * res$per_sample$nonchim / pmax(res$per_sample$reads_in, 1), 1)
                res$per_sample$gap_pp <- round(res$per_sample$predicted_merged_pct - res$per_sample$real_merged_pct, 1)
            }
            list(result = res, request_id = request_id)
        })
    })
    bind_task_button(validation_task, "optim_validate")

    observeEvent(input$optim_validate, {
        if (is.null(rv_exp$files) || is.null(rv_exp$files$fnFs)) {
            showNotification("Load quality profiles in Select before validating.", type = "warning"); return()
        }
        if (!requireNamespace("dada2", quietly = TRUE)) {
            rv_exp$validation_results <- list(
                summary = data.frame(
                    Note = "The 'dada2' package is not installed in this R environment, so empirical validation cannot run here.",
                    stringsAsFactors = FALSE),
                per_sample = NULL)
            add_console_msg("  X dada2 is not installed -- empirical validation cannot run here.", "error")
            showNotification("dada2 is not installed -- cannot run empirical validation here.", type = "error")
            return()
        }

        # Current parameters from the Quality Profiles sliders.
        tF <- input$exp_truncLenF; tR <- input$exp_truncLenR
        mF <- as.numeric(input$exp_maxEEF); mR <- as.numeric(input$exp_maxEER)
        paired_qd <- rv_exp$paired_quality_data

        # Surrogate prediction columns for the current params (so the validation
        # table's "Predicted merged %" = paired retention x mergeable fraction
        # still works). Computed across all loaded samples via the helper; the
        # mergeable fraction needs Visualize > Target Length Max.
        amp_p99 <- tryCatch({ r <- amplicon_length_range(); if (is.finite(r$max)) r$max else NA_real_ },
                            error = function(e) NA_real_)
        overall_ret <- NA_real_; merge_frac <- NA_real_
        if (!is.null(paired_qd)) {
            keep <- paired_qd[!vapply(paired_qd, is.null, logical(1))]
            if (length(keep) > 0) {
                per_sample <- lapply(keep, function(p)
                    simulate_paired_retention(
                        p$quality_matrix_fwd, p$quality_matrix_rev, tF, tR, mF, mR,
                        trunc_q = RETENTION_DEFAULT_TRUNCQ, amplicon_p99 = amp_p99,
                        min_overlap = DADA2_ABSOLUTE_MIN_OVERLAP))
                nsamp <- vapply(keep, function(p) nrow(p$quality_matrix_fwd), numeric(1))
                summ  <- tryCatch(summarize_paired_retention(per_sample, nsamp), error = function(e) NULL)
                if (!is.null(summ)) { overall_ret <- summ$read_weighted_retention; merge_frac <- summ$mergeable_fraction }
            }
        }
        cand <- data.frame(
            trunc_len_fwd = tF, trunc_len_rev = tR, max_ee_fwd = mF, max_ee_rev = mR,
            overall_paired_retention = overall_ret, mergeable_pair_fraction = merge_frac,
            stringsAsFactors = FALSE)

        # Select-all overrides the numeric field with the loaded sample count.
        n_total <- length(rv_exp$files$fnFs)
        n_val <- if (isTRUE(input$optim_validate_use_all)) n_total
                 else max(1, min(input$optim_validate_n_samples, n_total))
        add_console_msg("", "info")
        add_console_msg("=== Empirical DADA2 validation ===", "validation")
        add_console_msg(sprintf("  Validating truncLen %s/%s, maxEE %s/%s on %d sample(s). Real filterAndTrim/learnErrors/dada/mergePairs -- can take a while.",
                                tF, tR, mF, mR, n_val), "info")

        # Retention-stratified sample selection at the CURRENT params.
        retention_by_sample <- vapply(rv_exp$files$sample_names, function(sn) {
            p <- if (!is.null(paired_qd)) paired_qd[[sn]] else NULL
            if (is.null(p)) return(NA_real_)
            res_r <- tryCatch(simulate_paired_retention(
                p$quality_matrix_fwd, p$quality_matrix_rev, tF, tR, mF, mR,
                trunc_q = RETENTION_DEFAULT_TRUNCQ), error = function(e) NULL)
            if (is.null(res_r) || !is.finite(res_r$paired_retention)) return(NA_real_)
            round(res_r$paired_retention * 100, 1)
        }, numeric(1))
        val_selection  <- select_validation_samples(retention_by_sample, n_val)
        have_retention <- any(is.finite(val_selection$retention))
        if (have_retention) {
            add_console_msg("  Sample selection: retention-stratified across the cohort (ranked at the current truncLen/maxEE).", "info")
        } else {
            add_console_msg(paste0("  Sample selection: paired retention unavailable -- fell back to file order for the ", n_val, " sample(s)."), "info")
        }
        add_console_msg("  Samples used for validation:", "info")
        for (si in seq_len(nrow(val_selection))) {
            ret_txt <- if (is.finite(val_selection$retention[si]))
                sprintf("%.1f%% paired retention", val_selection$retention[si]) else "retention n/a"
            add_console_msg(sprintf("    - %s  (%s, %s)", val_selection$name[si],
                                    val_selection$role[si], ret_txt), "info")
        }

        rv_exp$validation_request_id <- rv_exp$validation_request_id + 1L
        request_id <- rv_exp$validation_request_id
        validation_task$invoke(
            cand, rv_exp$files$fnFs, rv_exp$files$fnRs, rv_exp$files$sample_names,
            n_val, val_selection$index, RETENTION_DEFAULT_TRUNCQ,
            paired_qd, amp_p99, tF, tR, mF, mR, request_id
        )
    })

    observe({
        completed <- validation_task$result()
        if (!identical(completed$request_id, rv_exp$validation_request_id)) return()
        res <- completed$result
        rv_exp$validation_results <- res
        if (!("Note" %in% names(res$summary))) {
            add_console_msg("  + DADA2 validation complete -- see the Validate tab.", "success")
        }
    })

    # Results are valid only for the exact filtering and amplicon settings used
    # to launch them. Changing any of those controls clears the table and causes
    # an older in-flight result to be ignored on completion.
    observeEvent(list(
        input$exp_truncLenF, input$exp_truncLenR,
        input$exp_maxEEF, input$exp_maxEER,
        input$vis_forward_primer_start, input$vis_reverse_primer_end,
        input$vis_forward_primer_length, input$vis_reverse_primer_length,
        input$vis_target_length_min, input$vis_target_length_max
    ), {
        rv_exp$validation_request_id <- rv_exp$validation_request_id + 1L
        rv_exp$validation_results <- NULL
    }, ignoreInit = TRUE)

    output$optim_validation_status <- renderUI({
        v <- rv_exp$validation_results
        if (is.null(v)) {
            return(tags$p(class = "text-muted mt-1 mb-0", style = "font-size:0.76rem;",
                          "Not yet run."))
        }
        if ("Note" %in% names(v$summary)) {
            return(div(class = "alert alert-warning py-1 px-2 mt-2 mb-0", style = "font-size:0.76rem;", v$summary$Note[1]))
        }
    })

    # Show a useful placeholder until a non-empty validation table exists.
    output$optim_validation_table_ui <- renderUI({
        v <- rv_exp$validation_results
        if (is.null(v) || is.null(v$per_sample) || nrow(v$per_sample) == 0) {
            return(tags$p(class = "text-muted mb-0", style = "font-size:0.82rem;",
                          "Not yet run. Set your sample count (or check \"Validate on all",
                          " loaded samples\") in the sidebar and click \"Validate with real",
                          " DADA2\" to see per-sample results here."))
        }
        DTOutput("optim_validation_table", fill = FALSE)
    })

    # Per-sample results preserve the selection order. Percentage and signed
    # difference columns use in-cell bars; Reads In remains an unscaled count.
    output$optim_validation_table <- renderDT({
        v <- rv_exp$validation_results
        if (is.null(v) || is.null(v$per_sample) || nrow(v$per_sample) == 0) return(NULL)
        ps <- v$per_sample

        display <- data.frame(
            Sample                 = ps$sample,
            `Reads In`             = format(ps$reads_in, big.mark = ","),
            `Filtered %`           = round(100 * ps$filtered  / pmax(ps$reads_in, 1), 1),
            `Denoised Fwd %`       = round(100 * ps$denoisedF / pmax(ps$reads_in, 1), 1),
            `Denoised Rev %`       = round(100 * ps$denoisedR / pmax(ps$reads_in, 1), 1),
            `Predicted Merged %`   = ps$predicted_merged_pct,
            `Real Merged %`        = ps$real_merged_pct,
            `Real Non-chim %`      = ps$real_nonchim_pct,
            `Difference (pp)`      = ps$gap_pp,
            check.names = FALSE, stringsAsFactors = FALSE
        )

        # Difference uses magnitude-dependent color tiers: <=5 pp teal,
        # 5-15 pp orange, and >15 pp red. The bar spans this run's observed
        # absolute range, with a floor of 1 to avoid division by zero.
        diff_max_abs <- max(1, max(abs(ps$gap_pp), na.rm = TRUE))
        diff_good_color     <- "#18bc9c"  # theme success teal --  <=5 pp: within noise
        diff_moderate_color <- "#f39c12"  # theme warning orange -- 5-15 pp: worth a look
        diff_high_color      <- "#e74c3c" # theme danger red --    >15 pp: surrogate meaningfully off
        diff_bar_js <- JS(sprintf(
            "isNaN(parseFloat(value)) ? '' : (function(v, maxAbs) {
                var absV = Math.abs(v);
                var pctTransparent = (maxAbs - Math.min(absV, maxAbs)) / maxAbs * 100;
                var tierColor = absV <= 5 ? '%s' : (absV <= 15 ? '%s' : '%s');
                return 'linear-gradient(270deg, transparent ' + pctTransparent + '%%, ' + tierColor + ' ' + pctTransparent + '%%)';
            })(value, %f)",
            diff_good_color, diff_moderate_color, diff_high_color, diff_max_abs
        ))

        datatable(
            display,
            # Show all rows without search chrome. Percentage widths plus the
            # fixed-layout CSS keep nine columns inside the main panel.
            options = list(pageLength = -1, dom = "t", ordering = TRUE,
                           autoWidth = FALSE,
                           columnDefs = list(
                               list(width = "13%", targets = 0),  # Sample
                               list(width = "9%",  targets = 1),  # Reads In
                               list(width = "10%", targets = 2),  # Filtered %
                               list(width = "11%", targets = 3),  # Denoised Fwd %
                               list(width = "11%", targets = 4),  # Denoised Rev %
                               list(width = "12%", targets = 5),  # Predicted Merged %
                               list(width = "11%", targets = 6),  # Real Merged %
                               list(width = "12%", targets = 7),  # Real Non-chim %
                               list(width = "11%", targets = 8)   # Difference (pp)
                           )),
            rownames = FALSE, class = "compact stripe hover"
        ) %>%
            # Bars are RIGHT-aligned (angle = 270 reverses the gradient so
            # the bar grows leftward from the right edge), same convention
            # as exp_sample_table's Forward/Reverse/Paired % columns.
            # Colors: Filtered % reuses PAIRED_COLOR (filtering is a
            # paired/combined step -- both directions must pass together);
            # Denoised Fwd/Rev % reuse this app's own FWD_COLOR/REV_COLOR
            # direction-identity colors; Predicted/Real Merged % use
            # blue/teal; Real Non-chim % uses a
            # neutral gray since chimera removal is not parameter-
            # controllable and is shown only for reference.
            formatStyle("Filtered %", background = styleColorBar(c(0, 100), PAIRED_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            formatStyle("Denoised Fwd %", background = styleColorBar(c(0, 100), FWD_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            formatStyle("Denoised Rev %", background = styleColorBar(c(0, 100), REV_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            formatStyle("Predicted Merged %", background = styleColorBar(c(0, 100), "#3498db", angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            formatStyle("Real Merged %", background = styleColorBar(c(0, 100), "#18bc9c", angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            formatStyle("Real Non-chim %", background = styleColorBar(c(0, 100), "#95a5a6", angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right") %>%
            # Signed, color-thresholded difference bar.
            formatStyle("Difference (pp)", background = diff_bar_js,
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right")
    })

    # -------------------------------------------------------------------------
    # 8.14 EXP: Export Generation
    # -------------------------------------------------------------------------
    # The same config and info reactives feed the live preview and saved
    # workbook. Parameters uses Step 5's exact variable names; Info stores
    # context and provenance. Retention fields remain blank until Select has
    # loaded quality profiles.

    # Human-readable platform value for the Info sheet.
    report_platform_note <- reactive({
        platform_label <- input$vis_platform
        if (is.null(platform_label)) {
            "not specified"
        } else {
            paste0(platform_label, " (", PLATFORM_READ_LENGTHS[[platform_label]], " bp nominal read length)")
        }
    })

    # Saving requires successfully loaded paired quality profiles as well as
    # physically possible overlap at the maximum expected target length.

    observe({
        ov <- tryCatch(exp_overlap_info(), error = function(e) NULL)
        shinyjs::toggleState("report_save_project",
                             condition = !is.null(rv_exp$files) &&
                                 !is.null(rv_exp$paired_quality_data) &&
                                 !is.null(ov) && ov$status$class != "overlap-fail")
    })

    report_config_table <- reactive({
        amp_range <- amplicon_length_range()

        data.frame(
            Parameter = c(
                "truncation_length_forward", "truncation_length_reverse",
                "max_expected_errors_forward", "max_expected_errors_reverse",
                "amplicon_min_length", "amplicon_max_length"
            ),
            # Explicit coercion protects the workbook's numeric Parameters
            # column if a UI control later returns a numeric-looking string.
            Value = c(
                input$exp_truncLenF, input$exp_truncLenR,
                as.numeric(input$exp_maxEEF), as.numeric(input$exp_maxEER),
                amp_range$min, amp_range$max
            ),
            Category = "config",
            Description = c(
                "DADA2 filterAndTrim() truncLen[1] (forward reads)",
                "DADA2 filterAndTrim() truncLen[2] (reverse reads)",
                "DADA2 filterAndTrim() maxEE[1] (forward reads)",
                "DADA2 filterAndTrim() maxEE[2] (reverse reads)",
                "Lower bound for the post-merge amplicon length filter from Visualize > Target Length; it applies to the primer-trimmed insert produced before DADA2",
                "Upper bound for the post-merge amplicon length filter from Visualize > Target Length; it applies to the primer-trimmed insert produced before DADA2"
            ),
            stringsAsFactors = FALSE
        )
    })

    report_info_table <- reactive({
        n_samples <- if (!is.null(rv_exp$files)) {
            length(rv_exp$files$sample_names)
        } else {
            rv_exp$n_detected_samples
        }

        primer_pair_label <- if (!is.null(input$vis_primer_preset_select) &&
                                  nchar(input$vis_primer_preset_select) > 0) {
            input$vis_primer_preset_select
        } else {
            "custom (see forward/reverse primer coordinates in Visualize)"
        }

        target_region_label <- paste(vis_covered_regions(), collapse = ", ")
        if (nchar(target_region_label) == 0) target_region_label <- "unknown"

        fwd_pct <- tryCatch(exp_retention_fwd(), error = function(e) NA)
        rev_pct <- tryCatch(exp_retention_rev(), error = function(e) NA)
        combined_pct <- tryCatch(exp_retention_combined(), error = function(e) NA)

        # Export the conservative overlap at the maximum expected target length,
        # matching the live status and save guard.
        overlap_bp <- tryCatch(exp_overlap_info()$overlap_at_p99, error = function(e) NA)

        # Paired-sampling metadata recorded when quality profiles are loaded.
        samp <- rv_exp$profile_sampling
        sampling_mode        <- if (!is.null(samp)) samp$mode else "(not recorded)"
        pairs_requested      <- if (!is.null(samp)) samp$requested else NA
        pairs_actual_median  <- if (!is.null(samp) && length(samp$actual) > 0) round(median(samp$actual, na.rm = TRUE)) else NA
        sampling_seed        <- if (!is.null(samp)) samp$base_seed else NA

        data.frame(
            Parameter = c(
                "generated_date", "sequencing_platform",
                "primer_pair", "target_region", "amplicon_size_bp",
                "estimated_overlap_bp",
                "input_directory", "file_naming_pattern", "n_samples",
                "sampling_mode", "sampling_pairs_requested",
                "sampling_pairs_actual_median", "sampling_random_seed",
                "truncation_quality",
                "retention_forward_pct",
                "retention_reverse_pct", "retention_combined_pct",
                "report_save_location"
            ),
            Value = c(
                format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
                report_platform_note(),
                primer_pair_label,
                target_region_label,
                vis_amplicon_coords()$size,
                if (!is.na(overlap_bp)) overlap_bp else NA,
                if (!is.null(rv_exp$selected_path)) {
                    to_project_relative_path(rv_exp$selected_path)
                } else {
                    "(not recorded)"
                },
                input$exp_pattern_choice,
                if (!is.null(n_samples)) n_samples else NA,
                sampling_mode,
                if (!is.na(pairs_requested)) pairs_requested else NA,
                if (!is.na(pairs_actual_median)) pairs_actual_median else NA,
                if (!is.na(sampling_seed)) sampling_seed else NA,
                RETENTION_DEFAULT_TRUNCQ,
                if (!is.na(fwd_pct)) fwd_pct else NA,
                if (!is.na(rev_pct)) rev_pct else NA,
                if (!is.na(combined_pct)) combined_pct else NA,
                REPORT_EXCEL_SAVE_LOCATION_LABEL
            ),
            Category = "info",
            Description = c(
                "Timestamp this table was generated",
                "Sequencing platform selected in Select > Read Length",
                "Primer pair selected in Visualize",
                "16S variable region(s) covered by the selected primers",
                "Expected amplicon size (bp), computed from the selected primer positions -- this is the full primer-to-primer span (includes both primers' own length); see amplicon_min_length/amplicon_max_length above for the primer-trimmed range actually used to filter merged reads",
                "Conservative forward/reverse overlap (bp) at the selected truncLen values, calculated against Visualize > Target Length Max; this is the value used for the GOOD/MODERATE/CRITICAL/LOW YIELD/FAIL status",
                "FASTQ directory loaded in Select > FASTQ Data, shown relative to the project root when possible",
                "Filename suffix selected in Select > FASTQ Data",
                "Number of samples detected in the loaded directory",
                "Paired-read sampling mode used when quality profiles were loaded",
                "Requested paired reads sampled per sample when quality profiles were loaded",
                "Median paired reads actually sampled per sample (may be lower than requested for small files)",
                "Base RNG seed for the paired sampler (each sample used base_seed + its index) -- makes the sampling reproducible",
                "Fixed DADA2 truncQ setting used by Step 4 and Step 5; recorded for provenance but not imported as a user-selectable Step 5 parameter",
                "Estimated forward-read retention (%) at the parameters above, from a per-file subsample",
                "Estimated reverse-read retention (%) at the parameters above, from a per-file subsample",
                "Estimated percentage of read pairs retained, requiring both the forward and reverse read to pass",
                "Path written when Save parameters is clicked; stored in the workbook so the original project location remains documented"
            ),
            stringsAsFactors = FALSE
        )
    })

    # Parameter-keyed numeric values used to restore numeric Excel cells in the
    # mixed-type Info sheet. The same cached reactives feed report_info_table().
    report_info_numeric_values <- reactive({
        n_samples <- if (!is.null(rv_exp$files)) {
            length(rv_exp$files$sample_names)
        } else {
            rv_exp$n_detected_samples
        }

        overlap_bp <- tryCatch(exp_overlap_info()$overlap_at_p99, error = function(e) NA)
        fwd_pct <- tryCatch(exp_retention_fwd(), error = function(e) NA)
        rev_pct <- tryCatch(exp_retention_rev(), error = function(e) NA)
        combined_pct <- tryCatch(exp_retention_combined(), error = function(e) NA)

        samp <- rv_exp$profile_sampling
        pairs_requested     <- if (!is.null(samp)) as.numeric(samp$requested) else NA_real_
        pairs_actual_median <- if (!is.null(samp) && length(samp$actual) > 0) as.numeric(round(median(samp$actual, na.rm = TRUE))) else NA_real_
        sampling_seed       <- if (!is.null(samp)) as.numeric(samp$base_seed) else NA_real_

        c(
            amplicon_size_bp = as.numeric(vis_amplicon_coords()$size),
            estimated_overlap_bp = as.numeric(overlap_bp),
            n_samples = as.numeric(if (!is.null(n_samples)) n_samples else NA),
            sampling_pairs_requested = pairs_requested,
            sampling_pairs_actual_median = pairs_actual_median,
            sampling_random_seed = sampling_seed,
            truncation_quality = as.numeric(RETENTION_DEFAULT_TRUNCQ),
            retention_forward_pct = as.numeric(fwd_pct),
            retention_reverse_pct = as.numeric(rev_pct),
            retention_combined_pct = as.numeric(combined_pct)
        )
    })

    # One display/export row order shared by the preview and both source tables.
    REPORT_ROW_ORDER <- c(
        "input_directory", "file_naming_pattern", "n_samples", "sequencing_platform",
        "target_region", "primer_pair", "amplicon_size_bp",
        "amplicon_min_length", "amplicon_max_length",
        "estimated_overlap_bp",
        "truncation_length_forward", "truncation_length_reverse",
        "max_expected_errors_forward", "max_expected_errors_reverse",
        "retention_forward_pct", "retention_reverse_pct",
        "retention_combined_pct",
        "sampling_mode", "sampling_pairs_requested",
        "sampling_pairs_actual_median", "sampling_random_seed",
        "truncation_quality",
        "report_save_location", "generated_date"
    )

    # Reorder any table containing Parameter and fail if a new row has not been
    # added to REPORT_ROW_ORDER.
    apply_report_row_order <- function(df) {
        missing_from_order <- setdiff(df$Parameter, REPORT_ROW_ORDER)
        if (length(missing_from_order) > 0) {
            stop("apply_report_row_order(): Parameter(s) not listed in REPORT_ROW_ORDER: ",
                 paste(missing_from_order, collapse = ", "))
        }
        df[order(match(df$Parameter, REPORT_ROW_ORDER)), , drop = FALSE]
    }

    report_table <- reactive({
        apply_report_row_order(rbind(report_config_table(), report_info_table()))
    })

    # Compact, unpaginated preview. The page handles overflow; the table does
    # not create a nested scroll container.
    output$report_table_preview <- renderDT({
        # Category becomes the preview's visible row-group label. The workbook
        # keeps config and info values on separate sheets.
        preview <- select(report_table(), -Category)
        preview$Section <- dplyr::case_when(
            preview$Parameter %in% c("input_directory", "file_naming_pattern", "n_samples",
                                     "sequencing_platform") ~ "Input data",
            preview$Parameter %in% c("target_region", "primer_pair", "amplicon_size_bp",
                                     "amplicon_min_length", "amplicon_max_length",
                                     "estimated_overlap_bp") ~ "Amplified target",
            preview$Parameter %in% c("truncation_length_forward", "truncation_length_reverse",
                                     "max_expected_errors_forward", "max_expected_errors_reverse",
                                     "truncation_quality") ~ "Filtering parameters",
            preview$Parameter %in% c("retention_forward_pct", "retention_reverse_pct",
                                     "retention_combined_pct") ~ "Estimated retention",
            TRUE ~ "Run record and reproducibility"
        )
        preview <- select(preview, Section, everything())

        datatable(
            preview,
            rownames = FALSE,
            extensions = "RowGroup",
            class = "compact stripe hover",
            options = list(
                dom = "t",
                paging = FALSE,
                ordering = FALSE,
                rowGroup = list(dataSrc = 0),
                columnDefs = list(list(visible = FALSE, targets = 0),
                                   list(width = "15%", targets = 1),
                                   list(width = "20%", targets = 2))
            )
        )
    })

    # Build a fresh complete workbook beside the destination, then publish it in
    # one replacement so obsolete sheets and partial saves cannot survive.
    write_report_excel_workbook <- function(config_table, info_table, info_numeric_values, path) {
        destination_dir <- dirname(path)
        temporary_path <- tempfile(".dada2-parameters-", tmpdir = destination_dir, fileext = ".xlsx")
        backup_path <- paste0(path, ".previous")
        on.exit(unlink(temporary_path, force = TRUE), add = TRUE)

        add_sheet_to_excel(temporary_path, REPORT_EXCEL_INFO_SHEET, info_table,
                            rownames = FALSE, overwrite = TRUE)
        fix_excel_numeric_typed_cells(temporary_path, REPORT_EXCEL_INFO_SHEET, info_table, info_numeric_values)

        add_sheet_to_excel(temporary_path, REPORT_EXCEL_PARAMETERS_SHEET, config_table,
                            rownames = FALSE, overwrite = TRUE)

        column_dictionary <- rbind(
            build_column_dictionary(
                sheet_name    = REPORT_EXCEL_INFO_SHEET,
                data          = info_table,
                descriptions  = REPORT_EXCEL_COLUMN_DESCRIPTIONS,
                workbook_path = temporary_path,
                rownames      = FALSE
            ),
            build_column_dictionary(
                sheet_name    = REPORT_EXCEL_PARAMETERS_SHEET,
                data          = config_table,
                descriptions  = REPORT_EXCEL_COLUMN_DESCRIPTIONS,
                workbook_path = temporary_path,
                rownames      = FALSE
            )
        )
        add_sheet_to_excel(temporary_path, "Column_Dictionary", column_dictionary,
                            rownames = FALSE, overwrite = TRUE)

        if (file.exists(backup_path)) unlink(backup_path, force = TRUE)
        if (file.exists(path) && !file.rename(path, backup_path)) {
            stop("Could not preserve the previous parameter workbook before replacement.")
        }
        if (!file.rename(temporary_path, path)) {
            if (file.exists(backup_path)) file.rename(backup_path, path)
            stop("Could not publish the completed parameter workbook; the previous workbook was restored.")
        }
        if (file.exists(backup_path)) unlink(backup_path, force = TRUE)

        invisible(path)
    }

    # Save parameters directly to Step 4's output folder so Step 5 can import
    # it. A reactive status remains visible after the notification fades.
    report_save_result <- reactiveVal(NULL)

    observeEvent(input$report_save_project, {
        if (is.null(rv_exp$files) || is.null(rv_exp$paired_quality_data)) {
            showNotification("Parameters were not saved: load the paired FASTQ quality profiles first.",
                             type = "error")
            return()
        }
        ov <- tryCatch(exp_overlap_info(), error = function(e) NULL)
        if (is.null(ov) || identical(ov$status$class, "overlap-fail")) {
            showNotification("Parameters were not saved: paired reads cannot overlap at the current truncLen settings.",
                             type = "error")
            return()
        }
        result <- tryCatch({
            dir.create(REPORT_EXCEL_OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
            output_path <- file.path(REPORT_EXCEL_OUTPUT_DIR, REPORT_EXCEL_FILENAME)
            write_report_excel_workbook(
                select(apply_report_row_order(report_config_table()), -Category),
                select(apply_report_row_order(report_info_table()), -Category),
                report_info_numeric_values(),
                output_path
            )
            list(success = TRUE, path = normalizePath(output_path, mustWork = FALSE))
        }, error = function(e) {
            list(success = FALSE, message = conditionMessage(e))
        })

        report_save_result(result)

        if (isTRUE(result$success)) {
            showNotification("Parameters saved to project.", type = "message")
        } else {
            showNotification(paste0("Could not save parameters: ", result$message), type = "error")
        }
    })

    output$report_save_status <- renderUI({
        result <- report_save_result()
        if (is.null(result)) return(NULL)

        if (isTRUE(result$success)) {
            div(class = "alert alert-success py-1 px-2 mt-2 mb-0",
                style = "font-size:0.74rem; line-height:1.35;",
                icon("check"), " Saved to ", tags$code(result$path))
        } else {
            div(class = "alert alert-danger py-1 px-2 mt-2 mb-0",
                style = "font-size:0.74rem; line-height:1.35;",
                icon("triangle-exclamation"), " ", result$message)
        }
    })

    # -------------------------------------------------------------------------
    # 8.15 EXP: Output Renderers
    # -------------------------------------------------------------------------
    
    output$exp_retention_F <- renderText({
        val <- tryCatch(exp_retention_fwd(), error = function(e) NA)
        if (is.na(val)) "--" else paste0(val, "%")
    })
    
    output$exp_retention_R <- renderText({
        val <- tryCatch(exp_retention_rev(), error = function(e) NA)
        if (is.na(val)) "--" else paste0(val, "%")
    })
    
    output$exp_retention_C <- renderText({
        val <- tryCatch(exp_retention_combined(), error = function(e) NA)
        if (is.na(val)) "--" else paste0(val, "%")
    })
    
    output$exp_n_samples <- renderText({
        if (is.null(rv_exp$files)) return("0 samples")
        paste(length(rv_exp$files$fnFs), "samples")
    })
    
    # Two stacked single-axis plots per direction:
    #   (1) make_quality_position_plot(): median Q vs position (truncLen view),
    #       with the truncLen line, shaded retained region + retained %.
    #   (2) make_retention_maxee_plot(): retained % vs maxEE 0-10 (maxEE view),
    #       with your-maxEE guide + marker. This is the ECDF of per-read EE at
    #       the current truncLen.

    # Median-Q-vs-position plot. The x-axis spans EXACTLY [0, x_max] where x_max
    # is the truncLen slider's own max (the selected platform's read length), so
    # the dashed truncLen line sits at the same horizontal fraction as the slider
    # handle below it, making the slider the visible x scale. The middle retained-
    # reads plot shows the aggregate percentage. fixedrange prevents zoom/pan from
    # breaking slider-to-line alignment.
    make_quality_position_plot <- function(stats, truncLen, x_max, qline_color) {
        plot_data <- stats[order(stats$Position), ]
        if (nrow(plot_data) == 0) return(NULL)
        if (!is.finite(x_max) || x_max <= 0) x_max <- max(plot_data$Position, na.rm = TRUE)
        plot_ly() %>%
            add_lines(x = plot_data$Position, y = plot_data$Median,
                      line = list(color = qline_color, width = 2.5), name = "Median Q",
                      hovertemplate = "Median Q<br>pos %{x} bp<br>Q %{y:.1f}<extra></extra>") %>%
            layout(
                # x-axis title + tick labels hidden (ticks = "" too), so the
                # truncLen slider directly below is the only visible x scale. Range
                # is EXACTLY [0, x_max] = the slider's own [min, max], and fixedrange
                # disables horizontal zoom so the plot range can never diverge from
                # the slider range.
                xaxis = list(title = NULL, range = c(0, x_max), zeroline = FALSE,
                             showticklabels = FALSE, ticks = "", fixedrange = TRUE),
                # Keep the y-axis, grid, curve, and hover visible.
                yaxis = list(title = "Quality Score (median)", range = c(0, 45),
                             titlefont = list(color = qline_color, size = 11),
                             showgrid = TRUE, gridcolor = "#e9ecef"),
                # The dashed truncLen line is a LAYOUT SHAPE (shapes[0]) rather than
                # a trace, so the truncLen observer can move it with a single
                # plotlyProxyInvoke("relayout", shapes[0].x0/x1) call -- updating
                # only the line without re-rendering the curve.
                shapes = list(list(type = "line", xref = "x", yref = "paper",
                                   x0 = truncLen, x1 = truncLen, y0 = 0, y1 = 1,
                                   line = list(color = TRUNC_LINE_COLOR, width = 2,
                                               dash = "dash"))),
                showlegend = FALSE,
                # Shared alignment margins (QP_MARGIN_*): the SAME numbers feed the
                # .axis-slider CSS variables, so the slider track and this plot's
                # data region begin/end at the same pixel.
                margin = list(l = QP_MARGIN_L, r = QP_MARGIN_R,
                              t = QP_MARGIN_T, b = QP_MARGIN_B)
            )
    }

    # Retention-vs-maxEE plot (ECDF of per-read total EE at the current truncLen).
    # x-axis spans EXACTLY [0, 10], matching the maxEE slider's own range, so the
    # dashed maxEE line sits at the same horizontal fraction as the slider handle
    # below it. The maxEE slider is the visible x scale; aggregate retention is
    # shown in the middle retained-reads plot.
    make_retention_maxee_plot <- function(ee, maxEE, curve_color = MAXEE_LINE_COLOR,
                                          direction = c("Fwd", "Rev")) {
        # Direction identifies the y-axis title.
        direction <- match.arg(direction)
        if (length(ee) == 0) return(NULL)
        me_grid   <- seq(0, 10, by = 0.1)
        ret_curve <- vapply(me_grid, function(m) 100 * mean(ee <= m), numeric(1))
        cur_ret   <- 100 * mean(ee <= maxEE)
        p <- plot_ly() %>%
            add_lines(x = me_grid, y = ret_curve,
                      line = list(color = curve_color, width = 2.5), name = "retained vs maxEE",
                      hovertemplate = "maxEE %{x:.1f}<br>retained %{y:.1f}%<extra></extra>")
        if (is.finite(maxEE)) {
            p <- p %>%
                add_segments(x = maxEE, xend = maxEE, y = 0, yend = 100,
                             line = list(color = "#e74c3c", width = 2, dash = "dash"),
                             showlegend = FALSE, hoverinfo = "none") %>%
                add_markers(x = maxEE, y = cur_ret,
                            marker = list(color = "#e74c3c", size = 9,
                                          line = list(color = "white", width = 1.5)),
                            showlegend = FALSE, hoverinfo = "text",
                            text = sprintf("maxEE %.1f<br>retained %.1f%%", maxEE, cur_ret))
        }
        p %>% layout(
            # x-axis hidden (the maxEE slider below is the x scale); range [0,10]
            # matches the slider, and the shared QP_MARGIN_* margins make the slider
            # track span the plot's data region (same alignment as the quality plot).
            xaxis = list(title = NULL, range = c(0, 10), zeroline = FALSE,
                         showticklabels = FALSE, ticks = "", fixedrange = TRUE),
            yaxis = list(title = sprintf("Retained %s reads %%", direction), range = c(0, 105),
                         titlefont = list(color = curve_color, size = 11)),
            showlegend = FALSE,
            margin = list(l = QP_MARGIN_L, r = QP_MARGIN_R,
                          t = QP_MARGIN_T, b = QP_MARGIN_B)
        )
    }

    # The quality plots read input$exp_truncLen* only via isolate(): the truncLen
    # value sets the dashed line's INITIAL position when the plot (re)renders on a
    # data change, but a subsequent slider move does NOT invalidate this render --
    # the line is moved instead by the lightweight plotlyProxy observers below, so
    # the curve is never recomputed just to reposition the line. exp_platform_max()
    # (the x-range = the slider's max) stays a live dependency so the plot re-fits
    # if the platform changes.
    output$exp_plot_fwd_quality <- renderPlotly({
        req(rv_exp$aggregated_fwd)
        make_quality_position_plot(rv_exp$aggregated_fwd$position_stats,
                                   isolate(input$exp_truncLenF),
                                   exp_platform_max(), FWD_COLOR)
    })
    # Debounced maxEE reads (Section 8.8): this is a full renderPlotly() of
    # the retention-vs-maxEE curve, so it should not re-render on every
    # intermediate maxEE tick while the user is still dragging.
    output$exp_plot_fwd_retmaxee <- renderPlotly({
        req(rv_exp$quality_data_fwd)
        make_retention_maxee_plot(exp_read_ee_fwd(), exp_maxEEF_db(), FWD_COLOR, "Fwd")
    })
    output$exp_plot_rev_quality <- renderPlotly({
        req(rv_exp$aggregated_rev)
        make_quality_position_plot(rv_exp$aggregated_rev$position_stats,
                                   isolate(input$exp_truncLenR),
                                   exp_platform_max(), REV_COLOR)
    })
    output$exp_plot_rev_retmaxee <- renderPlotly({
        req(rv_exp$quality_data_rev)
        make_retention_maxee_plot(exp_read_ee_rev(), exp_maxEER_db(), REV_COLOR, "Rev")
    })

    # Move only the dashed truncLen line while dragging; the debounced retention
    # reactives update after the user pauses.
    observeEvent(input$exp_truncLenF, {
        plotlyProxy("exp_plot_fwd_quality", session) %>%
            plotlyProxyInvoke("relayout",
                              list("shapes[0].x0" = input$exp_truncLenF,
                                   "shapes[0].x1" = input$exp_truncLenF))
    }, ignoreInit = TRUE)
    observeEvent(input$exp_truncLenR, {
        plotlyProxy("exp_plot_rev_quality", session) %>%
            plotlyProxyInvoke("relayout",
                              list("shapes[0].x0" = input$exp_truncLenR,
                                   "shapes[0].x1" = input$exp_truncLenR))
    }, ignoreInit = TRUE)

    # Retained-reads bar plot for the middle Quality Profiles pane: Forward,
    # Paired, and Reverse percentages at the current truncLen/maxEE. Forward and
    # Reverse are the marginal per-direction retention; Paired is the true joint
    # paired retention (both mates pass), so it is never higher than the smaller
    # marginal. Bars use the same direction identity colors as the rest of the
    # app (FWD teal, PAIRED violet, REV orange). The exact % sits above each bar.
    output$exp_retbar_plot <- renderPlotly({
        fwd    <- tryCatch(exp_retention_fwd(),      error = function(e) NA)
        paired <- tryCatch(exp_retention_combined(), error = function(e) NA)
        rev    <- tryCatch(exp_retention_rev(),      error = function(e) NA)
        cats   <- c("Fwd", "Paired", "Rev")
        vals   <- suppressWarnings(as.numeric(c(fwd, paired, rev)))
        cols   <- c(FWD_COLOR, PAIRED_COLOR, REV_COLOR)
        labs   <- ifelse(is.finite(vals), sprintf("%.1f%%", vals), "NA")
        # Explicit numeric centers and trace-level width keep all three bars
        # inside the narrow middle column. The percentage sits above each bar,
        # so the y-axis title and tick labels are hidden; x tick labels remain
        # horizontal. customdata supplies category names to hover text because
        # x is numeric.
        x_pos <- c(0, 0.7, 1.4)
        plot_ly() %>%
            add_trace(
                x = x_pos, y = vals, type = "bar",
                width = 0.45, marker = list(color = cols),
                text = labs, textposition = "outside", textfont = list(size = 11),
                customdata = cats, cliponaxis = FALSE,
                hovertemplate = "%{customdata}<br>retained %{y:.1f}%<extra></extra>"
            ) %>%
            layout(
                xaxis = list(title = "", tickfont = list(size = 12), tickangle = 0,
                             tickmode = "array", tickvals = x_pos, ticktext = cats,
                             range = c(min(x_pos) - 0.45, max(x_pos) + 0.45),
                             fixedrange = TRUE),
                yaxis = list(title = "", showticklabels = FALSE, range = c(0, 115),
                             showgrid = FALSE, zeroline = FALSE, fixedrange = TRUE),
                margin = list(l = 6, r = 6, t = 16, b = 26),
                showlegend = FALSE,
                paper_bgcolor = "white", plot_bgcolor = "white"
            ) %>%
            config(displayModeBar = FALSE)
    })

    # Per-sample table. Debounced truncLen/maxEE reads throughout (Section
    # 8.8): this recomputes calculate_retention_vectorized() and
    # simulate_paired_retention() for EVERY loaded sample, the
    # single most expensive reactive driven by these sliders, so it must
    # only fire once the user pauses rather than on every drag tick.
    output$exp_sample_table <- renderDT({
        req(rv_exp$quality_data_fwd, rv_exp$quality_data_rev, rv_exp$files)

        sample_retention <- data.frame(SampleID = rv_exp$files$sample_names,
                                       stringsAsFactors = FALSE)

        sample_retention$Fwd_Pct <- sapply(rv_exp$quality_data_fwd, function(qd) {
            if (is.null(qd)) return(NA)
            tryCatch(
                calculate_retention_vectorized(qd$qual_matrix, qd$cum_ee_matrix,
                                               exp_truncLenF_db(), exp_maxEEF_db(),
                                               RETENTION_DEFAULT_TRUNCQ),
                error = function(e) NA
            )
        })

        sample_retention$Rev_Pct <- sapply(rv_exp$quality_data_rev, function(qd) {
            if (is.null(qd)) return(NA)
            tryCatch(
                calculate_retention_vectorized(qd$qual_matrix, qd$cum_ee_matrix,
                                               exp_truncLenR_db(), exp_maxEER_db(),
                                               RETENTION_DEFAULT_TRUNCQ),
                error = function(e) NA
            )
        })

        # Combined % requires both reads in each row-aligned pair to pass.
        # Fwd % and Rev % remain marginal, single-direction estimates. The
        # product is used only if paired data is unavailable.
        paired <- rv_exp$paired_quality_data
        sample_retention$Combined <- vapply(rv_exp$files$sample_names, function(sn) {
            p <- if (!is.null(paired)) paired[[sn]] else NULL
            if (is.null(p)) {
                fwd <- sample_retention$Fwd_Pct[match(sn, rv_exp$files$sample_names)]
                rev <- sample_retention$Rev_Pct[match(sn, rv_exp$files$sample_names)]
                return(round((fwd / 100) * (rev / 100) * 100, 1))
            }
            res <- tryCatch(
                simulate_paired_retention(
                    p$quality_matrix_fwd, p$quality_matrix_rev,
                    exp_truncLenF_db(), exp_truncLenR_db(),
                    exp_maxEEF_db(), exp_maxEER_db(),
                    trunc_q = RETENTION_DEFAULT_TRUNCQ
                ),
                error = function(e) NULL
            )
            if (is.null(res) || !is.finite(res$paired_retention)) return(NA_real_)
            round(res$paired_retention * 100, 1)
        }, numeric(1))

        colnames(sample_retention) <- c("Sample ID", "Forward %", "Reverse %", "Paired %")

        datatable(
            sample_retention,
            # Show all rows at natural height with no internal scroll or search
            # chrome. Constrain Sample ID so the percentage columns use the width.
            options = list(pageLength = -1, dom = 't', ordering = TRUE,
                           autoWidth = TRUE,
                           columnDefs = list(list(width = "110px", targets = 0))),
            rownames = FALSE, class = "compact stripe hover"
        ) %>%
            # Keep value labels dark in every retention column. A bar can end
            # before the right-aligned number, so switching to white solely on
            # a high percentage can leave white text on the unfilled row
            # background. The dark label has sufficient contrast against all
            # three bar colors as well as the unfilled background.
            formatStyle("Forward %", background = styleColorBar(c(0, 100), FWD_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right",
                        color = "#212529",
                        fontWeight = "600") %>%
            formatStyle("Reverse %", background = styleColorBar(c(0, 100), REV_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right",
                        color = "#212529",
                        fontWeight = "600") %>%
            formatStyle("Paired %", background = styleColorBar(c(0, 100), PAIRED_COLOR, angle = 270),
                        backgroundSize = "98% 88%", backgroundRepeat = "no-repeat",
                        backgroundPosition = "right",
                        color = "#212529",
                        fontWeight = "600")
    })
}


# =============================================================================
# SECTION 9: RUN THE APPLICATION
# =============================================================================

# future::plan() is process-global. Activate the bounded multisession plan only
# for the app lifecycle and restore the caller's prior plan when the app stops.
shinyApp(
    ui = ui,
    server = server,
    onStart = function() {
        previous_future_plan <- future::plan()
        future::plan(future::multisession, workers = n_cores)
        shiny::onStop(function() {
            future::plan(previous_future_plan)
        })
    }
)
