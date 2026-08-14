################################################################################
# Script: add_sheet_to_excel_function.R
# Purpose:
#   - Define a reusable function `add_sheet_to_excel` for safely appending a new
#     worksheet containing a given dataset into an existing or new Excel workbook.
#   - Automatically apply styling (bold, centered headers), usability features
#     (freeze first row, auto-sized columns, filterable table), and a
#     colour-coded sheet tab to make the output suitable for publication-quality
#     reports, data analysis, or sharing with collaborators.
#
# NOTE: This script previously also saved an `add_sheet_to_excel_function.RData`
# copy of the function to the current working directory as a side effect of
# being source()'d. That was removed (2026-07-27): sourcing this .R file is
# already the reproducible way this function is loaded everywhere in this
# repository (nothing ever load()'d the .RData), and because RStudio's
# working directory during a knit is the currently-open .Rmd's own folder,
# that save() call was regenerating a stray, untracked .RData file inside
# R/notebooks/ every time any notebook sourcing this script was knitted.
#
# Change log (2026-07-29):
#   - Added `tab_colour` argument: colours each sheet tab automatically based
#     on keyword matches against `sheet_name` (e.g. "taxon" -> blue,
#     "differential"/"daa" -> orange, "biomarker"/"lefse" -> green,
#     "pathway"/"function" -> purple, "metadata"/"sample" -> grey), falling
#     back to a cycling default palette for sheet names that match no
#     keyword. Pass an explicit hex colour to override the automatic choice,
#     or NULL to disable tab colouring entirely.
#   - Added `overwrite` argument: previously any attempt to add a sheet name
#     that already existed in the workbook was an unconditional stop(); set
#     `overwrite = TRUE` to replace the existing sheet instead of erroring.
#   - Sheet names are now sanitised automatically: Excel forbids the
#     characters : \ / ? * [ ] in sheet names and truncates (or corrupts)
#     names longer than 31 characters, so both constraints are now enforced
#     here, with a warning() raised whenever a name had to be changed.
#   - Column-level documentation (per-header-cell popup comments) was added
#     here in an earlier revision, then moved out to
#     `build_column_dictionary_function.R` on the same day: that script
#     already owns the `descriptions` vector used to build each workbook's
#     trailing `Column_Dictionary` sheet, so reusing the same vector there to
#     also write the inline header popups keeps a single source of truth per
#     sheet instead of two description lists that could drift apart. This
#     function stays a generic "write and style one sheet" utility with no
#     knowledge of column-level documentation.
#   - Fixed a bug introduced earlier the same day: the existing-sheet check
#     was changed to call `getSheetNames(wb)` on the in-memory Workbook
#     object, but `getSheetNames()` only accepts a path to an already-saved
#     .xlsx file -- passing a Workbook object made every call fail
#     immediately with "invalid 'file' argument". Existing sheet names are
#     now read from `workbook_path` itself (only when that file already
#     exists), matching the original, working behaviour.
################################################################################

# ==============================================================================
# Dependencies
# ==============================================================================
# All openxlsx calls are namespace-qualified so sourcing this helper does not
# modify the caller's package search path.

# ==============================================================================
# Sheet-tab colour configuration
# ==============================================================================
# The keyword and fallback palettes are defined inside the function below so
# sourcing this script does not add configuration objects to the caller.
# ==============================================================================
# Define the function: add_sheet_to_excel
# ==============================================================================
add_sheet_to_excel <- function(workbook_path,
                                sheet_name,
                                data,
                                rownames = FALSE,
                                overwrite = FALSE,
                                tab_colour = "auto") {

  if (!is.character(workbook_path) || length(workbook_path) != 1L ||
      is.na(workbook_path) || !nzchar(trimws(workbook_path))) {
    stop("add_sheet_to_excel(): `workbook_path` must be one non-empty, non-missing character path.")
  }
  if (!is.character(sheet_name) || length(sheet_name) != 1L ||
      is.na(sheet_name) || !nzchar(trimws(sheet_name))) {
    stop("add_sheet_to_excel(): `sheet_name` must be one non-empty, non-missing character value.")
  }
  if (!is.logical(rownames) || length(rownames) != 1L || is.na(rownames)) {
    stop("add_sheet_to_excel(): `rownames` must be TRUE or FALSE.")
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("add_sheet_to_excel(): `overwrite` must be TRUE or FALSE.")
  }
  if (!is.null(tab_colour) &&
      (!is.character(tab_colour) || length(tab_colour) != 1L ||
       is.na(tab_colour) || !nzchar(tab_colour))) {
    stop("add_sheet_to_excel(): `tab_colour` must be NULL or one non-empty, non-missing character value.")
  }

  if (is.null(dim(data)) || ncol(data) < 1L) {
    stop("add_sheet_to_excel(): `data` must contain at least one column.")
  }

  # ---- Determine which sheet names already exist, then load/create the workbook ----
  # getSheetNames() only accepts a path to an already-saved .xlsx file, not an
  # in-memory Workbook object, so this must be evaluated from `workbook_path`
  # itself, before (or instead of) loading it into memory.
  existing_sheet_names <- if (file.exists(workbook_path)) openxlsx::getSheetNames(workbook_path) else character(0)

  if (file.exists(workbook_path)) {
    wb <- openxlsx::loadWorkbook(workbook_path)
  } else {
    wb <- openxlsx::createWorkbook()
  }

  # ---- Sanitise the sheet name --------------------------------------------------
  # Excel forbids the characters : \ / ? * [ ] in sheet names and caps the
  # length at 31 characters; exceeding either silently produces a corrupted or
  # truncated workbook rather than a clear error, so both are normalised here.
  forbidden_sheet_name_characters <- c(":", "\\", "/", "?", "*", "[", "]")
  sanitised_sheet_name <- sheet_name
  for (forbidden_character in forbidden_sheet_name_characters) {
    sanitised_sheet_name <- gsub(forbidden_character, "_", sanitised_sheet_name, fixed = TRUE)
  }
  if (nchar(sanitised_sheet_name) > 31) {
    sanitised_sheet_name <- substr(sanitised_sheet_name, 1, 31)
  }
  if (!identical(sanitised_sheet_name, sheet_name)) {
    warning(sprintf("Sheet name '%s' was not a valid Excel sheet name and was sanitised to '%s'.",
                     sheet_name, sanitised_sheet_name))
  }
  sheet_name <- sanitised_sheet_name

  # ---- Handle a sheet name that already exists in the workbook -------------------
  matching_sheet_index <- match(tolower(sheet_name), tolower(existing_sheet_names))
  if (!is.na(matching_sheet_index)) {
    existing_sheet_name <- existing_sheet_names[matching_sheet_index]
    if (isTRUE(overwrite)) {
      openxlsx::removeWorksheet(wb, existing_sheet_name)
    } else {
      stop(paste("Sheet with name", existing_sheet_name,
                 "already exists in the workbook. Set overwrite = TRUE to replace it."))
    }
  }

  # ---- Resolve the sheet tab colour -----------------------------------------------
  # tab_colour = "auto" (the default) looks up `sheet_name` against the keyword
  # palette above, falling back to a cycling default palette when no keyword
  # matches. Pass an explicit hex string (e.g. "#1F77B4") to force a specific
  # colour, or NULL to leave the tab uncoloured.
  resolved_tab_colour <- tab_colour
  if (identical(tab_colour, "auto")) {
    tab_colour_keyword_palette <- c(
      "taxon" = "#4E79A7", "differential" = "#F28E2B", "daa" = "#F28E2B",
      "biomarker" = "#59A14F", "lefse" = "#59A14F", "pathway" = "#B07AA1",
      "function" = "#B07AA1", "metadata" = "#BAB0AC", "sample" = "#BAB0AC"
    )
    tab_colour_fallback_cycle <- c("#76B7B2", "#EDC949", "#AF7AA1", "#FF9DA7", "#9C755F")
    keyword_matches <- vapply(
      names(tab_colour_keyword_palette),
      function(keyword) grepl(keyword, sheet_name, ignore.case = TRUE),
      logical(1)
    )
    if (any(keyword_matches)) {
      resolved_tab_colour <- tab_colour_keyword_palette[[names(keyword_matches)[keyword_matches][1]]]
    } else {
      # Same "invalid 'file' argument" pitfall as the existing-sheet check
      # above: getSheetNames() cannot take a Workbook object, so the sheet
      # count is derived from `existing_sheet_names` (read from disk before
      # this sheet was added) rather than re-querying `wb` here.
      next_sheet_position <- length(existing_sheet_names) + 1
      cycle_index <- ((next_sheet_position - 1) %% length(tab_colour_fallback_cycle)) + 1
      resolved_tab_colour <- tab_colour_fallback_cycle[cycle_index]
    }
  }

  # ---- Add the worksheet, applying the resolved tab colour (if any) --------------
  if (is.null(resolved_tab_colour)) {
    openxlsx::addWorksheet(wb, sheet_name)
  } else {
    openxlsx::addWorksheet(wb, sheet_name, tabColour = resolved_tab_colour)
  }

  # ---- Write the data as a filterable, styled Excel table -------------------------
  openxlsx::writeDataTable(wb = wb,
            sheet = sheet_name,
            x = data,
            colNames = TRUE,
            rowNames = rownames,
            tableStyle = "TableStyleLight9",
            withFilter = TRUE)

  # ---- Style the header row: bold, centered ----------------------------------------
  header_style <- openxlsx::createStyle(
    textDecoration = "bold",
    halign = "center"
  )
  written_column_count <- ncol(data) + as.integer(isTRUE(rownames))
  written_columns <- seq_len(written_column_count)
  openxlsx::addStyle(wb, sheet_name, header_style, rows = 1, cols = written_columns, gridExpand = TRUE)

  # ---- Usability: freeze the header row, auto-size columns --------------------------
  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet_name, cols = written_columns, widths = "auto")

  # ---- Save the workbook to disk -------------------------------------------------------
  destination_directory <- dirname(workbook_path)
  if (!dir.exists(destination_directory)) {
    stop("add_sheet_to_excel(): parent directory does not exist: ", destination_directory)
  }
  temporary_path <- tempfile(".workbook-", tmpdir = destination_directory, fileext = ".xlsx")
  on.exit(unlink(temporary_path), add = TRUE)
  openxlsx::saveWorkbook(wb, temporary_path, overwrite = TRUE)
  if (!file.rename(temporary_path, workbook_path)) {
    stop("add_sheet_to_excel(): could not atomically replace workbook: ", workbook_path)
  }

  message(paste("Sheet", sheet_name, "added successfully to the workbook:", workbook_path))
}

# ==============================================================================
# Example usage (uncomment to test)
# ==============================================================================
# add_sheet_to_excel(
#   workbook_path = "example_workbook.xlsx",
#   sheet_name    = "differential_abundance",
#   data          = my_dataframe,
#   rownames      = TRUE,
#   overwrite     = TRUE,
#   tab_colour    = "auto"
# )
