################################################################################
# Script: build_column_dictionary_function.R
# Purpose:
#   - Define a reusable function `build_column_dictionary()` that documents
#     every column of a data frame being exported to an Excel worksheet, as a
#     tidy Sheet / Column / Explanation table.
#   - Intended to be called once per sheet actually written with
#     `add_sheet_to_excel()`, with the resulting rows from every sheet in a
#     workbook row-bound together and written as a final trailing
#     "Column_Dictionary" sheet in that same workbook — so every workbook this
#     project produces documents its own columns, without a hand-maintained
#     description list that can silently drift out of sync with the real
#     exported data.
#   - Optionally (via `workbook_path`), also attaches each column's
#     explanation directly to that column's header cell, in that column's own
#     sheet, as a hover/click Excel comment -- a lightweight in-file data
#     dictionary alongside the trailing summary sheet. Both are driven by the
#     same `descriptions` vector, so there is exactly one place per sheet
#     where column explanations are written, rather than two lists that can
#     drift apart.
#
# Change log (2026-07-29):
#   - Added the optional `workbook_path` / `rownames` arguments and the
#     header-cell popup-comment step described above. This capability
#     previously lived in `add_sheet_to_excel_function.R` as a
#     `column_descriptions` argument; it was moved here because this script
#     already owns the validated `descriptions` vector for each sheet, so
#     reusing it here avoids maintaining the same column-explanation list in
#     two places. `add_sheet_to_excel()` remains a generic "write and style
#     one sheet" utility with no knowledge of column-level documentation.
#   - Popup comment boxes are now sized dynamically from each explanation's
#     length (see `estimate_comment_box_height_units()` below) instead of a
#     single small fixed size. Excel's legacy cell-comment boxes do not
#     auto-fit to their text, so a box sized for a short explanation clipped
#     longer ones instead of wrapping or scrolling to show the rest.
################################################################################

# ==============================================================================
# Load necessary library
# ==============================================================================
# Needed for the optional header-cell popup-comment step (loadWorkbook(),
# getSheetNames(), writeComment(), createComment(), saveWorkbook()). The
# dictionary-row-only code path (workbook_path = NULL) does not touch any
# Excel file itself and would not strictly need this, but it is loaded
# unconditionally for a single, predictable dependency list.
# All openxlsx calls are namespace-qualified so sourcing this helper does not
# modify the caller's package search path.

# ==============================================================================
# Comment box sizing
# ==============================================================================
# createComment()'s `width`/`height` arguments use the same units Excel
# itself uses for legacy cell-comment boxes: `width` in multiples of the
# workbook's default column width (~64 px each) and `height` in multiples of
# the default row height (~20 px each). Excel does not auto-fit these boxes
# to their text, so a box that is too small clips the explanation instead of
# wrapping or scrolling to reveal the rest.
#
# The box width is kept fixed so every popup
# in a workbook looks visually consistent; only the height is estimated per
# explanation, from its character count, so a one-sentence description gets
# a compact box while a long, multi-clause explanation gets a taller one
# instead of being clipped. The characters-per-line and lines-per-height-unit
# constants are approximations for Excel's default comment font (Tahoma,
# 8pt) -- if popups still feel too small or too large in practice, adjust
# the function-local sizing values.
# ==============================================================================
# Define the function: build_column_dictionary
# ==============================================================================
# Function arguments:
#   - sheet_name    : Character scalar, name of the worksheet `data` is (or
#                     will be) written to, exactly as passed to
#                     add_sheet_to_excel().
#   - data          : Data frame actually written to `sheet_name` — only its
#                     colnames(data) are used; the cell values themselves are
#                     never inspected.
#   - descriptions  : Named character vector. Names are column names as they
#                     appear in `data`; values are the explanation text for
#                     that column.
#   - workbook_path : Optional character scalar, path to the .xlsx workbook
#                     that already contains `sheet_name` (i.e. the same
#                     `workbook_path` already passed to add_sheet_to_excel()
#                     for this sheet). When supplied, each column's
#                     explanation is additionally written as a hover/click
#                     comment on that column's header cell. Defaults to NULL,
#                     which skips this step entirely and preserves this
#                     function's original, dictionary-rows-only behaviour.
#   - rownames      : Logical scalar, must match the `rownames` argument used
#                     in the matching add_sheet_to_excel() call for this
#                     sheet, so header columns line up correctly when
#                     attaching popup comments. Ignored when `workbook_path`
#                     is NULL. Defaults to FALSE.
#
# Workflow:
#   1. Read the actual column names from `data` (not from `descriptions`),
#      so the dictionary is always driven by what is really in the workbook.
#   2. Confirm every actual column has a matching entry in `descriptions`,
#      failing loudly and specifically if any column would otherwise be
#      exported to the workbook without documentation.
#   3. Assemble and return the Sheet / Column / Explanation rows for this
#      sheet, in the same left-to-right column order as the exported data.
#   4. If `workbook_path` is supplied, reopen that already-saved workbook and
#      attach each column's explanation as a popup comment on its header
#      cell in `sheet_name`, then resave.
#
# Build a Sheet / Column / Explanation Dictionary for One Exported Sheet
#
# Iterates over the actual column names of `data` and looks up a
# human-readable explanation for each one from a supplied lookup table.
# Deliberately driven by `colnames(data)` rather than by the lookup table's
# own names, so that:
#   - a column present in the workbook but missing from `descriptions` is
#     caught immediately (via `stop()`) instead of silently shipping an
#     undocumented column to collaborators; and
#   - a stale `descriptions` entry left over for a column that no longer
#     exists in `data` (e.g. after a later edit to the notebook) is simply
#     ignored, rather than appearing in the dictionary for a column that is
#     not actually present in the sheet.
#
build_column_dictionary <- function(sheet_name, data, descriptions, workbook_path = NULL, rownames = FALSE) {

  if (!is.character(sheet_name) || length(sheet_name) != 1L ||
      is.na(sheet_name) || !nzchar(trimws(sheet_name))) {
    stop("build_column_dictionary(): `sheet_name` must be one non-empty, non-missing character value.")
  }
  if (!is.character(descriptions) || is.null(names(descriptions)) ||
      anyNA(names(descriptions)) || any(!nzchar(names(descriptions)))) {
    stop("build_column_dictionary(): `descriptions` must be a named character vector with valid column names.")
  }
  if (!is.null(workbook_path) &&
      (!is.character(workbook_path) || length(workbook_path) != 1L ||
       is.na(workbook_path) || !nzchar(trimws(workbook_path)))) {
    stop("build_column_dictionary(): `workbook_path` must be NULL or one non-empty, non-missing character path.")
  }
  if (!is.logical(rownames) || length(rownames) != 1L || is.na(rownames)) {
    stop("build_column_dictionary(): `rownames` must be TRUE or FALSE.")
  }

  if (is.null(dim(data)) || ncol(data) < 1L) {
    stop("build_column_dictionary(): `data` must contain at least one column.")
  }

  # Match add_sheet_to_excel()'s Excel-name normalization so a caller can pass
  # the same original name to both helpers without the dictionary lookup failing
  # after the worksheet writer sanitizes it.
  resolved_sheet_name <- sheet_name
  for (forbidden_character in c(":", "\\", "/", "?", "*", "[", "]")) {
    resolved_sheet_name <- gsub(forbidden_character, "_", resolved_sheet_name, fixed = TRUE)
  }
  if (nchar(resolved_sheet_name) > 31L) {
    resolved_sheet_name <- substr(resolved_sheet_name, 1L, 31L)
  }

  # ---------------------------------------------------------------------------
  # Step 1: Pull the actual column names from the exported data, not from the
  # `descriptions` lookup table. This is what guarantees the dictionary can
  # never silently drift out of sync with what is really in the workbook.
  # ---------------------------------------------------------------------------
  actual_column_names <- colnames(data)

  # ---------------------------------------------------------------------------
  # Step 2: Confirm every actual column has a matching entry in `descriptions`.
  # Fail loudly and specifically (naming the sheet and the missing column(s))
  # rather than writing a blank/NA explanation into the workbook.
  # ---------------------------------------------------------------------------
  undocumented_columns <- setdiff(actual_column_names, names(descriptions))
  if (length(undocumented_columns) > 0) {
    stop(
      "build_column_dictionary(): sheet '", sheet_name, "' has column(s) with ",
      "no entry in `descriptions`: ", paste(undocumented_columns, collapse = ", "),
      ". Add a description for each before exporting, so the Column_Dictionary ",
      "sheet stays complete and accurate."
    )
  }

  # ---------------------------------------------------------------------------
  # Step 3: Assemble the Sheet / Column / Explanation rows for this sheet, in
  # the same column order as the actual exported data.
  # ---------------------------------------------------------------------------
  column_dictionary_rows <- data.frame(
    Sheet       = resolved_sheet_name,
    Column      = actual_column_names,
    Explanation = unname(descriptions[actual_column_names]),
    stringsAsFactors = FALSE
  )

  # ---------------------------------------------------------------------------
  # Step 4 (optional): Attach each column's explanation directly to that
  # column's own header cell, as a hover/click Excel comment. Only runs when
  # `workbook_path` is supplied, and only after `sheet_name` has already been
  # written to that workbook via add_sheet_to_excel() -- this function
  # documents sheets, it does not create them.
  # ---------------------------------------------------------------------------
  if (!is.null(workbook_path)) {

    if (!file.exists(workbook_path)) {
      stop(
        "build_column_dictionary(): workbook_path '", workbook_path, "' does not exist yet. ",
        "Call add_sheet_to_excel() to write sheet '", sheet_name, "' before documenting it here."
      )
    }

    # getSheetNames() only accepts a path to an already-saved .xlsx file, not
    # an in-memory Workbook object, so existing sheets are checked this way
    # before the workbook is loaded into memory.
    workbook_sheet_names <- openxlsx::getSheetNames(workbook_path)
    matching_sheet_index <- match(tolower(resolved_sheet_name), tolower(workbook_sheet_names))
    if (is.na(matching_sheet_index)) {
      stop(
        "build_column_dictionary(): sheet '", resolved_sheet_name, "' was not found in '", workbook_path,
        "'. Call add_sheet_to_excel() to write this sheet before documenting it here."
      )
    }

    resolved_sheet_name <- workbook_sheet_names[matching_sheet_index]
    wb <- openxlsx::loadWorkbook(workbook_path)

    comment_box_width_units <- 6
    estimate_comment_box_height_units <- function(comment_text) {
      estimated_wrapped_lines <- ceiling(nchar(comment_text) / 60)
      estimated_height_units <- ceiling(estimated_wrapped_lines / 1.6) + 1
      max(3, min(20, estimated_height_units))
    }

    # Column indices must shift by one when an extra row-names column was
    # written alongside `data` (i.e. rownames = TRUE in the matching
    # add_sheet_to_excel() call for this sheet).
    column_index_offset <- if (isTRUE(rownames)) 1L else 0L

    for (column_index in seq_along(actual_column_names)) {
      column_name <- actual_column_names[column_index]
      column_description <- descriptions[[column_name]]
      openxlsx::writeComment(
        wb, resolved_sheet_name,
        col = column_index + column_index_offset,
        row = 1,
        comment = openxlsx::createComment(
          comment = column_description,
          author = "",
          visible = FALSE,
          width = comment_box_width_units,
          height = estimate_comment_box_height_units(column_description)
        )
      )
    }

    destination_directory <- dirname(workbook_path)
    temporary_path <- tempfile(".workbook-", tmpdir = destination_directory, fileext = ".xlsx")
    on.exit(unlink(temporary_path), add = TRUE)
    openxlsx::saveWorkbook(wb, temporary_path, overwrite = TRUE)
    if (!file.rename(temporary_path, workbook_path)) {
      stop("build_column_dictionary(): could not atomically replace workbook: ", workbook_path)
    }
  }

  column_dictionary_rows
}
