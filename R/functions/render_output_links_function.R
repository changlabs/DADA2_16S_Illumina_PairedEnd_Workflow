################################################################################
# Script: render_output_links_function.R
# Purpose:
#   - Define a reusable function `render_output_links()` that prints one or
#     more file/folder paths as clickable Markdown links in a knitted R
#     Markdown report (HTML or GitHub-flavoured Markdown), so anybody running
#     this workflow can click straight through to a generated output file
#     instead of navigating the project's folder structure by hand.
#   - Intended to be called as the final statement of any chunk that writes
#     an output file (or files) -- an Excel workbook, a plot, a checkpoint
#     .RData file, an exported FASTA/Newick/CSV/TSV, an HTML widget, or an
#     entire output folder -- immediately documenting where that chunk's
#     result landed, for reproducibility across machines.
#   - Flags (rather than silently linking) any path that does not actually
#     exist on disk, since a clickable link to a file that was not created
#     (e.g. because an optional tool was unavailable, or a step failed) would
#     otherwise be misleading.
#   - Displays a short, project-relative path (starting at results/) as the
#     link text, and builds the link TARGET itself as a path relative to
#     R/notebooks/ (where every notebook's rendered .html always lives) --
#     so long absolute paths do not overflow the page margin in a knitted
#     HTML report, the link remains fully clickable, AND the link keeps
#     working if the whole project folder is copied, renamed, or shared with
#     someone else on a different machine (a machine-specific absolute path
#     as the href would silently break the moment the project folder moves).
################################################################################

# ==============================================================================
# Load necessary library
# ==============================================================================
# fs::file_exists() / fs::dir_exists() / fs::path_rel() are used to check for
# existing outputs and to shorten the displayed link text (see Step 2 below).
# knitr::asis_output() (called namespace-qualified below, so no
# library(knitr) is needed here) is what lets this function's return value be
# inserted into the rendered document as raw Markdown, without requiring the
# calling chunk to set the `results='asis'` chunk option. here::here() (also
# namespace-qualified, no library(here) needed) is used to find the project
# root that every displayed path is shortened relative to; every notebook in
# this project already depends on the `here` package for its own path
# construction, so it is always available by the time this function runs.
# Calls are namespace-qualified so sourcing this helper does not attach fs.

# ==============================================================================
# Define the function: render_output_links
# ==============================================================================
# Function arguments:
#   - paths  : Character vector (or fs_path vector) of one or more absolute
#              file or folder paths to link to. Typically the same path
#              variable(s) just used a few lines above to write the output
#              (e.g. `output_excel_path`, `tree_output_path`, a checkpoint
#              path, or an entire `output_folder`).
#   - labels : Optional character vector, the same length as `paths`, giving
#              a short human-readable description for each path (e.g.
#              "Alignment (FASTA)", "Phylogenetic tree (Newick)"). Defaults
#              to NULL, in which case each path is shown on its own with no
#              prefix. Use this whenever a single chunk writes more than one
#              distinct output, so each link is self-explanatory.
#
# Return value:
#   A `knit_asis` object (via `knitr::asis_output()`). Because knitr auto-
#   prints the final visible expression of a chunk, simply calling this
#   function as the last line of a chunk is enough to insert its Markdown
#   into the rendered report -- no `results='asis'` chunk option is needed.
#
# IMPORTANT -- calling this inside a `for` loop or `if`/`else` block:
#   A bare `for` loop always returns `invisible(NULL)`, regardless of its
#   body, so a call to this function as the last line inside a `for` loop
#   will NOT auto-print; wrap it in `print()` instead, e.g.
#   `print(render_output_links(path_for_this_iteration))`. A call inside an
#   `if () {...} else {...}` block does not need this, since the whole
#   `if`/`else` statement's value already auto-prints normally when it is
#   itself the chunk's last top-level expression (this is standard R
#   top-level auto-print behaviour, not specific to this function).
#
#' Render One or More File/Folder Paths as Clickable Markdown Links
#'
#' Prints a small Markdown bullet list -- one bullet per path -- where each
#' existing path is a clickable link. Both the link TEXT and the link TARGET
#' are portable, relative paths rather than machine-specific absolute paths:
#' the TEXT is shortened to a project-relative path starting at `results/`
#' (since every output this function links to lives somewhere under the
#' project's `results/` folder), and the TARGET (href) is expressed relative
#' to `R/notebooks/` -- the folder every notebook's rendered .html file
#' always lives in -- so the link actually resolves correctly from inside
#' that .html file. Using a relative href (instead of the full absolute
#' path) also means the link keeps working if the entire project folder is
#' copied, renamed, or shared with someone else on a different machine,
#' where the original absolute path would no longer exist. This keeps the
#' displayed text short enough to stay within a knitted HTML report's page
#' margin -- a full absolute path (e.g. `/Users/name/Downloads/workflows/...
#' /results/5_dada2_pipeline/file.pdf`) can otherwise overflow the page
#' width, especially on deeply nested project locations. A path that does
#' not exist on disk is flagged with a bolded "Not found" marker instead of
#' being linked, so a failed or skipped step is obvious rather than
#' producing a silently broken link.
#'
#' @param paths Character vector of one or more absolute file or folder
#'   paths.
#' @param labels Optional character vector, the same length as `paths`, of a
#'   short description to prefix each link with. Defaults to NULL, printing
#'   only the path itself.
#' @return A `knit_asis` object (see `knitr::asis_output()`), auto-printed as
#'   raw Markdown when this function is the last expression of a chunk (or
#'   explicitly `print()`-ed inside a `for` loop -- see above).
render_output_links <- function(paths, labels = NULL) {

  # ---------------------------------------------------------------------------
  # Step 1: Validate arguments
  # ---------------------------------------------------------------------------
  # as.character() strips any fs_path class before the string operations
  # below, since fs_path's own print/format methods are meant for console
  # display, not for embedding in a hand-built Markdown string.
  paths <- as.character(paths)

  if (length(paths) == 0L || anyNA(paths) || any(!nzchar(trimws(paths)))) {
    stop("render_output_links(): `paths` must contain one or more non-empty, non-missing paths.")
  }

  if (!is.null(labels) && length(labels) != length(paths)) {
    stop(
      "render_output_links(): `labels` (length ", length(labels), ") must be ",
      "the same length as `paths` (length ", length(paths), "), or NULL."
    )
  }
  if (!is.null(labels) && (!is.character(labels) || anyNA(labels) || any(!nzchar(trimws(labels))))) {
    stop("render_output_links(): `labels` must contain only non-empty, non-missing character values.")
  }

  # ---------------------------------------------------------------------------
  # Step 2: Build one Markdown bullet per path
  # ---------------------------------------------------------------------------
  # A path is treated as "found" if it exists as either a file or a
  # directory, since this function is used for both individual output files
  # and entire output folders.
  #
  # here::here() resolves to the project root (the folder containing the
  # .Rproj file, or wherever here::here() is otherwise anchored) -- the same
  # root every notebook already builds its output paths from (e.g.
  # here("results", "5_dada2_pipeline", ...) or here(output_folder, ...)).
  # fs::path_rel() re-expresses each absolute path relative to that root, so
  # a path like "/Users/name/Downloads/workflows/.../results/5_dada2_pipeline
  # /file.pdf" displays as "results/5_dada2_pipeline/file.pdf" -- short
  # enough to stay within the page margin of a knitted HTML report, while
  # still starting from (and clearly showing) the results/ folder the file
  # lives under. If a path happens to fall outside the project root (not
  # expected in normal use, but not fatal either), path_rel() simply returns
  # a "../"-prefixed relative path instead of failing.
  project_root <- here::here()

  # The link HREF, unlike the link TEXT above, must be expressed relative to
  # where the rendered .html file itself will live -- R/notebooks/, the same
  # folder every .Rmd notebook in this project is knit from -- not relative
  # to the project root. A browser resolves a relative href against the
  # current document's own location, so an href built relative to the
  # project root would be wrong by exactly the "R/notebooks/" prefix. Using
  # a relative href here (instead of the previous absolute path) is what
  # makes every link this function prints continue to work if the whole
  # project folder is copied, renamed, or shared with someone else.
  notebooks_root <- here::here("R", "notebooks")

  escape_markdown_text <- function(value) {
    value <- gsub("\\", "\\\\", as.character(value), fixed = TRUE)
    for (special in c("[", "]", "*", "_", "`")) {
      value <- gsub(special, paste0("\\", special), value, fixed = TRUE)
    }
    value
  }

  encode_markdown_href <- function(value) {
    encoded <- utils::URLencode(as.character(value), reserved = FALSE)
    encoded <- gsub("(", "%28", encoded, fixed = TRUE)
    encoded <- gsub(")", "%29", encoded, fixed = TRUE)
    encoded <- gsub("<", "%3C", encoded, fixed = TRUE)
    gsub(">", "%3E", encoded, fixed = TRUE)
  }

  markdown_bullets <- vapply(seq_along(paths), function(path_index) {

    current_path   <- paths[path_index]
    path_exists    <- fs::file_exists(current_path) || fs::dir_exists(current_path)
    display_path   <- as.character(fs::path_rel(current_path, start = project_root))
    href_path      <- as.character(fs::path_rel(current_path, start = notebooks_root))
    display_label  <- if (is.null(labels)) NULL else labels[path_index]
    link_text      <- if (is.null(display_label)) display_path else paste0(display_label, ": ", display_path)
    safe_link_text <- escape_markdown_text(link_text)
    safe_href      <- encode_markdown_href(href_path)

    if (path_exists) {
      paste0("- [", safe_link_text, "](<", safe_href, ">)")
    } else {
      paste0("- **Not found:** ", safe_link_text,
             " (expected but was not created -- check the chunk above for warnings or errors)")
    }
  }, character(1))

  # ---------------------------------------------------------------------------
  # Step 3: Return as raw Markdown, ready to auto-print
  # ---------------------------------------------------------------------------
  # A leading and trailing blank line ensure the bullet list is separated
  # from surrounding prose/output when inserted into the rendered document.
  knitr::asis_output(paste(c("", markdown_bullets, ""), collapse = "\n"))
}

# ==============================================================================
# Example usage (uncomment to test)
# ==============================================================================
# Single output file, no label -- call as the last line of the chunk that
# wrote output_excel_path:
# render_output_links(output_excel_path)
#
# Multiple distinct outputs from one chunk, each labelled:
# render_output_links(
#   c(alignment_output_path, tree_output_path),
#   labels = c("Multiple sequence alignment (FASTA)", "Phylogenetic tree (Newick)")
# )
#
# Inside a for loop -- must be wrapped in print():
# for (db_name in names(taxonomy_paths)) {
#   print(render_output_links(taxonomy_paths[[db_name]], labels = paste(db_name, "taxonomy table")))
# }
