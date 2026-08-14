################################################################################
# Script: render_output_tree_function.R
# Purpose:
#   - Define a reusable function `render_output_tree()` that prints a
#     notebook's own output folder as a clickable directory tree -- every
#     file and subfolder currently on disk, each with a working link, plus an
#     optional inline description after each file -- similar in spirit to
#     `render_output_links()` (same project, R/functions/), but shows the
#     *whole* output folder's structure at a glance (including nested
#     subfolders) instead of one bullet per individually-named path.
#   - The folder is scanned live, at knit time, via fs::dir_info(). The tree
#     therefore always reflects exactly what this notebook actually produced
#     on this run -- nothing here is a hand-maintained snapshot that could
#     drift out of sync with reality.
#   - IMPORTANT -- why raw HTML, not plain Markdown, and why this still works
#     in both of this project's knit targets (html_document AND
#     github_document): showing a tree's box-drawing connectors (|--, `--)
#     aligned requires a monospace, whitespace-preserving block. A fenced
#     code block (```) preserves that alignment, but Markdown link syntax
#     [text](url) is NOT rendered inside a fenced code block -- on GitHub, and
#     in a knitted HTML report, it shows up as literal, non-clickable text.
#     A raw <pre>...</pre> block preserves the same monospace alignment while
#     keeping <a href="..."> links fully clickable, because <pre> is ordinary
#     HTML, not a code fence. This function returns its output via
#     knitr::asis_output(), the same mechanism render_output_links() already
#     uses -- knitr inserts that raw text directly into the document's
#     Markdown source (never wrapped in a code fence). For html_document,
#     that raw HTML becomes part of the rendered page directly. For
#     github_document, rmarkdown knits to an intermediate Markdown file and
#     then lets pandoc write the final .md; pandoc's Markdown/GFM writer
#     preserves raw HTML blocks verbatim, and GitHub renders raw HTML in any
#     .md file exactly as a browser would -- the same technique this
#     project's root README.md project-structure tree already relies on.
#   - Intended to be called as the final statement of a chunk with
#     results='asis' (or simply as a chunk's last auto-printed expression),
#     placed after every output-producing chunk earlier in the notebook has
#     already finished writing to `output_folder`.
################################################################################

# ==============================================================================
# Load necessary library
# ==============================================================================
# fs::dir_info() / fs::path_rel() / fs::path() build the recursive file
# listing, the tree's nesting structure, and the relative paths used for
# links. here::here() (namespace-qualified below, so no library(here) needed)
# locates the R/notebooks/ folder every rendered notebook lives in, the same
# anchor render_output_links() already builds its hrefs from.
# Calls are namespace-qualified so sourcing this helper does not attach fs.

# ==============================================================================
# Internal helper 1: insert one file into a nested-list tree structure
# ==============================================================================
# Not intended to be called directly -- used internally by render_output_tree()
# only. `components` is a file's path relative to the output folder, already
# split on "/" (e.g. c("forward", "Sample1.fastq.gz")). Each intermediate
# component becomes (or reuses) a branch node with a `.children` list; the
# final component becomes a leaf node carrying that file's already-rendered
# HTML link (+ optional comment, appended later during alignment -- see
# render_output_tree() Step 5) and its raw comment text for later alignment.
.insert_into_output_tree <- function(node, components, rendered_link, comment_text) {

  first_component <- components[1]

  if (length(components) == 1) {
    # Final component: this is the file itself.
    node[[first_component]] <- list(
      .is_file  = TRUE,
      .rendered = rendered_link,
      .comment  = comment_text
    )
    return(node)
  }

  # Intermediate component: ensure a branch node exists here (creating one if
  # this is the first file encountered under this subfolder, or overwriting a
  # same-named leaf in the unlikely case a file and folder share a name --
  # not expected in practice, but keeps this helper from erroring either way).
  if (is.null(node[[first_component]]) || isTRUE(node[[first_component]]$.is_file)) {
    node[[first_component]] <- list(.is_file = FALSE, .children = list())
  }

  node[[first_component]]$.children <- .insert_into_output_tree(
    node[[first_component]]$.children, components[-1], rendered_link, comment_text
  )

  node
}

# ==============================================================================
# Internal helper 2: recursively print the tree's lines
# ==============================================================================
# Not intended to be called directly -- used internally by render_output_tree()
# only. Walks the nested list built by .insert_into_output_tree() and returns
# two parallel character vectors (line text so far, and that line's comment,
# NA if none) -- kept separate from the final printed lines so the caller can
# right-align every comment into a single column afterward (Step 5 below).
# Uses the classic recursive connector/continuation-prefix algorithm: every
# entry except the last child at a given level is prefixed with "|-- " and
# continues (for its own children) with "|   "; the last child at a given
# level is prefixed with "`-- " and continues with "    " (no vertical bar),
# since nothing remains below it to connect to.
.render_output_tree_lines <- function(node, output_folder, notebooks_root, box_prefix = "", path_prefix = "") {

  entry_names   <- names(node)
  n_entries     <- length(entry_names)
  line_text     <- character(0)
  line_comment  <- character(0)
  line_path     <- character(0)

  for (entry_index in seq_along(entry_names)) {

    entry_name  <- entry_names[entry_index]
    entry_value <- node[[entry_name]]
    is_last     <- entry_index == n_entries

    connector    <- if (is_last) "└── " else "├── "
    continuation <- if (is_last) "    " else "│   "

    if (isTRUE(entry_value$.is_file)) {

      # Leaf node: a file. Its link was already built when it was inserted.
      file_rel_path <- if (nzchar(path_prefix)) paste0(path_prefix, "/", entry_name) else entry_name
      line_text    <- c(line_text, paste0(box_prefix, connector, entry_value$.rendered))
      line_comment <- c(line_comment, entry_value$.comment)
      line_path    <- c(line_path, file_rel_path)

    } else {

      # Branch node: a folder. Build its own clickable label here (rather
      # than at insertion time), since only now do we know its full path.
      folder_rel_path <- if (nzchar(path_prefix)) paste0(path_prefix, "/", entry_name) else entry_name
      folder_href     <- as.character(fs::path_rel(fs::path(output_folder, folder_rel_path), start = notebooks_root))
      folder_href     <- utils::URLencode(paste0(folder_href, "/"), reserved = FALSE)
      folder_label    <- paste0(
        '<a href="', htmltools::htmlEscape(folder_href, attribute = TRUE), '">',
        htmltools::htmlEscape(entry_name), '/</a>'
      )

      line_text    <- c(line_text, paste0(box_prefix, connector, folder_label))
      line_comment <- c(line_comment, NA_character_)
      line_path    <- c(line_path, folder_rel_path)

      child_lines <- .render_output_tree_lines(
        entry_value$.children, output_folder, notebooks_root,
        paste0(box_prefix, continuation), folder_rel_path
      )
      line_text    <- c(line_text, child_lines$text)
      line_comment <- c(line_comment, child_lines$comment)
      line_path    <- c(line_path, child_lines$relative_path)
    }
  }

  list(text = line_text, comment = line_comment, relative_path = line_path)
}

# ==============================================================================
# Define the function: render_output_tree
# ==============================================================================
#' Render a Notebook's Output Folder as a Clickable Directory Tree
#'
#' Scans `output_folder` recursively and prints every file and subfolder it
#' currently contains as a directory tree, with a clickable link on every
#' entry (relative to `R/notebooks/`, so links resolve correctly from this
#' notebook's own rendered location) and an optional inline description after
#' each file, right-aligned into a single column -- matching the style of
#' this project's root README.md project-structure tree. Because the folder
#' is scanned live at knit time, the tree always reflects exactly what this
#' notebook actually produced on this run.
#'
#' @param output_folder Character string, the absolute path to the folder to
#'   render (typically this notebook's own `output_folder` variable).
#' @param descriptions Optional named character vector mapping a file's path
#'   *relative to output_folder* (e.g. "data_integrity_check.xlsx", or
#'   "forward/Sample1.fastq.gz" for a nested file) to a short human-readable
#'   description, appended as an inline comment after that file's link. A
#'   subfolder can also be given a description by using its own relative path
#'   (e.g. "checkpoints") as a key. Files/folders present on disk but missing
#'   from this vector are listed with no comment. Defaults to an empty
#'   vector (no comments anywhere).
#' @param root_label Character string used as the tree's own root line.
#'   Defaults to `output_folder`'s path relative to the project root (e.g.
#'   "results/1_data_integrity_check/"), matching the display convention
#'   `render_output_links()` already uses for its link text.
#' @return A `knit_asis` object (see `knitr::asis_output()`), auto-printed as
#'   raw HTML when this function is the last expression of a chunk.
render_output_tree <- function(output_folder, descriptions = character(0), root_label = NULL) {

  if (!is.character(output_folder) || length(output_folder) != 1L ||
      is.na(output_folder) || !nzchar(trimws(output_folder))) {
    stop("render_output_tree(): `output_folder` must be one non-empty, non-missing character path.")
  }
  if (!is.character(descriptions) || anyNA(descriptions)) {
    stop("render_output_tree(): `descriptions` must be a character vector without missing values.")
  }
  if (length(descriptions) > 0L && (is.null(names(descriptions)) ||
      anyNA(names(descriptions)) || any(!nzchar(names(descriptions))))) {
    stop("render_output_tree(): every `descriptions` value must have a non-empty, non-missing path name.")
  }
  if (!is.null(root_label) && (!is.character(root_label) || length(root_label) != 1L ||
      is.na(root_label) || !nzchar(root_label))) {
    stop("render_output_tree(): `root_label` must be NULL or one non-empty, non-missing character value.")
  }

  # ---------------------------------------------------------------------------
  # Step 1: Validate the folder and handle the "nothing written yet" case
  # ---------------------------------------------------------------------------
  # Mirrors render_output_links()'s "Not found" handling: a missing or empty
  # folder is flagged rather than silently producing a blank, misleading tree.
  if (!fs::dir_exists(output_folder)) {
    return(knitr::asis_output(paste0(
      "\n**Not found:** `", output_folder,
      "` (expected but was not created -- check the chunks above for warnings or errors)\n"
    )))
  }

  all_files <- fs::dir_info(output_folder, recurse = TRUE, type = "file")

  if (nrow(all_files) == 0) {
    return(knitr::asis_output(paste0(
      "\n*", output_folder, " exists but is currently empty -- run the chunks above first.*\n"
    )))
  }

  # ---------------------------------------------------------------------------
  # Step 2: Compute relative paths for display, linking, and tree structure
  # ---------------------------------------------------------------------------
  # notebooks_root is the folder every rendered notebook (.html and .md) lives
  # in -- the same anchor render_output_links() builds its hrefs from, so a
  # link built here resolves correctly from inside the knitted document.
  notebooks_root <- here::here("R", "notebooks")
  project_root   <- here::here()

  if (is.null(root_label)) {
    root_label <- paste0(as.character(fs::path_rel(output_folder, start = project_root)), "/")
  }

  # Path of each file relative to output_folder itself (e.g.
  # "forward/Sample1.fastq.gz") -- both the tree's own internal grouping key
  # (Step 3) and the lookup key into `descriptions`. Sorting alphabetically
  # here (rather than relying on dir_info()'s own listing order) keeps every
  # folder's contents displayed in a predictable, alphabetical order.
  rel_to_output <- as.character(fs::path_rel(all_files$path, start = output_folder))
  sort_order    <- order(rel_to_output)
  all_files     <- all_files[sort_order, ]
  rel_to_output <- rel_to_output[sort_order]

  # Path of each file relative to R/notebooks/ -- used as the clickable
  # link's href, identical in spirit to render_output_links()'s href_path.
  href_paths <- as.character(fs::path_rel(all_files$path, start = notebooks_root))

  # ---------------------------------------------------------------------------
  # Step 3: Build a nested list representing the folder's structure
  # ---------------------------------------------------------------------------
  tree <- list()

  for (file_index in seq_along(rel_to_output)) {

    path_components <- strsplit(rel_to_output[file_index], "/", fixed = TRUE)[[1]]
    file_basename    <- path_components[length(path_components)]

    safe_href <- htmltools::htmlEscape(
      utils::URLencode(href_paths[file_index], reserved = FALSE),
      attribute = TRUE
    )
    safe_basename <- htmltools::htmlEscape(file_basename)
    rendered_link <- paste0('<a href="', safe_href, '">', safe_basename, '</a>')

    comment_text <- if (rel_to_output[file_index] %in% names(descriptions)) {
      descriptions[[rel_to_output[file_index]]]
    } else {
      NA_character_
    }

    tree <- .insert_into_output_tree(tree, path_components, rendered_link, comment_text)
  }

  # ---------------------------------------------------------------------------
  # Step 4: Recursively render every line, with comments kept separate
  # ---------------------------------------------------------------------------
  rendered <- .render_output_tree_lines(tree, output_folder, notebooks_root)

  # Folders can also carry a caller-supplied comment. Each rendered line keeps
  # its exact path relative to output_folder, so nested folders with the same
  # basename remain unambiguous.
  is_folder_line <- is.na(rendered$comment) & grepl('/</a>$', rendered$text)
  if (any(is_folder_line) && length(descriptions) > 0) {
    for (line_index in which(is_folder_line)) {
      folder_key <- rendered$relative_path[line_index]
      if (folder_key %in% names(descriptions)) {
        rendered$comment[line_index] <- descriptions[[folder_key]]
      }
    }
  }

  # ---------------------------------------------------------------------------
  # Step 5: Right-align every comment into a single column, then assemble
  # ---------------------------------------------------------------------------
  # Alignment is computed on each line's *visible* length -- i.e. with HTML
  # anchor tags collapsed down to just their link text -- since the <a
  # href="..."> markup itself takes up no visual width once rendered.
  visible_length <- function(html_line) {
    nchar(gsub('<a href="[^"]*">([^<]*)</a>', '\\1', html_line))
  }

  has_comment <- !is.na(rendered$comment)
  target_col  <- if (any(has_comment)) {
    max(vapply(rendered$text[has_comment], visible_length, integer(1))) + 2L
  } else {
    0L
  }

  final_lines <- vapply(seq_along(rendered$text), function(i) {
    if (is.na(rendered$comment[i])) {
      rendered$text[i]
    } else {
      pad_width <- max(1L, target_col - visible_length(rendered$text[i]))
      safe_comment <- htmltools::htmlEscape(as.character(rendered$comment[i]))
      paste0(rendered$text[i], strrep(" ", pad_width), "# ", safe_comment)
    }
  }, character(1))

  # ---------------------------------------------------------------------------
  # Step 6: Wrap in a raw HTML <pre> block and return for auto-printing
  # ---------------------------------------------------------------------------
  tree_html <- paste0(
    "\n<pre>\n",
    htmltools::htmlEscape(root_label), "\n",
    paste(final_lines, collapse = "\n"),
    "\n</pre>\n"
  )

  knitr::asis_output(tree_html)
}

# Keep the recursive implementation helpers private to render_output_tree() so
# sourcing this script adds only its documented public function to the caller.
render_output_tree <- local({
  insert_into_output_tree <- .insert_into_output_tree
  render_output_tree_lines <- .render_output_tree_lines
  implementation <- render_output_tree
  .insert_into_output_tree <- insert_into_output_tree
  .render_output_tree_lines <- render_output_tree_lines
  environment(.insert_into_output_tree) <- environment()
  environment(.render_output_tree_lines) <- environment()
  environment(implementation) <- environment()

  function(output_folder, descriptions = character(0), root_label = NULL) {
    implementation(output_folder, descriptions, root_label)
  }
})
rm(.insert_into_output_tree, .render_output_tree_lines)

# ==============================================================================
# Example usage (uncomment to test)
# ==============================================================================
# Called as the last line of a chunk, after every output-writing chunk above
# it has already run, so the folder reflects this run's actual outputs:
#
# render_output_tree(
#   output_folder,
#   descriptions = c(
#     "data_integrity_check.xlsx" = "FASTQ pairing/read-count/validity report (Data_Integrity_Check + Column_Dictionary sheets)"
#   )
# )
