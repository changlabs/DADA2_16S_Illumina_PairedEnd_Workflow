################################################################################
# Script: download_reference_databases.R
# Purpose:
#   - Download the DADA2-formatted 16S rRNA gene reference databases for both
#     bacteria and archaea from two authoritative sources: the Genome Taxonomy
#     Database (GTDB release r220) and SILVA (version 138.2).
#   - Organise the downloaded files into a consistent directory hierarchy under
#     tools/trainsets/, with one subfolder per database (GTDB/ and SILVA/).
#   - Verify each downloaded file is non-empty and, where possible, confirm the
#     expected file size against the known Zenodo record sizes.
#   - Write a human-readable plain-text manifest (download_manifest.txt) inside
#     each subfolder recording the source URLs, download timestamp, file sizes,
#     and the exact R path strings that must be used in downstream DADA2
#     notebooks and scripts.
#   - Print a final copy-paste-ready R code block to the console so the user
#     can immediately configure their analysis environment.
#
# Background — why two databases?
#   DADA2 taxonomy assignment works in two passes that require separate files:
#
#   1. assignTaxonomy() — Uses a training FASTA where each sequence header
#      encodes the full taxonomic lineage (kingdom through genus). The function
#      trains a naive Bayes classifier on this file and classifies ASVs down to
#      genus level.
#      Files used: *_genus.fa.gz (GTDB)  |  *_toGenus_trainset.fa.gz (SILVA)
#
#   2. addSpecies() — Uses a species-assignment FASTA containing exact-match
#      reference sequences with genus+species labels. The function attempts to
#      assign species to ASVs that share 100% identity with a reference sequence.
#      Files used: *_species.fa.gz (GTDB)  |  *_assignSpecies.fa.gz (SILVA)
#
# Database versions downloaded by this script:
#   GTDB  r220   — Genome Taxonomy Database release 220 formatted for DADA2.
#                  Contains 46,891 bacterial and 2,812 archaeal 16S rRNA gene
#                  sequences derived from GTDB-curated genome trees.
#                  Zenodo record: https://zenodo.org/records/13984843
#
#   SILVA v138.2 — SILVA small-subunit (SSU) rRNA database version 138.2
#                  formatted for DADA2. Covers all three domains of life
#                  (Bacteria, Archaea, Eukarya) with comprehensive curation.
#                  Zenodo record: https://zenodo.org/records/14169026
#
# Note on database choice for your analysis:
#   - Use GTDB if your study focuses on bacterial/archaeal community composition
#     and you prefer the GTDB taxonomy, which is based on genome phylogeny and
#     is more consistent at higher ranks than SILVA's taxonomy.
#   - Use SILVA if you need eukaryotic sequences (e.g., from 18S amplicons or
#     mixed 16S/18S libraries), or if you need SILVA-compatible taxonomy strings
#     for downstream tools (e.g., PICRUSt2, SILVA-based phylogenetic analyses).
#   - It is common practice to run both and compare results for a subset of
#     samples before committing to one database for the full dataset.
#
# Directory layout created by this script:
#
#   <project_root>/
#   └── tools/
#       └── trainsets/
#           ├── GTDB/
#           │   ├── GTDB_bac120_arc53_ssu_r220_genus.fa.gz
#           │   ├── GTDB_bac120_arc53_ssu_r220_species.fa.gz
#           │   └── download_manifest.txt
#           └── SILVA/
#               ├── silva_nr99_v138.2_toGenus_trainset.fa.gz
#               ├── silva_v138.2_assignSpecies.fa.gz
#               └── download_manifest.txt
#
################################################################################


# ==============================================================================
# User-adjustable defaults
# ==============================================================================
# Modify these values to change the behaviour when no command-line arguments
# are provided. All settings can also be overridden via command-line flags.

# default_project_dir    : Absolute path to the root of your DADA2 project.
#                          The trainsets/ folder will be created at:
#                          <default_project_dir>/tools/trainsets/
default_project_dir <- getwd()

# default_trainsets_name : Name of the top-level reference database folder
#                          created inside tools/. "trainsets" is the standard
#                          name used throughout this workflow.
default_trainsets_name <- "trainsets"

# default_force_download : Set to TRUE to re-download files even if they
#                          already exist. FALSE (default) skips files that are
#                          already present and non-empty to save bandwidth.
default_force_download <- FALSE

# default_timeout_seconds: Maximum number of seconds to wait for each file
#                          download before aborting. The reference FASTA files
#                          are large (up to ~280 MB), so this is set generously
#                          to 3600 seconds (1 hour). Reduce if you are on a
#                          fast institutional network.
default_timeout_seconds <- 3600

# default_download_method: R download method passed to download.file().
#                          "auto" lets R pick the best available method
#                          (libcurl > curl > wget). Override to "curl" or
#                          "wget" if "auto" fails on your system.
default_download_method <- "auto"


# ==============================================================================
# Database definitions
# ==============================================================================
# These lists define every file this script will download. Each entry is a
# named list with:
#   url         : Direct download URL (Zenodo download=1 link)
#   filename    : Local filename to save as inside the database subfolder
#   expected_size_bytes : Exact byte size for the fixed Zenodo release
#   description : Human-readable description of the file's role in DADA2
#   dada2_role  : Either "assignTaxonomy" or "addSpecies" — the DADA2 function
#                 that uses this file
#
# To add a new file in the future, append another list entry here without
# touching the download logic below.

gtdb_files <- list(
  list(
    url         = "https://zenodo.org/records/13984843/files/GTDB_bac120_arc53_ssu_r220_genus.fa.gz?download=1",
    filename    = "GTDB_bac120_arc53_ssu_r220_genus.fa.gz",
    expected_size_bytes = 17453606,
    description = "GTDB r220 — genus-level training set for assignTaxonomy(). Contains full taxonomic lineages from domain to genus for 46,891 bacterial and 2,812 archaeal 16S rRNA sequences.",
    dada2_role  = "assignTaxonomy"
  ),
  list(
    url         = "https://zenodo.org/records/13984843/files/GTDB_bac120_arc53_ssu_r220_species.fa.gz?download=1",
    filename    = "GTDB_bac120_arc53_ssu_r220_species.fa.gz",
    expected_size_bytes = 17443214,
    description = "GTDB r220 — species-level assignment set for addSpecies(). Contains exact-match reference sequences with genus+species labels for 100%-identity species assignment.",
    dada2_role  = "addSpecies"
  )
)

silva_files <- list(
  list(
    url         = "https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1",
    filename    = "silva_nr99_v138.2_toGenus_trainset.fa.gz",
    expected_size_bytes = 139996892,
    description = "SILVA v138.2 — genus-level training set for assignTaxonomy(). Non-redundant (nr99) sequences from the SILVA SSU database covering Bacteria, Archaea, and Eukarya, annotated down to genus level.",
    dada2_role  = "assignTaxonomy"
  ),
  list(
    url         = "https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1",
    filename    = "silva_v138.2_assignSpecies.fa.gz",
    expected_size_bytes = 69899323,
    description = "SILVA v138.2 — species-level assignment set for addSpecies(). Contains exact-match reference sequences annotated with genus+species for 100%-identity species assignment.",
    dada2_role  = "addSpecies"
  )
)


# ==============================================================================
# Helper functions
# ==============================================================================


# ------------------------------------------------------------------------------
# message_line()
# ------------------------------------------------------------------------------
# Purpose : Print one line to stdout followed by a newline. Centralises all
#           console output so formatting is consistent throughout the script.
# Args    : ... — character scalars to concatenate and print.
# Returns : NULL invisibly.
message_line <- function(...) {
  cat(..., "\n", sep = "")
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# section_header()
# ------------------------------------------------------------------------------
# Purpose : Print a prominent section separator to visually divide console
#           output into labelled blocks, each corresponding to a major step.
# Args    : title — character scalar section label.
# Returns : NULL invisibly.
section_header <- function(title) {
  cat(
    "\n",
    paste(rep("=", 78), collapse = ""), "\n",
    title, "\n",
    paste(rep("=", 78), collapse = ""), "\n",
    sep = ""
  )
  invisible(NULL)
}


# ------------------------------------------------------------------------------
# normalize_path_safe()
# ------------------------------------------------------------------------------
# Purpose : Return a normalized absolute path without requiring the path to
#           already exist. The built-in normalizePath() emits a warning for
#           paths that do not yet exist unless mustWork = FALSE is set.
# Args    : path — character scalar path to normalize.
# Returns : Character scalar normalized path.
normalize_path_safe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

assert_project_root <- function(project_dir) {
  project_files <- list.files(project_dir, pattern = "[.]Rproj$", full.names = TRUE)
  required_paths <- file.path(project_dir, c("R", file.path("R", "notebooks"), "data"))
  if (length(project_files) == 0L || !all(dir.exists(required_paths))) {
    stop(
      "The selected project directory is not a DADA2 workflow root: ", project_dir, "\n",
      "Expected an .Rproj file plus R/, R/notebooks/, and data/. ",
      "Open the project in RStudio or pass --project-dir explicitly.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

project_relative_path <- function(path, project_dir) {
  normalized_path <- normalize_path_safe(path)
  normalized_root <- sub("/+$", "", normalize_path_safe(project_dir))
  prefix <- paste0(normalized_root, "/")
  if (!startsWith(normalized_path, prefix)) return(NA_character_)
  substring(normalized_path, nchar(prefix) + 1L)
}


# ------------------------------------------------------------------------------
# format_file_size()
# ------------------------------------------------------------------------------
# Purpose : Convert a raw byte count into a human-readable string (B, KB, MB,
#           or GB). Used when reporting file sizes in manifest files and console
#           output so the user can quickly judge whether a download looks right.
# Args    : bytes — numeric scalar, file size in bytes.
# Returns : Character scalar formatted size string (e.g., "278.4 MB").
format_file_size <- function(bytes) {
  if (is.na(bytes) || bytes < 0) return("unknown")
  if (bytes < 1024)              return(paste0(bytes, " B"))
  if (bytes < 1024^2)            return(paste0(round(bytes / 1024,    1), " KB"))
  if (bytes < 1024^3)            return(paste0(round(bytes / 1024^2,  1), " MB"))
  paste0(round(bytes / 1024^3, 2), " GB")
}


# ------------------------------------------------------------------------------
# parse_arguments()
# ------------------------------------------------------------------------------
# Purpose : Parse Rscript command-line arguments. All supported flags are
#           described in the Usage section at the top of the file. Unknown
#           flags produce a hard error to avoid silent misconfiguration.
# Returns : Named list with parsed argument values (or defaults).
parse_arguments <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)

  # Initialise all fields from the user-adjustable defaults.
  parsed <- list(
    project_dir     = default_project_dir,
    trainsets_name  = default_trainsets_name,
    force_download  = default_force_download,
    timeout_seconds = default_timeout_seconds,
    download_method = default_download_method
  )

  if (length(raw_args) == 0L) return(parsed)

  i <- 1L
  read_option_value <- function(option_name) {
    if (i >= length(raw_args)) {
      stop("Missing value after ", option_name, ".", call. = FALSE)
    }
    candidate <- raw_args[[i + 1L]]
    if (!nzchar(candidate) || startsWith(candidate, "--") || candidate %in% c("-h", "--help")) {
      stop("Missing value after ", option_name, ".", call. = FALSE)
    }
    i <<- i + 1L
    candidate
  }

  while (i <= length(raw_args)) {
    current_arg <- raw_args[[i]]

    if (current_arg == "--project-dir") {
      parsed$project_dir <- read_option_value(current_arg)

    } else if (current_arg == "--trainsets-name") {
      parsed$trainsets_name <- read_option_value(current_arg)

    } else if (current_arg == "--force") {
      # Re-download all files even if they already exist and are non-empty.
      parsed$force_download <- TRUE

    } else if (current_arg == "--timeout") {
      parsed$timeout_seconds <- suppressWarnings(as.integer(read_option_value(current_arg)))

    } else if (current_arg == "--method") {
      parsed$download_method <- read_option_value(current_arg)

    } else if (current_arg %in% c("-h", "--help")) {
      message_line(
        "Usage: Rscript download_reference_databases.R [options]\n",
        "\n",
        "Options:\n",
        "  --project-dir <path>     DADA2 project root (tools/trainsets/ created here)\n",
        "  --trainsets-name <name>  Override the trainsets folder name (default: trainsets)\n",
        "  --force                  Re-download files even if they already exist\n",
        "  --timeout <seconds>      Per-file download timeout (default: 3600)\n",
        "  --method <method>        download.file() method: auto, curl, wget, libcurl\n",
        "  -h, --help               Show this help message\n"
      )
      quit(save = "no", status = 0)

    } else {
      stop("Unknown argument: ", current_arg, call. = FALSE)
    }

    i <- i + 1L
  }

  if (is.na(parsed$timeout_seconds) || parsed$timeout_seconds <= 0L) {
    stop("--timeout must be a positive whole number of seconds.", call. = FALSE)
  }
  allowed_methods <- c("auto", "curl", "wget", "libcurl")
  if (!(parsed$download_method %in% allowed_methods)) {
    stop(
      "--method must be one of: ", paste(allowed_methods, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (parsed$trainsets_name %in% c(".", "..") ||
      basename(parsed$trainsets_name) != parsed$trainsets_name) {
    stop("--trainsets-name must be a single folder name, not a path.", call. = FALSE)
  }

  parsed
}


# ------------------------------------------------------------------------------
# verify_database_file()
# ------------------------------------------------------------------------------
# Validate the exact release size and stream through the complete gzip payload.
# Reading to EOF makes zlib check the compressed stream rather than accepting a
# non-empty partial download.
verify_database_file <- function(path, expected_size_bytes = NULL) {
  if (!file.exists(path)) {
    stop("Reference database file is missing: ", path, call. = FALSE)
  }

  actual_size <- unname(file.info(path)$size)
  if (!is.finite(actual_size) || actual_size <= 0) {
    stop("Reference database file is empty: ", path, call. = FALSE)
  }
  if (!is.null(expected_size_bytes) &&
      (!is.finite(expected_size_bytes) || actual_size != expected_size_bytes)) {
    stop(
      "Reference database size mismatch for ", basename(path), ": expected ",
      expected_size_bytes, " bytes, found ", actual_size, " bytes.",
      call. = FALSE
    )
  }

  if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    raw_connection <- file(path, open = "rb")
    gzip_magic <- readBin(raw_connection, what = "raw", n = 2L)
    close(raw_connection)
    if (!identical(gzip_magic, as.raw(c(0x1f, 0x8b)))) {
      stop(
        "Reference database does not have a gzip header: ", basename(path),
        call. = FALSE
      )
    }

    connection <- gzfile(path, open = "rb")
    on.exit(close(connection), add = TRUE)
    gzip_error <- tryCatch({
      withCallingHandlers({
        repeat {
          chunk <- readBin(connection, what = "raw", n = 1024L * 1024L)
          if (length(chunk) == 0L) break
        }
      }, warning = function(w) stop(conditionMessage(w), call. = FALSE))
      NULL
    }, error = function(e) conditionMessage(e))
    if (!is.null(gzip_error)) {
      stop(
        "Reference database is not a valid complete gzip stream: ",
        basename(path), " (", gzip_error, ")",
        call. = FALSE
      )
    }
  }

  invisible(actual_size)
}


# ------------------------------------------------------------------------------
# download_database_file()
# ------------------------------------------------------------------------------
# Purpose : Download a single reference database file from a URL to a local
#           path with robust error handling:
#             1. Verify and skip an existing destination unless
#                force_download = TRUE.
#             2. Attempt the download with a configurable timeout using
#                download.file(), which supports libcurl, curl, and wget.
#             3. Download to a temporary file in the destination folder.
#             4. Verify the exact release size and complete gzip stream.
#             5. Atomically replace the destination while retaining a backup
#                until the replacement succeeds.
#
# Args    : url             — character scalar, the full download URL.
#           dest_path       — character scalar, full path for the saved file.
#           description     — character scalar, human-readable description.
#           timeout_seconds — integer, download timeout.
#           download_method — character scalar, method for download.file().
#           force_download  — logical, if TRUE overwrite existing files.
#           expected_size_bytes — exact byte size recorded for this release.
# Returns : Named list: (url, dest_path, filename, file_size_bytes,
#           file_size_human, download_timestamp, verification, skipped).
download_database_file <- function(url,
                                   dest_path,
                                   description,
                                   timeout_seconds = 3600L,
                                   download_method = "auto",
                                   force_download  = FALSE,
                                   expected_size_bytes = NULL) {

  filename <- basename(dest_path)

  message_line("\n  File      : ", filename)
  message_line("  URL       : ", url)
  message_line("  Dest      : ", dest_path)

  # ---------------------------------------------------------------------------
  # Verify and skip an existing download
  # ---------------------------------------------------------------------------
  # Skipping avoids re-downloading large files on re-runs, which is important
  # for reference databases that are hundreds of megabytes. The user can
  # override this with --force if they suspect a previous download was corrupt.
  if (file.exists(dest_path) && !isTRUE(force_download)) {
    existing_size <- tryCatch(
      verify_database_file(dest_path, expected_size_bytes),
      error = function(e) {
        stop(
          conditionMessage(e),
          " Rerun with --force to replace the invalid file.",
          call. = FALSE
        )
      }
    )
    message_line("  Status    : SKIPPED (already exists, ",
                 format_file_size(existing_size), "; verified)")
    message_line("              Use --force to re-download.")

    return(list(
      url              = url,
      dest_path        = dest_path,
      filename         = filename,
      file_size_bytes  = existing_size,
      file_size_human  = format_file_size(existing_size),
      download_timestamp = format(file.info(dest_path)$mtime,
                                  tz = Sys.timezone(), usetz = TRUE),
      verification      = "exact size and gzip stream verified",
      skipped          = TRUE
    ))
  }

  # ---------------------------------------------------------------------------
  # Set the download timeout
  # ---------------------------------------------------------------------------
  # The global "timeout" option controls how long download.file() waits for a
  # response. We save the original value and restore it with on.exit() so that
  # this script does not permanently alter the user's R session options.
  original_timeout <- getOption("timeout")
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout = timeout_seconds)

  # ---------------------------------------------------------------------------
  # Perform the download
  # ---------------------------------------------------------------------------
  message_line("  Status    : Downloading...")
  download_start <- Sys.time()
  temp_extension <- if (grepl("\\.gz$", filename, ignore.case = TRUE)) ".gz" else ""
  temp_stem <- if (nzchar(temp_extension)) sub("\\.gz$", "", filename, ignore.case = TRUE) else filename
  download_temp <- tempfile(
    pattern = paste0(".", temp_stem, ".download-"),
    tmpdir = dirname(dest_path),
    fileext = temp_extension
  )
  on.exit(if (file.exists(download_temp)) unlink(download_temp, force = TRUE), add = TRUE)

  tryCatch({
    download.file(
      url      = url,
      destfile = download_temp,
      method   = download_method,
      mode     = "wb",   # Write binary — essential for .gz files
      quiet    = FALSE   # Show transfer progress in the console
    )
  }, error = function(e) {
    if (file.exists(download_temp)) unlink(download_temp, force = TRUE)
    stop("Download failed for: ", filename, "\n  Error: ", conditionMessage(e),
         call. = FALSE)
  })

  download_end <- Sys.time()
  elapsed_seconds <- round(as.numeric(difftime(download_end, download_start, units = "secs")), 1)

  # ---------------------------------------------------------------------------
  # Verify the downloaded file
  # ---------------------------------------------------------------------------
  # A successful download.file() call does not guarantee the file is complete;
  # the exact size and full gzip stream are checked before replacement.
  if (!file.exists(download_temp)) {
    stop("Download appeared to succeed but the temporary file is missing: ", download_temp,
         call. = FALSE)
  }

  final_size <- tryCatch(
    verify_database_file(download_temp, expected_size_bytes),
    error = function(e) {
      unlink(download_temp, force = TRUE)
      stop("Downloaded file failed verification: ", conditionMessage(e), call. = FALSE)
    }
  )

  # Preserve an existing destination until the verified replacement has been
  # installed successfully. Both temporary paths live in the destination
  # folder, keeping rename operations on the same filesystem.
  backup_path <- NULL
  if (file.exists(dest_path)) {
    backup_path <- tempfile(
      pattern = paste0(".", filename, ".backup-"),
      tmpdir = dirname(dest_path)
    )
    if (!file.rename(dest_path, backup_path)) {
      stop("Could not preserve the existing file before replacement: ", dest_path, call. = FALSE)
    }
  }
  if (!file.rename(download_temp, dest_path)) {
    restore_message <- ""
    if (!is.null(backup_path) && file.exists(backup_path)) {
      if (!file.rename(backup_path, dest_path)) {
        restore_message <- paste0(" The previous file remains at: ", backup_path)
      }
    }
    stop("Could not install the verified download at: ", dest_path, restore_message, call. = FALSE)
  }
  if (!is.null(backup_path) && file.exists(backup_path)) unlink(backup_path, force = TRUE)

  message_line("  Downloaded: ", format_file_size(final_size),
               " in ", elapsed_seconds, "s")

  list(
    url                = url,
    dest_path          = dest_path,
    filename           = filename,
    file_size_bytes    = final_size,
    file_size_human    = format_file_size(final_size),
    download_timestamp = format(Sys.time(), tz = Sys.timezone(), usetz = TRUE),
    verification       = "exact size and gzip stream verified",
    skipped            = FALSE
  )
}


# ------------------------------------------------------------------------------
# write_database_manifest()
# ------------------------------------------------------------------------------
# Purpose : Write a plain-text manifest file inside a database subfolder
#           recording everything a future user needs to understand, verify, and
#           reproduce the installation:
#             - Which database and version was installed
#             - The source Zenodo record URL
#             - Each downloaded file: URL, local path, file size, timestamp
#             - The exact R code to paste into DADA2 notebooks/scripts
#             - Guidance on which DADA2 function each file feeds
#
# Args    : manifest_path  — character scalar, full path for the manifest file.
#           db_name        — character scalar, database label ("GTDB" or "SILVA").
#           db_version     — character scalar, version label (e.g., "r220").
#           zenodo_record  — character scalar, Zenodo record URL.
#           file_results   — list of named lists returned by download_database_file().
#           r_var_genus    — character scalar, R variable name for the genus file path.
#           r_var_species  — character scalar, R variable name for the species file path.
# Returns : TRUE invisibly.
write_database_manifest <- function(manifest_path,
                                    db_name,
                                    db_version,
                                    zenodo_record,
                                    file_results,
                                    r_var_genus,
                                    r_var_species,
                                    project_dir) {

  # Identify genus and species file entries by their role indicator
  # (stored in the file definition lists at the top of the script).
  # We match against the expected filename patterns so the function is robust
  # even if the order of file_results changes.
  genus_result   <- Filter(function(x) grepl("genus|toGenus",   x$filename), file_results)[[1]]
  species_result <- Filter(function(x) grepl("species|assignSpecies", x$filename), file_results)[[1]]

  lines <- c(
    # -----------------------------------------------------------------------
    # Header block
    # -----------------------------------------------------------------------
    paste0(db_name, " REFERENCE DATABASE — DOWNLOAD MANIFEST"),
    paste0("Database version : ", db_version),
    paste0("Zenodo record    : ", zenodo_record),
    paste0("Script           : download_reference_databases.R"),
    paste0("Manifest written : ", format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
    "",

    # -----------------------------------------------------------------------
    # File records
    # -----------------------------------------------------------------------
    "FILES",
    paste0(rep("-", 70), collapse = ""),

    # Genus file
    paste0("Role        : assignTaxonomy() — genus-level training set"),
    paste0("Filename    : ", genus_result$filename),
    paste0("Source URL  : ", genus_result$url),
    paste0("Local path  : ", genus_result$dest_path),
    paste0("Project path: ", project_relative_path(genus_result$dest_path, project_dir)),
    paste0("File size   : ", genus_result$file_size_human,
           " (", genus_result$file_size_bytes, " bytes)"),
    paste0("Verified    : ", genus_result$verification),
    paste0("Downloaded  : ", genus_result$download_timestamp),
    "",

    # Species file
    paste0("Role        : addSpecies() — species-level assignment set"),
    paste0("Filename    : ", species_result$filename),
    paste0("Source URL  : ", species_result$url),
    paste0("Local path  : ", species_result$dest_path),
    paste0("Project path: ", project_relative_path(species_result$dest_path, project_dir)),
    paste0("File size   : ", species_result$file_size_human,
           " (", species_result$file_size_bytes, " bytes)"),
    paste0("Verified    : ", species_result$verification),
    paste0("Downloaded  : ", species_result$download_timestamp),
    "",

    # -----------------------------------------------------------------------
    # R path assignments
    # -----------------------------------------------------------------------
    "R PATH ASSIGNMENTS",
    paste0(rep("-", 70), collapse = ""),
    "Paste these lines into your DADA2 notebooks or config script:",
    "",
    paste0(r_var_genus,   ' <- "', genus_result$dest_path,   '"'),
    paste0(r_var_species, ' <- "', species_result$dest_path, '"'),
    "",

    # -----------------------------------------------------------------------
    # Usage guidance
    # -----------------------------------------------------------------------
    "HOW TO USE IN DADA2",
    paste0(rep("-", 70), collapse = ""),
    paste0("# Step 1 — Assign taxonomy to genus level"),
    paste0("taxa <- assignTaxonomy(seqtab_nochim, ", r_var_genus,   ", multithread = TRUE)"),
    "",
    paste0("# Step 2 — Add species-level assignments (exact 100% identity matches)"),
    paste0("taxa <- addSpecies(taxa, ", r_var_species, ")"),
    "",

    # -----------------------------------------------------------------------
    # Portability note
    # -----------------------------------------------------------------------
    "NOTES",
    paste0(rep("-", 70), collapse = ""),
    "The absolute paths above describe the installation location at download time.",
    "For portable workflow code, construct paths from the project root (for example with here::here()).",
    "If the project moves, the Project path entries remain valid relative to its new root.",
    "The .fa.gz files can be used directly — DADA2 reads gzipped FASTA natively."
  )

  writeLines(lines, con = manifest_path, useBytes = TRUE)
  invisible(TRUE)
}


# ==============================================================================
# Setup: parse arguments and resolve paths
# ==============================================================================
arguments       <- parse_arguments()
project_dir     <- normalize_path_safe(arguments$project_dir)
trainsets_name  <- arguments$trainsets_name
force_download  <- arguments$force_download
timeout_seconds <- arguments$timeout_seconds
download_method <- arguments$download_method

assert_project_root(project_dir)

# Resolve the two database subfolders.
trainsets_dir   <- file.path(project_dir, "tools", trainsets_name)
gtdb_dir        <- file.path(trainsets_dir, "GTDB")
silva_dir       <- file.path(trainsets_dir, "SILVA")


# ==============================================================================
# Print banner and configuration summary
# ==============================================================================
cat("
╔══════════════════════════════════════════════════════════════════════════════╗
║        DADA2 16S WORKFLOW - REFERENCE DATABASE DOWNLOAD                      ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")

section_header("Configuration")
message_line("Project directory  : ", project_dir)
message_line("Trainsets folder   : ", trainsets_dir)
message_line("GTDB subfolder     : ", gtdb_dir)
message_line("SILVA subfolder    : ", silva_dir)
message_line("Force re-download  : ", force_download)
message_line("Download timeout   : ", timeout_seconds, " seconds")
message_line("Download method    : ", download_method)
message_line("Files to download  : ", length(gtdb_files) + length(silva_files))


# ==============================================================================
# Create directory structure
# ==============================================================================
# Create the tools/trainsets/GTDB and tools/trainsets/SILVA directories.
# dir.create() with recursive = TRUE and showWarnings = FALSE is idempotent:
# it silently does nothing if the directories already exist.
section_header("Creating directory structure")

dir.create(gtdb_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(silva_dir, recursive = TRUE, showWarnings = FALSE)

message_line("  ✓ ", trainsets_dir)
message_line("  ✓ ", gtdb_dir)
message_line("  ✓ ", silva_dir)


# ==============================================================================
# Download GTDB r220 reference files
# ==============================================================================
# GTDB (Genome Taxonomy Database) release r220 provides a phylogenetically
# consistent taxonomy derived from concatenated protein trees of archaeal and
# bacterial genomes. The DADA2-formatted files are hosted on Zenodo record
# 13984843 and were prepared by the DADA2 community.
#
# Two files are downloaded:
#   1. Genus-level training FASTA  → used by assignTaxonomy()
#   2. Species-level FASTA         → used by addSpecies()
section_header("Downloading GTDB r220 reference files")
message_line("Source: https://zenodo.org/records/13984843")

gtdb_results <- list()

for (file_def in gtdb_files) {
  dest_path <- file.path(gtdb_dir, file_def$filename)

  # ---------------------------------------------------------------------------
  # Print the file's role to the console so the user understands what each
  # download is for without needing to read the documentation separately.
  # ---------------------------------------------------------------------------
  message_line("\n  DADA2 role: ", file_def$dada2_role, " — ", file_def$description)

  result <- download_database_file(
    url             = file_def$url,
    dest_path       = dest_path,
    description     = file_def$description,
    timeout_seconds = timeout_seconds,
    download_method = download_method,
    force_download  = force_download,
    expected_size_bytes = file_def$expected_size_bytes
  )

  gtdb_results[[length(gtdb_results) + 1L]] <- result
}


# ==============================================================================
# Download SILVA v138.2 reference files
# ==============================================================================
# SILVA (https://www.arb-silva.de/) is the most widely used rRNA reference
# database. Version 138.2 is the latest SILVA release with DADA2-formatted
# training files and covers all three domains of life. The DADA2-formatted
# files are hosted on Zenodo record 14169026.
#
# Two files are downloaded:
#   1. Genus-level training FASTA  → used by assignTaxonomy()
#   2. Species-level FASTA         → used by addSpecies()
section_header("Downloading SILVA v138.2 reference files")
message_line("Source: https://zenodo.org/records/14169026")

silva_results <- list()

for (file_def in silva_files) {
  dest_path <- file.path(silva_dir, file_def$filename)

  message_line("\n  DADA2 role: ", file_def$dada2_role, " — ", file_def$description)

  result <- download_database_file(
    url             = file_def$url,
    dest_path       = dest_path,
    description     = file_def$description,
    timeout_seconds = timeout_seconds,
    download_method = download_method,
    force_download  = force_download,
    expected_size_bytes = file_def$expected_size_bytes
  )

  silva_results[[length(silva_results) + 1L]] <- result
}


# ==============================================================================
# Write per-database manifests
# ==============================================================================
# Each database subfolder receives a download_manifest.txt file that acts as
# a permanent, human-readable record of what is installed and how to use it.
# The manifest is plain text so it can be opened in any editor or read with
# `cat` without requiring any software beyond a terminal.
section_header("Writing download manifests")

# ------------------------------------------------------------------------------
# GTDB manifest
# ------------------------------------------------------------------------------
gtdb_manifest_path <- file.path(gtdb_dir, "download_manifest.txt")

write_database_manifest(
  manifest_path  = gtdb_manifest_path,
  db_name        = "GTDB",
  db_version     = "r220 (bac120_arc53_ssu)",
  zenodo_record  = "https://zenodo.org/records/13984843",
  file_results   = gtdb_results,
  r_var_genus    = "gtdb_genus_db",
  r_var_species  = "gtdb_species_db",
  project_dir    = project_dir
)
message_line("  ✓ GTDB manifest written  : ", gtdb_manifest_path)

# ------------------------------------------------------------------------------
# SILVA manifest
# ------------------------------------------------------------------------------
silva_manifest_path <- file.path(silva_dir, "download_manifest.txt")

write_database_manifest(
  manifest_path  = silva_manifest_path,
  db_name        = "SILVA",
  db_version     = "v138.2 (nr99)",
  zenodo_record  = "https://zenodo.org/records/14169026",
  file_results   = silva_results,
  r_var_genus    = "silva_genus_db",
  r_var_species  = "silva_species_db",
  project_dir    = project_dir
)
message_line("  ✓ SILVA manifest written : ", silva_manifest_path)


# ==============================================================================
# Download summary
# ==============================================================================
# Collect results across both databases and report a pass/fail count.
section_header("Download summary")

all_results <- c(gtdb_results, silva_results)
n_total     <- length(all_results)
n_skipped   <- sum(sapply(all_results, function(x) isTRUE(x$skipped)))
n_fresh     <- n_total - n_skipped

message_line("Total files    : ", n_total)
message_line("Newly downloaded: ", n_fresh)
message_line("Skipped (existed): ", n_skipped)

for (res in all_results) {
  status_label <- if (isTRUE(res$skipped)) "SKIPPED" else "DOWNLOADED"
  message_line("  [", status_label, "]  ", res$filename,
               "  (", res$file_size_human, ")")
}


# ==============================================================================
# Final R path assignments — copy-paste into DADA2 notebooks
# ==============================================================================
# Print the complete set of R variable assignments the user needs for their
# analysis. These same lines are also written into the manifest files so they
# are available offline without re-running the script.
cat("\n
╔══════════════════════════════════════════════════════════════════════════════╗
║              R PATH ASSIGNMENTS — PASTE INTO YOUR NOTEBOOKS                  ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")

# Retrieve the exact local paths from the download results for display.
# Using the result objects (rather than constructing paths from variables)
# ensures the printed paths exactly match the files on disk.
gtdb_genus_result   <- Filter(function(x) grepl("genus",         x$filename), gtdb_results)[[1]]
gtdb_species_result <- Filter(function(x) grepl("species",       x$filename), gtdb_results)[[1]]
silva_genus_result  <- Filter(function(x) grepl("toGenus",       x$filename), silva_results)[[1]]
silva_species_result <- Filter(function(x) grepl("assignSpecies", x$filename), silva_results)[[1]]

message_line("# ── GTDB r220 ───────────────────────────────────────────────────────────────")
message_line('gtdb_genus_db   <- "', gtdb_genus_result$dest_path,    '"')
message_line('gtdb_species_db <- "', gtdb_species_result$dest_path,  '"')
message_line("")
message_line("# ── SILVA v138.2 ─────────────────────────────────────────────────────────────")
message_line('silva_genus_db   <- "', silva_genus_result$dest_path,   '"')
message_line('silva_species_db <- "', silva_species_result$dest_path, '"')
message_line("")
message_line("# ── Usage ────────────────────────────────────────────────────────────────────")
message_line("# Assign taxonomy (choose GTDB or SILVA genus file):")
message_line("# taxa <- assignTaxonomy(seqtab_nochim, gtdb_genus_db,   multithread = TRUE)")
message_line("# taxa <- assignTaxonomy(seqtab_nochim, silva_genus_db,  multithread = TRUE)")
message_line("")
message_line("# Add species (use the matching species file):")
message_line("# taxa <- addSpecies(taxa, gtdb_species_db)")
message_line("# taxa <- addSpecies(taxa, silva_species_db)")
message_line("")
message_line("The same path assignments are saved in each database's download_manifest.txt.")
