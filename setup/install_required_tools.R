# ============================================================================
# Script name: install_required_tools.R
# Description:
#   Installer for four tools used by Steps 2, 3, and 6 of the DADA2 16S
#   amplicon workflow, all placed inside a project-local tools/ directory:
#
#     1. Cutadapt  — Python, installed via pip into an isolated venv
#                    (used by Step 3, Primer Trimming)
#     2. MultiQC   — Python, installed via pip into an isolated venv
#                    (used by Step 2, FastQC Quality Reports)
#     3. FastQC    — Java, installed by downloading and unzipping the
#                    official release from the Babraham Bioinformatics site
#                    (used by Step 2, FastQC Quality Reports)
#     4. FastTree  — C, installed by downloading the source file from
#                    Morgan Price's FastTree page and compiling locally
#                    (used by Step 6, Phylogenetic Tree)
#
#   The script is designed for Linux and macOS. Python tools intentionally
#   avoid relying on conda, pixi, or system-wide PATH state. FastTree is
#   always built from source to ensure compatibility with the local CPU
#   architecture (x86-64 and ARM64 are both supported). FastQC requires
#   a Java Runtime Environment (JRE) to already be present on the system.
#
#   A fifth external tool, PICRUSt2, is required ONLY by the optional Step 7
#   (16S copy-number correction) notebook and is intentionally NOT installed
#   by this script. PICRUSt2 ships as a conda package with many compiled
#   phylogenetics dependencies (HMMER, EPA-ng, gappa, SEPP) that are not
#   practical to manage inside a plain pip venv, so it is installed into its
#   own dedicated conda environment by a separate script:
#
#     setup/install_picrust2.sh
#
#   Run that script once per machine (requires conda/miniconda already
#   installed) if you plan to use the optional Step 7 notebook. No other
#   step needs PICRUSt2 -- Step 8 (microbial load correction) in particular
#   depends only on R packages, not on PICRUSt2 or any tool installed by
#   this script.
#
#   Expected project layout after installation:
#
#   <project_root>/
#   ├── tools/
#   │   ├── cutadapt/
#   │   │   ├── venv/
#   │   │   ├── requirements.txt
#   │   │   └── install_manifest.txt
#   │   ├── multiqc/
#   │   │   ├── venv/
#   │   │   ├── requirements.txt
#   │   │   └── install_manifest.txt
#   │   ├── FastQC/
#   │   │   ├── fastqc               <- executable (chmod +x by this script)
#   │   │   ├── fastqc.jar (+ libs)
#   │   │   └── install_manifest.txt
#   │   └── fasttree/
#   │       ├── FastTree             <- compiled binary
#   │       ├── FastTree.c           <- downloaded source (kept for reference)
#   │       └── install_manifest.txt
#   ├── data/
#   ├── results/
#   └── R/
#
#   Each tool directory contains an install_manifest.txt that records
#   the installation timestamp, source URL, detected version, and the
#   exact R path variable to use in downstream notebooks.
#
# References:
#   1) Cutadapt — virtual-environment installation:
#      https://cutadapt.readthedocs.io/en/stable/installation.html
#   2) MultiQC — requires Python 3.9+, local venv recommended:
#      https://docs.seqera.io/multiqc/getting_started/installation
#   3) FastQC — Java application, distributed as a zip archive:
#      https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#   4) FastTree — C source distributed by Morgan Price; compile locally
#      for correct architecture support (x86-64 and ARM64):
#      https://morgannprice.github.io/fasttree/
#
# Path guidance for downstream R scripts and notebooks:
#
#   cutadapt_path  <- file.path(project_root, "tools", "cutadapt", "venv", "bin", "cutadapt")
#   multiqc_path   <- file.path(project_root, "tools", "multiqc",  "venv", "bin", "multiqc")
#   fastqc_path    <- file.path(project_root, "tools", "FastQC", "fastqc")
#   fasttree_path  <- file.path(project_root, "tools", "fasttree", "FastTree")
#
# System prerequisites (checked at runtime with clear error messages):
#   - Python 3.9+  : required for cutadapt and multiqc
#       Debian/Ubuntu note: python3-venv must also be installed separately —
#         sudo apt-get install python3-venv
#   - Java JRE     : required for FastQC
#       Linux  : sudo apt-get install default-jre
#       macOS  : https://adoptium.net/  or  brew install --cask temurin
#   - C compiler   : required to compile FastTree (gcc, cc, or clang)
#       Linux  : sudo apt-get install gcc
#       macOS  : xcode-select --install   (provides Apple clang)
#
# Notes:
#   - Python venvs are not copy-portable across machines or paths because
#     they embed absolute interpreter references. Recreate them on each
#     target machine instead of copying the installed directories.
#   - FastTree compiled binaries are architecture-specific. Recompile on
#     each target machine. On macOS with Apple clang (no Homebrew GCC),
#     the OpenMP multi-threaded build will fall back automatically to a
#     single-threaded build, which is still fully functional.
#   - This script uses base R only and can run on a fresh system before
#     any CRAN or Bioconductor packages are installed.
# ============================================================================

# -------------------------------
# User-adjustable defaults
# -------------------------------
default_project_dir        <- getwd()
default_cutadapt_version   <- "5.2"    # Tested workflow default
default_multiqc_version    <- "1.35"   # Tested workflow default
default_fastqc_version     <- "0.12.1" # FastQC release to download from Babraham
default_force_reinstall    <- FALSE
default_create_wrapper_bin <- FALSE     # Python tools only


# -------------------------------
# Helper functions — shared
# -------------------------------

# Print a formatted message line.
#
message_line <- function(...) {
  cat(..., "\n", sep = "")
  invisible(NULL)
}

# Print a section header.
#
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

# Normalize a path without requiring that it already exists.
#
normalize_path_safe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# Refuse to create a parallel workflow tree when the script was launched from
# the wrong working directory. Users can select another valid clone explicitly
# with --project-dir.
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

# Move an existing installation aside so a replacement can be built without
# destroying the last working copy. The caller restores backup_path on error.
backup_existing_directory <- function(path) {
  if (!dir.exists(path)) return(NULL)
  backup_path <- tempfile(pattern = paste0(".", basename(path), ".backup-"), tmpdir = dirname(path))
  if (!file.rename(path, backup_path)) {
    stop("Could not preserve the existing installation before replacement: ", path, call. = FALSE)
  }
  backup_path
}

restore_directory_backup <- function(path, backup_path) {
  if (is.null(backup_path) || !dir.exists(backup_path)) return(invisible(NULL))
  if (dir.exists(path)) unlink(path, recursive = TRUE, force = TRUE)
  if (!file.rename(backup_path, path)) {
    warning("Could not restore the previous installation. It remains at: ", backup_path)
  }
  invisible(NULL)
}

# Return TRUE if a path is an existing executable file.
#
file_is_executable <- function(path) {
  file.exists(path) && !dir.exists(path) && file.access(path, mode = 1) == 0
}

# Check that a path exists.
#
assert_exists <- function(path, description) {
  if (!file.exists(path)) {
    stop(description, " does not exist: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

# Run a command and capture combined stdout/stderr.
#
# Stops with an error if the exit status is non-zero.
#
run_command <- function(command, args = character(), wd = NULL) {
  output_file <- tempfile(pattern = "cmd_out_", fileext = ".log")
  original_wd <- getwd()

  on.exit({
    unlink(output_file, force = TRUE)
    if (!is.null(wd)) {
      setwd(original_wd)
    }
  }, add = FALSE)

  if (!is.null(wd)) {
    setwd(wd)
  }

  exit_status <- system2(
    command = command,
    args    = args,
    stdout  = output_file,
    stderr  = output_file,
    wait    = TRUE
  )

  command_output <- readLines(output_file, warn = FALSE)

  if (!identical(exit_status, 0L)) {
    stop(
      paste0(
        "Command failed with exit status ", exit_status, ":\n",
        paste(c(command, args), collapse = " "), "\n\n",
        paste(command_output, collapse = "\n")
      ),
      call. = FALSE
    )
  }

  command_output
}

# Run a command and capture output regardless of exit status.
#
# Unlike run_command(), this variant never stops on non-zero exit.
# Useful for programs that print version info to stderr and exit non-zero
# when called without required arguments (e.g., FastTree).
#
run_command_capture <- function(command, args = character()) {
  output_file <- tempfile(pattern = "cmd_capture_", fileext = ".log")
  on.exit(unlink(output_file, force = TRUE), add = TRUE)

  system2(
    command = command,
    args    = args,
    stdout  = output_file,
    stderr  = output_file,
    wait    = TRUE
  )

  readLines(output_file, warn = FALSE)
}

# Compare semantic version strings with base R only.
#
compare_versions <- function(version_a, version_b) {
  numeric_a <- package_version(version_a)
  numeric_b <- package_version(version_b)
  if (numeric_a < numeric_b) return(-1L)
  if (numeric_a > numeric_b) return(1L)
  0L
}

# Detect the operating system in a simple Linux/macOS-aware form.
#
detect_operating_system <- function() {
  sysname <- tolower(Sys.info()[["sysname"]])
  if (identical(sysname, "darwin")) return("mac")
  if (identical(sysname, "linux"))  return("linux")
  sysname
}


# -------------------------------
# Helper functions — Python tools
# -------------------------------

# Determine whether a Python path looks like a system/Homebrew or managed env.
#
classify_python_path <- function(path) {
  if (grepl("^/usr/bin/", path))                                       return("system")
  if (grepl("^/opt/homebrew/bin/|^/usr/local/bin/", path))            return("homebrew_or_local")
  if (grepl("conda|miniconda|anaconda|mamba|micromamba|pixi", path))   return("env_managed")
  return("unknown")
}

# Parse "python3 --version" output into a version string.
#
extract_python_version <- function(version_output) {
  first_line <- paste(version_output, collapse = " ")
  extracted  <- sub(".*Python\\s+([0-9]+\\.[0-9]+(?:\\.[0-9]+)?).*$", "\\1", first_line)
  if (!grepl("^[0-9]+\\.[0-9]+", extracted)) return(NA_character_)
  extracted
}

# Score a Python interpreter candidate. Lower scores are preferred.
#
score_python_candidate <- function(path, os_name) {
  path_class <- classify_python_path(path)

  if (identical(os_name, "mac")) {
    if (identical(path, "/opt/homebrew/bin/python3")) return(1L)
    if (identical(path, "/usr/local/bin/python3"))    return(2L)
    if (identical(path, "/usr/bin/python3"))          return(3L)
  }

  if (identical(os_name, "linux")) {
    if (identical(path, "/usr/bin/python3"))          return(1L)
    if (identical(path, "/usr/local/bin/python3"))    return(2L)
  }

  if (identical(path_class, "system"))            return(10L)
  if (identical(path_class, "homebrew_or_local")) return(11L)
  if (identical(path_class, "env_managed"))       return(40L)
  50L
}

# Validate a candidate Python interpreter.
#
validate_python_candidate <- function(path) {
  if (is.na(path) || identical(path, "") || !file_is_executable(path)) return(NULL)

  version_output <- tryCatch(
    run_command(command = path, args = "--version"),
    error = function(e) NULL
  )
  if (is.null(version_output)) return(NULL)

  version_string <- extract_python_version(version_output)
  if (is.na(version_string)) return(NULL)

  if (compare_versions(version_string, "3.9") < 0L) return(NULL)

  list(
    path           = normalize_path_safe(path),
    version        = version_string,
    source_class   = classify_python_path(path),
    version_output = paste(version_output, collapse = " ")
  )
}

# Auto-detect the best Python interpreter for tool installation.
#
detect_best_python <- function(user_python = NULL) {
  os_name <- detect_operating_system()

  if (!is.null(user_python)) {
    validated <- validate_python_candidate(user_python)
    if (is.null(validated)) {
      stop(
        "The user-specified Python interpreter is invalid, not executable, or older than Python 3.9: ",
        user_python, call. = FALSE
      )
    }
    validated$selection_reason <- "user_specified"
    return(validated)
  }

  candidates <- character()

  if (identical(os_name, "mac")) {
    candidates <- c("/opt/homebrew/bin/python3", "/usr/local/bin/python3",
                    "/usr/bin/python3", Sys.which("python3"))
  } else if (identical(os_name, "linux")) {
    candidates <- c("/usr/bin/python3", "/usr/local/bin/python3", Sys.which("python3"))
  } else {
    candidates <- c(Sys.which("python3"), "/usr/bin/python3")
  }

  unique_candidates    <- unique(candidates[nzchar(candidates)])
  validated_candidates <- Filter(Negate(is.null), lapply(unique_candidates, validate_python_candidate))

  if (length(validated_candidates) == 0L) {
    stop(
      "No suitable Python 3.9+ interpreter was detected. Install Python 3 and rerun, ",
      "or specify one explicitly with --python.", call. = FALSE
    )
  }

  candidate_scores <- vapply(
    validated_candidates,
    function(x) score_python_candidate(x$path, os_name),
    integer(1)
  )

  best_candidate <- validated_candidates[[which.min(candidate_scores)]]
  best_candidate$selection_reason <- "auto_detected"
  best_candidate
}

# Build a pip package specification.
#
build_package_spec <- function(package_name, version = NULL) {
  if (is.null(version) || identical(version, "")) return(package_name)
  paste0(package_name, "==", version)
}

# Write a per-tool installation manifest for Python venv-based tools.
#
write_tool_manifest <- function(manifest_path, tool_info, project_dir) {
  manifest_lines <- c(
    "TOOL INSTALLATION MANIFEST",
    paste0("Generated: ",                         format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
    paste0("Project directory: ",                 project_dir),
    paste0("Tool name: ",                         tool_info$tool_name),
    paste0("Tool directory: ",                    tool_info$tool_dir),
    paste0("Virtual environment directory: ",     tool_info$venv_dir),
    paste0("Python interpreter used: ",           tool_info$python_path),
    paste0("Python version used: ",               tool_info$python_version),
    paste0("Python source classification: ",      tool_info$python_source_class),
    paste0("Requested package specification: ",   tool_info$package_spec),
    paste0("Detected tool version output: ",      tool_info$version_output),
    paste0("Executable path to use in R: ",       tool_info$executable_path),
    paste0("Requirements file: ",                 tool_info$requirements_file),
    if (!is.na(tool_info$wrapper_executable_path)) {
      paste0("Optional wrapper path: ", tool_info$wrapper_executable_path)
    } else {
      "Optional wrapper path: not created"
    },
    "",
    "Recommended R code:",
    paste0(tool_info$r_path_assignment),
    "",
    "Important usage guidance:",
    "  Use the explicit executable path above in downstream scripts and notebooks.",
    "  Do not rely on PATH lookups such as just calling the tool name, because that",
    "  can silently switch to another installation from conda, pixi, Homebrew, or system.",
    "",
    "Portability note:",
    "  This installation should be recreated on each target machine rather than copied,",
    "  because Python virtual environments embed absolute interpreter references."
  )

  writeLines(manifest_lines, con = manifest_path, useBytes = TRUE)
  invisible(TRUE)
}

# Install one Python-based tool into its own project-local virtual environment.
#
install_python_tool <- function(tool_name,
                                package_name,
                                expected_executable,
                                project_dir,
                                python_info,
                                version             = NULL,
                                force_reinstall     = FALSE,
                                create_wrapper_bin  = FALSE) {
  tools_dir          <- file.path(project_dir, "tools")
  tool_dir           <- file.path(tools_dir, tool_name)
  venv_dir           <- file.path(tool_dir, "venv")
  venv_python        <- file.path(venv_dir, "bin", "python")
  venv_pip           <- file.path(venv_dir, "bin", "pip")
  tool_executable    <- file.path(venv_dir, "bin", expected_executable)
  wrapper_dir        <- file.path(tool_dir, "bin")
  wrapper_executable <- file.path(wrapper_dir, expected_executable)
  requirements_file  <- file.path(tool_dir, "requirements.txt")
  manifest_file      <- file.path(tool_dir, "install_manifest.txt")

  section_header(paste("Installing", tool_name))
  message_line("Target tool directory: ", tool_dir)

  backup_path <- NULL
  install_completed <- FALSE
  if (dir.exists(tool_dir) && isTRUE(force_reinstall)) {
    message_line("Preserving existing installation until its replacement is verified.")
    backup_path <- backup_existing_directory(tool_dir)
  }
  on.exit({
    if (!is.null(backup_path) && dir.exists(backup_path)) {
      if (isTRUE(install_completed)) unlink(backup_path, recursive = TRUE, force = TRUE)
      else restore_directory_backup(tool_dir, backup_path)
    }
  }, add = TRUE)

  dir.create(tool_dir, recursive = TRUE, showWarnings = FALSE)

  if (!dir.exists(venv_dir)) {
    message_line("Creating virtual environment with: ", python_info$path)
    run_command(command = python_info$path, args = c("-m", "venv", venv_dir))
  } else {
    message_line("Existing virtual environment detected; reusing it.")
  }

  message_line("Upgrading pip inside the virtual environment.")
  run_command(command = venv_python, args = c("-m", "pip", "install", "--upgrade", "pip"))

  package_spec <- build_package_spec(package_name = package_name, version = version)
  message_line("Installing package: ", package_spec)
  run_command(command = venv_python, args = c("-m", "pip", "install", package_spec))

  assert_exists(tool_executable, paste(tool_name, "executable"))

  if (isTRUE(create_wrapper_bin)) {
    dir.create(wrapper_dir, recursive = TRUE, showWarnings = FALSE)
    if (file.exists(wrapper_executable) || nzchar(Sys.readlink(wrapper_executable))) {
      unlink(wrapper_executable, force = TRUE)
    }
    file.symlink(
      from = file.path("..", "venv", "bin", expected_executable),
      to   = wrapper_executable
    )
    message_line("Created wrapper symlink: ", wrapper_executable)
  }

  detected_version <- paste(
    run_command(command = tool_executable, args = "--version"),
    collapse = " "
  )
  message_line("Detected version output: ", detected_version)

  freeze_output <- run_command(command = venv_pip, args = "freeze")
  writeLines(freeze_output, con = requirements_file, useBytes = TRUE)
  message_line("Wrote requirements file: ", requirements_file)

  r_assignment <- paste0(
    if (identical(tool_name, "multiqc")) "multiqc_path" else paste0(tool_name, "_path"),
    " <- \"", tool_executable, "\""
  )

  tool_info <- list(
    tool_name               = tool_name,
    tool_dir                = tool_dir,
    venv_dir                = venv_dir,
    executable_path         = tool_executable,
    wrapper_executable_path = if (isTRUE(create_wrapper_bin)) wrapper_executable else NA_character_,
    version_output          = detected_version,
    package_spec            = package_spec,
    requirements_file       = requirements_file,
    python_path             = python_info$path,
    python_version          = python_info$version,
    python_source_class     = python_info$source_class,
    r_path_assignment       = r_assignment,
    manifest_file           = manifest_file
  )

  write_tool_manifest(
    manifest_path = manifest_file,
    tool_info     = tool_info,
    project_dir   = project_dir
  )
  message_line("Wrote per-tool manifest: ", manifest_file)

  install_completed <- TRUE

  tool_info
}


# -------------------------------
# Helper functions — binary tools
# -------------------------------

# Write a per-tool installation manifest for non-Python binary tools.
#
# Used by install_fastqc() and install_fasttree() in place of the Python-
# specific write_tool_manifest(). Records the download URL, detected version,
# executable path, and the R code line to paste into notebooks.
#
write_binary_tool_manifest <- function(manifest_path,
                                       tool_name,
                                       tool_dir,
                                       executable_path,
                                       source_url,
                                       version_output,
                                       r_path_assignment,
                                       extra_notes   = character(),
                                       project_dir) {
  manifest_lines <- c(
    "TOOL INSTALLATION MANIFEST",
    paste0("Generated: ",                format(Sys.time(), tz = Sys.timezone(), usetz = TRUE)),
    paste0("Project directory: ",        project_dir),
    paste0("Tool name: ",               tool_name),
    paste0("Tool directory: ",          tool_dir),
    paste0("Executable path in R: ",    executable_path),
    paste0("Source URL: ",             source_url),
    paste0("Detected version output: ", version_output),
    "",
    "Recommended R code:",
    r_path_assignment,
    "",
    "Important usage guidance:",
    "  Use the explicit executable path above in downstream scripts and notebooks.",
    "  Do not rely on PATH lookups such as just calling the tool name.",
    "",
    "Portability note:",
    "  Binaries are architecture-specific. Recreate this installation on each",
    "  target machine rather than copying the compiled or unpacked binary."
  )

  if (length(extra_notes) > 0L) {
    manifest_lines <- c(manifest_lines, "", "Notes:", extra_notes)
  }

  writeLines(manifest_lines, con = manifest_path, useBytes = TRUE)
  invisible(TRUE)
}

# Check that a Java Runtime Environment is available on PATH.
#
# FastQC is a Java application and will not run without a JRE. This function
# performs a lightweight check before attempting the installation so the user
# gets a clear error message rather than a cryptic Java startup failure.
#
check_java <- function() {
  java_path <- Sys.which("java")
  if (!nzchar(java_path)) {
    stop(
      "Java was not found on PATH. FastQC requires a Java Runtime Environment (JRE).\n",
      "  On Linux  : sudo apt-get install default-jre\n",
      "  On macOS  : install from https://adoptium.net/ or via Homebrew (brew install --cask temurin)",
      call. = FALSE
    )
  }

  # `java -version` writes to stderr and always exits 0 on modern JREs.
  version_lines <- run_command_capture(command = java_path, args = "-version")
  version_string <- paste(version_lines, collapse = " ")
  message_line("  Java detected: ", version_string)
  version_string
}

# Install FastQC from the Babraham Bioinformatics zip release.
#
# Downloads fastqc_v{version}.zip, unpacks it into tools/FastQC/, and makes
# the fastqc shell-script launcher executable. The unpacked directory already
# contains all required JAR files so no further build step is needed.
#
install_fastqc <- function(project_dir,
                           version         = "0.12.1",
                           force_reinstall = FALSE) {

  # ---------------------------------------------------------------------------
  # Resolve paths
  # ---------------------------------------------------------------------------
  tools_dir     <- file.path(project_dir, "tools")
  tool_dir      <- file.path(tools_dir, "FastQC")
  fastqc_exe    <- file.path(tool_dir, "fastqc")
  manifest_file <- file.path(tool_dir, "install_manifest.txt")

  # URL pattern used by Babraham Bioinformatics for all FastQC releases.
  download_url  <- paste0(
    "https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v",
    version, ".zip"
  )

  section_header("Installing FastQC")
  message_line("Target directory : ", tool_dir)
  message_line("Version          : ", version)
  message_line("Download URL     : ", download_url)

  # ---------------------------------------------------------------------------
  # Skip if already installed and executable (unless force_reinstall)
  # ---------------------------------------------------------------------------
  if (file_is_executable(fastqc_exe) && !isTRUE(force_reinstall)) {
    detected <- paste(run_command_capture(fastqc_exe, "--version"), collapse = " ")
    message_line("FastQC already installed and executable — skipping download.")
    message_line("  Detected: ", detected)
    message_line("  Use --force-reinstall to replace the existing installation.")

    return(list(
      executable_path = fastqc_exe,
      manifest_file   = manifest_file,
      version_output  = detected,
      skipped         = TRUE
    ))
  }

  # ---------------------------------------------------------------------------
  # Check Java availability before downloading anything
  # ---------------------------------------------------------------------------
  message_line("Checking for Java Runtime Environment...")
  java_version <- check_java()

  # Preserve any incomplete, stale, or explicitly replaced directory until the
  # new archive has been unpacked and its executable verified.
  backup_path <- NULL
  install_completed <- FALSE
  if (dir.exists(tool_dir)) {
    message_line("Preserving existing FastQC directory until replacement is verified.")
    backup_path <- backup_existing_directory(tool_dir)
  }
  on.exit({
    if (!is.null(backup_path) && dir.exists(backup_path)) {
      if (isTRUE(install_completed)) unlink(backup_path, recursive = TRUE, force = TRUE)
      else restore_directory_backup(tool_dir, backup_path)
    }
  }, add = TRUE)

  # ---------------------------------------------------------------------------
  # Download the zip archive to a temporary file
  # ---------------------------------------------------------------------------
  zip_tmp <- tempfile(pattern = "fastqc_", fileext = ".zip")
  on.exit(unlink(zip_tmp, force = TRUE), add = TRUE)

  # Extend the download timeout: the FastQC zip is ~12 MB and can be slow on
  # institutional networks. Restore the original option on exit.
  original_timeout <- getOption("timeout")
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout = 600L)

  message_line("Downloading FastQC zip...")
  tryCatch({
    download.file(url = download_url, destfile = zip_tmp, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    stop("Failed to download FastQC: ", conditionMessage(e), call. = FALSE)
  })

  if (!file.exists(zip_tmp) || file.info(zip_tmp)$size == 0L) {
    stop("FastQC zip downloaded as empty file. Check the URL and network.", call. = FALSE)
  }
  message_line("  Downloaded: ", round(file.info(zip_tmp)$size / 1e6, 1), " MB")

  # ---------------------------------------------------------------------------
  # Unzip into tools/
  # ---------------------------------------------------------------------------
  # The zip archive contains a single top-level directory named "FastQC/".
  # Unpacking with exdir = tools_dir places all contents at tools/FastQC/.
  message_line("Unpacking into: ", tools_dir)
  dir.create(tools_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch({
    unzip(zipfile = zip_tmp, exdir = tools_dir)
  }, error = function(e) {
    stop("Failed to unzip FastQC archive: ", conditionMessage(e), call. = FALSE)
  })

  assert_exists(tool_dir,   "FastQC directory after unzip")
  assert_exists(fastqc_exe, "FastQC launcher script after unzip")

  # ---------------------------------------------------------------------------
  # Make the fastqc launcher script executable
  # ---------------------------------------------------------------------------
  # The zip ships the launcher without the executable bit set, so we must apply
  # chmod +x explicitly. Sys.chmod() takes an octal string; "0755" grants
  # rwxr-xr-x (owner can execute, group and others can read and execute).
  Sys.chmod(fastqc_exe, mode = "0755")

  if (!file_is_executable(fastqc_exe)) {
    stop("chmod +x failed for: ", fastqc_exe, call. = FALSE)
  }

  # ---------------------------------------------------------------------------
  # Verify the installation by running fastqc --version
  # ---------------------------------------------------------------------------
  detected_version <- paste(
    run_command(command = fastqc_exe, args = "--version"),
    collapse = " "
  )
  message_line("Detected version: ", detected_version)

  # ---------------------------------------------------------------------------
  # Write the manifest
  # ---------------------------------------------------------------------------
  r_assignment <- paste0('fastqc_path <- "', fastqc_exe, '"')

  write_binary_tool_manifest(
    manifest_path     = manifest_file,
    tool_name         = "FastQC",
    tool_dir          = tool_dir,
    executable_path   = fastqc_exe,
    source_url        = download_url,
    version_output    = detected_version,
    r_path_assignment = r_assignment,
    extra_notes       = c(
      paste0("Java version at install time: ", java_version),
      "FastQC requires Java to run. Ensure a JRE is on PATH when using this tool."
    ),
    project_dir       = project_dir
  )
  message_line("Wrote manifest: ", manifest_file)

  install_completed <- TRUE

  list(
    executable_path = fastqc_exe,
    manifest_file   = manifest_file,
    version_output  = detected_version,
    skipped         = FALSE
  )
}

# Find a working C compiler on the current system.
#
# Checks for gcc first (preferred because it supports -fopenmp for
# multi-threaded FastTree builds), then falls back to cc (which may
# be clang on macOS). Returns NULL if no compiler is found.
#
find_c_compiler <- function() {
  candidates <- c("gcc", "cc", "clang")

  for (compiler_name in candidates) {
    compiler_path <- Sys.which(compiler_name)
    if (nzchar(compiler_path)) {
      # Quick sanity check: ask the compiler for its version.
      version_out <- tryCatch(
        run_command_capture(command = compiler_path, args = "--version"),
        error = function(e) character()
      )
      if (length(version_out) > 0L) {
        message_line("  C compiler found: ", compiler_name, " (", compiler_path, ")")
        message_line("  Version: ", version_out[[1]])
        return(list(path = compiler_path, name = compiler_name))
      }
    }
  }

  NULL
}

# Install FastTree by downloading the C source and compiling locally.
#
# FastTree is distributed as a single C source file. Compiling from source
# ensures compatibility with the local CPU architecture (x86-64 or ARM64)
# without relying on pre-built binaries that may not match the host system.
#
# Compilation strategy:
#   1. First attempt: gcc with -fopenmp (enables multi-threaded FastTreeMP).
#      This is the preferred build because it can use all available CPU cores.
#   2. Fallback:      gcc/cc without -fopenmp (single-threaded but portable).
#      Used automatically if the OpenMP attempt fails (e.g., on macOS with
#      Apple clang which does not ship libomp by default).
#
# The resulting binary is named FastTree and placed at tools/fasttree/FastTree.
# The source file is kept alongside it for reference and reproducibility.
#
install_fasttree <- function(project_dir,
                             force_reinstall = FALSE) {

  # ---------------------------------------------------------------------------
  # Resolve paths
  # ---------------------------------------------------------------------------
  tools_dir     <- file.path(project_dir, "tools")
  tool_dir      <- file.path(tools_dir, "fasttree")
  fasttree_exe  <- file.path(tool_dir, "FastTree")
  source_file   <- file.path(tool_dir, "FastTree.c")
  manifest_file <- file.path(tool_dir, "install_manifest.txt")

  source_url <- "https://morgannprice.github.io/fasttree/FastTree.c"

  section_header("Installing FastTree")
  message_line("Target directory : ", tool_dir)
  message_line("Source URL       : ", source_url)
  message_line("Strategy         : compile from C source (architecture-independent)")

  # ---------------------------------------------------------------------------
  # Skip if already installed and executable (unless force_reinstall)
  # ---------------------------------------------------------------------------
  if (file_is_executable(fasttree_exe) && !isTRUE(force_reinstall)) {
    # FastTree only prints its version when stdin is a TTY (isatty check in main()).
    # When called from R via system2(), stdin is not a TTY, so bare invocation
    # blocks waiting for sequence data. Passing -help always prints the version
    # to stderr and exits immediately, regardless of TTY state.
    version_lines <- run_command_capture(fasttree_exe, args = "-help")
    detected <- grep("FastTree [0-9]", version_lines, value = TRUE)
    detected <- if (length(detected) > 0L) detected[[1]] else paste(version_lines[1:min(2, length(version_lines))], collapse = " ")
    message_line("FastTree already installed and executable — skipping compilation.")
    message_line("  Detected: ", detected)
    message_line("  Use --force-reinstall to recompile.")

    return(list(
      executable_path = fasttree_exe,
      manifest_file   = manifest_file,
      version_output  = detected,
      skipped         = TRUE
    ))
  }

  backup_path <- NULL
  install_completed <- FALSE
  if (dir.exists(tool_dir) && isTRUE(force_reinstall)) {
    message_line("Preserving existing FastTree directory until replacement is verified.")
    backup_path <- backup_existing_directory(tool_dir)
  }
  on.exit({
    if (!is.null(backup_path) && dir.exists(backup_path)) {
      if (isTRUE(install_completed)) unlink(backup_path, recursive = TRUE, force = TRUE)
      else restore_directory_backup(tool_dir, backup_path)
    }
  }, add = TRUE)

  # ---------------------------------------------------------------------------
  # Find a C compiler
  # ---------------------------------------------------------------------------
  message_line("Searching for a C compiler...")
  compiler <- find_c_compiler()

  if (is.null(compiler)) {
    stop(
      "No C compiler found (tried: gcc, cc, clang).\n",
      "  On Linux  : sudo apt-get install gcc\n",
      "  On macOS  : xcode-select --install   (installs Apple clang)\n",
      "  Or install gcc via Homebrew: brew install gcc",
      call. = FALSE
    )
  }

  # ---------------------------------------------------------------------------
  # Create tool directory and download the source file
  # ---------------------------------------------------------------------------
  dir.create(tool_dir, recursive = TRUE, showWarnings = FALSE)

  message_line("Downloading FastTree.c...")
  original_timeout <- getOption("timeout")
  on.exit(options(timeout = original_timeout), add = TRUE)
  options(timeout = 300L)

  tryCatch({
    download.file(url = source_url, destfile = source_file, mode = "wb", quiet = FALSE)
  }, error = function(e) {
    stop("Failed to download FastTree.c: ", conditionMessage(e), call. = FALSE)
  })

  if (!file.exists(source_file) || file.info(source_file)$size == 0L) {
    stop("FastTree.c downloaded as empty file. Check the URL and network.", call. = FALSE)
  }
  message_line("  Downloaded: ", round(file.info(source_file)$size / 1e3, 1), " KB")

  # ---------------------------------------------------------------------------
  # Attempt 1: compile with OpenMP for multi-threaded support
  # ---------------------------------------------------------------------------
  # -DOPEN_MP and -fopenmp enable the OpenMP parallelism that lets FastTree
  # use multiple CPU cores for tree building. This is the recommended build.
  # Apple clang on macOS does not include libomp, so this step may fail there.
  compile_flags_omp <- c(
    "-O3", "-finline-functions", "-funroll-loops",
    "-Wall", "-fomit-frame-pointer",
    "-DOPEN_MP", "-fopenmp",
    "-o", fasttree_exe,
    source_file,
    "-lm"
  )

  openmp_success <- FALSE
  message_line("Attempting OpenMP build (multi-threaded)...")

  omp_result <- tryCatch({
    run_command(command = compiler$path, args = compile_flags_omp, wd = tool_dir)
    TRUE
  }, error = function(e) {
    message_line("  OpenMP build failed (", conditionMessage(e), ")")
    FALSE
  })

  if (isTRUE(omp_result) && file_is_executable(fasttree_exe)) {
    openmp_success <- TRUE
    message_line("  OpenMP build succeeded.")
  } else {
    # Clean up any partial output from the failed attempt.
    if (file.exists(fasttree_exe)) unlink(fasttree_exe, force = TRUE)
  }

  # ---------------------------------------------------------------------------
  # Attempt 2: fallback — compile without OpenMP (single-threaded)
  # ---------------------------------------------------------------------------
  if (!openmp_success) {
    compile_flags_basic <- c(
      "-O3",
      "-o", fasttree_exe,
      source_file,
      "-lm"
    )

    message_line("Falling back to single-threaded build (no OpenMP)...")
    tryCatch({
      run_command(command = compiler$path, args = compile_flags_basic, wd = tool_dir)
    }, error = function(e) {
      stop(
        "FastTree compilation failed with both OpenMP and basic flags.\n",
        "  Compiler error: ", conditionMessage(e), "\n",
        "  Check that gcc is installed and the source file downloaded correctly.",
        call. = FALSE
      )
    })

    if (!file_is_executable(fasttree_exe)) {
      stop("Compilation appeared to succeed but the binary is missing or not executable: ",
           fasttree_exe, call. = FALSE)
    }
    message_line("  Single-threaded build succeeded.")
  }

  # ---------------------------------------------------------------------------
  # Detect the installed version
  # ---------------------------------------------------------------------------
  # FastTree's version string is only printed when stdin is a TTY (see the
  # isatty() guard in main()). system2() does not provide a TTY, so calling
  # FastTree with no arguments would block on stdin. Passing -help forces the
  # version/usage output to stderr immediately and exits non-zero — which
  # run_command_capture() handles correctly by ignoring the exit status.
  version_lines <- run_command_capture(fasttree_exe, args = "-help")
  detected_version <- grep("FastTree [0-9]", version_lines, value = TRUE)
  detected_version <- if (length(detected_version) > 0L) {
    detected_version[[1]]
  } else {
    paste(version_lines[1:min(2, length(version_lines))], collapse = " ")
  }
  message_line("Detected version: ", detected_version)

  build_type <- if (openmp_success) "OpenMP (multi-threaded)" else "basic (single-threaded, no OpenMP)"

  # ---------------------------------------------------------------------------
  # Write the manifest
  # ---------------------------------------------------------------------------
  r_assignment <- paste0('fasttree_path <- "', fasttree_exe, '"')

  write_binary_tool_manifest(
    manifest_path     = manifest_file,
    tool_name         = "FastTree",
    tool_dir          = tool_dir,
    executable_path   = fasttree_exe,
    source_url        = source_url,
    version_output    = detected_version,
    r_path_assignment = r_assignment,
    extra_notes       = c(
      paste0("Compiler used: ", compiler$name, " (", compiler$path, ")"),
      paste0("Build type: ", build_type),
      paste0("Source file retained at: ", source_file),
      "This binary is architecture-specific. Recompile on each new target machine."
    ),
    project_dir       = project_dir
  )
  message_line("Wrote manifest: ", manifest_file)

  install_completed <- TRUE

  list(
    executable_path = fasttree_exe,
    manifest_file   = manifest_file,
    version_output  = detected_version,
    build_type      = build_type,
    skipped         = FALSE
  )
}


# -------------------------------
# Argument parsing
# -------------------------------

# Parse command-line arguments.
#
parse_arguments <- function() {
  raw_args <- commandArgs(trailingOnly = TRUE)

  parsed <- list(
    project_dir        = default_project_dir,
    python             = NULL,
    cutadapt_version   = default_cutadapt_version,
    multiqc_version    = default_multiqc_version,
    fastqc_version     = default_fastqc_version,
    force_reinstall    = default_force_reinstall,
    create_wrapper_bin = default_create_wrapper_bin
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
    } else if (current_arg == "--base-dir") {
      parsed$project_dir <- dirname(read_option_value(current_arg))
    } else if (current_arg == "--python") {
      parsed$python <- read_option_value(current_arg)
    } else if (current_arg == "--cutadapt-version") {
      parsed$cutadapt_version <- read_option_value(current_arg)
    } else if (current_arg == "--multiqc-version") {
      parsed$multiqc_version <- read_option_value(current_arg)
    } else if (current_arg == "--fastqc-version") {
      parsed$fastqc_version <- read_option_value(current_arg)
    } else if (current_arg == "--force-reinstall") {
      parsed$force_reinstall <- TRUE
    } else if (current_arg == "--create-wrapper-bin") {
      parsed$create_wrapper_bin <- TRUE
    } else if (current_arg %in% c("-h", "--help")) {
      message_line(
        "Usage: Rscript install_required_tools.R [options]\n",
        "\n",
        "Options:\n",
        "  --project-dir <path>          DADA2 project root containing tools/\n",
        "  --base-dir <path>             Backward-compatible alias; project dir is the parent\n",
        "  --python <path>               Explicit Python 3.9+ interpreter to use\n",
        "  --cutadapt-version <ver>      Cutadapt version pin (default: 5.2)\n",
        "  --multiqc-version <ver>       MultiQC version pin (default: 1.35)\n",
        "  --fastqc-version <ver>        FastQC release to download (default: 0.12.1)\n",
        "  --force-reinstall             Replace existing tools, restoring them if installation fails\n",
        "  --create-wrapper-bin          Create tool-local bin/<tool> symlink wrappers (Python tools)\n",
        "  -h, --help                    Show this help message\n"
      )
      quit(save = "no", status = 0)
    } else {
      stop("Unknown argument: ", current_arg, call. = FALSE)
    }

    i <- i + 1L
  }

  parsed
}

# Resolve and normalise the project directory.
#
resolve_project_dir <- function(project_dir = NULL) {
  candidate <- if (is.null(project_dir) || identical(project_dir, "")) getwd() else project_dir
  normalize_path_safe(candidate)
}

# Create core DADA2 project directories if missing.
#
ensure_dada2_project_structure <- function(project_dir) {
  core_dirs <- c(
    "tools",
    file.path("data", "fastq"),
    "results",
    file.path("R", "notebooks")
  )

  for (relative_dir in core_dirs) {
    dir.create(file.path(project_dir, relative_dir), recursive = TRUE, showWarnings = FALSE)
  }

  invisible(TRUE)
}


# -------------------------------
# Setup: parse arguments
# -------------------------------
arguments          <- parse_arguments()
project_dir        <- resolve_project_dir(arguments$project_dir)
cutadapt_version   <- arguments$cutadapt_version
multiqc_version    <- arguments$multiqc_version
fastqc_version     <- arguments$fastqc_version
force_reinstall    <- arguments$force_reinstall
create_wrapper_bin <- arguments$create_wrapper_bin

assert_project_root(project_dir)

# Detect the best Python interpreter once; shared by both Python tools.
python_info <- detect_best_python(user_python = arguments$python)


# -------------------------------
# Configuration summary
# -------------------------------
section_header("Validating configuration")
message_line("Operating system:        ", detect_operating_system())
message_line("Project directory:       ", project_dir)
message_line("Tools directory:         ", file.path(project_dir, "tools"))
message_line("Selected Python path:    ", python_info$path)
message_line("Selected Python version: ", python_info$version)
message_line("Python source class:     ", python_info$source_class)
message_line("Selection mode:          ", python_info$selection_reason)
message_line("Cutadapt version:        ", cutadapt_version)
message_line("MultiQC version:         ", multiqc_version)
message_line("FastQC version:          ", fastqc_version)
message_line("FastTree:                compiled from source")
message_line("Force reinstall:         ", force_reinstall)
message_line("Create wrapper bin:      ", create_wrapper_bin)

ensure_dada2_project_structure(project_dir = project_dir)


# -------------------------------
# Install all four tools
# -------------------------------

cutadapt_info <- install_python_tool(
  tool_name          = "cutadapt",
  package_name       = "cutadapt",
  expected_executable = "cutadapt",
  project_dir        = project_dir,
  python_info        = python_info,
  version            = cutadapt_version,
  force_reinstall    = force_reinstall,
  create_wrapper_bin = create_wrapper_bin
)

multiqc_info <- install_python_tool(
  tool_name          = "multiqc",
  package_name       = "multiqc",
  expected_executable = "multiqc",
  project_dir        = project_dir,
  python_info        = python_info,
  version            = multiqc_version,
  force_reinstall    = force_reinstall,
  create_wrapper_bin = create_wrapper_bin
)

fastqc_info <- install_fastqc(
  project_dir     = project_dir,
  version         = fastqc_version,
  force_reinstall = force_reinstall
)

fasttree_info <- install_fasttree(
  project_dir     = project_dir,
  force_reinstall = force_reinstall
)


# -------------------------------
# Final summary
# -------------------------------
section_header("Recommended R paths for downstream scripts")
message_line('cutadapt_path  <- "', cutadapt_info$executable_path,  '"')
message_line('multiqc_path   <- "', multiqc_info$executable_path,   '"')
message_line('fastqc_path    <- "', fastqc_info$executable_path,    '"')
message_line('fasttree_path  <- "', fasttree_info$executable_path,  '"')

section_header("Manifest files written")
message_line("Cutadapt  : ", cutadapt_info$manifest_file)
message_line("MultiQC   : ", multiqc_info$manifest_file)
message_line("FastQC    : ", fastqc_info$manifest_file)
message_line("FastTree  : ", fasttree_info$manifest_file)

message_line(
  "\nAll tools installed successfully.\n",
  "correct project-local installations are always used."
)

section_header("Optional: PICRUSt2 (required only for the optional Step 7 notebook)")
message_line(
  "This script does not install PICRUSt2. If you plan to run the optional\n",
  "Step 7 (16S copy-number correction) notebook, install PICRUSt2 into its own\n",
  "conda environment by running, once per machine:\n\n",
  "  setup/install_picrust2.sh\n\n",
  "No other step needs PICRUSt2 -- Step 8 (microbial load correction) in\n",
  "particular depends only on R packages."
)
