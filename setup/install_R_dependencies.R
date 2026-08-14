################################################################################
# Script: install_R_dependencies.R
# Purpose:
#   - Install all R packages required to run the DADA2 16S amplicon sequencing
#     workflow (Steps 1-9, including the optional Step 7 copy-number-correction
#     and Step 8 microbial-load-correction notebooks), including core
#     data-manipulation libraries, Bioconductor packages (dada2, phyloseq,
#     ShortRead, Biostrings, DECIPHER), visualization tools, and the
#     Shiny-based interactive front-end packages.
#   - Provide clear progress feedback at every installation step so the user
#     can follow along and quickly identify failures.
#   - Verify that every package loaded correctly after installation and print a
#     human-readable summary showing how many packages succeeded or failed.
#   - Report the installed version of each critical Bioconductor package and the
#     active R / Bioconductor version for reproducibility documentation.
#
# Note on optional Steps 7 and 8:
#   Every R package used by the optional Step 7 (16S copy-number correction)
#   and Step 8 (microbial load / QMP correction) notebooks — Biostrings,
#   phyloseq, openxlsx, dplyr, readr, here — is already installed by this
#   script as part of the core dependency set below; no separate R
#   installation step is required for either optional notebook. Step 7 does,
#   however, require the external PICRUSt2 tool, which is NOT an R package
#   and is therefore installed separately — see setup/install_picrust2.sh.
#
# When to run this script:
#   Run this script ONCE after cloning or downloading the repository for the
#   first time, or whenever you set up a new R environment (e.g., a new conda
#   environment, a new RStudio project on a different machine, or after a major
#   R version upgrade). You do NOT need to re-run it on every analysis session.
#
# Prerequisites:
#   - R >= 4.1.0  (required by current DADA2 and phyloseq releases)
#   - An active internet connection (CRAN and Bioconductor servers must be reachable)
#   - Sufficient disk space (~2-3 GB for all packages and their compiled objects)
#
# Output:
#   - Packages are installed into your active R library path (.libPaths()[1]).
#   - Package files are written only to the active R library; project data and
#     analysis results are not modified.
#   - A final summary table is printed to the console showing success/failure
#     status and version strings for the key Bioconductor packages.
#
# Troubleshooting:
#   - If a CRAN package fails, try installing it manually:
#       install.packages("<package_name>", dependencies = TRUE)
#   - If a Bioconductor package fails, ensure BiocManager is up to date:
#       install.packages("BiocManager")
#       BiocManager::install("<package_name>", update = FALSE, ask = FALSE)
#   - On Linux, some packages (e.g., xml2, curl, openssl) require system
#     libraries. Install them with your package manager first, for example:
#       sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev
#   - On macOS with Apple Silicon, use an ARM-native R build from CRAN.
#
################################################################################


# ==============================================================================
# Print banner
# ==============================================================================
# A visual banner is printed at startup so the user can immediately confirm that
# the correct script is running. This is especially helpful when sourcing from a
# larger workflow or calling via Rscript in a batch job.
cat("
╔══════════════════════════════════════════════════════════════════════════════╗
║           DADA2 16S WORKFLOW - DEPENDENCY INSTALLATION                       ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")


# ==============================================================================
# Helper function: install_if_missing
# ==============================================================================
# Rather than calling install.packages() or BiocManager::install() directly
# everywhere in the script, we centralise the install-and-verify logic in this
# single helper. It checks whether the package is already present (skipping the
# download if it is), attempts the appropriate install command, and immediately
# confirms whether the install succeeded.
#
# Arguments:
#   packages  : Character vector of package names to check and install.
#   source    : Either "CRAN" (default) or "Bioconductor". Controls which
#               backend is used for packages not already present.
#
# Behaviour:
#   - Packages that are already available produce a "✓ already installed" line.
#   - Packages that need installing produce "Installing: <name> ..." followed
#     by either "✓ installed successfully" or "✗ Failed to install <name>".
#   - The function does NOT stop the script on failure; it merely reports the
#     failure so all packages are attempted before the user reviews the summary.
install_if_missing <- function(packages, source = "CRAN") {
    for (pkg in packages) {

        # -------------------------------------------------------------------
        # Check availability before attempting installation
        # -------------------------------------------------------------------
        # requireNamespace() loads the package namespace without attaching it,
        # making it the lightweight way to test whether a package is present.
        # quietly = TRUE suppresses the "there is no package called ..." message
        # so our own formatted output remains clean.
        if (!requireNamespace(pkg, quietly = TRUE)) {

            cat("Installing:", pkg, "...\n")

            # ---------------------------------------------------------------
            # Install from the appropriate source
            # ---------------------------------------------------------------
            install_error <- tryCatch({
                if (source == "CRAN") {
                    # dependencies = TRUE installs Suggests as well as Imports,
                    # which is important for packages like knitr and ggplot2 that
                    # depend on many optional components at runtime.
                    install.packages(pkg, dependencies = TRUE, quiet = TRUE)

                } else if (source == "Bioconductor") {
                    # update = FALSE prevents BiocManager from upgrading already-
                    # installed packages, which keeps the script idempotent.
                    # ask = FALSE allows unattended installation without user prompts.
                    BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
                } else {
                    stop("Unknown package source: ", source)
                }
                NULL
            }, error = function(e) conditionMessage(e))

            # ---------------------------------------------------------------
            # Verify that the install succeeded
            # ---------------------------------------------------------------
            # We call requireNamespace() again after installation. If it still
            # returns FALSE, the install silently failed (e.g., a missing system
            # library, network error, or compilation error).
            if (requireNamespace(pkg, quietly = TRUE)) {
                cat("  ✓", pkg, "installed successfully\n")
            } else {
                cat("  ✗ Failed to install", pkg)
                if (!is.null(install_error)) cat(": ", install_error, sep = "")
                cat("\n")
            }

        } else {
            # Package was already present; skip the download.
            cat("  ✓", pkg, "already installed\n")
        }
    }
}


# ==============================================================================
# Step 1 of 4 — Install BiocManager
# ==============================================================================
# BiocManager is the official gateway to Bioconductor packages. It must be
# present before we can call BiocManager::install() for dada2, phyloseq, etc.
# Installing it from CRAN first is safe even on systems that already have it;
# install.packages() will silently do nothing if the package is up to date.
cat("\n[1/4] Checking BiocManager...\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    cat("Installing BiocManager...\n")
    tryCatch(
        install.packages("BiocManager", quiet = TRUE),
        error = function(e) cat("  ✗ BiocManager installation failed: ", conditionMessage(e), "\n", sep = "")
    )
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    stop("BiocManager is unavailable; Bioconductor packages cannot be installed.", call. = FALSE)
}
cat("  ✓ BiocManager available\n")


# ==============================================================================
# Step 2 of 4 — CRAN packages
# ==============================================================================
# These packages are all available on the Comprehensive R Archive Network (CRAN)
# and cover four functional areas of the workflow:
#
#   Data manipulation  : data.table, tidyverse, dplyr, tidyr, stringr
#   File system / paths: here, fs
#   Excel & reporting  : openxlsx, knitr, rmarkdown, htmltools
#   Visualization      : ggplot2, plotly, RColorBrewer, htmlwidgets, png, grid
#   Interactive tables : DT
#   Parallel processing: parallel, pbapply
#   FastQC interface   : fastqcr
#
# All packages are passed as a single character vector to `install_if_missing`
# so that the loop handles each one consistently.
cat("\n[2/4] Installing CRAN packages...\n")

cran_packages <- c(
    # ------------------------------------------------------------------
    # Data manipulation
    # ------------------------------------------------------------------
    # data.table  : High-performance alternative to data.frame; used
    #               throughout the pipeline for fast tabular operations.
    "data.table",

    # tidyverse   : Meta-package that loads ggplot2, dplyr, tidyr, readr,
    #               purrr, tibble, stringr, and forcats in one call.
    "tidyverse",

    # dplyr       : Grammar of data manipulation; provides filter(), mutate(),
    #               summarise(), group_by(), etc.
    "dplyr",

    # tidyr       : Reshaping functions (pivot_wider, pivot_longer, etc.) for
    #               converting between wide and long formats.
    "tidyr",

    # stringr     : Consistent, pipe-friendly string manipulation functions
    #               built on top of stringi.
    "stringr",

    # ------------------------------------------------------------------
    # File system and paths
    # ------------------------------------------------------------------
    # here        : Constructs paths relative to the project root using
    #               .here or .Rproj anchors; prevents hardcoded absolute paths.
    "here",

    # fs          : Cross-platform file-system operations (path manipulation,
    #               directory creation, file copying) using a tidy interface.
    "fs",

    # ------------------------------------------------------------------
    # Excel and reporting
    # ------------------------------------------------------------------
    # openxlsx    : Read/write .xlsx files without Java; supports custom
    #               styling, multiple sheets, and formula writing.
    "openxlsx",

    # knitr       : Dynamic report generation that embeds R code inside
    #               Markdown, LaTeX, or HTML documents.
    "knitr",

    # rmarkdown   : R Markdown rendering engine used to knit the workflow's
    #               notebooks to HTML and GitHub Markdown.
    "rmarkdown",

    # htmltools   : HTML tag/dependency utilities used directly when notebooks
    #               assemble multiple widgets and report sections.
    "htmltools",

    # ------------------------------------------------------------------
    # Visualization
    # ------------------------------------------------------------------
    # ggplot2     : Grammar-of-graphics plotting library; used for quality
    #               control figures, alpha/beta diversity plots, etc.
    "ggplot2",

    # plotly      : Converts ggplot2 objects to interactive HTML widgets or
    #               builds standalone interactive plots.
    "plotly",

    # RColorBrewer: Provides perceptually uniform colour palettes suitable
    #               for categorical and sequential data.
    "RColorBrewer",

    # htmlwidgets : Framework for embedding JavaScript visualizations (e.g.,
    #               plotly, DT) in R Markdown reports.
    "htmlwidgets",

    # png         : Read and write PNG raster images; used when inserting
    #               figures into Excel or HTML outputs.
    "png",

    # ------------------------------------------------------------------
    # Interactive tables
    # ------------------------------------------------------------------
    # DT          : R interface to the DataTables JavaScript library; renders
    #               searchable, sortable, paginated tables in HTML reports.
    "DT",

    # ------------------------------------------------------------------
    # Parallel processing
    # ------------------------------------------------------------------
    # pbapply     : Drop-in replacement for apply-family functions that adds
    #               a progress bar, helpful for long DADA2 steps.
    "pbapply",

    # ------------------------------------------------------------------
    # FastQC R interface
    # ------------------------------------------------------------------
    # fastqcr     : Reads and summarises FastQC result files directly in R,
    #               enabling programmatic QC reporting without manual HTML review.
    "fastqcr",

    # readr       : Used directly for delimited input/output in Steps 7-9.
    "readr"
)

install_if_missing(cran_packages, source = "CRAN")


# ==============================================================================
# Step 3 of 4 — Bioconductor packages
# ==============================================================================
# These five packages are the analytical core of the DADA2 amplicon workflow.
# They must be installed via BiocManager because they are part of the
# Bioconductor project, which maintains coordinated release cycles and has
# stricter quality requirements than CRAN.
#
#   dada2       : The primary amplicon denoising engine. Implements the DADA2
#                 algorithm for high-resolution inference of exact amplicon
#                 sequence variants (ASVs) from Illumina short-read data.
#                 Reference: Callahan et al. (2016) Nature Methods.
#
#   ShortRead   : Input/output and quality assessment for short-read sequencing
#                 data. Used to read FASTQ files and generate per-base quality
#                 score summaries.
#
#   Biostrings  : Efficient manipulation of biological sequences (DNA, RNA,
#                 amino acid). Used for reverse-complement operations and
#                 primer trimming logic.
#
#   DECIPHER    : Sequence alignment toolkit used to align ASV sequences before
#                 phylogenetic-tree construction in Step 6.
#
#   phyloseq    : Container and analysis framework for microbiome data that
#                 integrates OTU/ASV tables, taxonomy tables, sample metadata,
#                 phylogenetic trees, and reference sequences into a single object.
#                 Also used by the optional Step 8 (microbial load / QMP
#                 correction) notebook, via phyloseq::rarefy_even_depth().
#
#   Biostrings is also used by the optional Step 7 (16S copy-number
#   correction) notebook to convert Step 5's ASV sequences into a FASTA file
#   for PICRUSt2's place_seqs.py. PICRUSt2 itself is an external (non-R) tool
#   and is installed separately — see setup/install_picrust2.sh.
cat("\n[3/4] Installing Bioconductor packages...\n")

bioc_packages <- c(
    "dada2",       # Core DADA2 denoising algorithm
    "ShortRead",   # FASTQ I/O and quality assessment
    "Biostrings",  # Biological sequence manipulation
    "DECIPHER",    # Sequence alignment and taxonomy
    "phyloseq"     # Microbiome data container and analysis
)

install_if_missing(bioc_packages, source = "Bioconductor")

# ------------------------------------------------------------------------------
# Optional: phylogenetics packages
# ------------------------------------------------------------------------------
# ape is required by Steps 6 and 9 for tree I/O and phyloseq integration.
# phangorn is soft-guarded in Step 6 and only adds optional tree statistics.
cat("\nInstalling phylogenetics packages...\n")
required_phylo_cran <- "ape"
optional_phylo_cran <- "phangorn"
install_if_missing(required_phylo_cran, source = "CRAN")
install_if_missing(optional_phylo_cran, source = "CRAN")


# ==============================================================================
# Step 4 of 4 — Shiny app packages
# ==============================================================================
# These packages power the interactive Shiny front-end that allows users to
# configure DADA2 filter/trim parameters through a graphical interface rather
# than editing R scripts directly.
#
#   shiny       : The core reactive web-application framework.
#   shinyFiles  : Adds server-side file system browsing widgets to Shiny UIs,
#                 allowing users to select FASTQ directories interactively.
#   shinyjs     : JavaScript utilities for Shiny (show/hide elements, disable
#                 buttons, run custom JS from R).
#   bslib       : Bootstrap-based theming for Shiny and R Markdown, providing
#                 a modern, responsive layout system.
#   future      : Asynchronous and parallel execution backend; allows Shiny
#                 apps to run long computations without blocking the UI.
#   furrr       : Parallel version of purrr's map functions, built on future;
#                 used to distribute per-sample DADA2 processing across cores.
#   Rcpp        : Seamlessly integrates C++ code into R; required as a compiled
#                 dependency by several packages in this list (e.g., dada2).
cat("\n[4/4] Installing Shiny app packages...\n")

shiny_packages <- c(
    "shiny",       # Core reactive web-application framework
    "shinyFiles",  # Server-side file/directory browser widgets
    "shinyjs",     # JavaScript utilities for Shiny
    "bslib",       # Bootstrap theming for Shiny and R Markdown
    "future",      # Asynchronous and parallel execution backend
    "furrr",       # Parallel purrr::map functions (built on future)
    "Rcpp"         # R/C++ interface (compiled dependency for dada2 etc.)
)

install_if_missing(shiny_packages, source = "CRAN")


# ==============================================================================
# Installation summary
# ==============================================================================
# After all installation attempts have completed, we compile the full list of
# packages and use requireNamespace() to test each one. The results are counted
# and displayed in a clear pass/fail format so the user can immediately see
# whether any packages need manual attention.
cat("\n
╔══════════════════════════════════════════════════════════════════════════════╗
║                         INSTALLATION SUMMARY                                 ║
╚══════════════════════════════════════════════════════════════════════════════╝
\n")

# Keep required and optional packages separate so an unavailable optional tree
# visualization dependency is not reported as a failed core installation.
required_packages <- unique(c(cran_packages, bioc_packages, shiny_packages, required_phylo_cran))
optional_packages <- unique(optional_phylo_cran)
all_packages <- c(required_packages, optional_packages)

# sapply returns a named logical vector: TRUE if the package loaded, FALSE if not.
installed <- sapply(all_packages, requireNamespace, quietly = TRUE)

cat("Required packages available:", sum(installed[required_packages]), "/", length(required_packages), "\n")
cat("Optional phylogenetics packages available:", sum(installed[optional_packages]), "/", length(optional_packages), "\n\n")

if (all(installed[required_packages])) {
    # All required packages are present; the core environment is ready.
    cat("✓ All required packages installed successfully!\n")
    cat("\nYou can now run the pipeline notebooks.\n")
} else {
    # At least one required package is missing. List each failed package so the user
    # knows exactly what to investigate or install manually.
    missing_required <- required_packages[!installed[required_packages]]
    cat("✗ Some required packages failed to install:\n")
    cat(paste("  -", missing_required, collapse = "\n"), "\n")
    cat("\nThe installer will not stop automatically. Install the packages listed above manually before running the affected notebooks.\n")
}

missing_optional <- optional_packages[!installed[optional_packages]]
if (length(missing_optional) > 0L) {
    cat("\nOptional phylogenetics packages not available:\n")
    cat(paste("  -", missing_optional, collapse = "\n"), "\n")
}


# ==============================================================================
# Installed version report
# ==============================================================================
# For reproducibility and collaborative troubleshooting it is important to know
# exactly which versions are active. We print version strings for the five key
# Bioconductor packages, plus the R and Bioconductor release versions.
# This output can be copied directly into a Methods section or a lab notebook.
cat("\n
══════════════════════════════════════════════════════════════════════════════
                          INSTALLED VERSIONS
══════════════════════════════════════════════════════════════════════════════
\n")

# Check and report version for each key package.
# packageVersion() returns a package_version object; we coerce it to character
# for clean printing.
key_packages <- c("dada2", "phyloseq", "ShortRead", "DECIPHER", "Biostrings")
for (pkg in key_packages) {
    if (requireNamespace(pkg, quietly = TRUE)) {
        ver <- as.character(packageVersion(pkg))
        # sprintf with %-15s left-aligns the package name in a 15-character field
        # so all version numbers line up in a tidy column.
        cat(sprintf("%-15s %s\n", pkg, ver))
    }
}

# Report the R and Bioconductor release versions.
# These are essential context for anyone trying to reproduce the analysis on a
# different machine or at a later date.
cat("\nR version:", R.version.string, "\n")
cat("Bioconductor version:", as.character(BiocManager::version()), "\n")
