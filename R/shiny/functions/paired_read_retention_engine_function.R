# =============================================================================
# Script:  paired_read_retention_engine_function.R
# Purpose: Computational helpers used by dada2_parameter_selection_app.R for
#          paired FASTQ detection, sampling, retention estimation, and validation.
#
#          The filtering and summary functions use base R only, which keeps the
#          calculations independently testable. extract_paired_quality_profile()
#          is the sole FASTQ-content I/O function and requires ShortRead at runtime.
#
# Conventions:
#   - Quality matrices are numeric [n_reads x n_positions] Phred scores, padded
#     with NA beyond each read's length. NA positions are absent bases.
#   - Expected error (EE) for a base of Phred quality Q is 10^(-Q/10).
#   - truncQ defaults to 2 and remains fixed in the app.
#   - minLen defaults to 20, matching DADA2.
#
# The simulator follows DADA2's filtering order: truncate to truncLen, truncate
# before the first base with Q <= truncQ, reject reads shorter than truncLen,
# reject effective lengths below minLen, and require EE <= maxEE. It is a fast
# estimator rather than a bit-identical replacement for dada2::filterAndTrim();
# the Validate tab runs the real DADA2 pipeline on representative samples.
# =============================================================================


# =============================================================================
# SECTION 1: ENGINE-LEVEL CONSTANTS
# =============================================================================

# Absolute minimum forward/reverse overlap accepted by DADA2 mergePairs().
DADA2_ABSOLUTE_MIN_OVERLAP <- 12

# Default minimum retained effective read length; matches DADA2 minLen.
RETENTION_DEFAULT_MIN_LEN <- 20

# Fixed truncQ used by the app and simulation helpers.
RETENTION_DEFAULT_TRUNCQ <- 2

# FASTQ extensions accepted by the shared pairing helpers.
PAIRED_FASTQ_EXTENSION_PATTERN <- "(?i)\\.(fastq|fq)(\\.gz)?$"


# =============================================================================
# SECTION 2: PAIRED FASTQ FILE DETECTION
# =============================================================================

# Pair forward and reverse FASTQ files by their complete filename stem. The
# direction-specific suffix and FASTQ extension are removed first; files are
# paired only when those remaining stems match exactly. Sample IDs are then
# derived from the text before the first underscore, as required by the
# workflow, and must be unique so named result objects and temporary output
# files cannot collide downstream.
detect_paired_files <- function(path,
                                fwd_pattern,
                                rev_pattern,
                                extension_pattern = PAIRED_FASTQ_EXTENSION_PATTERN) {
  empty_result <- function(message) {
    list(
      fnFs = character(0),
      fnRs = character(0),
      sample_names = character(0),
      file_stems = character(0),
      error = message
    )
  }

  if (!is.character(path) || length(path) != 1L || is.na(path) || !nzchar(trimws(path))) {
    return(empty_result("The selected FASTQ directory path is invalid"))
  }
  if (!is.character(fwd_pattern) || length(fwd_pattern) != 1L ||
      is.na(fwd_pattern) || !nzchar(fwd_pattern) ||
      !is.character(rev_pattern) || length(rev_pattern) != 1L ||
      is.na(rev_pattern) || !nzchar(rev_pattern)) {
    return(empty_result("Forward and reverse filename patterns must be non-empty character values"))
  }
  if (!is.character(extension_pattern) || length(extension_pattern) != 1L ||
      is.na(extension_pattern) || !nzchar(extension_pattern)) {
    return(empty_result("The FASTQ extension pattern is invalid"))
  }

  if (!dir.exists(path)) {
    return(empty_result("The selected FASTQ directory does not exist"))
  }

  all_files <- list.files(
    path,
    pattern = extension_pattern,
    full.names = TRUE
  )
  if (length(all_files) == 0L) {
    return(empty_result("No FASTQ files found in directory"))
  }

  escape_regex <- function(value) {
    gsub("([.^$*+?{}\\[\\]\\\\|()])", "\\\\\\1", value)
  }
  fwd_regex <- paste0(escape_regex(fwd_pattern), extension_pattern)
  rev_regex <- paste0(escape_regex(rev_pattern), extension_pattern)

  fnFs <- all_files[grepl(fwd_regex, basename(all_files))]
  fnRs <- all_files[grepl(rev_regex, basename(all_files))]
  if (length(fnFs) == 0L) {
    return(empty_result(paste0("No forward files found matching pattern '", fwd_pattern, "'")))
  }
  if (length(fnRs) == 0L) {
    return(empty_result(paste0("No reverse files found matching pattern '", rev_pattern, "'")))
  }

  fwd_stems <- sub(fwd_regex, "", basename(fnFs))
  rev_stems <- sub(rev_regex, "", basename(fnRs))

  duplicate_fwd <- unique(fwd_stems[duplicated(fwd_stems)])
  duplicate_rev <- unique(rev_stems[duplicated(rev_stems)])
  if (length(duplicate_fwd) > 0L || length(duplicate_rev) > 0L) {
    duplicate_text <- unique(c(duplicate_fwd, duplicate_rev))
    return(empty_result(paste0(
      "Multiple FASTQ files resolve to the same paired stem: ",
      paste(duplicate_text, collapse = ", ")
    )))
  }

  missing_reverse <- setdiff(fwd_stems, rev_stems)
  missing_forward <- setdiff(rev_stems, fwd_stems)
  if (length(missing_reverse) > 0L || length(missing_forward) > 0L) {
    details <- c(
      if (length(missing_reverse) > 0L) {
        paste0("missing reverse mate(s): ", paste(missing_reverse, collapse = ", "))
      },
      if (length(missing_forward) > 0L) {
        paste0("missing forward mate(s): ", paste(missing_forward, collapse = ", "))
      }
    )
    return(empty_result(paste0(
      "Forward and reverse files do not pair by complete filename stem (",
      paste(details, collapse = "; "), ")"
    )))
  }

  ordered_stems <- sort(fwd_stems)
  fnFs <- fnFs[match(ordered_stems, fwd_stems)]
  fnRs <- fnRs[match(ordered_stems, rev_stems)]
  sample_names <- sub("_.*$", "", ordered_stems)

  if (any(!nzchar(sample_names))) {
    return(empty_result("At least one paired FASTQ stem has an empty sample ID"))
  }
  duplicated_samples <- unique(sample_names[duplicated(sample_names)])
  if (length(duplicated_samples) > 0L) {
    return(empty_result(paste0(
      "Sample IDs must be unique before the first underscore; duplicate ID(s): ",
      paste(duplicated_samples, collapse = ", ")
    )))
  }

  list(
    fnFs = fnFs,
    fnRs = fnRs,
    sample_names = unname(sample_names),
    file_stems = unname(ordered_stems),
    error = NULL
  )
}


# Return an exact paired-sample count for one selected filename convention.
# Invalid, incomplete, or duplicate pair sets report zero rather than claiming
# that an empty directory contains one sample.
estimate_paired_sample_count <- function(path, fwd_pattern, rev_pattern) {
  paired <- detect_paired_files(path, fwd_pattern, rev_pattern)
  if (!is.null(paired$error)) 0L else length(paired$fnFs)
}


# =============================================================================
# SECTION 3: SINGLE-DIRECTION DADA2-LIKE FILTERING SIMULATION
# =============================================================================

# Simulate DADA2-style filtering of one read direction from a quality matrix.
#
# Arguments:
#   quality_matrix : numeric [n_reads x n_positions] Phred scores, NA-padded.
#   trunc_len      : truncation length (bp) for this direction.
#   max_ee         : maxEE threshold for this direction.
#   trunc_q        : truncQ (truncate at first base with Q <= trunc_q).
#   min_len        : minimum retained effective length (bp).
#
# Returns a list with per-read vectors (length n_reads):
#   original_length  : each read's own (non-NA) length.
#   effective_length : retained length after truncLen + truncQ (0 if discarded
#                      for being shorter than trunc_len).
#   expected_errors  : summed expected errors over the retained bases.
#   pass             : logical, TRUE if the read passes this direction's filter.
simulate_read_filtering <- function(quality_matrix,
                                    trunc_len,
                                    max_ee,
                                    trunc_q = RETENTION_DEFAULT_TRUNCQ,
                                    min_len = RETENTION_DEFAULT_MIN_LEN) {

  if (!is.numeric(quality_matrix) || any(!is.finite(quality_matrix[!is.na(quality_matrix)])) ||
      any(quality_matrix[!is.na(quality_matrix)] < 0)) {
    stop("simulate_read_filtering(): `quality_matrix` must contain non-negative numeric Phred scores or NA.")
  }
  scalar_nonnegative <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) && x >= 0
  scalar_positive_integer <- function(x) scalar_nonnegative(x) && x >= 1 && x == as.integer(x)
  if (!scalar_positive_integer(trunc_len)) {
    stop("simulate_read_filtering(): `trunc_len` must be one positive integer.")
  }
  if (!scalar_nonnegative(max_ee)) {
    stop("simulate_read_filtering(): `max_ee` must be one finite, non-negative number.")
  }
  if (!scalar_nonnegative(trunc_q)) {
    stop("simulate_read_filtering(): `trunc_q` must be one finite, non-negative number.")
  }
  if (!scalar_positive_integer(min_len)) {
    stop("simulate_read_filtering(): `min_len` must be one positive integer.")
  }

  # Coerce a single-read vector to a 1-row matrix so the code below is uniform.
  if (is.null(dim(quality_matrix))) {
    quality_matrix <- matrix(quality_matrix, nrow = 1)
  }

  n_reads       <- nrow(quality_matrix)
  n_positions   <- ncol(quality_matrix)

  if (n_reads == 0L) {
    return(list(original_length = integer(0), effective_length = integer(0),
                expected_errors = numeric(0), pass = logical(0)))
  }

  # Each read's own length is the count of non-NA positions (NA is "no base",
  # never a Q40 base).
  original_length <- rowSums(!is.na(quality_matrix))

  # Reads whose ORIGINAL length is shorter than the requested trunc_len are
  # discarded outright (DADA2 removes reads shorter than truncLen).
  long_enough <- original_length >= trunc_len

  # Nothing can pass if trunc_len exceeds the matrix width (no read is that long).
  if (trunc_len > n_positions || trunc_len < 1) {
    return(list(
      original_length  = original_length,
      effective_length = rep(0, n_reads),
      expected_errors  = rep(NA_real_, n_reads),
      pass             = rep(FALSE, n_reads)
    ))
  }

  # Work within the truncation window [1 .. trunc_len].
  window <- quality_matrix[, seq_len(trunc_len), drop = FALSE]

  # truncQ: locate the FIRST position (within the window) whose quality is
  # <= trunc_q. max.col on the 0/1 low-quality mask returns the first such
  # column when ties.method = "first"; rows with no low-quality base are handled
  # via has_low_quality below.
  low_quality_mask <- window <= trunc_q
  low_quality_mask[is.na(low_quality_mask)] <- FALSE      # NA positions are not "low quality"
  has_low_quality  <- rowSums(low_quality_mask) > 0
  first_low_col    <- max.col(low_quality_mask, ties.method = "first")

  # Effective length = bases kept before the first low-quality base, capped at
  # trunc_len. Reads with no low-quality base in-window keep the full trunc_len.
  effective_length <- ifelse(has_low_quality, first_low_col - 1L, trunc_len)
  # A read discarded for being too short contributes no effective length.
  effective_length[!long_enough] <- 0L

  # Expected errors per position (NA -> 0 so they contribute nothing and never
  # propagate NA through the cumulative sum).
  error_prob <- 10^(-window / 10)
  error_prob[is.na(error_prob)] <- 0
  cumulative_error <- t(apply(error_prob, 1, cumsum))
  # Guard the 1-row case: apply() drops to a vector and t() would mis-shape it.
  if (n_reads == 1) {
    cumulative_error <- matrix(cumsum(error_prob[1, ]), nrow = 1)
  }

  # EE summed over the retained bases [1 .. effective_length].
  clamped_len     <- pmax(effective_length, 1L)
  expected_errors <- cumulative_error[cbind(seq_len(n_reads), clamped_len)]
  expected_errors[effective_length == 0] <- 0

  # Final per-direction pass: long enough, effective length >= min_len, and the
  # expected-error budget is met.
  pass <- long_enough &
    (effective_length >= min_len) &
    (expected_errors <= max_ee)

  list(
    original_length  = original_length,
    effective_length = effective_length,
    expected_errors  = expected_errors,
    pass             = pass
  )
}


# =============================================================================
# SECTION 4: PAIRED RETENTION SIMULATION (one sample)
# =============================================================================

# Evaluate one parameter set against a sample's row-aligned forward and reverse
# quality matrices.
#
# quality_matrix_fwd / quality_matrix_rev MUST be row-aligned: row k of each is
# the forward / reverse read of the same physical pair. amplicon_p99 is the
# conservative maximum expected amplicon length from the app; min_overlap is
# the feasibility floor used for the mergeable-fraction estimate.
simulate_paired_retention <- function(quality_matrix_fwd,
                                      quality_matrix_rev,
                                      trunc_len_fwd,
                                      trunc_len_rev,
                                      max_ee_fwd,
                                      max_ee_rev,
                                      trunc_q      = RETENTION_DEFAULT_TRUNCQ,
                                      min_len      = RETENTION_DEFAULT_MIN_LEN,
                                      amplicon_p99 = NA_real_,
                                      min_overlap  = DADA2_ABSOLUTE_MIN_OVERLAP) {

  if (!is.numeric(amplicon_p99) || length(amplicon_p99) != 1L ||
      (!is.na(amplicon_p99) && (!is.finite(amplicon_p99) || amplicon_p99 <= 0))) {
    stop("simulate_paired_retention(): `amplicon_p99` must be NA or one positive finite number.")
  }
  if (!is.numeric(min_overlap) || length(min_overlap) != 1L || is.na(min_overlap) ||
      !is.finite(min_overlap) || min_overlap < 0) {
    stop("simulate_paired_retention(): `min_overlap` must be one finite, non-negative number.")
  }

  n_fwd <- if (is.null(dim(quality_matrix_fwd))) 1L else nrow(quality_matrix_fwd)
  n_rev <- if (is.null(dim(quality_matrix_rev))) 1L else nrow(quality_matrix_rev)
  if (n_fwd != n_rev) {
    stop("simulate_paired_retention(): forward and reverse matrices must contain the same number of reads.")
  }

  fwd <- simulate_read_filtering(quality_matrix_fwd, trunc_len_fwd, max_ee_fwd, trunc_q, min_len)
  rev <- simulate_read_filtering(quality_matrix_rev, trunc_len_rev, max_ee_rev, trunc_q, min_len)

  # A pair is retained only when both directions pass.
  paired_pass      <- fwd$pass & rev$pass
  n_pairs          <- length(paired_pass)
  paired_retention <- if (n_pairs > 0) mean(paired_pass) else NA_real_

  # Per-pair realized overlap uses each read's EFFECTIVE length (which varies
  # once truncQ trims low-quality tails), measured against the conservative
  # (longest) expected amplicon length p99. Used only for the mergeable
  # fraction, and only over pairs that were actually retained.
  retained_idx <- which(paired_pass)
  if (length(retained_idx) > 0 && is.finite(amplicon_p99)) {
    per_pair_overlap  <- fwd$effective_length[retained_idx] + rev$effective_length[retained_idx] - amplicon_p99
    mergeable_fraction <- mean(per_pair_overlap >= min_overlap)
  } else {
    mergeable_fraction <- NA_real_
  }

  list(
    paired_retention   = paired_retention,
    n_pairs            = n_pairs,
    n_retained         = length(retained_idx),
    mergeable_fraction = mergeable_fraction
  )
}


# =============================================================================
# SECTION 5: CROSS-SAMPLE SUMMARY
# =============================================================================

# Combine per-sample paired-retention results into cohort-level estimates.
# n_sampled contains the sampled-pair count for each sample.
summarize_paired_retention <- function(per_sample_results, n_sampled) {

  if (!is.list(per_sample_results) || !is.numeric(n_sampled) || anyNA(n_sampled) ||
      any(!is.finite(n_sampled)) || any(n_sampled < 0)) {
    stop("summarize_paired_retention(): supply a result list and finite, non-negative sampling counts.")
  }

  if (length(per_sample_results) != length(n_sampled)) {
    stop("summarize_paired_retention(): n_sampled must contain one value per sample.")
  }

  retentions <- vapply(per_sample_results, function(x) x$paired_retention, numeric(1))
  weights    <- as.numeric(n_sampled)

  # Only samples with a defined retention contribute.
  valid <- is.finite(retentions) & is.finite(weights) & weights > 0
  retentions_v <- retentions[valid]
  weights_v    <- weights[valid]

  if (length(retentions_v) == 0) {
    return(list(
      read_weighted_retention = NA_real_,
      mergeable_fraction      = NA_real_,
      n_samples               = 0L
    ))
  }

  # Read-weighted retention equals total retained pairs / total sampled pairs.
  read_weighted_retention <- sum(retentions_v * weights_v) / sum(weights_v)

  # Mergeability is defined among retained pairs, so samples are weighted by
  # their retained-pair counts rather than their original sampling depths.
  mergeable <- vapply(per_sample_results, function(x) x$mergeable_fraction, numeric(1))
  n_retained <- vapply(per_sample_results, function(x) x$n_retained, numeric(1))
  merge_ok <- is.finite(mergeable) & is.finite(n_retained) & n_retained > 0
  mergeable_fraction <- if (any(merge_ok)) {
    sum(mergeable[merge_ok] * n_retained[merge_ok]) / sum(n_retained[merge_ok])
  } else {
    NA_real_
  }

  list(
    read_weighted_retention = read_weighted_retention,
    mergeable_fraction      = mergeable_fraction,
    n_samples               = length(retentions_v)
  )
}




# =============================================================================
# SECTION 6: SAMPLING DEPTH GUIDANCE
# =============================================================================

# Adaptive automatic policy for how many paired reads to sample per sample,
# balancing statistical precision against runtime and memory.
recommended_pairs_per_sample <- function(n_samples) {
  if (!is.numeric(n_samples) || length(n_samples) != 1L || is.na(n_samples)) return(15000L)
  if (!is.finite(n_samples) || n_samples <= 0) return(15000L)
  if (n_samples <= 20)  return(20000L)
  if (n_samples <= 50)  return(15000L)
  if (n_samples <= 100) return(9000L)
  # More than 100 samples: keep per-sample sampling modest to bound total
  # runtime/memory (precision is already high with many samples).
  4000L
}


# =============================================================================
# SECTION 7: PAIRED FASTQ SAMPLER (I/O)
# =============================================================================

# This is the only function in this file that reads FASTQ contents or touches
# Bioconductor. It
# streams a forward/reverse FASTQ pair in lockstep and reservoir-samples matched
# read PAIRS across the WHOLE file (not just the head), so row k of the returned
# forward and reverse quality matrices is the same physical pair. Reproducible
# via `seed`. Requires ShortRead (FastqStreamer/quality) at runtime.
#
# It is defined here so sampling and retention calculations share one helper
# file. The regression tests exercise both this I/O path and the pure
# calculation functions above.

extract_paired_quality_profile <- function(fwd_file,
                                           rev_file,
                                           n_pairs    = 15000,
                                           seed       = 1L,
                                           chunk_size = 50000,
                                           validate_ids = TRUE) {

  valid_file <- function(x) is.character(x) && length(x) == 1L && !is.na(x) && nzchar(trimws(x)) && file.exists(x)
  positive_integer <- function(x) is.numeric(x) && length(x) == 1L && !is.na(x) &&
    is.finite(x) && x >= 1 && x == as.integer(x)
  if (!valid_file(fwd_file) || !valid_file(rev_file)) {
    stop("extract_paired_quality_profile(): forward and reverse FASTQ files must both exist.")
  }
  if (!positive_integer(n_pairs) || !positive_integer(chunk_size)) {
    stop("extract_paired_quality_profile(): `n_pairs` and `chunk_size` must be positive integers.")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) || !is.finite(seed)) {
    stop("extract_paired_quality_profile(): `seed` must be one finite number.")
  }
  if (!is.logical(validate_ids) || length(validate_ids) != 1L || is.na(validate_ids)) {
    stop("extract_paired_quality_profile(): `validate_ids` must be TRUE or FALSE.")
  }

  if (!requireNamespace("ShortRead", quietly = TRUE)) {
    stop("extract_paired_quality_profile(): ShortRead is required but not installed.")
  }
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("extract_paired_quality_profile(): Biostrings is required but not installed.")
  }

  had_random_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_random_seed) previous_random_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_random_seed) {
      assign(".Random.seed", previous_random_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(as.integer(seed))

  # Namespace-qualified calls keep direct use and tests independent of whether
  # ShortRead has been attached to the caller's search path.
  fwd_stream <- ShortRead::FastqStreamer(fwd_file, n = chunk_size)
  rev_stream <- ShortRead::FastqStreamer(rev_file, n = chunk_size)
  on.exit({ close(fwd_stream); close(rev_stream) }, add = TRUE)

  # Reservoir of up to n_pairs sampled pairs; each slot holds one read's numeric
  # Phred vector (integer, NA-padded), matrixized once at the end.
  reservoir_fwd <- vector("list", n_pairs)
  reservoir_rev <- vector("list", n_pairs)
  filled       <- 0L      # how many reservoir slots are populated
  seen         <- 0L      # how many pairs streamed so far (for reservoir math)
  chunk_index  <- 0L

  repeat {
    fwd_chunk <- ShortRead::yield(fwd_stream)
    rev_chunk <- ShortRead::yield(rev_stream)
    chunk_index <- chunk_index + 1L
    n_fwd <- length(fwd_chunk); n_rev <- length(rev_chunk)
    if (n_fwd == 0 && n_rev == 0) break
    # Forward and reverse must stay in lockstep; unequal chunk sizes mean the
    # files are not properly paired.
    if (n_fwd != n_rev) {
      stop("extract_paired_quality_profile(): forward/reverse read counts differ ",
           "in a chunk (", n_fwd, " vs ", n_rev, ") -- files are not properly paired.")
    }
    n_chunk <- n_fwd

    # Validate every chunk. Mate suffixes /1 /2 and the Casava mate field are
    # stripped before comparison; a mismatch invalidates row-aligned pairing
    # and therefore stops processing rather than producing an estimate.
    if (validate_ids) {
      strip_mate <- function(ids) sub("[ /][12].*$", "", ids)
      fwd_ids <- strip_mate(as.character(ShortRead::id(fwd_chunk)))
      rev_ids <- strip_mate(as.character(ShortRead::id(rev_chunk)))
      if (!identical(fwd_ids, rev_ids)) {
        first_mismatch <- which(fwd_ids != rev_ids)[1]
        stop(
          "extract_paired_quality_profile(): forward/reverse read IDs do not match ",
          "in chunk ", chunk_index, " at row ", first_mismatch,
          " -- files are not properly paired."
        )
      }
    }

    # Numeric Phred matrices for the whole chunk (rows = reads, NA-padded to the
    # chunk's longest read). Row k of qf and qr is the same physical read pair.
    qf <- as(Biostrings::quality(fwd_chunk), "matrix")
    qr <- as(Biostrings::quality(rev_chunk), "matrix")

    # Reservoir sampling over pairs: fill first, then replace with decreasing
    # probability so every pair in the file has equal inclusion odds.
    for (k in seq_len(n_chunk)) {
      seen <- seen + 1L
      if (filled < n_pairs) {
        filled <- filled + 1L
        reservoir_fwd[[filled]] <- qf[k, ]
        reservoir_rev[[filled]] <- qr[k, ]
      } else {
        j <- sample.int(seen, 1L)
        if (j <= n_pairs) {
          reservoir_fwd[[j]] <- qf[k, ]
          reservoir_rev[[j]] <- qr[k, ]
        }
      }
    }
  }

  if (filled == 0L) return(NULL)
  reservoir_fwd <- reservoir_fwd[seq_len(filled)]
  reservoir_rev <- reservoir_rev[seq_len(filled)]

  # Combine the sampled per-read Phred vectors into one integer matrix, padding
  # every row to the longest sampled read (trailing NA = no base). Integer
  # storage roughly halves memory vs. double, with no precision loss since Phred
  # scores are whole numbers.
  rows_to_matrix <- function(rows) {
    max_len <- max(vapply(rows, length, integer(1)))
    mat <- matrix(NA_integer_, nrow = length(rows), ncol = max_len)
    for (i in seq_along(rows)) {
      v <- as.integer(rows[[i]])
      if (length(v) > 0) mat[i, seq_along(v)] <- v
    }
    mat
  }

  quality_matrix_fwd <- rows_to_matrix(reservoir_fwd)
  quality_matrix_rev <- rows_to_matrix(reservoir_rev)

  list(
    quality_matrix_fwd = quality_matrix_fwd,
    quality_matrix_rev = quality_matrix_rev,
    n_sampled          = filled,
    n_total            = seen,
    max_len_fwd        = ncol(quality_matrix_fwd),
    max_len_rev        = ncol(quality_matrix_rev),
    seed               = seed
  )
}
