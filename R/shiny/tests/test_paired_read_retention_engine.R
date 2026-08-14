# =============================================================================
# Script:  test_paired_read_retention_engine.R
# Purpose: Regression tests for paired-file detection, paired-read sampling,
#          and retention calculations used by dada2_parameter_selection_app.R.
#
# Run from the repository root:
#   Rscript R/shiny/tests/test_paired_read_retention_engine.R
# =============================================================================

# Locate the engine when invoked from the repository root or this test folder.
engine_candidates <- c(
  file.path("R", "shiny", "functions", "paired_read_retention_engine_function.R"),
  file.path("..", "functions", "paired_read_retention_engine_function.R"),
  "paired_read_retention_engine_function.R"
)
engine_path <- engine_candidates[file.exists(engine_candidates)][1]
if (is.na(engine_path)) stop("Could not locate paired_read_retention_engine_function.R")
source(engine_path)

.tests_run <- 0L
.tests_failed <- 0L

check <- function(condition, message) {
  .tests_run <<- .tests_run + 1L
  if (isTRUE(condition)) {
    cat(sprintf("  [ok]   %s\n", message))
  } else {
    .tests_failed <<- .tests_failed + 1L
    cat(sprintf("  [FAIL] %s\n", message))
  }
}

approx <- function(a, b, tol = 1e-8) {
  is.finite(a) && is.finite(b) && abs(a - b) <= tol
}

raises_error <- function(expr) {
  inherits(try(force(expr), silent = TRUE), "try-error")
}

qmat_const <- function(n, len, q) {
  matrix(q, nrow = n, ncol = len)
}

qread <- function(q_vec, pad_to = length(q_vec)) {
  out <- matrix(NA_real_, nrow = 1, ncol = pad_to)
  out[1, seq_along(q_vec)] <- q_vec
  out
}


cat("== Engine surface ==\n")
engine_functions <- ls(pattern = "^(detect_paired_files$|estimate_paired_sample_count$|recommended_pairs_per_sample$|simulate_read_filtering$|simulate_paired_retention$|summarize_paired_retention$|extract_paired_quality_profile$)")
expected_functions <- c(
  "detect_paired_files",
  "estimate_paired_sample_count",
  "extract_paired_quality_profile",
  "recommended_pairs_per_sample",
  "simulate_read_filtering",
  "simulate_paired_retention",
  "summarize_paired_retention"
)
check(setequal(engine_functions, expected_functions),
      "engine exposes only helpers used by the app")

cat("== Public argument validation ==\n")
check(raises_error(simulate_read_filtering(qmat_const(1, 10, 30), NA, 2)),
      "missing truncation lengths are rejected clearly")
check(raises_error(simulate_read_filtering(qmat_const(1, 10, 30), 10, -1)),
      "negative maxEE values are rejected clearly")
check(raises_error(simulate_read_filtering(matrix("Q", nrow = 1), 1, 2)),
      "non-numeric quality matrices are rejected clearly")
invalid_pattern <- detect_paired_files(tempdir(), NA_character_, "_R2")
check(!is.null(invalid_pattern$error),
      "invalid filename patterns return a structured detection error")


cat("== Paired FASTQ detection ==\n")
pairing_root <- tempfile("paired_fastq_detection_")
dir.create(pairing_root)

good_dir <- file.path(pairing_root, "good")
dir.create(good_dir)
invisible(file.create(file.path(good_dir, c(
  "Alpha_S1_L001_R1_001.fastq", "Alpha_S1_L001_R2_001.fastq",
  "Beta_S2_L001_R1_001.fastq", "Beta_S2_L001_R2_001.fastq"
))))
good_pairs <- detect_paired_files(good_dir, "_L001_R1_001", "_L001_R2_001")
check(is.null(good_pairs$error) && identical(good_pairs$sample_names, c("Alpha", "Beta")),
      "complete filename stems pair in deterministic order")
check(identical(good_pairs$file_stems, c("Alpha_S1", "Beta_S2")),
      "full paired stems are retained for validation")

mismatch_dir <- file.path(pairing_root, "mismatch")
dir.create(mismatch_dir)
invisible(file.create(file.path(mismatch_dir, c(
  "F3D0_A_R1.fastq", "F3D0_A_R2.fastq",
  "F3D0_B_R1.fastq", "F3D0_C_R2.fastq"
))))
mismatched_pairs <- detect_paired_files(mismatch_dir, "_R1", "_R2")
check(!is.null(mismatched_pairs$error) && grepl("complete filename stem", mismatched_pairs$error),
      "shared short sample IDs cannot hide mismatched file stems")

duplicate_dir <- file.path(pairing_root, "duplicate_sample")
dir.create(duplicate_dir)
invisible(file.create(file.path(duplicate_dir, c(
  "F3D0_A_R1.fastq", "F3D0_A_R2.fastq",
  "F3D0_B_R1.fastq", "F3D0_B_R2.fastq"
))))
duplicate_pairs <- detect_paired_files(duplicate_dir, "_R1", "_R2")
check(!is.null(duplicate_pairs$error) && grepl("duplicate ID", duplicate_pairs$error),
      "duplicate IDs before the first underscore are rejected")

empty_dir <- file.path(pairing_root, "empty")
dir.create(empty_dir)
check(estimate_paired_sample_count(empty_dir, "_R1", "_R2") == 0L,
      "an empty folder reports zero paired samples")

unlink(pairing_root, recursive = TRUE, force = TRUE)


cat("== Single-direction filtering ==\n")
q20 <- simulate_read_filtering(
  qmat_const(3, 100, 20),
  trunc_len = 100,
  max_ee = 2,
  trunc_q = 0,
  min_len = 20
)
check(approx(q20$expected_errors[1], 1.0),
      "100 Q20 bases have one expected error")
check(all(q20$pass), "Q20 reads pass maxEE 2")
check(all(q20$effective_length == 100),
      "reads without a truncQ hit retain truncLen bases")

q20_strict <- simulate_read_filtering(
  qmat_const(3, 100, 20),
  trunc_len = 100,
  max_ee = 0.5,
  trunc_q = 0
)
check(all(!q20_strict$pass), "Q20 reads fail maxEE 0.5")

short <- simulate_read_filtering(
  qread(rep(35, 50), pad_to = 120),
  trunc_len = 100,
  max_ee = 2
)
check(isTRUE(!short$pass[1]), "reads shorter than truncLen are discarded")
check(short$effective_length[1] == 0,
      "discarded short reads have zero effective length")

q_early_low <- rep(35, 100)
q_early_low[5] <- 2
early <- simulate_read_filtering(
  qread(q_early_low),
  trunc_len = 100,
  max_ee = 5,
  trunc_q = 2,
  min_len = 20
)
check(early$effective_length[1] == 4,
      "truncQ removes the triggering base and all following bases")
check(isTRUE(!early$pass[1]),
      "effective lengths below minLen fail")

q_late_low <- rep(35, 100)
q_late_low[90] <- 2
late <- simulate_read_filtering(
  qread(q_late_low),
  trunc_len = 100,
  max_ee = 5,
  trunc_q = 2,
  min_len = 20
)
check(late$effective_length[1] == 89,
      "a late truncQ hit preserves the preceding bases")
check(isTRUE(late$pass[1]), "a retained read above minLen can pass")

na_padded <- simulate_read_filtering(
  qread(rep(20, 30), pad_to = 100),
  trunc_len = 30,
  max_ee = 2,
  trunc_q = 0
)
check(approx(na_padded$expected_errors[1], 0.30),
      "NA padding does not contribute expected errors")


cat("== Paired retention and overlap ==\n")
good <- qmat_const(10, 120, 35)
bad <- qmat_const(10, 120, 5)

poor_reverse <- simulate_paired_retention(
  good, bad, 100, 100, 2, 2,
  trunc_q = 2,
  amplicon_p99 = 180
)
check(approx(poor_reverse$paired_retention, 0),
      "a pair fails when its reverse read fails")
check(poor_reverse$n_retained == 0,
      "no failed pairs are counted as retained")

both_good <- simulate_paired_retention(
  good, good, 100, 100, 2, 2,
  trunc_q = 2,
  amplicon_p99 = 180,
  min_overlap = 12
)
check(approx(both_good$paired_retention, 1),
      "pairs pass when both directions pass")
check(both_good$n_pairs == 10 && both_good$n_retained == 10,
      "paired counts are reported explicitly")
check(approx(both_good$mergeable_fraction, 1),
      "20 bp effective overlap clears a 12 bp minimum")

insufficient_overlap <- simulate_paired_retention(
  good, good, 100, 100, 2, 2,
  trunc_q = 2,
  amplicon_p99 = 190,
  min_overlap = 12
)
check(approx(insufficient_overlap$mergeable_fraction, 0),
      "10 bp effective overlap fails a 12 bp minimum")

correlated_fwd <- rbind(qmat_const(5, 120, 35), qmat_const(5, 120, 5))
correlated_rev <- rbind(qmat_const(5, 120, 35), qmat_const(5, 120, 5))
correlated <- simulate_paired_retention(
  correlated_fwd, correlated_rev, 100, 100, 2, 2,
  trunc_q = 2
)
check(approx(correlated$paired_retention, 0.5),
      "paired retention uses row-aligned pair outcomes")

check(
  raises_error(simulate_paired_retention(
    qmat_const(2, 100, 35),
    qmat_const(3, 100, 35),
    90, 90, 2, 2
  )),
  "mismatched forward and reverse row counts fail clearly"
)


cat("== Cross-sample summaries ==\n")
per_sample <- list(
  list(
    paired_retention = 0.50,
    n_retained = 50,
    mergeable_fraction = 0.20
  ),
  list(
    paired_retention = 1.00,
    n_retained = 900,
    mergeable_fraction = 0.80
  )
)
summary_result <- summarize_paired_retention(
  per_sample,
  n_sampled = c(100, 900)
)
check(approx(summary_result$read_weighted_retention, 0.95),
      "retention is weighted by sampled-pair counts")
check(
  approx(summary_result$mergeable_fraction, (0.20 * 50 + 0.80 * 900) / 950),
  "mergeability is weighted by retained-pair counts"
)
check(summary_result$n_samples == 2,
      "summary reports the number of contributing samples")

empty_summary <- summarize_paired_retention(list(), numeric())
check(is.na(empty_summary$read_weighted_retention) &&
        is.na(empty_summary$mergeable_fraction) &&
        empty_summary$n_samples == 0,
      "an empty cohort returns an explicit empty summary")

check(
  raises_error(summarize_paired_retention(per_sample, n_sampled = 100)),
  "sample results and sampling counts must have matching lengths"
)


cat("== Sampling-depth guidance ==\n")
check(recommended_pairs_per_sample(0) == 15000L,
      "invalid or empty cohorts use the fallback depth")
check(recommended_pairs_per_sample(20) == 20000L,
      "up to 20 samples use 20,000 pairs per sample")
check(recommended_pairs_per_sample(21) == 15000L,
      "21-50 samples use 15,000 pairs per sample")
check(recommended_pairs_per_sample(51) == 9000L,
      "51-100 samples use 9,000 pairs per sample")
check(recommended_pairs_per_sample(101) == 4000L,
      "larger cohorts use 4,000 pairs per sample")


cat("== Paired FASTQ streaming ==\n")
if (requireNamespace("ShortRead", quietly = TRUE) &&
    requireNamespace("Biostrings", quietly = TRUE)) {
  write_test_fastq <- function(path, ids, mate) {
    records <- unlist(lapply(ids, function(read_id) {
      c(paste0("@", read_id, "/", mate), "ACGTACGT", "+", "IIIIIIII")
    }), use.names = FALSE)
    writeLines(records, path, useBytes = TRUE)
  }

  io_root <- tempfile("paired_fastq_io_")
  dir.create(io_root)
  fwd_path <- file.path(io_root, "sample_R1.fastq")
  rev_path <- file.path(io_root, "sample_R2.fastq")
  ids <- paste0("read", 1:4)
  write_test_fastq(fwd_path, ids, 1L)
  write_test_fastq(rev_path, ids, 2L)

  sampled <- extract_paired_quality_profile(
    fwd_path, rev_path, n_pairs = 4L, chunk_size = 2L
  )
  check(sampled$n_total == 4L && sampled$n_sampled == 4L,
        "the sampler works without attaching ShortRead to the search path")
  check(identical(dim(sampled$quality_matrix_fwd), c(4L, 8L)) &&
          identical(dim(sampled$quality_matrix_rev), c(4L, 8L)),
        "paired quality matrices remain row-aligned")

  set.seed(987L)
  expected_next_random <- runif(1)
  set.seed(987L)
  invisible(extract_paired_quality_profile(
    fwd_path, rev_path, n_pairs = 2L, seed = 44L, chunk_size = 2L
  ))
  check(identical(runif(1), expected_next_random),
        "the paired sampler restores the caller's RNG state")
  check(raises_error(extract_paired_quality_profile(
    fwd_path, rev_path, n_pairs = 0L, chunk_size = 2L
  )), "invalid sampling sizes are rejected clearly")

  mismatched_ids <- c("read1", "read2", "different3", "read4")
  write_test_fastq(rev_path, mismatched_ids, 2L)
  check(
    raises_error(extract_paired_quality_profile(
      fwd_path, rev_path, n_pairs = 4L, chunk_size = 2L
    )),
    "a read-ID mismatch after the first streamed chunk stops processing"
  )
  unlink(io_root, recursive = TRUE, force = TRUE)
} else {
  cat("  [skip] ShortRead/Biostrings unavailable; FASTQ streaming tests skipped\n")
}


cat(sprintf("\n%d checks run, %d failed.\n", .tests_run, .tests_failed))
if (.tests_failed > 0) {
  quit(status = 1)
}
cat("ALL TESTS PASSED\n")
