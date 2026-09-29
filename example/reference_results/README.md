# Reference results

These are curated outputs from `Rscript example/run_example.R --include-qmp`, generated on 29 September 2026 with R 4.5.2, DADA2 1.36.0, Cutadapt 5.2, PICRUSt2 2.6.3, and GTDB r220.

The `reports/` HTML files provide the easiest way to browse each executed step. The step folders contain compact analytical deliverables: integrity and primer-trimming workbooks; ASV, taxonomy, processing-summary, quality-profile, error-model, and amplicon-length outputs; PICRUSt2 copy-number predictions and the corrected abundance table; the synthetic-cell-count QMP table; and final GTDB phyloseq objects for raw, copy-number-corrected, and microbial-load-corrected abundances.

The cell-count inputs are synthetic test values, not measurements from the source study. They must not be used for biological interpretation. See [`../data/cell_count/README.md`](../data/cell_count/README.md) for their range, assumptions, mock-community placeholder, and supporting references.

Large reproducible intermediates are intentionally omitted: primer-trimmed and filtered FASTQs, checkpoints, PICRUSt2 placement files, console logs, and interactive barplot dependency folders. Rerunning the example writes those files to the ignored `example/run_results/` tree and does not change these reference results automatically.
