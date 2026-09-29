# Reference results

These are curated outputs from `Rscript example/run_example.R`, generated on 29 September 2026 with R 4.5.2, DADA2 1.36.0, Cutadapt 5.2, and GTDB r220.

The `reports/` HTML files provide the easiest way to browse each executed step. The step folders contain compact analytical deliverables: integrity and primer-trimming workbooks; ASV, taxonomy, processing-summary, quality-profile, error-model, and amplicon-length outputs; and the final GTDB phyloseq object and summary workbook.

Large reproducible intermediates are intentionally omitted: primer-trimmed and filtered FASTQs, DADA2 checkpoints, and console logs. Rerunning the example writes those files to the ignored `example/run_results/` tree and does not change these reference results automatically.
