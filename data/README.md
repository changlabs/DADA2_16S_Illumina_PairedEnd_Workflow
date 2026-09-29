# Data Directory Guide

The `data/` directory holds the user-provided inputs for this workflow:

- paired-end sequencing reads in `data/fastq/`;
- optional sample metadata at `data/metadata.tsv` for Step 9;
- optional microbial-load measurements at `data/cell_count/cell_count.tsv` for Step 8.

FASTQ datasets placed under `data/fastq/` are ignored by Git. The separately documented tutorial fixture under `example/data/fastq/` is intentionally tracked. The metadata and cell-count files in this normal `data/` directory are templates and must be reviewed or replaced before analysis.

------------------------------------------------------------------------

## Bundled Example Dataset

A clone-ready test dataset is maintained separately under [`example/`](../example/) so it cannot be discovered together with files placed in `data/fastq/`. Its generated reruns also go to the ignored `example/run_results/` directory rather than the normal `results/` tree. This separation prevents an example run from overwriting a user's results and prevents newly added user FASTQs from changing the example analysis.

The example contains 20 paired samples (19 longitudinal mouse-feces samples and one defined mock community) of **Illumina MiSeq paired-end 2 x 250 bp reads targeting the V4 region**, from the official DADA2 paired-end tutorial / mothur MiSeq SOP dataset. Filenames use descriptive IDs such as `Mouse-F3-Day-141` and `Mock-Community`; [`example/data/metadata.tsv`](../example/data/metadata.tsv) retains each source code in `OriginalSampleID` and otherwise follows the metadata layout documented below.

The example inputs also include the validated Step 4 export [`example/data/dada2_filter_parameters.xlsx`](../example/data/dada2_filter_parameters.xlsx). Step 5 reads this file directly during the example run, using V4/2 x 250 settings (`truncLen` 240/160, `maxEE` 2/2, and a 250–256 bp merged-amplicon range) instead of its generic fallback values.

A [`example/data/cell_count/cell_count.tsv`](../example/data/cell_count/cell_count.tsv) fixture provides synthetic, biologically plausible-scale values so the complete example can exercise Steps 7 and 8. These are not measurements from the source study and must not be used for biological interpretation. The mock-community value is a computational placeholder because a fecal cell density is not defined for that sample. The normal `Rscript example/run_example.R` command uses this fixture.

The source tutorial reads were distributed after primer removal. To test this workflow's Step 3, the bundled `fastq.gz` records physically contain a synthetic concrete 515F prefix on every R1 sequence and a synthetic concrete 806R prefix on every R2 sequence. Matching Phred-40 characters were added to the quality strings. The example runner configures the corresponding degenerate 515F/806R definitions, runs Cutadapt, and recovers the original primer-free reads. These added primer bases are a transparent test-fixture transformation, not original sequencer observations.

After installing dependencies and reference databases, run the isolated example from the repository root:

``` bash
Rscript example/run_example.R
```

Browsable outputs are included under [`example/reference_results/`](../example/reference_results/). Full provenance, run scope, filename mapping, and third-party data notice are in the [example guide](../example/README.md).

Dataset reference:

> Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013). Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Applied and Environmental Microbiology*, 79(17), 5112–5120. https://doi.org/10.1128/AEM.01043-13

The source archive is linked from both the [DADA2 tutorial](https://benjjneb.github.io/dada2/tutorial) and the [mothur MiSeq SOP](https://mothur.org/wiki/miseq_sop/).

------------------------------------------------------------------------

## Paired-End FASTQ Files

Place raw paired-end FASTQ files directly in [`data/fastq/`](fastq/) before running the pipeline.

### File requirements

- **Format:** Standard FASTQ using `.fastq`, `.fq`, or either gzip-compressed variant. Extensions are matched case-insensitively and formats may be mixed in one run.
- **Read type:** Paired-end, with one forward and one reverse file per sample.
- **Location:** `data/fastq/`, not directly under `data/`.

### Default naming convention

Use the default Illumina tokens below for an unedited, pipeline-wide run. Apart from the read-direction token, the complete forward and reverse filename stems must be identical.

| Forward                              | Reverse                              |
|------------------------------------|------------------------------------|
| `Treatment1_S1_L001_R1_001.fastq.gz` | `Treatment1_S1_L001_R2_001.fastq.gz` |
| `SampleA_L001_R1_001.fastq.gz`       | `SampleA_L001_R2_001.fastq.gz`       |
| `SampleA_L001_R1_001.fq`             | `SampleA_L001_R2_001.fq`             |

The default discovery patterns match `_L001_R1_001` and `_L001_R2_001`, followed by `.fastq` or `.fq` and an optional `.gz` suffix, case-insensitively.

### Sample ID extraction

The text before the first underscore (`_`) is used as the sample ID throughout the workflow.

| Filename                             | Sample ID    |
|--------------------------------------|--------------|
| `SampleA_L001_R1_001.fastq.gz`       | `SampleA`    |
| `Treatment1_S1_L001_R1_001.fastq.gz` | `Treatment1` |
| `WT-Control_L001_R1_001.fastq.gz`    | `WT-Control` |

Use hyphens rather than underscores inside a sample ID. The same IDs must be used in the optional metadata and cell-count files.

### Alternative filename tokens

If your files do not use `_L001_R1_001` and `_L001_R2_001`, review these settings together:

| Component | Setting to review |
|------------------------------------|------------------------------------|
| [Step 1 — Data Integrity Check](../R/notebooks/1_data_integrity_check.md) | `forward_token`, `reverse_token` |
| [Step 3 — Primer Trimming](../R/notebooks/3_cutadapt_primer_trimming.md) | `forward_pattern`, `reverse_pattern` |
| [Step 4 — Parameter Explorer](../R/shiny/dada2_parameter_selection_app.R) | Select the matching pattern under **Input Data**; the app supports `_L001_R1_001`/`_L001_R2_001`, `_R1`/`_R2`, and `_1`/`_2` |
| [Step 5 — DADA2 Pipeline](../R/notebooks/5_dada2_pipeline.md) | `forward_pattern`, `reverse_pattern` |

[Step 2 — FastQC Quality Reports](../R/notebooks/2_fastqc_quality_reports.md) analyzes every discovered FASTQ file, but its direction-specific summary and paired report table identify reads using `_R1_` and `_R2_` within filenames. Plain `_R1`/`_R2` immediately before the extension and `_1`/`_2` therefore require changes to Step 2's parsing logic for full compatibility.

No configuration change is needed solely to switch among supported FASTQ extensions or compression states.

------------------------------------------------------------------------

## Sample Metadata (Optional)

Sample metadata is optional and is used by [Step 9 — Phyloseq Object Construction](../R/notebooks/9_phyloseq_object.md). A phyloseq object can be built without metadata, but grouping, faceting, and statistical comparisons by variables such as treatment, time point, subject, or sample type require it.

### File requirements

- **Default path:** [`data/metadata.tsv`](metadata.tsv).
- **Format:** Tab-separated values with a header row. A comma-separated `.csv` is also accepted; Step 9 selects the delimiter from the extension.
- **Required column:** `SampleID`.
- **Optional columns:** Any number of sample-level variables. You may rename, remove, or add these columns to match your experimental design; only `SampleID` is required by the workflow.

| Column | Description |
|------------------------------------|------------------------------------|
| `SampleID` | Must exactly match the ID extracted from the corresponding FASTQ filename—the text before its first underscore. |
| `SubjectID` | Optional subject, participant, animal, site, or other experimental-unit identifier. Useful for paired or repeated-measures designs. |
| `Treatment` | Optional experimental group, exposure, intervention, or control assignment. |
| `Timepoint` | Optional collection time, visit, or experimental stage. It may be numeric or categorical. |
| `SampleType` | Optional specimen, material, habitat, or control type. |
| `Batch` | Optional extraction, library-preparation, sequencing, or processing batch. |
| Any additional column | Allowed. Add any sample-level covariates needed for your study, such as `Site`, `Sex`, `Age`, `Diet`, `Condition`, or technical variables. Step 9 passes them to phyloseq `sample_data()`. Use `NA` when a value is unknown or does not apply. |

The included [`metadata.tsv`](metadata.tsv) is a ten-sample template using the same generic `S01`–`S10` identifiers as the cell-count template. Every value is illustrative. Replace the sample IDs and metadata values, add or remove rows to match your samples, and adapt the optional columns before analysis.

Do not rename or remove `SampleID`. All other template columns are optional: retain only those useful for your study and add as many additional columns as needed. Column names should be unique and descriptive.

### How Step 9 uses metadata

1.  `metadata_path` defaults to `here("data", "metadata.tsv")`. Set it to `NULL` to build objects without metadata, or change it to use another file.
2.  Step 9 matches `SampleID` values to the abundance table.
3.  If an abundance-table sample lacks metadata, Step 9 warns and fills its metadata fields with `NA`, preserving the sample in every phyloseq object.
4.  Metadata rows absent from the abundance table are ignored with a warning.
5.  Attached metadata columns become available for downstream grouping, faceting, and statistical analysis.

------------------------------------------------------------------------

## Microbial-Load Cell Counts (Optional)

Cell-count data are only required for [Step 8 — Microbial Load Correction](../R/notebooks/8_microbial_load_correction.md). Use independently measured microbial loads from flow cytometry, qPCR, or an equivalent method.

The included [`data/cell_count/cell_count.tsv`](cell_count/cell_count.tsv) is a template. Its `S01`–`S10` sample IDs and `100,000`–`1,000,000` counts are generic placeholders, not experimental measurements.

### File requirements

- **Default path:** `data/cell_count/cell_count.tsv`.
- **Format:** Tab-separated values with a header row.
- **Required columns:** `SampleID` and `Cell_Count`. Do not rename them.

| Column | Description |
|------------------------------------|------------------------------------|
| `SampleID` | Must exactly match the sample IDs in the abundance table and those derived from FASTQ filenames. |
| `Cell_Count` | Independently measured, finite, positive microbial load, such as cells per gram or cells per millilitre. Units may vary by study but must be consistent across all samples. |

Before running Step 8:

1.  Replace `S01`–`S10` with the exact sample identifiers in your abundance table.
2.  Replace every placeholder count with an independently measured microbial load.
3.  Add or remove rows so the cell-count file and abundance table contain exactly the same samples.
4.  Keep the filename, tab-separated format, and column names unchanged to use Step 8's default configuration. If the file lives elsewhere, update `cell_counts_path`.

Step 8 uses this file for Quantitative Microbiome Profiling following [Vandeputte et al. (2017), *Nature* 551:507–511](https://doi.org/10.1038/nature24460). The workflow adapts the original [`raeslab/QMP`](https://github.com/raeslab/QMP) reference implementation.

------------------------------------------------------------------------

## Input Consistency Checklist

Before processing a new dataset, confirm that:

- every FASTQ sample has exactly one forward and one reverse file;
- sample IDs are unique before the first underscore;
- metadata `SampleID` values use those same IDs;
- if Step 8 will be run, the cell-count file contains exactly the same sample set as its input abundance table;
- all template metadata and cell-count values have been replaced with real study data.
