# Bundled DADA2 tutorial example

This directory contains a self-contained example profile for the workflow. It is deliberately separate from [`data/`](../data/) and [`results/`](../results/): running the example cannot mix tutorial reads with a user's FASTQ files, and adding another dataset to `data/fastq/` cannot change the example run.

## Contents

- `data/fastq/`: 20 paired samples (40 gzip-compressed FASTQ files) with documented synthetic V4 primer prefixes.
- `data/metadata.tsv`: the repository's normal metadata layout, plus `TimePeriod` and `OriginalSampleID` provenance columns.
- `config/dada2_filter_parameters.xlsx`: tutorial-appropriate DADA2 settings (`truncLen = 240/160`, `maxEE = 2/2`, retained amplicon length 250–256 bp).
- `reference_results/`: curated outputs from Steps 1, 3, 5, and 9, generated with GTDB r220, for inspection without running the workflow.
- `run_example.R`: recreates the core results under the ignored `example/run_results/` directory.

Run from the repository root after installing the R dependencies and downloading the reference databases:

``` bash
Rscript example/run_example.R
```

The runner uses `DADA2_DATA_DIR`, `DADA2_RESULTS_DIR`, and `DADA2_TAXONOMY_DATABASE` only inside its R process. Normal notebook runs remain unchanged and continue to use `data/`, `results/`, and both taxonomy databases by default.

The official DADA2 tutorial files are already demultiplexed and have had barcodes, adapters, and primers removed. To make the bundled data exercise this workflow's trimming stage, the clone-facing copies have a concrete 515F sequence (`GTGCCAGCAGCCGCGGTAA`) prefixed to every R1 sequence and a concrete 806R sequence (`GGACTACAAGGGTATCTAAT`) prefixed to every R2 sequence. Matching high-quality (`I`, Phred 40) characters were prefixed to the quality strings. Step 3 is configured with the corresponding degenerate 515F/806R definitions and one Cutadapt removal round, so it removes those 5' prefixes and reproduces the original primer-free tutorial reads exactly. Opposite-primer read-through trimming is disabled for this synthetic fixture because no 3' primer sequence was added. This is a transparent test fixture transformation, not original sequencing data at the added primer positions.

Steps 2, 6, 7, and 8 are optional and are not part of the compact reference run. Step 7 additionally needs PICRUSt2, and Step 8 needs independently measured microbial-load values that this dataset does not provide.

## Renamed samples

The original sample codes remain in `OriginalSampleID`. Clone-facing filenames use descriptive identifiers while retaining the Illumina lane/read suffix expected by the workflow:

| Original | Bundled `SampleID` | Meaning |
|---|---|---|
| `F3D0` | `Mouse-F3-Day-000` | Female mouse 3, day 0 post-weaning |
| `F3D141` | `Mouse-F3-Day-141` | Female mouse 3, day 141 post-weaning |
| `Mock` | `Mock-Community` | Defined bacterial mock community |

The same rule is applied to all mouse samples. `Early` covers days 0–9 represented in the archive; `Late` covers days 141–150.

## Data provenance and citation

Downloaded on 29 September 2026 from the official [mothur MiSeq SOP archive](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip), which is also the example dataset linked by the [DADA2 paired-end tutorial](https://benjjneb.github.io/dada2/tutorial). The bundled FASTQs differ only by descriptive filenames, deterministic gzip compression, and the synthetic primer prefixes described above.

Please cite the dataset-generating study:

> Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013). Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Applied and Environmental Microbiology*, 79(17), 5112–5120. https://doi.org/10.1128/AEM.01043-13

When describing the analysis, also cite DADA2:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. https://doi.org/10.1038/nmeth.3869

## Third-party data notice

The FASTQ reads are third-party research data from the source above. Their inclusion does not place them under this repository's MIT software license; attribution and any rights remain with the original authors and source. No human participant data are included in this mouse/microbial mock-community example.
