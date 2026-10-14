# Example run

This directory provides saved results from a 20-sample **Illumina MiSeq paired-end 2 x 250 bp V4** example run. View the [HTML reports on the example landing page](https://changlabs.github.io/DADA2_16S_Illumina_PairedEnd_Workflow/) or browse the [example results directory](reference_results/) without installing dependencies or running the workflow.

**Input FASTQ files are not included in the repository.** The metadata, synthetic cell counts, and parameter workbook document the saved run; they do not supply its sequencing reads. Use the main workflow with your own FASTQs as described in [Setup](../README.md#setup).

## Contents

- `data/metadata.tsv`: the repository's normal metadata layout, plus `TimePeriod` and `OriginalSampleID` provenance columns.
- `data/cell_count/cell_count.tsv`: clearly labelled synthetic microbial loads matched to all 20 samples, used to generate the saved QMP results.
- `data/dada2_filter_parameters.xlsx`: the prepared Step 4 export with tutorial-appropriate DADA2 settings (`truncLen = 240/160`, `maxEE = 2/2`, retained amplicon length 250–256 bp).
- `reference_results/`: curated outputs from the complete example run, including parallel SILVA and GTDB taxonomy-dependent results, for inspection without running the workflow.

## Preparation used for the saved run

The official DADA2 tutorial files are already demultiplexed and have had barcodes, adapters, and primers removed. For the saved run's trimming demonstration, the input reads had a concrete 515F sequence (`GTGCCAGCAGCCGCGGTAA`) prefixed to every R1 sequence and a concrete 806R sequence (`GGACTACAAGGGTATCTAAT`) prefixed to every R2 sequence. Matching high-quality (`I`, Phred 40) characters were prefixed to the quality strings. Step 3 is configured with the corresponding degenerate 515F/806R definitions and one Cutadapt removal round, so it removes those 5' prefixes and reproduces the original primer-free tutorial reads exactly. Opposite-primer read-through trimming is disabled for this synthetic fixture because no 3' primer sequence was added. This is a transparent test fixture transformation, not original sequencing data at the added primer positions.

## Synthetic microbial loads

The source tutorial dataset does not contain cell counts. Values under [`data/cell_count/`](data/cell_count/) were therefore created only to test the workflow and are not observations from the Kozich et al. study. The 19 mouse-feces values range from `8.4e10` to `1.38e11` cells per gram of wet feces, a plausible conventional-mouse range informed by published murine intestinal/fecal measurements. They include modest day-to-day variation but no imposed early-versus-late trend.

The mock community does not have a biologically meaningful fecal cell density. Its value is explicitly marked `SYNTHETIC_COMPUTATIONAL_PLACEHOLDER` and exists only because Step 8 requires a positive value for every abundance-table sample. Do not use any of these values for biological inference or method validation. See the [cell-count fixture notes](data/cell_count/README.md) for the assumptions and supporting references.

## Renamed samples

The original sample codes remain in `OriginalSampleID`. The saved results use descriptive sample identifiers mapped to the original sample codes:

| Original | Example `SampleID` | Meaning |
|---|---|---|
| `F3D0` | `Mouse-F3-Day-000` | Female mouse 3, day 0 post-weaning |
| `F3D141` | `Mouse-F3-Day-141` | Female mouse 3, day 141 post-weaning |
| `Mock` | `Mock-Community` | Defined bacterial mock community |

The same rule is applied to all mouse samples. `Early` covers days 0–9 represented in the archive; `Late` covers days 141–150.

## Data provenance and citation

The reads used for the saved run were downloaded on 29 September 2026 from the official [mothur MiSeq SOP archive](https://mothur.s3.us-east-2.amazonaws.com/wiki/miseqsopdata.zip), which is also the example dataset linked by the [DADA2 paired-end tutorial](https://benjjneb.github.io/dada2/tutorial). The input reads were prepared with descriptive filenames, deterministic gzip compression, and the synthetic primer prefixes described above.

Please cite the dataset-generating study:

> Kozich JJ, Westcott SL, Baxter NT, Highlander SK, Schloss PD. (2013). Development of a dual-index sequencing strategy and curation pipeline for analyzing amplicon sequence data on the MiSeq Illumina sequencing platform. *Applied and Environmental Microbiology*, 79(17), 5112–5120. https://doi.org/10.1128/AEM.01043-13

When describing the analysis, also cite DADA2:

> Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJA, Holmes SP. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13, 581–583. https://doi.org/10.1038/nmeth.3869

## Third-party data notice

The input FASTQ reads are third-party research data. Input reads are not included in the repository, which provides saved analysis results and supporting configuration. The MIT software license does not grant rights to the original sequencing data. Source attribution is preserved above. No human participant data are involved in this mouse/microbial mock-community example.
