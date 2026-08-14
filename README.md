# DADA2 Workflow for Illumina Paired-End 16S rRNA Amplicon Sequencing

[![R Version](https://img.shields.io/badge/R-%3E%3D4.1-blue)](https://www.r-project.org/) [![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE) [![DADA2](https://img.shields.io/badge/DADA2-Bioconductor-orange)](https://benjjneb.github.io/dada2/) [![Platform](https://img.shields.io/badge/platform-Linux-lightgrey)](#4-install-external-tools)

A reproducible R-based pipeline for processing **Illumina paired-end** 16S rRNA amplicon sequencing data using the [DADA2](https://benjjneb.github.io/dada2/) algorithm. The pipeline is organized as RMarkdown notebooks covering quality control, primer trimming, denoising, taxonomy assignment, phylogenetic tree construction, 16S copy-number correction, microbial load correction, and phyloseq object generation.

------------------------------------------------------------------------

## Table of Contents

- [Key Features](#key-features)
- [Pipeline Overview](#pipeline-overview)
- [Quick Start](#quick-start)
- [Setup](#setup)
  - [1. Clone the Repository](#1-clone-the-repository)
  - [2. Open the R Project in RStudio](#2-open-the-r-project-in-rstudio)
  - [3. Install R Dependencies](#3-install-r-dependencies)
  - [4. Install External Tools](#4-install-external-tools)
  - [5. Install PICRUSt2 (Optional — for Step 7)](#5-install-picrust2-optional--for-step-7)
  - [6. Download Reference Databases](#6-download-reference-databases)
  - [7. Copy Your FASTQ Files](#7-copy-your-fastq-files)
  - [8. Copy Your Cell Count File (Optional)](#8-copy-your-cell-count-file-optional)
- [Running the Pipeline](#running-the-pipeline)
  - [Step 1 — Data Integrity Check](#step-1--data-integrity-check-optional)
  - [Step 2 — FastQC Quality Reports](#step-2--fastqc-quality-reports-optional)
  - [Step 3 — Primer Trimming](#step-3--primer-trimming)
  - [Step 4 — DADA2 Parameter Explorer](#step-4--dada2-parameter-explorer)
  - [Step 5 — DADA2 Pipeline](#step-5--dada2-pipeline)
  - [Step 6 — Phylogenetic Tree](#step-6--phylogenetic-tree-optional)
  - [Step 7 — 16S Copy Number Correction](#step-7--16s-copy-number-correction-optional)
  - [Step 8 — Microbial Load Correction](#step-8--microbial-load-correction-optional)
  - [Step 9 — Phyloseq Object](#step-9--phyloseq-object-optional)
- [Column Dictionaries](#column-dictionaries)
- [Project Structure](#project-structure)
- [References](#references)
- [License](#license)
- [Acknowledgments](#acknowledgments)

------------------------------------------------------------------------

## Key Features

- **Modular design** — Self-contained RMarkdown notebooks that can be run independently or as a complete end-to-end workflow
- **Dual taxonomy databases** — Assign taxonomy using [SILVA](https://www.arb-silva.de/), [GTDB](https://gtdb.ecogenomic.org/), or both simultaneously
- **Interactive quality-filter selection** — Built-in Shiny app for visualizing the amplified amplicon target and its coordinates, then inspecting data while choosing optimal DADA2 filtering parameters
- **16S copy-number correction** — [PICRUSt2](https://github.com/picrust/picrust2/wiki)-based hidden-state prediction of each ASV's genomic 16S copy number, applied to correct raw abundances toward relative cell (genome) abundance
- **Optional microbial load correction** — Converts relative ASV abundances into absolute, cell-count-scaled Quantitative Microbiome Profiles following [Vandeputte et al. (2017)](https://doi.org/10.1038/nature24460)
- **Multi-object phyloseq generation** — The final step builds one object per taxonomy database and available abundance source
- **Project-local tool installation** — [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [FastTree](https://morgannprice.github.io/fasttree/), [Cutadapt](https://cutadapt.readthedocs.io/), and [MultiQC](https://multiqc.info/) are installed inside the project; the tools themselves do not require a system-wide installation. [PICRUSt2](https://github.com/picrust/picrust2/wiki) is installed separately into its own conda environment.
- **Comprehensive outputs** — Per-step Excel summaries (each ending in a `Column_Dictionary` sheet documenting every column of every sheet), [phyloseq](https://joey711.github.io/phyloseq/) objects, Newick trees, and interactive HTML composition plots
- **Reproducible paths** — All notebooks use project-relative paths, making the workflow portable and reproducible

------------------------------------------------------------------------

## Pipeline Overview

```         
       Raw FASTQ Files
              │
              ▼
    ┌───────────────────┐
    │ Step 1 (Optional) │  Verify file integrity, pair matching,
    │  Integrity Check  │  read counts, and FASTQ validation
    └─────────┬─────────┘
              │
              ▼
    ┌───────────────────┐
    │ Step 2 (Optional) │  Per-sample FastQC reports aggregated
    │  FastQC Reports   │  into a MultiQC summary
    └─────────┬─────────┘
              │
              ▼
    ┌───────────────────┐
    │ Step 3 (Required) │  Remove PCR primer sequences using
    │  Primer Trimming  │  Cutadapt, with before/after verification
    └─────────┬─────────┘
              │
              ▼
    ┌───────────────────┐
    │       Step 4      │
    │    (Recommended)  │  Interactive Shiny app for visualizing
    │   DADA2 Parameter │  target/primer coordinates and selecting
    │      Explorer     │  truncLen/maxEE before the pipeline
    └─────────┬─────────┘
              │
              ▼
    ┌───────────────────┐
    │ Step 5 (Required) │  Quality filtering, error learning,
    │  DADA2 Pipeline   │  denoising, merging, chimera removal,
    │                   │  ASV table, and taxonomy assignment
    └─────────┬─────────┘
              │
              ├───────────────────────────┬───────────────────────────┐
              ▼                           ▼                           │
    ┌───────────────────┐       ┌───────────────────┐                 │
    │ Step 6 (Optional) │       │ Step 7 (Optional) │                 │
    │ Phylogenetic Tree │       │    Copy Number    │                 │
    │                   │       │    Correction     │                 │
    └─────────┬─────────┘       └─────────┬─────────┘                 │
              │                           │                           │
              │                           ▼                           │
              │                 ┌───────────────────┐                 │
              │                 │ Step 8 (Optional) │                 │
              │                 │   Microbial Load  │                 │
              │                 │    Correction     │                 │
              │                 └─────────┬─────────┘                 │
              │                           │                           │
              └─────────────┬─────────────┴───────────────────────────┘
                            ▼
                   ┌───────────────────┐
                   │ Step 9 (Optional) │  Auto-detects every abundance
                   │  Phyloseq Object  │  source you ran (raw / copy-
                   │                   │  number- / microbial-load-
                   │                   │  corrected) and builds one
                   │                   │  phyloseq object per
                   │                   │  database x source
                   └─────────┬─────────┘
                             │
                             ▼
                    Analysis-Ready Data
                (phyloseq .RData, corrected
                 abundance tables, & plots)
```

Step 4 is optional but recommended. Its DADA2 Parameter Explorer Shiny app is run interactively, while the accompanying RMarkdown guide can be used as a reference document. The app visualizes the amplified target and primer coordinates, helps you inspect quality and retention while choosing `truncLen`/`maxEE`, and keeps `truncQ` fixed at 2 before Step 5. Steps 6 and 7 are independent optional branches off Step 5. Step 8 follows Step 7 because microbial load correction requires the copy-number-corrected table; the tree from Step 6 is independent of both. Step 9 runs last, validates the provenance of optional corrected inputs, and builds every valid taxonomy-database x abundance-source combination.

------------------------------------------------------------------------

## Setup

### 1. Clone the Repository

Go to the repository page, click the green **`<> Code`** button, and select **Download ZIP**. Once downloaded, extract the `.zip` file and move the folder to your desired location.

or from the terminal via:

``` bash
git clone https://github.com/changlabs/DADA2_16S_Illumina_PairedEnd_Workflow.git
```

### 2. Open the R Project in RStudio

Open the R-project file `DADA2_16S_Illumina_PairedEnd_Workflow.Rproj` by double-clicking it, or from inside [RStudio](https://posit.co/products/open-source/rstudio/) via **File → Open Project**. All notebook paths use `here::here()` and resolve relative to this project root — always work from within the `.Rproj` session.

### 3. Install R Dependencies

Open [`setup/install_R_dependencies.R`](setup/install_R_dependencies.R) in RStudio and run it with **Source**. This installs all required CRAN and [Bioconductor](https://bioconductor.org/) packages, including [`dada2`](https://benjjneb.github.io/dada2/), [`DECIPHER`](http://www2.decipher.codes/), [`phyloseq`](https://joey711.github.io/phyloseq/), [`Biostrings`](https://bioconductor.org/packages/release/bioc/html/Biostrings.html), [`ShortRead`](https://bioconductor.org/packages/release/bioc/html/ShortRead.html), [`ape`](https://cran.r-project.org/package=ape), [`phangorn`](https://cran.r-project.org/package=phangorn), [`data.table`](https://cran.r-project.org/package=data.table), [`openxlsx`](https://cran.r-project.org/package=openxlsx), and the Shiny packages needed for the DADA2 Parameter Explorer.

This package set covers all notebooks, including optional Steps 6–9. Step 7 additionally requires the external PICRUSt2 installation described below, but none of the optional notebooks needs a separate R-package installation step.

On Linux you may need system libraries before running the script:

``` bash
sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev libfontconfig1-dev
```

### 4. Install External Tools

Before running the installer, make sure the following system prerequisites are available:

- Python 3.9 or newer; Debian/Ubuntu also requires `python3-venv`
- A Java Runtime Environment (JRE) for FastQC
- A C compiler (`gcc`, `cc`, or `clang`) for FastTree

Installing these prerequisites may require administrator access. On Debian/Ubuntu, for example, use `sudo apt-get install python3-venv default-jre gcc`. On macOS, install a JRE and use `xcode-select --install` for Apple clang.

Then open [`setup/install_required_tools.R`](setup/install_required_tools.R) in RStudio and run it with **Source**. The script downloads and installs the following tools into the project-local [`tools/`](tools/) directory without requiring conda or a system-wide installation of the tools themselves:

| Tool | Installed location | Notes |
|------------------------|------------------------|------------------------|
| [Cutadapt](https://cutadapt.readthedocs.io/) | [`tools/cutadapt/venv/`](tools/cutadapt/venv/) | Python virtual environment |
| [MultiQC](https://multiqc.info/) | [`tools/multiqc/venv/`](tools/multiqc/venv/) | Python virtual environment |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | [`tools/FastQC/`](tools/FastQC/) | Java-based, requires a JRE on PATH |
| [FastTree](https://morgannprice.github.io/fasttree/) | [`tools/fasttree/FastTree`](tools/fasttree/FastTree) | Compiled from source (C); requires `gcc` or `clang` |

After installation, verify that the executables have execute permissions:

``` bash
ls -l tools/FastQC/fastqc tools/fasttree/FastTree
```

If either shows `-rw-r--r--` instead of `-rwxr-xr-x`, restore the permission with:

``` bash
chmod +x tools/FastQC/fastqc tools/fasttree/FastTree
```

[FastTree](https://morgannprice.github.io/fasttree/) is compiled for your specific architecture (x86-64 or ARM64) and will attempt an OpenMP multi-threaded build first, falling back to single-threaded if OpenMP is not available (which is normal on macOS with Apple clang).

### 5. Install PICRUSt2 (Optional — for Step 7)

**Skip this step if you do not plan to run the optional [Step 7](#step-7--16s-copy-number-correction-optional) notebook for 16S copy number correction.**

[PICRUSt2](https://github.com/picrust/picrust2/wiki) is a conda package with several compiled phylogenetics dependencies (HMMER, EPA-ng, gappa, SEPP) that are not practical to manage inside the plain pip virtual environments used above, so it is installed separately, into its own dedicated conda environment. This requires a working conda/miniconda installation already present on your machine.

``` bash
bash setup/install_picrust2.sh
```

This installer supports **Linux only** and exits without making changes on macOS or other operating systems. By default it creates a conda environment named `picrust2`, which is the name [Step 7 (16S Copy Number Correction)](R/notebooks/7_copy_number_correction.md) expects. Run it once per Linux machine.

### 6. Download Reference Databases

Open [`setup/download_reference_databases.R`](setup/download_reference_databases.R) in RStudio and run it with **Source**. This downloads DADA2-formatted reference databases into [`tools/trainsets/`](tools/trainsets/):

A `download_manifest.txt` file is written to each subfolder recording the source URLs, download timestamps, file sizes, exact file paths, and the exact-size/gzip-stream verification result.

[SILVA](https://www.arb-silva.de/) is appropriate for most studies. [GTDB](https://gtdb.ecogenomic.org/) uses a rank-normalized, genome-based taxonomy and is better suited for prokaryote-focused analyses where consistent genus/species nomenclature matters.

### 7. Copy Your FASTQ Files

Copy your paired-end raw FASTQ files into the [`data/fastq/`](data/fastq/) directory.

Default file naming requirements:

- Forward reads: `{SampleID}_..._L001_R1_001.fastq.gz` (or `.fastq`, uncompressed)
- Reverse reads: `{SampleID}_..._L001_R2_001.fastq.gz` (or `.fastq`, uncompressed)
- The complete forward and reverse filename stems must match apart from the configured read token.
- The text before the first underscore is used as the **sample identifier**; use hyphens rather than underscores inside the identifier itself.
- Both gzip-compressed and uncompressed files are detected automatically.

Use the default Illumina tokens for an unedited, pipeline-wide run. Other naming conventions require coordinated configuration changes and may not be recognized by Step 2's forward/reverse summary. See [`data/README.md`](data/README.md) for the full requirements.

### 8. Copy Your Cell Count File (Optional)

Only needed for Step 8 (Microbial Load Correction) if you have independent microbial load (cell count) measurements per sample -- e.g. from flow cytometry or qPCR. Skip this if you don't; every other step works on relative abundances without it.

The cell-count template is provided at [`data/cell_count/cell_count.tsv`](data/cell_count/cell_count.tsv), which is also this pipeline's default expected input path. Replace its generic sample IDs and placeholder counts with your measurements. See [`data/README.md`](data/README.md) for the full file format.

### 9. Add Your Metadata File (Optional)

Only needed for Step 9 if you want experimental sample information included in the generated phyloseq objects. Replace the template at [`data/metadata.tsv`](data/metadata.tsv) with metadata whose sample identifiers match your data. If no metadata file is provided, Step 9 still runs and generates the phyloseq objects without metadata. - See [`data/README.md`](data/README.md) for the full file format.

------------------------------------------------------------------------

## Running the Pipeline

Open each required notebook in RStudio and run **chunkwise** or click **Run All** (Ctrl+Alt+R / Cmd+Alt+R). Use **Knit** when you want an HTML report. Steps 1 and 2 are optional preflight checks; Step 3 precedes the recommended Step 4 app and required Step 5 pipeline. After Step 5, Steps 6 and 7 are independent optional branches. Step 8 is optional but requires Step 7. Run Step 9 last if you want phyloseq objects assembled from whichever upstream outputs are available.

### Step 1 — [Data Integrity Check](R/notebooks/1_data_integrity_check.md) *(optional)*

Open [`1_data_integrity_check.Rmd`](R/notebooks/1_data_integrity_check.Rmd) and run as-is. Verifies that all FASTQ files are intact, correctly paired, and consistently named before proceeding.

### Step 2 — [FastQC Quality Reports](R/notebooks/2_fastqc_quality_reports.md) *(optional)*

Open [`2_fastqc_quality_reports.Rmd`](R/notebooks/2_fastqc_quality_reports.Rmd) and run as-is. Generates per-sample [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) reports and an aggregated [MultiQC](https://multiqc.info/) summary. Use the quality profiles to inform the filter parameters in [Step 4](#step-4--dada2-parameter-explorer) and [Step 5](#step-5--dada2-pipeline).

### Step 3 — [Primer Trimming](R/notebooks/3_cutadapt_primer_trimming.md) *(required)*

Open [`3_cutadapt_primer_trimming.Rmd`](R/notebooks/3_cutadapt_primer_trimming.Rmd). **Edit the primer sequences** in the configuration section to match your library preparation before running.

### Step 4 — [DADA2 Parameter Explorer](R/notebooks/4_dada2_parameter_selection.md) *(recommended)*

After primer trimming, open the interactive Shiny app ([`R/shiny/dada2_parameter_selection_app.R`](R/shiny/dada2_parameter_selection_app.R)) in RStudio and click **Run App**, or run the following command in the R console:

``` r
shiny::runApp("R/shiny/dada2_parameter_selection_app.R")
```

The app follows a four-tab workflow: **Visualize** the amplified target and primer coordinates, **Select** `truncLen` and `maxEE` from the loaded FASTQ quality/retention profiles, optionally **Validate** representative samples, and **Export** the selected values. `truncQ` remains fixed at 2. Clicking **Save parameters** writes `results/4_dada2_parameter_selection/dada2_filter_parameters.xlsx` with `Info`, `Parameters`, and `Column_Dictionary` sheets; [Step 5](#step-5--dada2-pipeline) imports the six selectable values from `Parameters` automatically. See the [DADA2 Parameter Explorer User Guide](R/notebooks/4_dada2_parameter_selection.md) for the full walkthrough and screenshots.

### Step 5 — [DADA2 Pipeline](R/notebooks/5_dada2_pipeline.md) *(required)*

Open [`5_dada2_pipeline.Rmd`](R/notebooks/5_dada2_pipeline.Rmd). At the beginning of its configuration, Step 5 checks for the report saved by the [Step 4 app](#step-4--dada2-parameter-explorer). If it exists, six saved parameters are imported automatically: two `truncLen` values, two `maxEE` values, and two amplicon-length bounds. `truncQ` remains fixed at 2 in Step 5. If the report does not exist, Step 5 uses the clearly marked `default_dada2_parameters` values instead; review that block before running with a different dataset. An active-parameter table shows both the final value and its source. Step 5 also stops early with a clear message if a `truncLen` exceeds the observed primer-trimmed read lengths.

Select the taxonomy database(s) in the configuration section (`"SILVA"`, `"GTDB"`, or `"BOTH"`).

### Step 6 — [Phylogenetic Tree](R/notebooks/6_phylogenetic_tree.md) *(optional)*

Open [`6_phylogenetic_tree.Rmd`](R/notebooks/6_phylogenetic_tree.Rmd) and run as-is. Performs multiple sequence alignment with [DECIPHER](http://www2.decipher.codes/) and constructs a maximum-likelihood phylogenetic tree using [FastTree](https://morgannprice.github.io/fasttree/). Runs directly off Step 5's output — does not depend on, and is not blocked by, Steps 7 or 8.

### Step 7 — [16S Copy Number Correction](R/notebooks/7_copy_number_correction.md) *(optional)*

Open [`7_copy_number_correction.Rmd`](R/notebooks/7_copy_number_correction.Rmd) once PICRUSt2 is installed. It phylogenetically places Step 5's ASVs into PICRUSt2's reference tree (`place_seqs.py`), predicts each ASV's 16S rRNA gene copy number via hidden-state prediction (`hsp.py`), and divides raw ASV counts by their predicted copy number. The provided [`setup/install_picrust2.sh`](setup/install_picrust2.sh) installer supports Linux only; users on another operating system must supply a compatible `picrust2` conda environment themselves. If you need neither Step 7 nor Step 8, skip directly to [Step 6](#step-6--phylogenetic-tree-optional) or [Step 9](#step-9--phyloseq-object-optional).

### Step 8 — [Microbial Load Correction](R/notebooks/8_microbial_load_correction.md) *(optional)*

Open [`8_microbial_load_correction.Rmd`](R/notebooks/8_microbial_load_correction.Rmd). Requires [Step 7](#step-7--16s-copy-number-correction-optional) to have been run first, plus an independent microbial-load (cell count) measurement per sample — typically from flow cytometry or qPCR — that this pipeline cannot generate for you. Prepare a `SampleID`/`Cell_Count` TSV file (see [`data/README.md`](data/README.md) and the current file at [`data/cell_count/cell_count.tsv`](data/cell_count/cell_count.tsv)), then run the notebook as-is. This step rarefies each sample to a common sampling depth per cell using [phyloseq](https://joey711.github.io/phyloseq/)'s `rarefy_even_depth()`, then rescales by each sample's cell count to produce an absolute-abundance Quantitative Microbiome Profile, following [Vandeputte et al. (2017)](https://doi.org/10.1038/nature24460). If you do not have microbial-load data, skip this notebook — every other step works on relative abundances without it.

### Step 9 — [Phyloseq Object](R/notebooks/9_phyloseq_object.md) *(optional)*

Open [`9_phyloseq_object.Rmd`](R/notebooks/9_phyloseq_object.Rmd) after whichever of Steps 5–8 you plan to use. Before running it, replace the generic [`data/metadata.tsv`](data/metadata.tsv) template with matching experimental metadata, change `metadata_path` to another supported TSV/CSV file, or set `metadata_path <- NULL`. Step 9 integrates the available ASV table(s), taxonomy assignments, optional metadata, and optional Step 6 tree into analysis-ready [phyloseq](https://joey711.github.io/phyloseq/) objects. It builds one object per taxonomy database for each available abundance source. Optional Step 7 and Step 8 outputs are accepted only when their provenance workbooks and checksums match the current upstream inputs; otherwise rerun the stale upstream notebook. Step 9 stages and verifies a complete new deliverable set before replacing its previous output.

------------------------------------------------------------------------

## Column Dictionaries

Every Excel workbook produced in this pipeline—including the Step 4 app export—ends with a trailing **`Column_Dictionary`** sheet that documents each data column in plain language. This is generated automatically by [`R/functions/build_column_dictionary_function.R`](R/functions/build_column_dictionary_function.R). If you add or rename an exported column, update the corresponding `descriptions` vector so the dictionary stays accurate.

------------------------------------------------------------------------

## Project Structure

The tree below summarizes the maintained workflow files plus the project-local tools and reference databases created by Setup. Notebook HTML and Markdown companions are retained for GitHub viewing.

```         
DADA2_16S_Illumina_PairedEnd_Workflow/
├── DADA2_16S_Illumina_PairedEnd_Workflow.Rproj
├── README.md
├── LICENSE
├── .gitignore
├── .gitattributes
│
├── .github/
│   └── workflows/
│       └── pages.yml
│
├── pages/
│   └── index.html
│
├── setup/
│   ├── install_R_dependencies.R
│   ├── install_required_tools.R
│   ├── download_reference_databases.R
│   └── install_picrust2.sh
│
├── R/
│   ├── notebooks/
│   │   ├── 1_data_integrity_check.Rmd
│   │   ├── 2_fastqc_quality_reports.Rmd
│   │   ├── 3_cutadapt_primer_trimming.Rmd
│   │   ├── 4_dada2_parameter_selection.Rmd
│   │   ├── 5_dada2_pipeline.Rmd
│   │   ├── 6_phylogenetic_tree.Rmd
│   │   ├── 7_copy_number_correction.Rmd
│   │   ├── 8_microbial_load_correction.Rmd
│   │   ├── 9_phyloseq_object.Rmd
│   │   ├── *.md
│   │   └── *.html
│   │
│   ├── functions/
│   │   ├── add_sheet_to_excel_function.R
│   │   ├── build_column_dictionary_function.R
│   │   ├── render_output_links_function.R
│   │   └── render_output_tree_function.R
│   │
│   ├── shiny/
│   │   ├── dada2_parameter_selection_app.R
│   │   ├── functions/
│   │   │   └── paired_read_retention_engine_function.R
│   │   └── tests/
│   │       ├── test_dada2_parameter_selection_app.R
│   │       └── test_paired_read_retention_engine.R
│   │
│   └── images/
│       ├── 4_dada2_parameter_selection_visualize.png
│       ├── 4_dada2_parameter_selection_select_quality_profiles.png
│       ├── 4_dada2_parameter_selection_select_retained_reads.png
│       ├── 4_dada2_parameter_selection_validate.png
│       ├── 4_dada2_parameter_selection_export.png
│       ├── 5_error_rate_plots.png
│       └── 5_read_quality_profiles.png
│
├── data/
│   ├── README.md
│   ├── metadata.tsv
│   ├── fastq/
│   │   └── .gitkeep
│   └── cell_count/
│       └── cell_count.tsv
│
└── tools/
    ├── cutadapt/
    │   ├── requirements.txt
    │   ├── install_manifest.txt
    │   └── venv/
    ├── multiqc/
    │   ├── requirements.txt
    │   ├── install_manifest.txt
    │   └── venv/
    ├── FastQC/
    │   ├── fastqc
    │   └── install_manifest.txt
    ├── fasttree/
    │   ├── FastTree
    │   ├── FastTree.c
    │   └── install_manifest.txt
    └── trainsets/
        ├── SILVA/
        │   ├── silva_nr99_v138.2_toGenus_trainset.fa.gz
        │   ├── silva_v138.2_assignSpecies.fa.gz
        │   └── download_manifest.txt
        └── GTDB/
            ├── GTDB_bac120_arc53_ssu_r220_genus.fa.gz
            ├── GTDB_bac120_arc53_ssu_r220_species.fa.gz
            └── download_manifest.txt
```

------------------------------------------------------------------------

## References

### Methods

- Callahan BJ, et al. (2016). DADA2: High-resolution sample inference from Illumina amplicon data. *Nature Methods*, 13(7):581-583. [DOI:10.1038/nmeth.3869](https://doi.org/10.1038/nmeth.3869)
- Price MN, et al. (2010). FastTree 2 – Approximately Maximum-Likelihood Trees for Large Alignments. *PLoS ONE*, 5(3):e9490. [DOI:10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)
- McMurdie PJ, Holmes S (2013). phyloseq: An R Package for Reproducible Interactive Analysis and Graphics of Microbiome Census Data. *PLoS ONE*, 8(4):e61217. [DOI:10.1371/journal.pone.0061217](https://doi.org/10.1371/journal.pone.0061217)

### 16S Copy Number Correction

- Douglas GM, Maffei VJ, Zaneveld JR, et al. (2020). PICRUSt2 for prediction of metagenome functions. *Nat Biotechnol* 38, 685-688. [DOI:10.1038/s41587-020-0548-6](https://doi.org/10.1038/s41587-020-0548-6)
- Langille MGI, Zaneveld J, Caporaso JG, et al. (2013). Predictive functional profiling of microbial communities using 16S rRNA marker gene sequences. *Nat Biotechnol* 31, 814-821. [DOI:10.1038/nbt.2676](https://doi.org/10.1038/nbt.2676)
- Kembel SW, Wu M, Eisen JA, Green JL (2012). Incorporating 16S Gene Copy Number Information Improves Estimates of Microbial Diversity and Abundance. *PLoS Comput Biol* 8(10):e1002743. [DOI:10.1371/journal.pcbi.1002743](https://doi.org/10.1371/journal.pcbi.1002743)
- [PICRUSt2 GitHub repository](https://github.com/picrust/picrust2) / [PICRUSt2 wiki](https://github.com/picrust/picrust2/wiki)

### Microbial Load Correction

- Vandeputte D, Kathagen G, D'hoe K, et al. (2017). Quantitative microbiome profiling links gut community variation to microbial load. *Nature* 551, 507-511. [DOI:10.1038/nature24460](https://doi.org/10.1038/nature24460)
- [raeslab/QMP GitHub repository](https://github.com/raeslab/QMP) — original reference implementation (`QMP.R`) this workflow's Step 8 algorithm follows (original concept: Gwen Falony; original script contributors: Doris Vandeputte, Gunter Kathagen, Kevin d'Hoe, Joao Sabino, Mireia Valles-Colomer, Sara Vieira-Silva).

### Databases

- Quast C, et al. (2013). The SILVA ribosomal RNA gene database project: improved data processing and web-based tools. *Nucleic Acids Research*, 41(D1):D590-D596. [DADA2-formatted SILVA reference database](https://zenodo.org/records/14169026)
- Parks DH, et al. (2022). GTDB: an ongoing census of bacterial and archaeal diversity through a phylogenetically consistent, rank normalized and complete genome-based taxonomy. *Nucleic Acids Research*, 50(D1):D199-D207. [DADA2-formatted GTDB reference database](https://zenodo.org/records/13984843)

------------------------------------------------------------------------

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

------------------------------------------------------------------------

## Acknowledgments

- [DADA2 Team](https://benjjneb.github.io/dada2/) for the core algorithm
- [PICRUSt2 Team](https://github.com/picrust/picrust2) for the phylogenetic-placement / hidden-state-prediction methodology adapted in Step 7
- [raeslab/QMP](https://github.com/raeslab/QMP) authors for the original Quantitative Microbiome Profiling reference implementation adapted in Step 8
- [Bioconductor](https://bioconductor.org/) community for Bioconductor packages
- [CRAN - R Project](https://cran.r-project.org/) community for R packages
- [Claude AI](https://claude.ai/) for assistance in code and documentation

------------------------------------------------------------------------

**Author**: Amro Abbas - Generated with Claude AI assistance\
**Last Updated**: August 2026
