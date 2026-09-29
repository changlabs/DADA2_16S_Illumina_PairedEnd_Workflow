# Synthetic cell-count fixture

The values in [`cell_count.tsv`](cell_count.tsv) are **synthetic test data**. They were not measured in the Kozich et al. samples and must not be used for biological conclusions.

The 19 mouse-feces values represent plausible total microbial loads in cells per gram of wet feces. They span `8.4e10`–`1.38e11`, vary from day to day, and deliberately contain no imposed early-versus-late trend. This order of magnitude is consistent with published mouse intestinal and fecal bacterial-density measurements. Sarma-Rupavtarm et al. summarize approximately `1e10`–`1e12` bacteria per gram in the murine large intestine, while Nair et al. report conventional-mouse fecal loads on the order of `1e11` 16S copies per gram. These sources support the scale of the fixture, not the individual sample values; gene-copy and cell counts are not interchangeable measurements.

The mock-community sample has no biologically meaningful fecal density. Its `1.10e11` value is solely a computational placeholder on the same numeric scale, included because Step 8 requires one positive value for every sample in the abundance table. It is marked separately in `Data_Status` and uses the unit `cells_per_g_input_equivalent`.

Only `SampleID` and `Cell_Count` are consumed by Step 8. `Unit` and `Data_Status` are audit columns that make the synthetic status explicit.

References used to calibrate the plausible range:

- Sarma-Rupavtarm RB, Ge Z, Schauer DB, Fox JG, Polz MF. (2004). Spatial distribution and stability of the eight microbial species of the altered Schaedler flora in the mouse gastrointestinal tract. *Applied and Environmental Microbiology*, 70(5), 2791–2800. https://doi.org/10.1128/AEM.70.5.2791-2800.2004
- Nair AB, Jacob S, et al. (2022). Jejunoileal mucosal growth in mice with a limited microbiome. *PLOS ONE*, 17(4), e0266251. https://doi.org/10.1371/journal.pone.0266251

These synthetic values exist only to exercise the copy-number/microbial-load branch of the workflow. Replace them with independently measured values from flow cytometry, qPCR with an appropriate conversion, or another validated absolute-quantification method for real analyses.
