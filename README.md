# NF1 mutant melanoma Xenium spatial transcriptomics — manuscript reproduction code
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22258628.svg)](https://doi.org/10.5281/zenodo.22258628)
This repository contains the cleaned Xenium/spatial-transcriptomics code used to reproduce the manuscript analyses for **“Genotype-Driven Tumor Ecosystems Drive Immune Evasion and Immunotherapy Resistance in Melanoma.”**

## What this clean version preserves

The scientific settings and figure logic are kept as validated during the manuscript reproduction. The cleanup removes obsolete/diagnostic copies, gives one ordered script flow, and uses a single generated-output root (`results/`). It does **not** substitute new statistical methods for the versions that reproduced the manuscript.

Key publication-specific behaviors retained include:

- the saved 0–25 Seurat clusters and saved UMAP as the downstream starting point;
- exact 40-µm / 12-meta-niche neighborhood reconstruction followed by the documented publication-reference boundary corrections;
- Figure 3A contact effect `log2(mean distance WT / mean distance Mut)`, where positive means closer in NF1Mut;
- core-level CD8 pseudobulk DESeq2 for Figure 3F-G;
- cell-level melanoma Wilcoxon testing for Figure 4E;
- historical meta-niche limma-voom/Hallmark workflow for Figure 5A;
- raw-count EGFR representation for Figure 5B and historical response representation for Figure 5C/S6A;
- frozen historical Hallmark modules for Figure 5D;
- SCT/core-level EGFR correlations for Figure 5F;
- the recovered historical Figure S5B heatmap logic: raw RNA counts → `log1p` → mean per tissue core → `pheatmap(scale="row")`.

## Required local inputs

Place these large processed objects in `data/processed/`:

```text
merged_obj_filtered.rds
xenium_list.rds
```

They are intentionally excluded from GitHub because of size.

For exact publication-reference analyses, also add the two **small public reference files** to `data/reference/`:

```text
published_meta_niche_boundary_corrections.csv
figure5D_historical_pathway_modules.rds
```

These two reference files **should be committed to GitHub** because they are required for exact manuscript reproduction and are not the giant raw/processed data objects.

For Figure 5C / Figure S6A, the historical response annotation is expected locally at:

```text
private_provenance/patient_treatment_response_annotation.xlsx
```

This file is private by default and is excluded by `.gitignore`; do **not** commit it publicly unless a de-identified release has been approved.

## Run

From the repository root:

```r
source("scripts/00_setup/00_check_inputs_and_packages.R")
source("run_xenium_pipeline.R")
```

CellChat is separate because it is memory-intensive and version-sensitive:

```r
source("run_cellchat.R")
```

All generated files go under:

```text
results/
```

## Repository organization

```text
R/                         shared helpers/constants
config/                    project paths and analysis constants
scripts/00_setup/          preflight
scripts/01_preprocessing/  optional upstream provenance
scripts/02_annotation/     cluster annotation + Figure 1
scripts/03_meta_niches/    meta-niches + Figure 2
scripts/04_cellchat/       CellChat (separate runner)
scripts/05_spatial_distance/ Figure 3 + S4
scripts/06_antigen_presentation/ Figure 4
scripts/07_egfr/           Figure 5 + S6
scripts/08_supplement/     Figure S5
scripts/09_reproducibility/ session information
data/reference/            small frozen public reproduction inputs
private_provenance/        local-only private metadata
results/                   generated outputs
```

See `docs/PIPELINE.md` and `docs/FIGURE_TO_CODE.md` for the ordered flow and manuscript panel mapping.
