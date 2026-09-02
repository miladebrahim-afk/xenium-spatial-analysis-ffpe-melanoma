# Clean Xenium reproduction flow

## Analysis entry point

The validated downstream analysis starts from two processed files:

- `data/processed/merged_obj_filtered.rds`
- `data/processed/xenium_list.rds`

The first is the final filtered/integrated Seurat object containing the published 0–25 clusters and saved UMAP. The second is the 39-core spatially intact Xenium list.

The upstream scripts under `scripts/01_preprocessing/` are retained for provenance but are **not** required when those two starting objects are available.

## Main order

1. `scripts/00_setup/00_check_inputs_and_packages.R`
2. `scripts/02_annotation/`
3. `scripts/03_meta_niches/`
4. `scripts/05_spatial_distance/`
5. `scripts/06_antigen_presentation/`
6. `scripts/07_egfr/`
7. `scripts/08_supplement/`
8. `scripts/09_reproducibility/`

Run the ordered workflow with:

```r
source("run_xenium_pipeline.R")
```

CellChat is intentionally separate:

```r
source("run_cellchat.R")
```

## Publication-reference steps

The de novo meta-niche reconstruction is preserved. Exact manuscript reproduction then applies the documented public 140-cell correction table to create the publication-reference objects. Downstream manuscript panels use the publication-reference objects.

## Analysis unit

In the reconstructed core-level analyses, each of the 39 Xenium tissue cores/specimens is treated as an independent sample. The metadata field `patient_id` functions operationally as the tissue-core/specimen identifier in these objects; paired-looking core IDs are not merged.

## Outputs

All generated outputs live under `results/`. No production script in the clean flow should create a new top-level `figures/` directory.
