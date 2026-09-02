# Data inputs

## Required for the main pipeline

Place in `data/processed/`:

- `merged_obj_filtered.rds` — filtered/integrated Seurat object with `RNA`, `patient_id`, `NF1`, and final numeric `seurat_clusters`; retain the final UMAP if available.
- `xenium_list.rds` — list of per-core spatial Seurat objects with intact Xenium coordinates/FOVs and `patient_id`/`NF1` metadata.

These files are large and intentionally gitignored.

## Treatment-response metadata

For Figure 2H-I and response-dependent EGFR panels, provide a de-identified table at:

`data/metadata/patient_treatment_response_annotation.csv`

Required columns: `patient_id`, `NF1`, `treatment`, `response`, `Site`. Use the template and do not commit PHI.

## Core-coordinate provenance

TMA cores were manually outlined in Xenium Explorer using its annotation tools. The region/core boundary coordinates were exported as CSV files and used to define the individual cores. De-identified coordinate CSVs may be stored under `data/core_coordinates/` if approved for sharing.
